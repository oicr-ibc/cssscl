#!/usr/bin/env python

# Copyright 2015(c) The Ontario Institute for Cancer Research. All rights reserved.
#
# This program and the accompanying materials are made available under the
# terms of the GNU Public License v3.0.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division 
import pymongo 
import os.path, os, sys
import subprocess
import math
from database import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
import glob
import pickle
from functools import partial
import multiprocessing
from subroutines import *


def load_taxonomy_names(name_file):  
    '''
    Loads NCBI taxonomy names into a couple of dicts that are returned, 
    as a scientific name collection and a common name collection.
    '''

    logger.info("Loading taxonomy names: {0}".format(name_file))

    NamesFile = open(name_file)

    SciNameDict  = {}
    ComNameDict  = {}
    for line in NamesFile:
        taxID, Name, Uname, NameClass = line.split('\t|\t')
        taxID = int(taxID)
        if NameClass[0:3] == 'sci': # Process scientific name
            if Name == 'environmental samples':
                # recode as ENV-xxx with xxx being phylogrouping
                # example: environmental samples <Bacillariophyta> becomes ENV-Bacillariophyta
                groups = Uname.split('<')
                if len(groups) == 2:
                    Name = 'ENV-' + groups[1][:-1]
                else:
                    Name = 'ENV-Undefined'
            SciNameDict[taxID] = Name
        elif NameClass[0:11] == 'genbank com': # Process Genbank common name
            ComNameDict[taxID] = Name  

    return SciNameDict, ComNameDict



def load_taxonomy_data(name_file, node_file):
    '''
    Loads NCBI taxonomy node data into MongoDB in the taxa collection.
    '''

    # drop indexes from the taxa collection 
    #db.taxa.dropIndex( 'parent' )

    # drop the taxa collection
    db.taxa.drop()

    SciNameDict, ComNameDict = load_taxonomy_names(name_file)

    logger.info("Loading taxonomy nodes: {0}".format(node_file))

    NodesFile = open(node_file)

    for line in NodesFile:
        fields  = line.split('\t|\t')
        taxID   = int(fields[0])
        if int(fields[1]) == 1 and fields[2] == 'no rank':
            first = 0
            #first = None
        else:
            first = int(fields[1])
        SciName = SciNameDict.get(taxID, "<%d>" % taxID)
        NCBIdata = {
            "_id": taxID,  
            "parent"  : first,  # Parent taxID
            "rank"    : fields[2], # Taxonomic rank name
            "sciName" : SciName,                # Scientific name
            "comName" : ComNameDict.get(taxID, SciName) # Genbank common name
            }
        db.taxa.update({'_id': taxID},NCBIdata,True)

    logger.info("Adding parent index")
    db.taxa.ensure_index('parent')



def load_taxonomy_genomes(fasta_file, genome_file):

    logger.info("Get gi identifiers for the training set")
    
    gi_ids=[]

    handle_fasta = open(fasta_file, "rU")

    gi_file = os.path.dirname(os.path.abspath(fasta_file)) + '/gis.txt'

    handle_out_gi = open(gi_file,'w')

    for record in SeqIO.parse(handle_fasta, "fasta"):
        if record.id.split("|")[0] == 'gi':
            handle_out_gi.write("%s\n" % (int(record.id.split("|")[1])))   
        else:
            logger.error('There is no gi for this entry please check the identifier of the sequence in the training FASTA file provided')
            sys.exit(1)
 
    handle_out_gi.close()
    handle_fasta.close()  

    temp = os.path.dirname(os.path.abspath(fasta_file))

    outfile_gi_sorted = open(gi_file+'.sorted','w')
    p = subprocess.Popen(['sort', gi_file, '-t\t', '-d', '--key=1,1', '-T', temp], stdout=subprocess.PIPE)
    m = subprocess.Popen(['uniq'], stdin=p.stdout, stdout=outfile_gi_sorted)
    wait = p.wait()
    outfile_gi_sorted.close()
    if wait != 0:
        logger.error('Was not able to sort the gi file. Exiting')
        sys.exit(1)


    if subprocess.call(['ls', genome_file+'.sorted'], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
        logger.info("Sorting {0}".format(genome_file))
        if subprocess.call(['sort', genome_file, '-t\t', '-d', '--key=1,1', '-T', temp, '-o', genome_file+'.sorted']) != 0:
            logger.error("Error sorting {0}".format(genome_file))
            sys.exit(1)
    

    proc = subprocess.Popen(['join', '-t\t', '-1', '1', '-2', '1', gi_file+'.sorted', genome_file+'.sorted'], stdout=subprocess.PIPE)
    if not proc:
        logger.error('Could not perform the join. Exiting')
    out_join=proc.communicate()[0]
    
    gi_taxons = out_join.split('\n')[:-1]

    logger.info("Loading genome identifiers")

    ## drop the gitaxid collection first only if --btax True
    if build_tax_db:
        db.gitaxid.drop()

    gi_in_db = []    

    for gi_tax in gi_taxons:
        gi = int(gi_tax.split('\t')[0])
        gi_in_db.append(gi)
        taxid = int(gi_tax.split('\t')[1])
        db.gitaxid.update({'_id': gi}, {'$set' : {'taxid': taxid}}, True) 

    logger.info("Finished Loading genome identifiers")

    logger.info("Outputting sequences in the training set with gi entries that have taxon ids in the db.")

    cont_gi_fasta = 0
    count_gi_fasta_left = 0

    handle_out_fasta=open(os.path.dirname(os.path.abspath(fasta_file)) + '/train.fasta', "w")
    handle = open(fasta_file, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        cont_gi_fasta += 1
        if int(record.id.split("|")[1]) in gi_in_db:
            count_gi_fasta_left += 1 
            SeqIO.write(record, handle_out_fasta, "fasta")
    handle.close()
    handle_out_fasta.close()

    subprocess.call(['rm', gi_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.call(['rm', gi_file+'.sorted'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    logger.info("Total number of sequences in the training set: {0} total number of sequences kept that have taxon ids: {1}".format(str(cont_gi_fasta), str(count_gi_fasta_left)))



def load_taxonomy(fasta_file, directory):
    '''
    Loads NCBI taxonomy data files into MongoDB in the taxa collection.
    '''

    logger.info("Build the taxonomy db")

    if build_tax_db:
        name_file = os.path.join(directory, 'names.dmp')
        node_file = os.path.join(directory, 'nodes.dmp')
        load_taxonomy_data(name_file, node_file)
    else:
        logger.info("Skip loading taxon data.")

    genome_file = os.path.join(directory, 'gi_taxid_nucl.dmp')

    load_taxonomy_genomes(fasta_file, genome_file)




def build_training_blast_db(fasta_file, fasta_type):
    '''
    Build BLAST database 
    '''

    if not dbs_directory:
        blast_dir = os.path.dirname(os.path.abspath(fasta_file)) + '/blastdb/'
    else:
        blast_dir = os.path.abspath(dbs_directory) + '/blastdb/'

    logger.info("Building BLAST db in {0}".format(str(blast_dir)))

    subprocess.call(['rm','-rf',blast_dir]) 

    if subprocess.call(['mkdir', '-p', blast_dir]) != 0:
       logger.error("Failed to create the blastdb directory {0}".format(str(blast_dir)))  
       sys.exit(1)

    blast_input_file = blast_dir + 'train.fasta'

    train_fasta = os.path.dirname(os.path.abspath(fasta_file)) + '/train.fasta'

    if subprocess.call(['ln', '-s', train_fasta, blast_input_file]) != 0:
       logger.error("Failed to create a simlink to {0}".format(str(blast_input_file)))   

    if subprocess.call(['makeblastdb', '-dbtype', fasta_type, '-parse_seqids', '-in', blast_input_file]) != 0:
        logger.error("Error building BLAST DB for {0}".format(fasta_file))
        sys.exit(1)
  



def build_compression_db(fasta_file):
    '''
    Build the compression database 
    '''

    if not dbs_directory:
        compress_db_dir = os.path.dirname(os.path.abspath(fasta_file)) + '/compressdb/'
    else:
        compress_db_dir = os.path.abspath(dbs_directory) + '/compressdb/'

    train_fasta_cb = os.path.dirname(os.path.abspath(fasta_file)) + '/train.fasta'

    max_cpus = multiprocessing.cpu_count()

    plzip_threads = max_cpus
    
    subprocess.call(['rm','-rf', compress_db_dir]) 

    if subprocess.call(['mkdir', '-p', compress_db_dir]) != 0:
       logger.error("Failed to create the compression db directory {0}".format(str(blast_dir)))  
       sys.exit(1)

    subprocess.call(['ln', '-s', train_fasta_cb, compress_db_dir + 'train.fasta'])

    # compress data in the training set 
    handle = open(train_fasta_cb, "rU")
    seqs = SeqIO.parse(handle, "fasta")
    partial_compress_data_train_mp = partial(compress_data_train_mp, dir_out=compress_db_dir, plzip_threads=plzip_threads)
    pool = multiprocessing.Pool(num_threads)
    train_compressed = dict(pool.map(partial_compress_data_train_mp, seqs))
    pool.close()
    pool.join()
    handle.close()
    pickle.dump(train_compressed, open(compress_db_dir + 'train_compressed.p', "wb" ))




def get_kmer_counts_jellyfish_mp(rec, kmer, counts_dir, blast_db_dir, file_size, seq_len_hash):
    '''
    Calculate kmer counts using Jellyfish
    '''
    max_cpus = multiprocessing.cpu_count()
    jellyfish_threads = max_cpus - num_threads
    if jellyfish_threads <= 0:
       jellyfish_threads = int(1)
    if file_size == 'small':
        gi = int(rec.id.split("|")[1])
    elif file_size == 'large':
        gi = rec
    jellyfish_out = counts_dir + '/jellyfish/' + str(gi) + '/'
    jellyfish_temp_out = jellyfish_out + 'temp/'
    jellyfish_exec = os.popen("which %s" % "jellyfish").read().strip()
    if len(jellyfish_exec) == 0:
        logger.error("Could not find jellyfish. Exiting")
        sys.exit(1)
    if subprocess.call(['mkdir', '-p', jellyfish_out]) != 0:
        logger.error("Failed to create temp counts directory for jellyfish")
        sys.exit(1)        
    if subprocess.call(['mkdir', '-p', jellyfish_temp_out]) != 0:
        logger.error("Failed to create temp counts directory for jellyfish")
        sys.exit(1) 
    # make fasta file for jellyfish 
    if file_size == 'small':
        seq_len = len(rec.seq)
        handle_out_fasta=open(jellyfish_temp_out + '/fasta_out.fa', "w")
        SeqIO.write(rec, handle_out_fasta, "fasta")
        handle_out_fasta.close()
    elif file_size == 'large':
        # using BLAST
        seq_len = seq_len_hash[rec]
        jelly_file_out = jellyfish_temp_out + '/fasta_out.fa'
        train_blast_db = blast_db_dir + 'train.fasta'
        if subprocess.call(['blastdbcmd','-db', train_blast_db, '-entry', str(rec), '-out', jelly_file_out]) != 0:
            logger.error("Failed to output sequences from the blastdb. Exiting.")  
            sys.exit(1)           
    # run jellyfish 
    if subprocess.call(['jellyfish', 'count', '-C', '-t', str(jellyfish_threads), '-m', str(kmer), '-o', jellyfish_out + '/output', '-c', str(3), '-s',  str(seq_len), jellyfish_temp_out + '/fasta_out.fa' ]) != 0:
        logger.error("Failed to run jellyfish")  
        sys.exit(1)
    proc = glob.glob(jellyfish_out + '/output_*')
    nmb_files = len(proc)
    if nmb_files > 1:
        merge_out = os.system(jellyfish_exec + ' merge -o ' + jellyfish_out + 'output.jf' + ' ' + jellyfish_out + 'output\_*')
        if merge_out != 0:
            logger.error("Failed to run jellyfish merge {0}".format(str(merge_out)))
            sys.exit(1)
        proc = subprocess.Popen([jellyfish_exec, 'dump', jellyfish_out  + '/' + 'output.jf'], stdout=subprocess.PIPE)
        out=proc.communicate()
        os.system('rm -rf ' + jellyfish_out)
    elif nmb_files == 1:
        proc = subprocess.Popen([jellyfish_exec, 'dump', jellyfish_out + 'output_0'], stdout=subprocess.PIPE)
        out=proc.communicate()
        os.system('rm -rf ' + jellyfish_out)
    res_jellyfish = out[0].split("\n")
    if res_jellyfish is not None:
       dummy_dict_counts={}
       counts_kmers = res_jellyfish
       counts = [int(counts_kmers[i].split(">")[1]) for i in range(0,(len(counts_kmers) - 1),2)]
       kmers = [str(counts_kmers[i]) for i in range(1,(len(counts_kmers) - 1),2)]
       for i,kmer in enumerate(kmers):
            dummy_dict_counts[kmer]=counts[i]
    else:
       pass
    return (gi, dummy_dict_counts)





def build_counts_mp(fasta_file, fasta_type, kmer_min):
    '''
    Build count data using jellyfish
    '''

    logger.info("Building training count data using jellyfish for {0}".format(str(fasta_file)))

    if not dbs_directory:
        counts_dir = os.path.dirname(os.path.abspath(fasta_file)) + '/countsdb/'
        blast_db_dir = os.path.dirname(os.path.abspath(fasta_file)) + '/blastdb/'
    else:
        counts_dir = os.path.abspath(dbs_directory) + '/countsdb/'
        blast_db_dir = os.path.abspath(dbs_directory) + '/blastdb/'                

    subprocess.call(['rm','-rf',counts_dir]) 

    if subprocess.call(['mkdir', '-p', counts_dir]) != 0:
       logger.error("Failed to create the counts directory {0}".format(str(counts_dir)))  
       sys.exit(1)
       
    count_input_file = counts_dir + 'train.fasta'  

    train_fasta = os.path.dirname(os.path.abspath(fasta_file)) + '/train.fasta'

    if subprocess.call(['ln', '-s', train_fasta, count_input_file]) != 0:
       logger.error("Failed to create a simlink to {0}".format(str(count_input_file)))   

    file_size_bytes = os.path.getsize(train_fasta)   

    if file_size_bytes/1000000000 > 5.0:
        file_size = 'large'
        if not build_blast_db:
            logger.info("The size of the training set exceeds 5GB. Will need to build a BLAST DB for this datsets for indexing the genomes in the training set")
            build_training_blast_db(train_fasta, fasta_type)
    else:
        file_size = 'small' 
    
       
    # get gi identifiers in the fasta file
    gi_ids=[]
    seq_len=[]
    seq_len_hash={}
    handle = open(train_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        if record.id.split("|")[0] == 'gi':
            gi_ids.append(int(record.id.split("|")[1]))
            seq_len.append(int(len(record.seq)))
            seq_len_hash[int(record.id.split("|")[1])] = int(len(record.seq))
        else:
            logger.error('There is no gi for this entry please check the identifier of the sequence in the training FASTA file provided')
            sys.exit(1)
    handle.close()  

    # set the upper bound for k-mer size to 6 
    if fasta_type == 'nucl':
        max_kmer=round(math.log(min(seq_len),4))
        if max_kmer > 6:
            max_kmer = 6
    if fasta_type == 'prot':
        max_kmer=round(math.log(min(seq_len),20))
        if max_kmer > 6:
            max_kmer = 6
     
    pool=multiprocessing.Pool(num_threads)

    dict_kmer_counts = {}

    # calculate the k-mer counts for each sequence           
    for kmer_size in range(int(kmer_min),int(max_kmer)+1):
        logger.info("Getting counts for k-mer: {0}".format(str(kmer_size)))
        dict_kmer_counts[kmer_size] = {}
        if file_size == 'small':
            handle = open(train_fasta, "rU")
            records = SeqIO.parse(handle, "fasta")
            partial_get_kmer_counts_jellyfish_mp = partial(get_kmer_counts_jellyfish_mp, kmer = kmer_size, counts_dir = counts_dir, blast_db_dir = blast_db_dir, file_size = file_size, seq_len_hash = seq_len_hash)
            dict_kmer_counts[kmer_size] = dict(pool.map(partial_get_kmer_counts_jellyfish_mp, records))
            handle.close()
        elif file_size == 'large':
            records = gi_ids
            partial_get_kmer_counts_jellyfish_mp = partial(get_kmer_counts_jellyfish_mp, kmer = kmer_size, counts_dir = counts_dir, blast_db_dir = blast_db_dir, file_size = file_size, seq_len_hash = seq_len_hash)
            dict_kmer_counts[kmer_size] = dict(pool.map(partial_get_kmer_counts_jellyfish_mp, records))
 
    pool.close()
    pool.join()

    for kmer_size in range(int(kmer_min),int(max_kmer)+1):     
        pickle.dump(dict_kmer_counts[kmer_size], open( counts_dir + "/dict_kmer_" + str(kmer_size) + "_counts.p", "wb" ))

    subprocess.call(['rm', '-r', counts_dir + '/jellyfish/'])    
    logger.info("Finished building training count data.")



def check_program_in_path(program):
    #if subprocess.call(['which', program], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
    #       logger.error("{0} not in the path. Exiting.".format(program))
    #       sys.exit(1)
    executable = os.popen("type %s" % str(program)).read().strip('\n')
    if executable == program + ': not found':
           logger.error("{0} not in the path. Exiting.".format(program))
           sys.exit(1)
    

def main(args):
    '''Build dbs (blast, taxonomy and counts) for cssscl'''

    global db, logger, build_tax_db, dbs_directory, num_threads, build_blast_db, TAXON_RANK     

    logger = args.logging.getLogger(__name__)

    db = connect(args)
  
    cnt_gi = db.gitaxid.count()
    cnt_taxa = db.taxa.count()

    if cnt_gi == 0 and cnt_taxa == 0:
        db_present = False
    else:
        db_present = True


    if args.build_taxonomy_data:
        build_tax_db = True
    else:
        build_tax_db = False


    if args.disable_kmer_count:
        use_counts = False
    else:
        use_counts = True


    if not use_counts and not args.use_blast:
       logger.error("kmer count based measures can only be disabled if BLAST is specified.")
       sys.exit(1) 


    dbs_directory = args.dbs_directory     


    if dbs_directory:
        logger.info("dbs will be outtputed to {0}".format(str(dbs_directory)))  
    else:
        logger.info("dbs will be outtputed to {0}".format(str(os.path.dirname(os.path.abspath(args.file_fasta)))))  

    #if args.number_threads:   
    #    num_threads = int(args.number_threads)

    num_threads = int(args.number_threads)

    max_cpus = multiprocessing.cpu_count()

    if num_threads > max_cpus:
        logger.info("Warning. You have specified {0} processes and the machine has a maximum of {1} processors.".format(num_threads, max_cpus))
    else:
        logger.info("Number of CPUs specified: {0}".format(num_threads))

    if args.use_blast:
        check_program_in_path('blastn')
    if args.use_compression:
        check_program_in_path('plzip')
    if use_counts:
        check_program_in_path('jellyfish')

    # check if blastdb will need to be built     
    build_blast_db = args.use_blast    

    if not args.build_taxonomy_data and not db_present:
        logger.error("No entries detected in the db. You first need to build the databases for the training set using the option -btax. Exiting.")
        sys.exit(1)
    elif db_present and not args.build_taxonomy_data:
        logger.info("Detected a total of {0} gi entries in the db, total number of taxa entries in the db: {1}.".format(cnt_gi, cnt_taxa))  
        logger.info("Will reload the gi entries from the training set into the db.")
        load_taxonomy(args.file_fasta, args.taxon_directory)
        # REMOVE ALL THE cl_train*.p that might exist
        subprocess.call(['rm', os.path.abspath(dbs_directory) + '/cl_train_*.p'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #sys.exit(1)
    elif args.build_taxonomy_data:
        logger.info("Building taxon db and loading entries from the training set into db.")
        load_taxonomy(args.file_fasta, args.taxon_directory)
        # REMOVE ALL THE cl_train*.p that might exist
        subprocess.call(['rm', os.path.abspath(dbs_directory) + '/cl_train_*.p'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
 
    if args.use_blast: 
        build_training_blast_db(args.file_fasta, args.fastatype)


    if args.use_compression: 
        logger.info("Creating the compression db.")
        build_compression_db(args.file_fasta)
        logger.info("Finished creating the compression db.")


    if not use_counts and args.use_blast and not args.use_compression:
        logger.info("CSSSCL will be run with BLAST only, kmer and compresson based measures have been turned off.")
    elif not use_counts and args.use_blast and args.use_compression:
        logger.info("CSSSCL will be run with BLAST and compression measures only, kmer measure has been turned off.")
    else:
        build_counts_mp(args.file_fasta, args.fastatype, args.kmer_min)



if __name__ == '__main__':
    print 'This program should be run as part of the cssscl package:\n\t$ cssscl build_dbs -h\n\tor\n\t$ /path/to/cssscl/bin/cssscl configure -h'
