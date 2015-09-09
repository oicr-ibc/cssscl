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
from Bio import SearchIO
from Bio.Seq import Seq
import glob
import pickle
import random
import numpy as np
from numpy import *
from subroutines import *
from sklearn import cross_validation
import multiprocessing
from multiprocessing import Value, Lock
from multiprocessing import Manager
import multiprocessing.pool
from functools import partial
from random import randint
from collections import defaultdict
import timeit
from datetime import datetime, timedelta

# It turns out that certain Python modules 
# (numpy, scipy, tables, pandas, skimage...) 
# mess with core affinity on import. 
# A workaround is to reset the task affinity 
# after importing these modules see below  
# this problem seems to be specifically caused by them linking against multithreaded OpenBLAS libraries.
# reset the task affinity 
mask = hex(int(''.join(multiprocessing.cpu_count() * ['1']),2))
os.system("taskset -p " + mask + " %d" % os.getpid())



class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False
    def _set_daemon(self, value):
        pass
    daemon = property(_get_daemon, _set_daemon)



class MyPool(multiprocessing.pool.Pool):
    Process = NoDaemonProcess


def get_procs():
    if num_threads < 2:
        num_procs_c = int(1)
        num_procs_o = int(1)
    elif num_threads >= 2:
        num_procs_c = int(1)
        num_procs_o = int(num_threads)
    return (num_procs_c, num_procs_o)



def blast_data(fasta_file, dbdirectory, fasta_type, which):
    '''
    BLAST file against a BLAST database 
    '''  

    dir_basename_test = os.path.dirname(os.path.abspath(fasta_file))

    blastdb = dbdirectory + '/blastdb/'
    
    blast_out_train = dbdirectory + '/blast_train/'

    blast_out_test = dir_basename_test + '/blast_test/'
   
    if which == 'test':
        logger.info("Found test set: {0}".format(str(fasta_file)))
    elif which == 'train': 
        logger.info("Found train set: {0}".format(str(fasta_file)))

    if which == 'test':
        logger.info("BLAST sequences in the test set against the BLAST database in {0}".format(blastdb))
    elif which == 'train':
        logger.info("BLAST sequences in the training set against the BLAST database in {0}".format(blastdb))

    if subprocess.call(['mkdir', '-p', blast_out_train]) != 0:
       logger.error("Failed to create the output directory for BLAST {0}. Exiting.".format(str(blast_out_train)))  
       sys.exit(1)

    if subprocess.call(['mkdir', '-p', blast_out_test]) != 0:
       logger.error("Failed to create the output directory for BLAST {0}. Exiting.".format(str(blast_out_test)))  
       sys.exit(1)
    
    if fasta_type == 'nucl':
        try:
            #blast_db_file = os.path.basename(glob.glob(blastdb + '*.nhr')[0]).split('.nhr')[0]
            os.path.basename(glob.glob(blastdb + '*.nhr')[0]).split('.nhr')[0]
            blast_db_file = 'train.fasta'
        except IndexError:
            logger.error("Nucleotide BLAST database not found. Exiting.") 
            sys.exit(1)
    elif fasta_type == 'prot':
        try:
            #blast_db_file = os.path.basename(glob.glob(blastdb + '*.phr')[0]).split('.phr')[0]
            os.path.basename(glob.glob(blastdb + '*.phr')[0]).split('.phr')[0]
            blast_db_file = 'train.fasta'
        except IndexError:   
            logger.error("Protein BLAST database not found. Exiting.") 
            sys.exit(1)
     
    if which == 'test':
        if blast_type == 'blastn':
            if subprocess.call(['blastn', '-task', 'blastn', '-num_alignments', str(num_alignments), '-word_size', str(11), '-num_threads', str(num_threads), '-outfmt', str(6), '-db', blastdb + blast_db_file, '-query', fasta_file, '-out', blast_out_test + 'blast_test.txt']) != 0:
                logger.error("Error running BLAST for {0}. Exiting.".format(fasta_file))
                sys.exit(1)
            else:
                # index the blast output file 
                idx_filename = blast_out_test + 'blast_test.idx'

                if subprocess.call(['ls', idx_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0:
                    subprocess.call(['rm', idx_filename])

                logger.debug("Indexing the BLAST aligned output file {0}".format(idx_filename))
                SearchIO.index_db(idx_filename, blast_out_test + 'blast_test.txt', 'blast-tab')
               
        elif blast_type == 'megablast':
            if subprocess.call(['blastn', '-task', 'megablast', '-num_alignments', str(num_alignments), '-num_threads', str(num_threads), '-outfmt', str(6), '-db', blastdb + blast_db_file, '-query', fasta_file, '-out', blast_out_test + 'blast_test.txt']) != 0:
                logger.error("Error running BLAST for {0}. Exiting.".format(fasta_file))
                sys.exit(1)
            else:
                # index the blast output file 
                idx_filename = blast_out_test + 'blast_test.idx'           
            
                if subprocess.call(['ls', idx_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0:
                    subprocess.call(['rm', idx_filename])

                logger.debug("Indexing the BLAST aligned output file {0}".format(idx_filename))
                SearchIO.index_db(idx_filename, blast_out_test + 'blast_test.txt', 'blast-tab')
            
    if which == 'train':
        if blast_type == 'blastn':
            if subprocess.call(['blastn', '-task', 'blastn', '-num_alignments', str(10), '-word_size', str(11), '-num_threads', str(num_threads), '-outfmt', str(6), '-db', blastdb + blast_db_file, '-query', fasta_file, '-out', blast_out_train + 'blast_train.txt']) != 0:
                logger.error("Error running BLAST for {0}. Exiting.".format(fasta_file))
                sys.exit(1)
        elif blast_type == 'megablast':
            if subprocess.call(['blastn', '-task', 'megablast', '-num_alignments', str(10), '-num_threads', str(num_threads), '-outfmt', str(6), '-db', blastdb + blast_db_file, '-query', fasta_file, '-out', blast_out_train + 'blast_train.txt']) != 0:
                logger.error("Error running BLAST for {0}. Exiting.".format(fasta_file))
                sys.exit(1)
            


        
def build_blast_db(fasta_file, fasta_type):
    '''
    Build BLAST database 
    '''

    blast_dir = os.path.dirname(os.path.abspath(fasta_file)) + '/blastdb/'

    subprocess.call(['rm','-rf',blast_dir]) 

    if subprocess.call(['mkdir', '-p', blast_dir]) != 0:
       logger.error("Failed to create the blastdb directory {0}".format(str(blast_dir)))  
       sys.exit(1)

    blast_input_file = blast_dir + 'train.fasta'

    train_fasta = os.path.dirname(os.path.abspath(fasta_file)) + '/optimized_train.fasta'

    if subprocess.call(['ln', '-s', train_fasta, blast_input_file]) != 0:
       logger.error("Failed to create a simlink to {0}".format(str(blast_input_file)))   

   
    if subprocess.call(['makeblastdb', '-dbtype', fasta_type, '-parse_seqids', '-in', blast_input_file]) != 0:
        logger.error("Error building BLAST DB for {0}".format(fasta_file))
        sys.exit(1)
  



def cb_loocv_mp(loo_iter, dist_cb, classes_array):
    train_index, test_index = list(loo_iter)
    dist_cb_sub = dist_cb[test_index,:][:,train_index]
    accuracy, ftest_score = get_accuracy_scores_tune_knn_and_predict_loo(dist_cb_sub,classes_array[train_index],classes_array[test_index],K=1)
    return accuracy, ftest_score



def evaluate_compression_train_data(dbdirectory):
    '''
    Compress data using plzip and evaluate its classification performance 
    '''

    # split process/threads equally between plzip and the dist calculation 
    if num_threads > 16:
            plzip_threads = int(16)
    else:
            plzip_threads = int(16)
    
    num_procs_c = int(num_threads)

    optimized_kmer_dir = os.path.abspath(dbdirectory) + '/optimum_kmer/'

    train_fasta_cb = dbdirectory + '/compressdb/train.fasta'

    train_fasta_cb_opt = optimized_kmer_dir + 'optimized_train.fasta'

    #cl_train = pickle.load(open(os.path.abspath(dbdirectory) + '/cl_train.p', 'rb'))

    cl_train = pickle.load(open(optimized_kmer_dir + '/cl_train_opt_' + TAXON_RANK + '.p', 'rb'))

    compress_db_dir = os.path.abspath(dbdirectory) + '/compressdb/'

    train_data_compressed = pickle.load(open(compress_db_dir + 'train_compressed.p', 'rb'))

    if subprocess.call(['ls', compress_db_dir], stdout=subprocess.PIPE) != 0:
       logger.error("Can not find the compression db {0}. Exiting.".format(compress_db_dir))
       sys.exit(1)

    if subprocess.call(['ls', train_fasta_cb], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
       logger.error("Train fasta file in {0} not found. Exiting.".format(compress_db_dir))
       sys.exit(1)
   
    handle = open(train_fasta_cb, "rU")
    # get the size of the fasta file
    nseqs = sum(1 for _ in SeqIO.parse(handle, 'fasta'))
    handle.close()

    dist_cb = np.zeros((nseqs,nseqs))

    gis_cb = []

    pool=multiprocessing.Pool(num_procs_c)

    # find compression distances 
    # open the simulated test set fasta file 
    handle = open(train_fasta_cb_opt, 'rU')
    for i, recseq in enumerate(SeqIO.parse(handle, 'fasta')):
        gis_cb.append(int(recseq.id.split("|")[1]))
        # open the train set fasta file
        handle2 = open(train_fasta_cb, 'rU') 
        seqs = SeqIO.parse(handle2, 'fasta')
        partial_get_ncd_distance_mp_train = partial(get_ncd_distance_mp_train, seq1=recseq, plzip_threads=plzip_threads, dict_train = train_data_compressed) 
        dist_compression = pool.map(partial_get_ncd_distance_mp_train, seqs)
        dist_cb[i, ] = dist_compression
        handle2.close()

    pool.close()
    pool.join()
    handle.close()

    order_gis = [item for sublist in cl_train.values() for item in sublist] 

    indx = []
    for gi_ord in order_gis:
        indx.append(gis_cb.index(gi_ord))

    indx = array(indx)

    dist_cb = dist_cb[indx,:][:,indx]

    pickle.dump(indx, open(compress_db_dir + 'indx_train_opt_' +  TAXON_RANK + '.p', 'wb' ))
    pickle.dump(dist_cb, open(compress_db_dir + 'dist_cb_train_opt_' + TAXON_RANK + '.p', 'wb' ))      

    classes_array = []
    for taxon_name in cl_train.keys():
        for gi in cl_train[taxon_name]:
            classes_array.append(taxon_name)

    classes_array = array(classes_array)

    dist_accuracies_train = {}
    
    LOO = cross_validation.LeaveOneOut(len(classes_array)) 
    results_cb_train = {}
    results_cb_train_accu_loo = [] 

    pool=multiprocessing.Pool(num_procs_c)

    partial_cb_loocv_mp = partial(cb_loocv_mp, dist_cb = dist_cb, classes_array = classes_array)
    results_cb = pool.map(partial_cb_loocv_mp, LOO)
    pool.close()
    pool.join()

    for res in results_cb:
        results_cb_train_accu_loo.append(res[0])

    dist_accuracies_train['d_cb'] = [array(results_cb_train_accu_loo).mean(), TAXON_RANK]    
    pickle.dump(dist_accuracies_train, open(optimized_kmer_dir + 'cb_train_results_' + TAXON_RANK + '.p', 'wb' ))    
    return dist_accuracies_train




def get_compression_sequences_sequential_mp(record, dbdirectory, num_procs_c):
    '''
    Get compression dist to the training set  
    '''
    
    dist_compression = []
    dist = []
    dist_cb = []

    if num_procs_c > 16: 
        plzip_threads = int(16)
    else:
        plzip_threads = int(16)

    pool = multiprocessing.Pool(num_procs_c)

    #optimized_kmer_dir = os.path.abspath(dbdirectory) + '/optimum_kmer/'
    #train_kmer_dir = os.path.abspath(dbdirectory) + '/countsdb/' 

    cl_train = pickle.load(open(os.path.abspath(dbdirectory) + '/cl_train_' + TAXON_RANK + '.p', 'rb'))

    compress_db_dir = dbdirectory + '/compressdb/'

    train_fasta_cb = compress_db_dir + 'train.fasta'

    if subprocess.call(['ls', compress_db_dir + 'indx_train_' + TAXON_RANK + '.p'], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
        gis_cb = []
        handle = open(train_fasta_cb, "rU")
        for seq in SeqIO.parse(handle, "fasta"):
            gis_cb.append(int(seq.id.split("|")[1]))
        handle.close()
        order_gis = [item for sublist in cl_train.values() for item in sublist]
        indx = []
        for gi_ord in order_gis:
            indx.append(gis_cb.index(gi_ord))      
        pickle.dump(indx, open(compress_db_dir + 'indx_train_' + TAXON_RANK + '.p', 'wb' ))
    else:
        indx = pickle.load(open(compress_db_dir + 'indx_train_' + TAXON_RANK + '.p', 'rb'))    
        
    # load dict with the training set 

    train_data_compressed = pickle.load(open(compress_db_dir + 'train_compressed.p', 'rb'))

    if subprocess.call(['ls', compress_db_dir], stdout=subprocess.PIPE) != 0:
       logger.error("Can not find the compression db {0}. Exiting.".format(compress_db_dir))
       sys.exit(1)
  
    if subprocess.call(['ls', train_fasta_cb], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
       logger.error("Train fasta file in {0} not found. Exiting.".format(compress_db_dir))
       sys.exit(1)

    handle = open(train_fasta_cb, "rU") 
    seqs = SeqIO.parse(handle, "fasta")
    partial_get_ncd_distance_mp_test = partial(get_ncd_distance_mp_test, seq1=record, plzip_threads = plzip_threads, dict_train = train_data_compressed)
    dist_compression = pool.map(partial_get_ncd_distance_mp_test, seqs)
    pool.close()
    pool.join()
    handle.close()
    dist = array(dist_compression)[array(indx)]
    dist_cb = array([dist])
    return dist_cb



def get_blast_hits_sequential(seqid, dbdirectory, gis_blast):
    '''
    '''

    blast_idx = dbdirectory + '/blast_test/blast_test.idx'
    
    if subprocess.call(['ls', blast_idx], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
       logger.error("BLAST index for the BLAST output file not found {0}".format(blast_idx))
       sys.exit(1)

    try:
        db_idx = SearchIO.index_db(blast_idx)
        query_results = db_idx[seqid]
        blast_scores = np.zeros((1,len(gis_blast)))
        for qresult in query_results:
            for hit in qresult:
                if int(hit.hit_id.split('|')[1]) in gis_blast and blast_scores[0,gis_blast.index(int(hit.hit_id.split('|')[1]))] <= float(hit.bitscore):
                    blast_scores[0,gis_blast.index(int(hit.hit_id.split('|')[1]))] = float(hit.bitscore)
        max_val = np.apply_along_axis(return_max_in_array, axis=1, arr=blast_scores)
        # for those entries that have 0.0 max_val
        max_val[max_val == 0.0] = 1.0
        dist_blast = 1 - blast_scores/max_val[:,None] # this is like a transpose of max_val
        db_idx.close()
    except KeyError:
        dist_blast = None

    return dist_blast 



def get_kmer_counts_jellyfish_mp_2(rec, kmer, counts_dir):
    '''
    Calculate kmer counts using Jellyfish
    '''
    max_cpus = multiprocessing.cpu_count()
    jellyfish_threads = max_cpus - num_threads
    if jellyfish_threads <= 0:
       jellyfish_threads = int(1) 
    gi = int(rec.id.split("|")[1])
    jellyfish_out = counts_dir + '/jellyfish/' + str(gi) + '/'
    jellyfish_temp_out = jellyfish_out + 'temp/'
    #jellyfish_exec = os.popen("which %s" % "jellyfish").read().strip()
    jellyfish_exec = os.popen("type %s" % "jellyfish").read().strip('jellyfish is ').strip('\n')
    #if len(jellyfish_exec) == 0:
    if jellyfish_exec == 'jellyfish: not found':
        logger.error("Could not find jellyfish. Exiting")
        sys.exit(1)
    if subprocess.call(['mkdir', '-p', jellyfish_out]) != 0:
        logger.error("Failed to create temp counts directory for jellyfish")
        sys.exit(1)        
    if subprocess.call(['mkdir', '-p', jellyfish_temp_out]) != 0:
        logger.error("Failed to create temp counts directory for jellyfish")
        sys.exit(1) 
    seq_len = len(rec.seq)
    # make fasta file for jellyfish    
    handle_out_fasta=open(jellyfish_temp_out + '/fasta_out.fa', "w")
    SeqIO.write(rec, handle_out_fasta, "fasta")
    handle_out_fasta.close()
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





def get_kmer_counts_jellyfish_mp_2_blast(rec, kmer, counts_dir, seq_len_hash):
    '''
    Calculate kmer counts using Jellyfish
    '''
    max_cpus = multiprocessing.cpu_count()
    jellyfish_threads = max_cpus - num_threads
    if jellyfish_threads <= 0:
       jellyfish_threads = int(1) 
    gi = rec
    jellyfish_out = counts_dir + '/jellyfish/' + str(gi) + '/'
    jellyfish_temp_out = jellyfish_out + 'temp/'
    #jellyfish_exec = os.popen("which %s" % "jellyfish").read().strip()
    jellyfish_exec = os.popen("type %s" % "jellyfish").read().strip('jellyfish is ').strip('\n')
    #if len(jellyfish_exec) == 0:
    if jellyfish_exec == 'jellyfish: not found':
        logger.error("Could not find jellyfish. Exiting")
        sys.exit(1)
    if subprocess.call(['mkdir', '-p', jellyfish_out]) != 0:
        logger.error("Failed to create temp counts directory for jellyfish")
        sys.exit(1)        
    if subprocess.call(['mkdir', '-p', jellyfish_temp_out]) != 0:
        logger.error("Failed to create temp counts directory for jellyfish")
        sys.exit(1) 
    seq_len = seq_len_hash[rec]
    # make fasta file for jellyfish
    blast_db_dir = counts_dir + '/blastdb/'    
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







def get_kmer_counts_jellyfish_mp(rec, kmer, counts_dir):
    '''
    Calculate kmer counts using Jellyfish
    '''
    while True:
        rand = randint(0,1e+10)
        jellyfish_out = counts_dir + '/jellyfish/' + str(rand) + '/' 
        jellyfish_temp_out = jellyfish_out + 'temp/'
        if not os.path.isdir(jellyfish_out):
            break 
    max_cpus = multiprocessing.cpu_count()
    jellyfish_threads = max_cpus - num_threads
    if jellyfish_threads <= 0:
       jellyfish_threads = int(1)     
    #jellyfish_exec = os.popen("which %s" % "jellyfish").read().strip()
    jellyfish_exec = os.popen("type %s" % "jellyfish").read().strip('jellyfish is ').strip('\n')
    #if len(jellyfish_exec) == 0:
    if jellyfish_exec == 'jellyfish: not found':
        logger.error("Could not find jellyfish. Exiting")
        sys.exit(1)
    if subprocess.call(['mkdir', '-p', jellyfish_out]) != 0:
        logger.error("Failed to create temp counts directory for jellyfish")
        sys.exit(1)        
    if subprocess.call(['mkdir', '-p', jellyfish_temp_out]) != 0:
        logger.error("Failed to create temp counts directory for jellyfish")
        sys.exit(1) 
    seq_len = len(rec.seq)
    # make fasta file for jellyfish    
    handle_out_fasta=open(jellyfish_temp_out + '/fasta_out.fa', "w")
    SeqIO.write(rec, handle_out_fasta, "fasta")
    handle_out_fasta.close()
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
        ##if subprocess.call(['jellyfish', 'merge', '-o', jellyfish_out + 'output.jf', jellyfish_out + 'output\_*' ]) != 0:
        ##    logger.error("Failed to run jellyfish merge")  
        ##sys.exit(1)
        proc = subprocess.Popen([jellyfish_exec, 'dump', jellyfish_out  + '/' + 'output.jf'], stdout=subprocess.PIPE)
        out=proc.communicate()
        os.system('rm -rf ' + jellyfish_out)
    elif nmb_files == 1:
        proc = subprocess.Popen([jellyfish_exec, 'dump', jellyfish_out + 'output_0'], stdout=subprocess.PIPE)
        out=proc.communicate()
        os.system('rm -rf ' + jellyfish_out)
    return out[0].split("\n")



def get_stats_on_sequences_test_set(fasta_file, fasta_type, call_type):    
    seq_len=[]
    gis = []
    mean_seq_test = 0
    std_seq_test = 0
    handle = open(fasta_file, "rU")
    if call_type == 'test':
        count_test_seq = 0
        for record in SeqIO.parse(handle, "fasta"):
            seq_len.append(int(len(record.seq)))
            count_test_seq += 1
            if count_test_seq > 9999:
                break 
    if call_type == 'train':
        for record in SeqIO.parse(handle, "fasta"):
           gis.append(int(record.id.split("|")[1]))
           seq_len.append(int(len(record.seq)))
    handle.close()
    if call_type == 'train':        
        if fasta_type == 'nucl':
            max_kmer=round(math.log(min(seq_len),4))
            if max_kmer > 6:
                max_kmer = 6
        if fasta_type == 'prot':
            max_kmer=round(math.log(min(seq_len),20))
            if max_kmer > 6:
                max_kmer = 6
    # get the stats on the sequences in the test set
    if call_type == 'test':
        max_seq_len = max(seq_len)
        min_seq_len = min(seq_len)
        mean_seq_test = mean(seq_len)   
        std_seq_test = std(seq_len) 
        # set the upper bound for k-mer size to 6
    if call_type == 'train':
        return gis, max_kmer
    elif call_type == 'test':
        return mean_seq_test, std_seq_test, max_seq_len, min_seq_len 




def get_taxon_ids(gis):
    '''
    Get taxons from mongo db 
    '''
    gi_taxon={}
    for gi in gis:
        gi = int(gi)
        #print gi
        rank_taxid = []
        sciName = []
        taxon_id=[]
        try: 
            taxid = db.gitaxid.find({'_id': gi})[0]
        except IndexError:
            logger.error("The taxon id for the gi id: {0} is not available in the db, check or update your training fasta file or your taxon file".format(str(gi)))
            sys.exit(1)
        if db.taxa.find({"_id" : taxid['taxid']})[0]['rank'] != 'no rank':
            rank_taxid.append(db.taxa.find({"_id" : taxid['taxid']})[0]['rank'])
            sciName.append(db.taxa.find({"_id" : taxid['taxid']})[0]['sciName'])
            taxon_id.append(str(taxid['taxid']))
        else:
            rank_taxid.append('NO_VALUE')
            sciName.append('NO_VALUE') 
            taxon_id.append('NO_VALUE')
        parent = db.taxa.find({"_id" : taxid['taxid']})[0]['parent']
        while True:
            tax_info = db.taxa.find({"_id" : parent})[0]
            if tax_info['rank'] != 'no rank':
                rank_taxid.append(tax_info['rank'])   
                sciName.append(tax_info['sciName'])
                taxon_id.append(str(parent))
            else:
                rank_taxid.append('NO_VALUE')
                sciName.append('NO_VALUE')
                taxon_id.append('NO_VALUE')
            parent = tax_info['parent']
            if (parent == 1 or parent == 0 or parent == 2157 or parent == 2 or  parent == 629395 or parent == 2759 or parent == 12884 or parent == 10239 or parent == 28384 or parent == 12908):
                break
        gi_taxon[gi] = {'rank': rank_taxid, 'taxid': taxon_id}      
    return gi_taxon     
    



def format_data_test_count(file_fasta_test, data_kmer_test_counts, kmer_size):
    '''
    Format test data for classification
    '''
    data_set= {}
    counts_arrayS = []
    freq_array = []
    kmer_keys = None
    counts = pickle.load(open(data_kmer_test_counts + "/dict_kmer_" + str(kmer_size) + "_counts.p", "rb" ))
    # get gis in the same order as in the fasta test file
    handle = open(file_fasta_test, "rU")
    seq_id = []
    for record in SeqIO.parse(handle, "fasta"):
        seq_id.append(record.id)
    handle.close()
    for gi in seq_l:
             if kmer_keys is None:
                kmer_keys = counts[gi].keys()
                dummy_list = []
                for val in kmer_keys:
                    dummy_list.append(counts[gi][val])
                counts_array.append(dummy_list)
             else:
                dummy_list=[] 
                if kmer_keys == counts[gi].keys():
                    for val in kmer_keys:
                        dummy_list.append(counts[gi][val])
                else:
                    kmer_add = set(counts[gi].keys()) - set(kmer_keys)
                    if len(kmer_add) > 0:
                        kmer_keys = kmer_keys + list(kmer_add)
                    for val in kmer_keys:
                            if val in counts[gi]:
                                dummy_list.append(counts[gi][val])
                            else:
                                dummy_list.append(0)
                counts_array.append(dummy_list)
    for i in range(0,len(counts_array)):
        difference = len(kmer_keys) - len(counts_array[i])
        if difference > 0:
            counts_array[i] = counts_array[i] + difference*[0]
            #print i
        else:
            break
    for val in counts_array:
        freq_array.append(list(val/sum(val)))
    data_set['counts'] = counts_array
    data_set['freq'] = freq_array
    data_set['features'] = kmer_keys
    data_set['names'] = seq_id
    return data_set
        




def format_data_train_count(data_kmer_dir, cl_train, kmer_size, train_kmers = False):
    '''
    Format train data for classification
    '''
    data_set = {}
    counts_array = []
    classes_array = []
    kmer_keys = None
    for taxon_name in cl_train.keys():
        counts = pickle.load(open(data_kmer_dir + "/dict_kmer_" + str(kmer_size) + "_counts.p", "rb" ))
        for gi in cl_train[taxon_name]:
             if kmer_keys is None:
                kmer_keys = counts[gi].keys()
                dummy_list = []
                for val in kmer_keys:
                    dummy_list.append(counts[gi][val])
                counts_array.append(dummy_list)
                classes_array.append(taxon_name)
             else:
                dummy_list=[] 
                if kmer_keys == counts[gi].keys():
                    for val in kmer_keys:
                        dummy_list.append(counts[gi][val])
                else:
                    kmer_add = set(counts[gi].keys()) - set(kmer_keys)
                    if len(kmer_add) > 0:
                        kmer_keys = kmer_keys + list(kmer_add)
                    for val in kmer_keys:
                            if val in counts[gi]:
                                dummy_list.append(counts[gi][val])
                            else:
                                dummy_list.append(0)
                counts_array.append(dummy_list)
                classes_array.append(taxon_name)
    if not train_kmers:
        for i in range(0,len(counts_array)):
            difference = len(kmer_keys) - len(counts_array[i])
            if difference > 0:
                counts_array[i] = counts_array[i] + difference*[0]
    elif train_kmers:
        for i in range(0,len(counts_array)):
            difference = len(train_kmers) - len(counts_array[i])
            if difference > 0:
                counts_array[i] = counts_array[i] + difference*[0]
        counts_array_ord=[]
        for r in counts_array:
            dummy=[0]*(len(r))
            for i,k in enumerate(kmer_keys):
                indx = train_kmers.index(k)
                dummy[indx] = r[i]
            counts_array_ord.append(dummy)
    if not train_kmers:
        data_set['counts'] = counts_array
    elif train_kmers:
        data_set['counts'] = counts_array_ord
    data_set['target'] = classes_array
    data_set['features'] = kmer_keys    
    return data_set



def get_classes_in_train(gis, taxon_rank):
    ''' Get Taxons in the training set
    '''
    gi_taxon = get_taxon_ids(gis)
    cl_train = {}
    for gi in gi_taxon.keys():
        try:
            taxon_id = gi_taxon[gi]['taxid'][gi_taxon[gi]['rank'].index(taxon_rank)]
            cl_train[taxon_id] = []
        except ValueError:
            # gi does not contain taxon info skip to the next 
            continue 
    for gi in gi_taxon.keys():
        try:
            taxon_id = gi_taxon[gi]['taxid'][gi_taxon[gi]['rank'].index(taxon_rank)]
            cl_train[taxon_id].append(gi)
        except ValueError:
            # gi does not contain taxon info skip to the next 
            continue 
    for key in cl_train.keys():
        if len(cl_train[key]) < 2:
           cl_train.pop(key)
    if len(cl_train.keys()) < 1:
        logger.error("The training set does not contain enough data for classification, at least 2 data points per class are required")
        sys.exit(1)
    return cl_train



def kmer_loocv_mp(loo_iter, data_set_train, data_set_test):
    train_index, test_index = list(loo_iter)
    if sum(data_set_test['counts'][array(test_index)]) != 0:
        #accuracy, ftest_score = training_kmer_performance_loo(data_set_train, array(data_set_train['target']), data_set_test, array(train_index), array(test_index))
        accuracy, ftest_score = training_kmer_performance_loo(data_set_train, array(data_set_test['target']), data_set_test, array(train_index), array(test_index))
        return accuracy, ftest_score



def perform_kmer_optimization(dbdirectory, gis, kmer_range):
    '''
    Perform kmer optimization and distance (based on counts) selection
    '''
    pool=multiprocessing.Pool(num_threads)
    taxon_rank = TAXON_RANK
    test_kmer_dir = os.path.abspath(dbdirectory) + '/optimum_kmer/'
    train_kmer_dir = os.path.abspath(dbdirectory) + '/countsdb/' 
    results_kmer_optimum={}
    optimized_kmer_results={}
    cl_train_opt = pickle.load(open(test_kmer_dir + 'cl_train_opt_' + TAXON_RANK + '.p', 'rb' ))
    cl_train = pickle.load(open(os.path.abspath(dbdirectory) + '/cl_train_' + TAXON_RANK + '.p', 'rb' ))
    results_kmer_optimum = {}
    accuracy_kmer = []
    individual_accuracy = defaultdict(list)
    individual_accuracy_opt = {}
    for kmer_size in kmer_range:  
    	kmer = 'kmer' + str(kmer_size)
        logger.info("kmer value: {0}".format(str(kmer_size)))
        logger.info("Loading and pre-processing count data for the training set")
        data_set_train = format_data_train_count(train_kmer_dir, cl_train_opt, kmer_size)
        data_set_test = format_data_train_count(test_kmer_dir, cl_train_opt, kmer_size, train_kmers = data_set_train['features'])
        logger.debug("Opt: Performing cross validation.")

        LOO = cross_validation.LeaveOneOut(array(data_set_test['target']).shape[0])

        results_kmer_optimum[kmer]={}

        partial_kmer_loocv_mp = partial(kmer_loocv_mp, data_set_train = data_set_train, data_set_test = data_set_test)
        results_kmer_optimum[kmer] = pool.map(partial_kmer_loocv_mp, LOO)

    pool.close()
    pool.join()

    for kmer_size in kmer_range:
        kmer = 'kmer' + str(kmer_size)
        logger.debug("Opt: outputting results kmer: {0}".format(str(kmer)))
        # overall accuracy for each k-mer
        measures = results_kmer_optimum[kmer][0][0].keys()
        accuracy = []
        for res in results_kmer_optimum[kmer]:
            try:
                accuracy.append(mean(res[0].values()))
            except TypeError:
                continue
        accuracy_kmer.append(mean(accuracy))     
        accuracy = []
        for dist in measures:
            for res in results_kmer_optimum[kmer]:
                try:
                    accuracy.append(res[0][dist])
                except TypeError:
                    continue
            individual_accuracy[dist].append(mean(accuracy))
    optimum_kmer = kmer_range[accuracy_kmer.index(max(accuracy_kmer))]
    for d in measures:
        individual_accuracy_opt[d] = individual_accuracy[d][accuracy_kmer.index(max(accuracy_kmer))]
    optimized_kmer_results[taxon_rank] = [optimum_kmer, individual_accuracy_opt]    
    pickle.dump(optimized_kmer_results, open(test_kmer_dir + 'optimized_kmer_results_' + TAXON_RANK + '.p', "wb" ))    
    return optimized_kmer_results
 




def build_optimization_set(fasta_file, dbdirectory, fasta_type, kmin, find_kmer_optimum, kmer_specified):
            '''
            Build test count data using jellyfish
            '''
            proc = glob.glob(os.path.abspath(dbdirectory) + '/countsdb/' + '*kmer*')

            nmb_kmer_files = len(proc)

            if nmb_kmer_files < 1:
                logger.error("Failed to find counts for the training set, need to run buildDB first. Exiting") 
                sys.exit(1)

            optimized_kmer_dir = os.path.abspath(dbdirectory) + '/optimum_kmer/'

            subprocess.call(['rm','-rf', optimized_kmer_dir])

            if subprocess.call(['mkdir', '-p', optimized_kmer_dir]) != 0:
                logger.error("Failed to create the directory {0}. Exiting.".format(str(optimized_kmer_dir)))  
                sys.exit(1)

            optim_train_fasta = optimized_kmer_dir + 'optimized_train.fasta'

            # get sequences from the training set
            train_fasta = os.path.abspath(dbdirectory) + '/train.fasta'

            if subprocess.call(['ls', train_fasta], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
               logger.error("Failed to find the training fasta file in {0}. Exiting.".format(os.path.abspath(dbdirectory)))
               sys.exit(1)

            call_type='test'               

            mu, std, max_seq_len, min_seq_len = get_stats_on_sequences_test_set(fasta_file, fasta_type, call_type)

            if min_seq_len < 200:
               min_seq_len = 200 

        
            call_type='train'         
            
            gis_in_train, max_kmer = get_stats_on_sequences_test_set(train_fasta, fasta_type, call_type)

            if find_kmer_optimum:
                kmer_range=range(int(kmin),int(max_kmer)+1)
            elif kmer_specified:
                kmer_range=[int(kmer_specified)]

            if std > 0.0:
                seq_len_test_training = abs(array(map(round,np.random.normal(mu, std, len(gis_in_train)))))
                seq_len_test_training[seq_len_test_training < min_seq_len] = min_seq_len
                if max(seq_len_test_training) < max_seq_len:
                   seq_len_test_training[seq_len_test_training == max(seq_len_test_training)] = max_seq_len
            elif std == 0.0:
                seq_len_test_training = array(mu)

            #print seq_len_test_training
            logger.info("Found mean sequence length mu: {0} with std: {1} from sequences in the test set".format(str(mu),str(std)))

            if find_kmer_optimum:
               logger.info("Generate count data for the pre-selection of sequence similarity measures and k-mer optimization, counts outputted to {0}".format(str(optimized_kmer_dir)))  
            elif kmer_specified:
               logger.info("Generate count data for the pre-selection of sequence similarity measures for the k-mer value of: {0} , counts outputted to {1}".format(str(kmer_specified),str(optimized_kmer_dir)))
 
            seq_optimization_count = 0 
            seq_in_train = 0

            handle_optim_fasta_out = open(optim_train_fasta, "w")
            # output the training/test fasta file to be parallelized  
            for record in SeqIO.parse(open(train_fasta, "rU"), "fasta"):
                seq_in_train += 1
                #print seq_in_train
                if 'N' in record.seq.upper():
                    record.seq = Seq(str(record.seq).upper().replace('N',''))
                diff_seq_len = len(record.seq) - seq_len_test_training
                try:
                  if std > 0.0:  
                      diff_seq_value = int(seq_len_test_training[list(diff_seq_len).index(min(diff_seq_len[diff_seq_len > 0]))])
                      seq_len_test_training[list(diff_seq_len).index(min(diff_seq_len[diff_seq_len > 0]))] = 0.0
                  elif std == 0.0:
                      if diff_seq_len > 0:
                          diff_seq_value =  mu
                      else:
                          diff_seq_value =  ValueError
                except ValueError:
                    diff_seq_value = 0
                if diff_seq_value != 0:  
                    start = random.uniform(0, len(record.seq) - int(diff_seq_value))
                    end = start + diff_seq_value                   
                else:
                    start = 0.0
                    end = len(record.seq)
                if len(record.seq) >= min_seq_len and len(set(record.seq)) > 1:
                    record_count = SeqRecord(record.seq[int(start):int(end)],id=record.id,name=record.name,description=record.description) 
                    SeqIO.write(record_count,handle_optim_fasta_out,"fasta")
                    seq_optimization_count += 1
                elif len(set(record.seq)) > 1:
                    record_count = SeqRecord(record.seq,id=record.id,name=record.name,description=record.description) 
                    SeqIO.write(record_count,handle_optim_fasta_out,"fasta") 
                    seq_optimization_count += 1

            handle_optim_fasta_out.close()
            logger.debug("Created the FASTA file for the selection of similarity measures with mu: {0}, std: {1}, files outputted to {2}".format(str(mu),str(std),str(optimized_kmer_dir)))
 
            percent_simulated = seq_optimization_count/seq_in_train * 100
            logger.debug("The simulated fasta file used in the optimization contains {0} % of sequences compared to the training set".format(str(percent_simulated)))

            logger.info("Get taxon classes for the optimization training set.")
            get_train_classes_opt(optimized_kmer_dir)


            # output k-mer counts
            logger.debug("Build kmer counts for the the selection of similarity measures in {0}".format(str(optimized_kmer_dir)))
            
            pool=multiprocessing.Pool(num_threads)

            dict_kmer_counts = {}

            # get the size of the optim_train_fasta 
            file_size_bytes = os.path.getsize(optim_train_fasta)   

            if file_size_bytes/1000000000 <= 5.0:
                for kmer_size in kmer_range:
                    dict_kmer_counts[kmer_size] = {}
                    handle = open(optim_train_fasta, "rU")
                    records = SeqIO.parse(handle, "fasta")     
                    partial_get_kmer_counts_jellyfish_mp_2 = partial(get_kmer_counts_jellyfish_mp_2, kmer = kmer_size, counts_dir = optimized_kmer_dir)
                    dict_kmer_counts[kmer_size] = dict(pool.map(partial_get_kmer_counts_jellyfish_mp_2, records))
                    handle.close()
                pool.close()
                handle.close()
            else:
                logger.info("The size of the training set for optimization exceeds 5GB. Will need to build a BLAST DB for this datsets for indexing the genomes.")
                build_blast_db(optim_train_fasta, fasta_type='nucl')
                gi_ids=[]
                seq_len_hash={}
                handle = open(optim_train_fasta, "rU")
                for record in SeqIO.parse(handle, "fasta"):
                   gi_ids.append(int(record.id.split("|")[1]))
                   seq_len_hash[int(record.id.split("|")[1])] = int(len(record.seq))
                handle.close()
                records = gi_ids
                for kmer_size in kmer_range: 
                   dict_kmer_counts[kmer_size] = {}
                   partial_get_kmer_counts_jellyfish_mp_2_blast = partial(get_kmer_counts_jellyfish_mp_2_blast, kmer = kmer_size, counts_dir = optimized_kmer_dir, seq_len_hash = seq_len_hash)
                   dict_kmer_counts[kmer_size] = dict(pool.map(partial_get_kmer_counts_jellyfish_mp_2_blast, records))
                pool.close()    


            for kmer_size in kmer_range:
                pickle.dump(dict_kmer_counts[kmer_size], open( optimized_kmer_dir + "/dict_kmer_" + str(kmer_size) + "_counts.p", "wb" ))

            if find_kmer_optimum:    
                logger.info("Performing kmer optimization and measure selection")
                results_opt = perform_kmer_optimization(dbdirectory, gis_in_train, kmer_range)
                logger.info("Found optimum k-mer size and completed measure selection")
                return results_opt
            elif kmer_specified:
                logger.info("Performing measure selection for k-mer: {0}".format(str(kmer_specified)))
                results_opt = perform_kmer_optimization(dbdirectory, gis_in_train, kmer_range)
                logger.info("Finished measure selection for k-mer: {0}".format(str(kmer_specified)))
                return results_opt



def build_optimization_set_blast(fasta_file, dbdirectory, fasta_type):
            '''
            Build the simulated test set for estimating the BLAST accuracy 
            '''
            
            optimized_kmer_dir = os.path.abspath(dbdirectory) + '/optimum_kmer/'

            subprocess.call(['rm','-rf', optimized_kmer_dir])

            if subprocess.call(['mkdir', '-p', optimized_kmer_dir]) != 0:
                logger.error("Failed to create the directory {0}. Exiting.".format(str(optimized_kmer_dir)))  
                sys.exit(1)

            optim_train_fasta = optimized_kmer_dir + 'optimized_train.fasta'

            # get sequences from the training set
            train_fasta = os.path.abspath(dbdirectory) + '/train.fasta'

            if subprocess.call(['ls', train_fasta], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
               logger.error("Failed to find the training fasta file in {0}. Exiting.".format(os.path.abspath(dbdirectory)))
               sys.exit(1)

            call_type='test'               

            mu, std, max_seq_len, min_seq_len = get_stats_on_sequences_test_set(fasta_file, fasta_type, call_type)
            
            if min_seq_len < 200:
               min_seq_len = 200 

            call_type='train'         
            
            gis_in_train, max_kmer = get_stats_on_sequences_test_set(train_fasta, fasta_type, call_type)


            if std > 0.0:    
                seq_len_test_training = abs(array(map(round,np.random.normal(mu, std, len(gis_in_train)))))
                seq_len_test_training[seq_len_test_training < min_seq_len] = min_seq_len
                if max(seq_len_test_training) < max_seq_len:
                   seq_len_test_training[seq_len_test_training == max(seq_len_test_training)] = max_seq_len
            elif std == 0.0:
                seq_len_test_training = array(mu)



            #print seq_len_test_training
            logger.info("Found mean sequence length mu: {0} with std: {1} from sequences in the test set".format(str(mu),str(std)))

            logger.info("Generate sequences for estimating the BLAST accuracy, sequences outputted to {0}".format(str(optimized_kmer_dir)))

            
            seq_optimization_count = 0 
            seq_in_train = 0

            handle_optim_fasta_out = open(optim_train_fasta, "w")
            # output the training/test fasta file to be parallelized  
            for record in SeqIO.parse(open(train_fasta, "rU"), "fasta"):
                seq_in_train += 1
                if 'N' in record.seq.upper():
                    record.seq = Seq(str(record.seq).upper().replace('N',''))
                diff_seq_len = len(record.seq) - seq_len_test_training
                try:
                  if std > 0.0:  
                      diff_seq_value = int(seq_len_test_training[list(diff_seq_len).index(min(diff_seq_len[diff_seq_len > 0]))])
                      seq_len_test_training[list(diff_seq_len).index(min(diff_seq_len[diff_seq_len > 0]))] = 0.0
                  elif std == 0.0:
                      if diff_seq_len > 0:
                          diff_seq_value =  mu
                      else:
                          diff_seq_value =  ValueError
                except ValueError:
                    diff_seq_value = 0
                if diff_seq_value != 0:  
                    start = random.uniform(0, len(record.seq) - int(diff_seq_value))
                    end = start + diff_seq_value                   
                else:
                    start = 0.0
                    end = len(record.seq)
                if len(record.seq) >= min_seq_len and len(set(record.seq)) > 1:
                    record_count = SeqRecord(record.seq[int(start):int(end)],id=record.id,name=record.name,description=record.description) 
                    SeqIO.write(record_count,handle_optim_fasta_out,"fasta")
                    seq_optimization_count += 1
                elif len(set(record.seq)) > 1:
                    record_count = SeqRecord(record.seq,id=record.id,name=record.name,description=record.description) 
                    SeqIO.write(record_count,handle_optim_fasta_out,"fasta")    
                    seq_optimization_count += 1

            handle_optim_fasta_out.close()
            logger.debug("Created the FASTA file for the selection of similarity measures with mu: {0}, std: {1}, files outputted to {2}".format(str(mu),str(std),str(optimized_kmer_dir)))

            percent_simulated = seq_optimization_count/seq_in_train * 100

            logger.debug("The simulated fasta file used in the optimization contains {0} % of sequences compared to the training set ".format(str(percent_simulated)))
            logger.info("Get taxon classes for the optimization training set.")
            get_train_classes_opt(optimized_kmer_dir)





def blast_loocv_mp(loo_iter, dist_blast, classes_array):
    train_index, test_index = list(loo_iter)
    dist_blast_sub = dist_blast[test_index,:][:,train_index]
    accuracy, ftest_score = get_accuracy_scores_tune_knn_and_predict_loo(dist_blast_sub,classes_array[train_index],classes_array[test_index],K=1)
    return accuracy, ftest_score



def evaluate_blast_train_data(dbdirectory, fasta_type):
    '''
    Format blast data for classification according to classes in cl_train 
    '''
    pool=multiprocessing.Pool(num_threads)
    #train_kmer_dir = os.path.abspath(dbdirectory) + '/countsdb/' 
    optimized_kmer_dir = os.path.abspath(dbdirectory) + '/optimum_kmer/'
    optim_train_fasta = optimized_kmer_dir + 'optimized_train.fasta'
    call_type='train'         
    #cl_train = pickle.load(open(os.path.abspath(dbdirectory) + '/cl_train.p', 'rb'))
    cl_train = pickle.load(open(optimized_kmer_dir + '/cl_train_opt_' + TAXON_RANK + '.p', 'rb'))
    classes_array = []
    for taxon_name in cl_train.keys():
        for gi in cl_train[taxon_name]:
            classes_array.append(taxon_name)
    classes_array = array(classes_array)
    dist_accuracies_train = {}
    gis_blast = [item for sublist in cl_train.values() for item in sublist]
    blast_output_train_file_fasta =os.path.abspath(dbdirectory) + '/blast_train/blast_train.txt'
    n = len(gis_blast)
    blast_scores = np.zeros((n,n))
    with open(blast_output_train_file_fasta) as f:
        for line in f:
            if (int(line.split('\t')[1].split('|')[1]) in gis_blast and int(line.split('\t')[0].split('|')[1]) in gis_blast) and (blast_scores[gis_blast.index(int(line.split('\t')[0].split('|')[1])),gis_blast.index(int(line.split('\t')[1].split('|')[1]))] <= float(line.split('\t')[-1].split('\n')[0])):
                blast_scores[gis_blast.index(int(line.split('\t')[0].split('|')[1])),gis_blast.index(int(line.split('\t')[1].split('|')[1]))] = float(line.split('\t')[-1].split('\n')[0])
    max_val = np.apply_along_axis(return_max_in_array, axis=1, arr=blast_scores)
    # for those entries that have 0.0 max_val
    max_val[max_val == 0.0] = 1.0
    dist_blast = 1 - blast_scores/max_val[:,None] # this is like a transpose of max_val
    #replicate_runs = 0
    # LOO
    LOO = cross_validation.LeaveOneOut(len(classes_array)) 
    results_blast_train = {}
    results_blast_train_accu_loo = [] 

    partial_blast_loocv_mp = partial(blast_loocv_mp, dist_blast = dist_blast, classes_array = classes_array)
    results_blast = pool.map(partial_blast_loocv_mp, LOO)
    pool.close()
    pool.join()
       
    for res in results_blast:
        results_blast_train_accu_loo.append(res[0])

    dist_accuracies_train['d_blast'] = [array(results_blast_train_accu_loo).mean(), TAXON_RANK]
    pickle.dump(dist_accuracies_train, open(optimized_kmer_dir + 'BLAST_train_results_' + TAXON_RANK + '.p', "wb" ))
    return dist_accuracies_train




def format_jellyfish_counts(res_jellyfish, seqid):
    '''
    '''
    dict_kmer_counts = {}
    if res_jellyfish is not None:
        dummy_dict_counts={}
        counts_kmers = res_jellyfish
        counts = [int(counts_kmers[i].split(">")[1]) for i in range(0,(len(counts_kmers) - 1),2)]
        kmers = [str(counts_kmers[i]) for i in range(1,(len(counts_kmers) - 1),2)]
        for i,kmer in enumerate(kmers):
            dummy_dict_counts[kmer]=counts[i]
        dict_kmer_counts[seqid] = dummy_dict_counts
    else:
        dict_kmer_counts = None
    return dict_kmer_counts




def format_data_test_count_sequential(dict_kmer_counts, kmers_train_features):
    '''
    Format test data for classification
    '''
    data_set= {}
    counts_array = []
    freq_array = []
    kmer_keys = None
    counts = dict_kmer_counts
    # get gis in the same order as in the fasta test file
    for gi in dict_kmer_counts.keys():
             if kmer_keys is None:
                kmer_keys = counts[gi].keys()
                dummy_list = []
                for val in kmer_keys:
                    dummy_list.append(counts[gi][val])
                counts_array.append(dummy_list)
             else:
                dummy_list=[] 
                if kmer_keys == counts[gi].keys():
                    for val in kmer_keys:
                        dummy_list.append(counts[gi][val])
                else:
                    kmer_add = set(counts[gi].keys()) - set(kmer_keys)
                    if len(kmer_add) > 0:
                        kmer_keys = kmer_keys + list(kmer_add)
                    for val in kmer_keys:
                            if val in counts[gi]:
                                dummy_list.append(counts[gi][val])
                            else:
                                dummy_list.append(0)
                counts_array.append(dummy_list)

    for i in range(0,len(counts_array)):
        difference = len(kmers_train_features) - len(counts_array[i])
        if difference > 0:
            counts_array[i] = counts_array[i] + difference*[0]

    counts_array_ord=[]
    for r in counts_array:
        dummy=[0]*(len(r))
        for i,k in enumerate(kmer_keys):
            indx = kmers_train_features.index(k)
            dummy[indx] = r[i]
    counts_array_ord.append(dummy)

    data_set['counts'] = counts_array_ord
    data_set['features'] = kmers_train_features
    data_set['names'] = gi
    return data_set



def estimate_accuracy_loocv_mp(loo_iter, data_set_train, data_set_test, dist_blast, dist_cb, classes_array, dist_accuracies_train):
    train_index, test_index = list(loo_iter)
    if data_set_train is not None:
        if sum(data_set_test['counts'][test_index]) != 0:
                        y_train = classes_array[train_index]
                        distances_test = {}
                        true_classe = classes_array[test_index]
                        if 'd_js' in dist_accuracies_train.keys():                  
                            feature = 'counts'
                            distances_test['d_js'] = pairwise_distances(array(data_set_test[feature])[test_index], array(data_set_train[feature])[train_index], metric=js_pairwise)
                        else:
                            distances_test['d_js'] = None
                        if 'd_eucl' in dist_accuracies_train.keys():
                            feature = 'counts'
                            distances_test['d_eucl'] = pairwise_distances(array(data_set_test[feature])[test_index], array(data_set_train[feature])[train_index], metric=norm_euclidean)
                        else:
                            distances_test['d_eucl'] = None
                        if 'd_blast' in dist_accuracies_train.keys():
                            distances_test['d_blast'] = dist_blast[test_index,:][:,train_index]
                        else:
                            distances_test['d_blast'] = None
                        if 'd_cb' in dist_accuracies_train.keys():
                            dist_cb = dist_cb
                            distances_test['d_cb'] = dist_cb[test_index,:][:,train_index]                
                        else:
                            distances_test['d_cb'] = None                            
                        kn = 1
                        f_stat_p = partial(f_stat,Y = y_train)
                        f_sum = 0
                        d_combined = 0
                        d_count = 0
                        for d in distances_test.keys():
                            if distances_test[d] is not None:
                                d_count += 1
                                f = np.apply_along_axis(f_stat_p,axis=1,arr=distances_test[d]) 
                                f_sum = f_sum + f
                                d_combined = d_combined + distances_test[d]*f[:,None]
                        if d_count > 0: 
                            f_sum[f_sum == 0.0] = 1.0
                            d_combined = d_combined/f_sum[:,None]
                            pred = run_knn_and_predict(d_combined, y_train, K=kn)
                        elif d_count == 0:
                            pred = 'NA'
                        if true_classe == pred[0]:
                            accu = 1
                        else:
                            accu = 0
    else:
                        y_train = classes_array[train_index]
                        distances_test = {}
                        true_classe = classes_array[test_index]
                        if 'd_js' in dist_accuracies_train.keys():                  
                            feature = 'counts'
                            distances_test['d_js'] = pairwise_distances(array(data_set_test[feature])[test_index], array(data_set_train[feature])[train_index], metric=js_pairwise)
                        else:
                            distances_test['d_js'] = None
                        if 'd_eucl' in dist_accuracies_train.keys():
                            feature = 'counts'
                            distances_test['d_eucl'] = pairwise_distances(array(data_set_test[feature])[test_index], array(data_set_train[feature])[train_index], metric=norm_euclidean)
                        else:
                            distances_test['d_eucl'] = None
                        if 'd_blast' in dist_accuracies_train.keys():
                            distances_test['d_blast'] = dist_blast[test_index,:][:,train_index]
                        else:
                            distances_test['d_blast'] = None
                        if 'd_cb' in dist_accuracies_train.keys():
                            dist_cb = dist_cb
                            distances_test['d_cb'] = dist_cb[test_index,:][:,train_index]                
                        else:
                            distances_test['d_cb'] = None                            
                        kn = 1
                        f_stat_p = partial(f_stat,Y = y_train)
                        f_sum = 0
                        d_combined = 0
                        d_count = 0
                        for d in distances_test.keys():
                            if distances_test[d] is not None:
                                d_count += 1
                                f = np.apply_along_axis(f_stat_p,axis=1,arr=distances_test[d]) 
                                f_sum = f_sum + f
                                d_combined = d_combined + distances_test[d]*f[:,None]
                        if d_count > 0: 
                            f_sum[f_sum == 0.0] = 1.0
                            d_combined = d_combined/f_sum[:,None]
                            pred = run_knn_and_predict(d_combined, y_train, K=kn)
                        elif d_count == 0:
                            pred = 'NA'
                        if true_classe == pred[0]:
                            accu = 1
                        else:
                            accu = 0
    return (accu, float(f_sum))
 





def estimate_accuracy(dbdirectory):
    '''
    '''

    pool=multiprocessing.Pool(num_threads)

    test_kmer_dir = os.path.abspath(dbdirectory) + '/optimum_kmer/'
 
    train_kmer_dir = os.path.abspath(dbdirectory) + '/countsdb/' 

    train_compress_dir = os.path.abspath(dbdirectory) + '/compressdb/'  

    # BLAST                 
    if use_blast:
        try:
            results_blast_train = pickle.load(open(os.path.abspath(dbdirectory) + '/optimum_kmer/BLAST_train_results_' + TAXON_RANK + '.p', 'rb'))
        except IOError:
            logger.error("Optimized BLAST data not available in {0}".format(open(os.path.abspath(dbdirectory))))
            sys.exit(1)
    else:
        results_blast_train = None

    # Compression                
    if use_compress:
        try:
            results_cb_train = pickle.load(open(os.path.abspath(dbdirectory) + '/optimum_kmer/cb_train_results_'  + TAXON_RANK + '.p', 'rb'))
        except IOError:
            logger.error("Compression data not available in {0}".format(open(os.path.abspath(dbdirectory))))
            sys.exit(1)
    else:
        results_cb_train = None

    # Counts
    if use_counts:
        try:
            results_opt_counts = pickle.load(open(os.path.abspath(dbdirectory) + '/optimum_kmer/optimized_kmer_results_' + TAXON_RANK + '.p', 'rb'))
        except IOError:
            logger.error("Optimized counts data not available in {0}".format(open(os.path.abspath(dbdirectory))))
            sys.exit(1)
    else:
       results_opt_counts = None 
 
    # Get classes in the training set 
    try:
         cl_train = pickle.load(open(os.path.abspath(dbdirectory) + '/optimum_kmer/cl_train_opt_' + TAXON_RANK + '.p', 'rb'))
    except IOError:
        logger.error("Train data classes not available in {0}".format(open(os.path.abspath(dbdirectory))))
        sys.exit(1)

        
    for taxon in [TAXON_RANK]:
            dist_accuracies_train = {}
            if use_counts:
                opt_kmer = results_opt_counts[taxon][0] 
                kmer_size = opt_kmer
            else:
                opt_kmer = None
                kmer_size = opt_kmer                
            logger.info("Estimate the model classification accuracy with kmer={0}, cl acc threshold = {1} at taxon: {2}".format(opt_kmer, min_accu_threshold, taxon))
            if use_counts:
                for d in results_opt_counts[taxon][1]:
                    dist_accuracies_train[d] = results_opt_counts[taxon][1][d]
            if use_blast:
                dist_accuracies_train['d_blast'] =  results_blast_train['d_blast'][0]
            else:
                dist_accuracies_train['d_blast'] = None 
            if use_compress:
                dist_accuracies_train['d_cb'] =  results_cb_train['d_cb'][0]
            else:
                dist_accuracies_train['d_cb'] = None 
            for smeasure in dist_accuracies_train.keys():
                if dist_accuracies_train[smeasure] < min_accu_threshold: 
                    dist_accuracies_train.pop(smeasure)
                    logger.info("Eliminated measure from the model : {0}".format(str(smeasure)))
                if len(dist_accuracies_train.keys()) == 0:
                    logger.info("No similarity measures selected for this dataset will not attempt the classification. Exiting")
                    sys.exit(1)
            if len(dist_accuracies_train.keys()) == 1:
                parameters = {}
                parameters['measures'] = [] 
                if 'd_blast' in dist_accuracies_train.keys():
                    parameters['accu'] = results_blast_train['d_blast'][0]
                elif 'd_cb' in  dist_accuracies_train.keys():
                    parameters['accu'] = results_cb_train['d_cb'][0]
                parameters['thr'] =  min_accu_threshold
                parameters['measures'] = dist_accuracies_train.keys()
                parameters['taxon'] = taxon
                parameters['kmer'] = None
            elif len(dist_accuracies_train.keys()) == 2 and ('d_js' in dist_accuracies_train.keys() and 'd_eucl' in dist_accuracies_train.keys()):
                parameters = {}
                parameters['measures'] = [] 
                parameters['accu'] = mean([results_opt_counts[taxon][1]['d_eucl'],results_opt_counts[taxon][1]['d_js']])  
                parameters['thr'] =  min_accu_threshold
                parameters['measures'] = [m for m in dist_accuracies_train.keys()]
                parameters['taxon'] = taxon
                parameters['kmer'] = opt_kmer                 
            else: 
                if use_counts:
                    data_set_train = format_data_train_count(train_kmer_dir, cl_train, kmer_size)
                    data_set_test = format_data_train_count(test_kmer_dir, cl_train, kmer_size, train_kmers = data_set_train['features'])
                else:
                    data_set_train = None
                    data_set_test = None                    
                classes_array = []
                for taxon_name in cl_train.keys():
                    for gi in cl_train[taxon_name]:
                        classes_array.append(taxon_name)
                classes_array = array(classes_array)
                # get blast data from the training set
                if use_blast:    
                    gis_blast = [item for sublist in cl_train.values() for item in sublist]
                    blast_output_train_file_fasta = os.path.abspath(dbdirectory) + '/blast_train/blast_train.txt'
                    n = len(gis_blast)
                    blast_scores = np.zeros((n,n))
                    with open(blast_output_train_file_fasta) as f:
                        for line in f:
                            if (int(line.split('\t')[1].split('|')[1]) in gis_blast and int(line.split('\t')[0].split('|')[1]) in gis_blast) and (blast_scores[gis_blast.index(int(line.split('\t')[0].split('|')[1])),gis_blast.index(int(line.split('\t')[1].split('|')[1]))] <= float(line.split('\t')[-1].split('\n')[0])):
                                blast_scores[gis_blast.index(int(line.split('\t')[0].split('|')[1])),gis_blast.index(int(line.split('\t')[1].split('|')[1]))] = float(line.split('\t')[-1].split('\n')[0])
                    max_val = np.apply_along_axis(return_max_in_array, axis=1, arr=blast_scores)
                    # for those entries that have 0.0 max_val
                    max_val[max_val == 0.0] = 1.0
                    dist_blast = 1 - blast_scores/max_val[:,None] # this is like a transpose of max_val
                    # get the cross validation
                else:
                    dist_blast = None
                if use_compress:
                   dist_cb = pickle.load(open(train_compress_dir + 'dist_cb_train_opt_' + TAXON_RANK + '.p', 'rb'))
                else:
                   dist_cb = None 

                LOO = cross_validation.LeaveOneOut(len(classes_array))

                partial_est_acc_loocv_mp = partial(estimate_accuracy_loocv_mp, data_set_train = data_set_train, data_set_test = data_set_test, dist_blast = dist_blast, dist_cb = dist_cb, classes_array = classes_array, dist_accuracies_train = dist_accuracies_train)
                results_acc = pool.map(partial_est_acc_loocv_mp, LOO)
                pool.close()
                pool.join()

                acc = []

                for res in results_acc:
                    acc.append(res[0])

                est_accu = sum(acc)/len(acc)

                parameters = {}
                parameters['measures'] = [] 
                parameters['accu'] = est_accu
                parameters['thr'] =  min_accu_threshold
                parameters['measures'] = [m for m in dist_accuracies_train.keys()]
                parameters['taxon'] = taxon
                if 'd_eucl' in dist_accuracies_train.keys() or 'd_js' in dist_accuracies_train.keys():
                    parameters['kmer'] = opt_kmer
                else:
                    parameters['kmer'] = None
            pickle.dump(parameters, open( test_kmer_dir + 'estimated_model_accuracy_' + TAXON_RANK + '.p', 'wb'))
            return parameters




def perform_classification_mp(record, opt_kmer, data_set_train, data_set_test_counts_dir, dbdirectory, train_taxon_cl, gis_blast, out_file, dist_accuracies_train, num_procs_c):
    '''
    '''
    global count_seq_processed

    test_data_dir = os.path.dirname(os.path.abspath(out_file))
    
    if data_set_train != None:
        res_jellyfish = get_kmer_counts_jellyfish_mp(record, opt_kmer, data_set_test_counts_dir)
        #print record.id
        dict_kmer_counts = format_jellyfish_counts(res_jellyfish, record.id)
        #print dict_kmer_counts
        data_count_test_set = format_data_test_count_sequential(dict_kmer_counts, data_set_train['features'])
    else:
        data_count_test_set = None
    if 'd_blast' in dist_accuracies_train.keys():
        dist_blast_test_set = get_blast_hits_sequential(record.id, test_data_dir, gis_blast)
    else:
        dist_blast_test_set = None
    if 'd_cb' in dist_accuracies_train.keys():
        dist_cb_test_set = get_compression_sequences_sequential_mp(record, dbdirectory, num_procs_c)
    else:
        dist_cb_test_set = None
    res_knn = knn_predict(data_set_train, array(train_taxon_cl), data_count_test_set, dist_blast_test_set, dist_cb_test_set, dist_accuracies_train)
    with lock:
        count_seq_processed.value += 1
    if count_seq_processed.value % 1000000 == 0:
        logger.info("Processed {0} sequences".format(count_seq_processed.value))
    FILE_OUT = open(out_file, "a")
    FILE_OUT.write("%s\t%s\n" % (record.id, res_knn[0][0]))
    FILE_OUT.close()
    




def perform_classification(file_fasta_test, dbdirectory, fasta_type, results_opt_counts, results_blast_train, results_cb_train, est_accu):
    '''
    Perform classsification of sequences in the test set using the training set with the optimized kmer value and selected similarity measures 
    '''
    global count_seq_processed


    dir_basename_test = os.path.dirname(file_fasta_test)

    optimized_train_fasta = os.path.abspath(dbdirectory) + '/optimum_kmer/' + 'optimized_train.fasta'
 
    optimized_kmer_dir = os.path.abspath(dbdirectory) + '/optimum_kmer/'
   
    data_set_train_counts_dir = os.path.abspath(dbdirectory) + '/countsdb/'

    data_set_test_counts_dir = os.path.abspath(dbdirectory) + '/counts_test/'

    for taxon in [TAXON_RANK]:
        dist_accuracies_train = {}
        if use_counts:
            opt_kmer = results_opt_counts[taxon][0] 
        cl_train = pickle.load(open(os.path.abspath(dbdirectory) + '/cl_train_' + TAXON_RANK + '.p', 'rb'))
        gis_blast = [item for sublist in cl_train.values() for item in sublist]
        if use_counts:
            for d in results_opt_counts[taxon][1]:
                dist_accuracies_train[d] = results_opt_counts[taxon][1][d]
        if use_blast:
            dist_accuracies_train['d_blast'] =  results_blast_train['d_blast'][0]
        if use_compress:
            dist_accuracies_train['d_cb'] =  results_cb_train['d_cb'][0]
        train_taxon_cl = []
        for taxon_name in cl_train.keys():
            for gi in cl_train[taxon_name]:
                train_taxon_cl.append(taxon_name)
        for smeasure in dist_accuracies_train.keys():
            if dist_accuracies_train[smeasure] < min_accu_threshold:
               dist_accuracies_train.pop(smeasure)
        if min_accu_threshold != est_accu['thr'] or taxon != est_accu['taxon'] or set(dist_accuracies_train.keys()) != set(est_accu['measures']):
               logger.info("Warning! The parameters used in this run are different from those used for estimating the model accuracy. The estimated accuracy of the model might not be correct.") 
        if len(dist_accuracies_train.keys()) == 0:
             logger.info("No similarity measures selected for this dataset will not attempt the classification. Exiting")
             sys.exit(1)
        else:
            if ('d_eucl' in dist_accuracies_train.keys()) or ('d_js' in dist_accuracies_train.keys()):     
                # get training set
                logger.info("Loading counts for the training set from {0}".format(data_set_train_counts_dir))
                data_set_train = format_data_train_count(data_set_train_counts_dir, cl_train, opt_kmer)
            else:
                data_set_train = None
                opt_kmer = None
            if output_file_name:
                out_file = os.path.abspath(dir_basename_test) + '/' + output_file_name
            else:
                out_file = os.path.abspath(dir_basename_test) + '/' + 'cssscl_results_' + taxon + '.txt'
            if os.path.isfile(out_file):
                subprocess.call(['rm', out_file]) 
            handle = open(file_fasta_test, "rU")
            logger.info("Start the classification of sequences with the following set of parameters: \nkmer opt = {0} \nmin cl acc threshold = {1} \ntaxon: {2} \nselected similarity measures: {3}".format(opt_kmer, min_accu_threshold, taxon, ",".join(dist_accuracies_train.keys())))
            records = SeqIO.parse(handle, "fasta")
            if use_compress and 'd_cb' in dist_accuracies_train.keys():
                num_procs_c, num_procs_o = get_procs()
            else:
                num_procs_c, num_procs_o = 0, num_threads
            partial_perform_classification_mp = partial(perform_classification_mp, opt_kmer = opt_kmer, data_set_train = data_set_train, data_set_test_counts_dir = data_set_test_counts_dir, dbdirectory = dbdirectory, train_taxon_cl = train_taxon_cl, gis_blast = gis_blast, out_file = out_file, dist_accuracies_train = dist_accuracies_train, num_procs_c = num_procs_c)
            pool = MyPool(num_procs_o)
            pool.map(partial_perform_classification_mp, records)
            pool.close()
            pool.join()
            FILE_OUT = open(out_file, "a")
            FILE_OUT.write("%s %s\n" % ('Overall estimated accuracy :', str(est_accu['accu'])))
            FILE_OUT.close()
            handle.close()
            logger.info("Finished classification")

       

def get_train_classes(dbdirectory):
    train_fasta = os.path.abspath(dbdirectory) + '/train.fasta'  
    gi_ids=[]
    handle = open(train_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        if record.id.split("|")[0] == 'gi':
            gi_ids.append(int(record.id.split("|")[1]))
    handle.close()  
    cl_train = get_classes_in_train(gi_ids, TAXON_RANK)   
    pickle.dump(cl_train, open(os.path.abspath(dbdirectory) + '/cl_train_' + TAXON_RANK + '.p', 'wb' ))



def get_train_classes_opt(dbdirectory):
    train_fasta = os.path.abspath(dbdirectory) + '/optimized_train.fasta'  
    gi_ids=[]
    handle = open(train_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        if record.id.split("|")[0] == 'gi':
            gi_ids.append(int(record.id.split("|")[1]))
    handle.close()  
    cl_train_opt = get_classes_in_train(gi_ids, TAXON_RANK)   
    pickle.dump(cl_train_opt, open(os.path.abspath(dbdirectory) + '/cl_train_opt_' + TAXON_RANK + '.p', 'wb' ))



def check_program_in_path(program):
    #if subprocess.call(['which', program], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
    executable = os.popen("type %s" % str(program)).read().strip('\n')
    if executable == program + ': not found':
           logger.error("{0} not in the path. Exiting.".format(program))
           sys.exit(1)



def main(args):

    '''Estimate parameters and perform classification'''

    global db, logger, num_threads, TAXON_RANK, use_blast, use_compress, use_counts, count_seq_processed, min_accu_threshold, output_file_name, lock, blast_type, num_alignments 

    random.seed(12345)

    logger = args.logging.getLogger(__name__)
    
    db = connect(args)

    cnt_gi = db.gitaxid.count()
    cnt_taxa = db.taxa.count()

    if cnt_gi == 0 and cnt_taxa == 0:
        db_present = False
    else:
        db_present = True

    if not db_present:   
        logger.error("No entries detected in the db. You first need to build the databases for the training set using build_dbs with the option -btax. Exiting.")
        sys.exit(1)

    # shared object count
    count_seq_processed = Value('i', 0)
    lock = Lock()

    num_threads = int(args.number_threads)

    max_cpus = multiprocessing.cpu_count()

    output_file_name = args.output_file
   
    min_accu_threshold = float(args.min_accu_threshold)

    if args.use_blast:
        use_blast = True
    else:
        use_blast = False

    blast_type = args.use_blast

    num_alignments = args.num_alignments

    if args.use_compression:
        use_compress = True
    else:
        use_compress = False


    if args.disable_kmer_count:
        use_counts = False
    else:
        use_counts = True


    if not use_counts and not use_blast:
       logger.error("kmer count based measures can only be disabled if -blast is specified.")
       sys.exit(1) 

    if not use_counts and args.specify_kmer:
       logger.error("-kmeroff and -kmer can not be specified together. Either specify -kmeroff or -kmer")
       sys.exit(1) 


    if num_alignments != 250 and not use_blast:
       logger.error("-num_alignments can only be specified if -blast is specified.")
       sys.exit(1) 
        

    #check that programs are in the path
    if use_blast:
        check_program_in_path('blastn')
    if use_compress:
        check_program_in_path('plzip')
    if use_counts:
        check_program_in_path('jellyfish')


    TAXON_RANK = args.taxonRank

    if args.optimize and use_counts:
        logger.info("The classifier will be run with the following set of parameters: \nmin cl acc threshold: {0} \nuse counts: {6} \nuse blast: {1} \nuse compression: {2} \ntaxon: {3} \nperform kmer optimization: {4} \nk-mer specified: {5}".format(args.min_accu_threshold, args.use_blast, args.use_compression, args.taxonRank, args.optimize, args.specify_kmer, use_counts))
    elif args.optimize and not use_counts:
        logger.info("The classifier will be run with the following set of parameters: \nmin cl acc threshold: {0} \nuse counts: {6} \nuse blast: {1} \nuse compression: {2} \ntaxon: {3} \nperform optimization: {4} \nk-mer specified: {5}".format(args.min_accu_threshold, args.use_blast, args.use_compression, args.taxonRank, args.optimize, args.specify_kmer, use_counts))


    if num_threads > max_cpus:
        logger.info("Warning. You have specified {0} processes and the machine has a maximum of {1} processors.".format(num_threads, max_cpus))
    else:
        logger.info("Number of CPUs specified: {0}".format(num_threads))

    if subprocess.call(['ls', os.path.abspath(args.dbdirectory) + '/cl_train_' + TAXON_RANK + '.p'], stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
       #subprocess.call(['rm', os.path.abspath(args.dbdirectory) + '/cl_train_' + TAXON_RANK + '.p'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
       logger.info("Get taxon classes for the training set.") 
       get_train_classes(args.dbdirectory)


    # IB    
    if use_blast:
        if num_alignments != 250:
            logger.info("-blast parameter -num_alignments specified: {0}".format(num_alignments))
        blast_data(args.filefasta, args.dbdirectory, args.fastatype, which = 'test')


    if not args.optimize and not args.specify_kmer:
         est_accu = pickle.load(open(os.path.abspath(args.dbdirectory) + '/optimum_kmer/estimated_model_accuracy_' + TAXON_RANK + '.p', 'rb'))
         logger.info("The CSSS model classification accuracy: {0} was estimated using the following set of parameters: \nkmer opt = {1} \nmin cl acc threshold = {2} \ntaxon: {3} \nselected similarity measures: {4}".format(est_accu['accu'], est_accu['kmer'], est_accu['thr'], est_accu['taxon'], ",".join(est_accu['measures'])))
         if args.taxonRank != est_accu['taxon']:
                     logger.error("Taxon rank provided {0} is different from the one used for estimating the model accuracy {1}. Exiting.".format(args.taxonRank, est_accu['taxon']))
                     sys.exit(1)             
         if use_counts:
             try:
                 results_opt_counts = pickle.load(open(os.path.abspath(args.dbdirectory) + '/optimum_kmer/optimized_kmer_results_' + TAXON_RANK + '.p', 'rb'))
                 logger.info("Using the pre-calculated value for the optimal k-mer size from {0}".format(os.path.abspath(args.dbdirectory) + '/optimum_kmer/optimized_kmer_results_' + TAXON_RANK + '.p'))
                 if args.taxonRank != results_opt_counts.keys()[0]:
                     logger.error("Taxon rank provided {0} is different from the one used in the optimization {1}. Exiting.".format(args.taxonRank, results_opt_counts.keys()[0]))
                     sys.exit(1)
             except IOError:
                 logger.error("Optimized train data for counts are not available - will need to run model parameter optimization first - please see the -help page")
                 sys.exit(1)
         else:
            results_opt_counts = None
            logger.info("Similarity measures based on kmer counts will not be run!")
            if ('d_eucl' in est_accu['measures'] or 'd_js' in est_accu['measures']):
               logger.error("Please rerun the optimization without kmer counts to estimate the correct model accuracy") 
         if use_blast:
             try:
                 results_blast_train = pickle.load(open(os.path.abspath(args.dbdirectory) + '/optimum_kmer/BLAST_train_results_' + TAXON_RANK + '.p', 'rb')) 
             except IOError:
                 logger.error("Estimates for the model accuracy with BLAST are not available you will need to run model parameter optimization first - please see the -help page")
                 sys.exit(1)
         else:
             results_blast_train = None
             logger.info("BLAST will not be run!")
             if 'd_blast' in est_accu['measures']:
                 logger.error("Please rerun the optimization without BLAST to estimate the correct model accuracy")
         if use_compress:
             try:
                 results_cb_train = pickle.load(open(os.path.abspath(args.dbdirectory) + '/optimum_kmer/cb_train_results_' + TAXON_RANK + '.p', 'rb')) 
             except IOError:
                 logger.error("Compression train data for optimization not available will need to run the optimization of the model parameters first - please see the -help page")
                 sys.exit(1)
         else:
             results_cb_train = None
             logger.info("Compression based measure will not be run!")
             if 'd_cb' in est_accu['measures']:
                logger.error("Please rerun the optimization without CB (compression) to estimate the correct model accuracy")             
    else:
        if not use_counts and not use_compress and use_blast:
            logger.info("Estimated accuracy will be based on BLAST only, kmer and compresson based measures have been turned off.")
            # IB 
            build_optimization_set_blast(args.filefasta, args.dbdirectory, args.fastatype)
            results_opt_counts = None
        elif not use_counts and use_compress and use_blast:
            logger.info("Estimated accuracy will be based on the BLAST and compression based measures only, kmer based measure has been turned off.")
            # IB 
            build_optimization_set_blast(args.filefasta, args.dbdirectory, args.fastatype)
            results_opt_counts = None
        elif use_counts:
            # IB 
            results_opt_counts = build_optimization_set(args.filefasta, args.dbdirectory, args.fastatype, args.kmer_min, args.optimize, args.specify_kmer)
            logger.info("Estimated the kmer (kmer = {0}) classification accuracy: {1}".format(results_opt_counts[TAXON_RANK][0], mean(results_opt_counts[TAXON_RANK][1].values())))
            # IB 
            #pass
        if use_blast:
            logger.info("Estimating the BLAST accuracy using the training set.")
            logger.debug("BLAST {0} against the BLAST database".format(os.path.abspath(args.dbdirectory) + '/optimum_kmer/optimized_train.fasta'))
            # IB 
            blast_data(os.path.abspath(args.dbdirectory) + '/optimum_kmer/optimized_train.fasta', args.dbdirectory, args.fastatype, which = 'train')
            results_blast_train = evaluate_blast_train_data(args.dbdirectory, args.fastatype)
            logger.info("Estimated the BLAST classification accuracy: {0}".format(results_blast_train['d_blast'][0]))
            if not use_counts and not use_compress:
                opt_kmer_dir = os.path.abspath(args.dbdirectory) + '/optimum_kmer/'
                subprocess.call(['rm', opt_kmer_dir + 'estimated_model_accuracy_' + TAXON_RANK + '.p'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                parameters = {}
                parameters['accu'] = results_blast_train['d_blast']
                parameters['thr'] =  min_accu_threshold
                parameters['measures'] = ['d_blast']
                parameters['taxon'] = TAXON_RANK
                parameters['kmer'] = None
                pickle.dump(parameters, open(opt_kmer_dir  + 'estimated_model_accuracy_' + TAXON_RANK + '.p', 'wb'))
        else:
            results_blast_train = None
        if use_compress:
            logger.info("Estimating the accuracy of the compression based measure using sequences in the training set.")
            logger.debug("Using {0} against compressed sequence in {1}".format(os.path.abspath(args.dbdirectory) + '/optimum_kmer/optimized_train.fasta', os.path.abspath(args.dbdirectory) + '/compressdb/'))
            # IB 
            results_cb_train = evaluate_compression_train_data(os.path.abspath(args.dbdirectory)) 
            logger.info("Estimated CB classification accuracy: {0}".format(results_cb_train['d_cb'][0]))
            if not use_counts and not use_blast:
                opt_kmer_dir = os.path.abspath(args.dbdirectory) + '/optimum_kmer/'
                subprocess.call(['rm', opt_kmer_dir + 'estimated_model_accuracy_' + TAXON_RANK + '.p'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                parameters = {}
                parameters['accu'] = results_cb_train['d_cb']
                parameters['thr'] =  min_accu_threshold
                parameters['measures'] = ['d_cb']
                parameters['taxon'] = TAXON_RANK
                parameters['kmer'] = None
                pickle.dump(parameters, open(opt_kmer_dir  + 'estimated_model_accuracy_' + TAXON_RANK + '.p', 'wb'))               
        else:
            results_cb_train = None       

            
        if (not use_counts and not use_compress and use_blast) or (not use_counts and use_compress and not use_blast):
            # case 1: use blast only, or case 2: use compression only 
            est_accu = pickle.load(open(os.path.abspath(args.dbdirectory) + '/optimum_kmer/' + 'estimated_model_accuracy_' + TAXON_RANK + '.p', 'rb'))
            logger.info("The CSSS model classification accuracy: {0} was estimated using the following set of parameters: \nkmer opt = {1} \nmin cl acc threshold = {2} \ntaxon: {3} \nselected similarity measures: {4}".format(est_accu['accu'], est_accu['kmer'], est_accu['thr'], est_accu['taxon'], est_accu['measures']))
        else: 
            # in all the other cases estimate the overall accuracy (with at least two measures)
            est_accu = estimate_accuracy(args.dbdirectory)     
            logger.info("The CSSS model classification accuracy: {0} was estimated using the following set of parameters: \nkmer opt = {1} \nmin cl acc threshold = {2} \ntaxon: {3} \nselected similarity measures: {4}".format(est_accu['accu'], est_accu['kmer'], est_accu['thr'], est_accu['taxon'], ",".join(est_accu['measures']))) 


    file_fasta_test = os.path.abspath(args.filefasta)
    
    # perform classification 
    s_t = timeit.default_timer()
    perform_classification(file_fasta_test, args.dbdirectory, args.fastatype, results_opt_counts, results_blast_train, results_cb_train, est_accu)
    elapsed_time = timedelta(seconds=int(timeit.default_timer() - s_t))
    dtime = datetime(1,1,1) + elapsed_time


    logger.info("Total number of sequences processed: {0} in {1}d {2}h {3}m {4}s".format(count_seq_processed.value, dtime.day-1, dtime.hour, dtime.minute, dtime.second))


if __name__ == '__main__':
    print 'This program should be run as part of the cssscl package:\n\t$ cssscl classify -h\n\tor\n\t$ /path/to/cssscl/bin/cssscl -h'
