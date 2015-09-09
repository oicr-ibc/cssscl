import pymongo
import getpass
import os
import base64
import ConfigParser
import sys

from database import *
from pymongo.errors import DuplicateKeyError

db, logger = None, None

def setup_config(args):
    '''Saves MongoDB settings to a configuration file'''

    config = ConfigParser.SafeConfigParser()

    config.add_section('MongoDB')

    print 'Please enter the settings for your MongoDB server:'
    config.set('MongoDB', 'host', args.host or raw_input('Host [localhost]: ') or 'localhost')
    config.set('MongoDB', 'port', args.port or raw_input('Port [27017]: ') or '27017')
    config.set('MongoDB', 'database', args.database or raw_input('Database [cssscl]: ') or 'cssscl')
    #config.set('MongoDB', 'username', args.username or raw_input('Username [none]: '))
    #config.set('MongoDB', 'password', args.password or getpass.getpass('Password [none]: '))

    # Writing our configuration file
    with open(os.path.expanduser('~/.cssscl/cssscl.cfg'), 'wb') as configfile:
        config.write(configfile)


def main(args):
    '''Setup MongoDB for use by cssscl'''

    global db, logger

    logger = args.logging.getLogger(__name__)

    # Setup config files
    setup_config(args)

    db = connect(args)

    logger.info('Done!')


if __name__ == '__main__':
    print 'This program should be run as part of the cssscl package:\n\t$ cssscl configure -h\n\tor\n\t$ /path/to/cssscl/bin/cssscl configure -h'
