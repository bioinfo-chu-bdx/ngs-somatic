#!/usr/bin/python
import logging
import time

# set up logging to file
logging.basicConfig(level=logging.INFO,format='%(levelname)-7s %(message)s',filename='/DATA/work/scripts/test.log',filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.INFO)
# add the handler to the root logger
logging.getLogger('').addHandler(console)

logging.info('Start time= %s \n'% (time.strftime("%m/%d/%Y %H:%M:%S")))

logging.info('%s une info'% (time.strftime("[%H:%M:%S]")))
logging.warning('%s un avertissement'% (time.strftime("[%H:%M:%S]")))
logging.error('%s une erreur'% (time.strftime("[%H:%M:%S]")))
