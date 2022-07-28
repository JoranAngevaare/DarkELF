import sys, os, glob
from darkelf import logger


def list_all(path=os.path.dirname(__file__)+"/../data/"):
    """
    List available target materials in data directory
    """

    for file in glob.glob(path+"*"):
        logger.debug('\t',file.split("/")[-1])

def files(target,path=os.path.dirname(__file__)+"/../data/"):
    """
    List available config files, epsilon grids, and density of states files in data directory

    target: string
    """
    logger.debug('Available configuration files: ')
    for file in glob.glob(path +str(target)+"/*yaml"):
            logger.debug('\t',file.split("/")[-1])
    logger.debug(" ")
    logger.debug('Available data for epsilon: ')
    for file in glob.glob(path +str(target)+"/*dat"):
        if ('DoS' not in file) and ('pDoS' not in file) and ('Fn' not in file) and ('Zion' not in file):
            logger.debug('\t',file.split("/")[-1])
    logger.debug(" ")
    logger.debug('Available data for phonon (partial) density of states: ')
    for file in glob.glob(path +str(target)+"/*DoS.dat"):
            logger.debug('\t',file.split("/")[-1])
    logger.debug(" ")
    logger.debug('Available data for Fn(omega) functions: ')
    for file in glob.glob(path +str(target)+"/*Fn.dat"):
            logger.debug('\t',file.split("/")[-1])
    logger.debug(" ")
    logger.debug('Available data for form factors: ')
    for file in glob.glob(path +str(target)+"/*Zion.dat"):
            logger.debug('\t',file.split("/")[-1])
