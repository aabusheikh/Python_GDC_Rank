# Copyright (C) 2018  Ahmad A. A. (https://github.com/bbpgrs/)

import os
import errno
import logging

"""
Common code (global variables, functions) that will be used by multiple modules of this program
for configuration or functionality
"""

# excluded cancer types
SKIP_CTYPE = ["Nervous System"]

# excluded gender
SKIP_GENDER = ["unknown"]

# PTEN code
PTEN_GENE_CODE = "ENSG00000171862.8"

# correlation table file name template/pattern
CORR_FNAME = "%s_%s_%s_%s_correlation.txt"

# combined PTEN correlation file name
PTEN_CORR_FNAME = "PTEN_correlation.txt"

# combined PTEN correlation file name
PTEN_COMB_FNAME = "PTEN_combined.txt"

# combined rankings file
COL_RANK_FNAME = "column_ranks.txt"

# avg correlation file
AVG_CORR_FNAME = "PTEN_combined.fixed.txt"

# ensemble to gene name dictionary
GENE_ID_DICT_FNAME = "mart_export_Gene_Name_Id.txt"

# other paper's target sequence analysis results
OTHER_PAPER_FNAME = "PTEN_ceRNAs.txt"

# global ranking file
RANK_FNAME = "ranks.txt"

# absolute directory path to download files to (recommended if project files are stored in
# a space-sensitive location such as a cloud drive directory (Dropbox, OneDrive, etc...)
ABS_DL_DIR = os.path.abspath(os.path.join(os.path.expanduser('~'), 'Downloads', 'GDC_Downloads'))

# relative directory path to download files to
# should be in .gitignore
REL_DL_DIR = 'GDC_Downloads'

# chosen downloads directory
DL_DIR = ABS_DL_DIR

# gene/mirna list file
RNA_LIST_PATH = os.path.join(DL_DIR, "rna_list.txt")
MIRNA_LIST_PATH = os.path.join(DL_DIR, "mirna_list.txt")

# log file directory
LOG_DIR = 'logs'

# logging level
LOG_LEVEL = logging.DEBUG

# logging message format
LOG_MSG_FORMAT = '%(asctime)s - %(levelname)s:%(message)s'

# logging time format
LOG_TIME_FORMAT = '%H:%M:%S'

# log file name template
LOG_FILE_NAME = 'logfile_%s.log'

# log file name date/time format
LOG_FILE_NAME_TIME = '%m%d-%H%M%S'


def make_dir(directory):
    """
    Recursively create directories so that the directory path provided is created

    :param directory: String directory path to create
    """
    if not os.path.exists(directory):
        try:
            os.makedirs(directory)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


def list_dir(directory):
    """
    wrapper for os.listdir()

    :param directory:
    :return:
    """
    if os.path.isdir(directory):
        return os.listdir(directory)
    else:
        return []
