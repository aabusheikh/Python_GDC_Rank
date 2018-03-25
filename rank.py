# Copyright (C) 2018  Ahmad A. A. (https://github.com/bbpgrs/)

import pandas as pd
import common as cmn
import os
import logging
import operator
import math


def col_rank_genes():
    """
    """
    combined_file_path = os.path.join(cmn.DL_DIR, cmn.PTEN_COMB_FNAME)
    if os.path.isfile(combined_file_path):
        logging.info("Found combined PTEN correlation file, generating rankings ...\n")

        df = pd.read_csv(combined_file_path, sep="\t", index_col=0)

        column_ranks_file = os.path.join(cmn.DL_DIR, cmn.COL_RANK_FNAME)
        col_ranks = pd.DataFrame(index=df.index, columns=df.columns)
        #i = 0

        for col in df.columns:
            logging.info("Generating rankings for '%s' ..." % col)

            col_rank_dict = {}
            for row in df.index:
                if not math.isnan(float(df.loc[row][col])):
                    col_rank_dict[row] = abs(float(df.loc[row][col]))
                else:
                    col_rank_dict[row] = 0

            col_rank = enumerate(sorted(col_rank_dict.items(), key=operator.itemgetter(1), reverse=True))
            #col_ranks.insert(i, col, [t[0] for t in col_rank])
            # i += 1
            for r, t in col_rank:
                col_ranks.loc[t[0]][col] = r

            logging.info("Done.\n")

        logging.info("Writing column rankings to file ...")
        col_ranks.to_csv(column_ranks_file, sep="\t")
        logging.info("Done.\n")

    else:
        logging.warning("Combined PTEN correlation file not found.\n")


def global_rank_genes():
    """

    :return:
    """
    column_ranks_file = os.path.join(cmn.DL_DIR, cmn.COL_RANK_FNAME)

    if os.path.isfile(column_ranks_file):
        logging.info("Creating global ranking ...")

        df = pd.read_csv(column_ranks_file, sep="\t", index_col=0)

        rank_dict = {}
        n = len(df.index)
        i = 1

        for row in df.index:
            logging.info("Processing global ranking for [Gene %s out of %s] '%s' ..." % (i, n, row))

            r_vals = df.loc[row].values

            ''' # Less efficient, but looks cleaner. also, O(3n) instead of O(n), so still O(n)
            r_sum = sum(r_vals)
            r_min = min(r_vals)
            r_max = max(r_vals)
            '''

            # Might be slightly more efficient, but looks less 'clean' -> O(n) instead of O(3n)
            r_sum, r_min, r_max = 0, r_vals[0], r_vals[0]
            for val in r_vals:
                r_sum += val
                if val > r_max:
                    r_max = val
                if val < r_min:
                    r_min = val


            rank_dict[row] = (r_sum//len(df.loc[row].values), r_min, r_max)
            i += 1

            logging.info("Done.\n")

        rank = enumerate(sorted(rank_dict.items(), key=operator.itemgetter(1)))

        logging.info("Writing global ranking to file ...")

        out_lines = ["\t".join(("Rank", "Gene_Code", "Avg_Col_Rank", "Min_Col_Rank", "Max_Col_Rank"))]+\
                    ["\t".join((str(r), str(t[0]), str(t[1][0]), str(t[1][1]), str(t[1][2]))) for r, t in rank]
        global_rank_file = os.path.join(cmn.DL_DIR, cmn.RANK_FNAME)
        with open(global_rank_file, 'w') as out_file:
            out_file.write("\n".join(out_lines))

        logging.info("Done.\n")

    else:
        logging.warning("Column rankings file not found.\n")


def run(g=False):
    """
    """
    if not g:
        col_rank_genes()
    global_rank_genes()
