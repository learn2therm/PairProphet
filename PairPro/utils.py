"""
The following is importable random utilites.
You will find:
- logger function
- pairwise sequence builder
"""
import pandas as pd
import itertools

# import distributed
import logging
logger = logging.getLogger(__name__)


def start_logger_if_necessary(logger_name: str, log_file: str, log_level,
                              filemode: str = 'a', worker: bool = False):
    """Quickly configure and return a logger that respects parallel processes.

    Parameters
    ----------
    logger_name : str
        name of logger to start or retrieve
    log_file : str
        path to file to log to
    log_level
        log level to respect
    worker: str
        name of worker using this logger
    filemode : str
        mode to apply to log file eg "a" for append
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(log_level)
    fh = logging.FileHandler(log_file, mode=filemode)
    if worker:
        # worker_name = distributed.get_worker().name
        # fh.setFormatter(logging.Formatter('%(filename)s %(worker)s - %(asctime)s %(levelname)-8s %(message)s'))
        # if len(logger.handlers) == 0:
        #     logger.addHandler(fh)
        # else:
        #     logger.handlers[-1] = fh
        # logger = logging.LoggerAdapter(logger, extra={'worker': worker_name})
        pass
    else:
        fh.setFormatter(logging.Formatter('%(filename)s - %(asctime)s %(levelname)-8s %(message)s'))
        if len(logger.handlers) == 0:
            logger.addHandler(fh)
        else:
            logger.handlers[-1] = fh
    return logger


def make_pairs(seq1_list, seq2_list, seq1_name='seq1', seq2_name='seq2',
               csv_path='./paired_seqs.csv', save=True):
    '''
    Function for building a combinatorial set of sequences from two lists.

    Args:
        seq1_list (list): List of protein sequence strings
        seq2_list (list): List of protein sequence strings
        seq1_name (str): Column name for first sequence column
        seq2_name (str): Column name for second sequence column
        csv_path (str): Path for saved .csv file
        save (bool): Saves paired sequences as .csv when True

    Returns:
        combined_df (pd.DataFrame): A dataframe with rows as all possible
        sequence pairs.
    '''
    combined = list(itertools.product(seq1_list, seq2_list))
    combined_df = pd.DataFrame(combined, columns=[seq1_name, seq2_name])

    if save is True:

        combined_df.to_csv(csv_path, index=False)

    return combined_df
