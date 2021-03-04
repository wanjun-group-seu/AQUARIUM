# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import logging
import os
import time
import datetime

__doc__ = '''
'''
__author__ = 'zerodel'

format_default = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')


def default_logger(logger_name):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)

    to_stdout = logging.StreamHandler()
    to_stdout.setFormatter(format_default)

    logger.addHandler(to_stdout)

    return logger


def log_to_file(logger, log_file_path):
    to_file = logging.FileHandler(log_file_path)
    to_file.setLevel(logging.DEBUG)
    logger.addHandler(to_file)
    return logger


def set_logger_file(logger, log_file, file_default="top_level.log"):
    path_abs = os.path.abspath(os.curdir)
    log_file_path = log_file if log_file else os.path.join(path_abs, file_default)
    return log_to_file(logger, log_file_path)


def raw_timestamp():
    x = datetime.datetime.now()
    return "{year}_{month}_{day}_{hour}_{minute}_{second}_{misec}".format(year=x.year,
                                                                          month=x.month,
                                                                          day=x.day,
                                                                          hour=x.hour,
                                                                          minute=x.minute,
                                                                          second=x.second,
                                                                          misec=x.microsecond)
