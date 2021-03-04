# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
__author__ = 'zerodel'

import argparse

if __name__ == "__main__":
    your_parser = argparse.ArgumentParser()
    your_parser.add_argument("-o", default="", help="some output")
    your_parser.add_argument("files", nargs="+", help="list of files you need")

    args = your_parser.parse_args()
    print(args.files)