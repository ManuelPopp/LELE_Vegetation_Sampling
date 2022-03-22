#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 21:46:55 2021

@author: manuel
"""
import argparse

def parseArguments():
    parser = argparse.ArgumentParser()
    # Positional mandatory arguments
    parser.add_argument("inp", help = "Input file",\
                        type = str)
    parser.add_argument("outp", help = "Output file",\
                        type = str)
    # Parse arguments
    args = parser.parse_args()
    return args
if __name__ == "__main__":
    # Parse the arguments
    args = parseArguments()

in_path = args.inp
out_path = args.out

import pandas as pd
data = pd.read_html(in_path)
data.to_csv(out_path)
