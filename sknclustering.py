#!/usr/bin/python

import sys
import os

import pandas as pd



if __name__=="__main__":
	for elem in sys.argv[1:]:
		print(elem)
		x=pd.read_csv(elem,sep="\t")
		
		print(x.sum(axis=0,skipna = True)["word_n"] )
