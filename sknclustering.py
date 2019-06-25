#!/usr/bin/python

import sys
import os
import pandas as pd
from sknfeatures import HEADERS


if __name__=="__main__":
	murteita=dict()
	for elem in sys.argv[1:]:
		feats=dict()
		print(elem)
		x=pd.read_csv(elem,sep="\t")
		
		data=x.sum(axis=0,skipna = True)
		lpho=data["longueurpho"]
		
		for e in HEADERS:
			if e.endswith("ass"):
				if e[0]=="t":
					feats[e]=(data[e]/data["word_t"])
				elif e[0] == "n":
					feats[e]=data[e]/data["word_n"]
				elif e.startswith("glott"):
					feats[e]=data[e]/data["word_glott"]
								
			
			elif e.startswith("syn_")
				feats[e]=data[e]/data["props"]		
			
			elif e.startswith("var_"):
				feats[e]=(data[e]/lpho)
		
		
				
