#!/usr/bin/python

import sys
import os
import pandas as pd
from sknfeatures import HEADERS


if __name__=="__main__":
	murteita=dict()
	finalfeats=set()
	for elem in sys.argv[1:]:
		feats=dict()
		nimi=" ".join(os.path.basename(elem).split(" ")[0:2])
		print(nimi)
		
		x=pd.read_csv(elem,sep="\t")
		
		data=x.sum(axis=0,skipna = True)
		lpho=data["longueurpho"]
		
		for e in x.columns:
			if e.endswith("ass"):
				finalfeats.add(e)
				if e[0]=="t":
					feats[e]=(data[e]/data["word_t"])
				elif e[0] == "n":
					feats[e]=data[e]/data["word_n"]
				elif e.startswith("glott"):
					feats[e]=data[e]/data["word_glott"]
				
				elif e[0]== "o":
					feats[e]=data[e]/data['longueur']
				
				elif e.startswith("mp"):
					feats[e]=data[e]/data["word_n"]
			
			elif e.startswith("syn_"):
				finalfeats.add(e)
				feats[e]=data[e]/data["props"]		
			
			elif e.startswith("var_"):
				finalfeats.add(e)
				feats[e]=(data[e]/lpho)
		
		murteita[nimi]=feats
	
	f=list(sorted(finalfeats))
	
	with open("sortie.arff","w") as outfile:
		outfile.write("@relation murteita\n")
		
		for e in f:
			outfile.write("@attribute \""+e+"\" "+"numeric\n")
		
		outfile.write("@attribute nimi { \""+"\",\"".join(murteita.keys())+"\"}\n")
		
		outfile.write("@data\n")
		for elem in murteita:
			outfile.write(",".join((str(10000*murteita[elem].get(e, 0)) for e in f)) +",\""+elem+"\"\n")
			
	
	
		
		
		
		
	
	
		
		
				
		
		
		
				
