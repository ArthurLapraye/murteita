#!/usr/bin/python

import sys
import os
import pandas as pd
import folium
import csv
from folium.plugins import HeatMap


	

if __name__=="__main__":
	
	murteita=dict()
	finalfeats=set()
	kyl=dict()
	
	with open("parishes.csv") as fi:
		for (ville,n,lat,lo) in csv.reader(fi):
			if ville != "Nom":
				print(ville)
				kyl[n]=(ville,float(lat),float(lo))
	
	for elem in sys.argv[1:]:
		feats=dict()
		nimi=" ".join(os.path.basename(elem).split(" ")[0:2])
		nimi2=nimi.split(" ")[0]
		print(nimi2)
		
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
				
			elif e.startswith("var2_"):
				finalfeats.add(e)
				feats[e]=(data[e]/(2*lpho))
		
		murteita[nimi2]=feats
	
	features=list(sorted(finalfeats))
	
	with open("sortie.arff","w") as outfile:
		outfile.write("@relation murteita\n")
		
		for f in features:
			outfile.write("@attribute \""+f+"\" "+"numeric\n")
		
		outfile.write("@attribute nimi { \""+"\",\"".join(murteita.keys())+"\"}\n")
		
		outfile.write("@data\n")
		for elem in murteita:
			outfile.write(",".join((str(10000*murteita[elem].get(e, 0)) for f in features)) +",\""+elem+"\"\n")
	
	
	with open("sortie.csv","w") as outfile:
		outfile.write("nimi,lat,lo,"+",".join(features)+"\n")
		
		for elem in murteita:
			outfile.write(",".join(str(kyl[elem]))+",".join( ( str(murteita[elem].get(e, 0)) for f in features ) )+"\n")
	
	
		
		
		
		
		
	
	
		
		
		
		
	
	
		
		
				
		
		
		
				
