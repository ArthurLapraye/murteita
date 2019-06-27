#!/usr/bin/python

import sys
import os
import pandas as pd
import folium
import csv

def regeneratemap(dic):
	carte = folium.Map(location=[65,25],zoom_start=5,control_scale=True)
	for elem in dic:
		folium.Marker(location=dic[elem],
		popup=elem,
		icon=folium.Icon(color="red",icon="ok-sign")).add_to(dialectloc)
	
	return carte
	

if __name__=="__main__":
	
	murteita=dict()
	finalfeats=set()
	kyl=dict()
	
	with open("parishes") as fi:
		for (ville,n,lat,lo) in csv.reader(fi):
			if ville != "Nom":
				print(ville)
				kyl[n]=(ville,float(lat),float(lo))
				
	
	dialectloc=regeneratemap(kyl)
	
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
	
	f=list(sorted(finalfeats))
	
	with open("sortie.arff","w") as outfile:
		outfile.write("@relation murteita\n")
		
		for e in f:
			outfile.write("@attribute \""+e+"\" "+"numeric\n")
		
		outfile.write("@attribute nimi { \""+"\",\"".join(murteita.keys())+"\"}\n")
		
		outfile.write("@data\n")
		for elem in murteita:
			outfile.write(",".join((str(10000*murteita[elem].get(e, 0)) for e in f)) +",\""+elem+"\"\n")
	
	for trait in f:
		maximum=max(murteita[elem].get(f,0) for elem in murteita)
		
		locations=list()
		
		for ville in kyl:
			lat=kyl[ville][1]
			lo=kyl[ville][2]
			featvalue=murteita[ville][f]
			locations.append( (lat,lo,featvalue))
			
		HeatMap(locations,max_val=maximum).add_to(dialectloc)
		dialectloc.save("/tmp/sortie_"+f)
		dialectloc=regeneratemap(kyl)
		
		
		
		
	
	
		
		
		
		
	
	
		
		
				
		
		
		
				
