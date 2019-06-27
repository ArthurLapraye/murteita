#!/usr/bin/python

import folium

import csv

dialectloc = folium.Map(location=[65,25],zoom_start=5,control_scale=True)

kyl=dict()

with open("parishes") as fi:
	for (ville,_,lat,lo) in csv.reader(fi):
		if ville != "Nom":
			print(ville)
			kyl[ville]=(float(lat),float(lo))

for elem in kyl:
	folium.Marker(location=kyl[elem],
	popup=elem,
	icon=folium.Icon(color="red",icon="ok-sign")).add_to(dialectloc)
	
dialectloc.save("/tmp/sortie.html")
		
