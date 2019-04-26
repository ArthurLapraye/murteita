#!/usr/bin/python

import json
#from Levenshtein import distance as lev
import sys
import logging
from Bio.pairwise2 import format_alignment, align
import functools

from collections import defaultdict

from sknCONLL import load_skn

alphabet={'ᴱ', 'h', '’', 'ᴉ', 'ϑ', 'Ä', 'c', 'ᵏ', 'J', '˘', 'ˉ', 'E', 'B', 'a', 'ᵉ', '₍', '0', 'Ǹ', 'ŧ', '”', '“', '1', '»', 'ᵃ', 'ĕ', 'u', 'i', '%', 'q', '̬', 'd', 'ò', 'ê', 'P', '‿', 'R', 'ᶠ', 'ʸ', '"', 'ᴋ', 'Z', '–', 'ᵑ', 'Y', 'ʷ', '?', 'M', ' ', 'û', '͔', 'ẅ', 'ʜ', '̮', ']', 'ʲ', 'ŏ', 'Ŏ', '^', 'ᴹ', 'Ì', '\uf21d', ')', 'á', '̳', 'ẁ', '\u0ff1', 'w', '+', 'ᴍ', 'Ĺ', 'å', 'ᴡ', 'χ', 'A', ',', '.', 'é', '*', '§', 'ŕ', 'N', '_', 'I', 'ŋ', 'b', 'ᵐ', 'ỳ', 'e', 'ŷ', 'k', 'ɦ', 'p', 'Ŭ', 'ə', 'ʟ', '[', 'ɢ', 'V', 'H', 'r', '°', 'W', 'ä', '̇', 'ḙ', 'ᴰ', 'ɪ', '̰', 'ʏ', 'ð', '\u1cff', 'ᴅ', 'ᴇ', '̹', 'ᵗ', 'ȯ', 'ì', '̲', 'ʀ', 'ᴺ', 'C', 's', 'U', 'ḭ', 'ᴵ', 'ᴠ', 'ɔ', 'ʾ', 'ʼ', 'ᴬ', 'y', 'ᴜ', '͇', 'Ń', 'ᵎ', 'ʙ', 'ᵻ', 'ö', '!', '5', '́', 'è', 'f', 'î', 'И', 't', '̀', 'Ö', 'ù', 'ị', '̄', '=', 'ṋ', '∼', 'β', 'δ', '-', 'à', 'D', 'ᴼ', 'Ȯ', 'Î', 'ʳ', 'v', 'ĺ', '(', 'ⁱ', 'ᴛ', 'ė', '~', '̥', 'ᵒ', '̭', 'ᵛ', 'ᵊ', 'ᵁ', '̆', 'È', 'ⁿ', '̜', 'ˈ', 'ú', 'ʰ', 'Å', '̤', 'ᴊ', '\xa0', '̈', 'ʿ', 'ś', 'n', 'ń', 'ḱ', 'ɐ', 'Ò', 'F', 'ǹ', 'ˢ', 'í', 'š', 'ǰ', '̗', 'ᵘ', 'O', 'ŭ', 'ᴮ', 'T', 'ᴏ', 'o', 'ĭ', 'ḛ', 'ˡ', 'ó', 'ᴶ', 'À', '\u0fff', 'ᴎ', 'l', 'z', "'", 'g', 'ₒ', '#', 'j', 'ṵ', 'ô', '<', 'm', 'L', '´', 'K', '`', 'ɴ', 'G', '̂', '͕', 'ɛ', 'ˀ', 'S', 'ᵖ', 'â', '3', 'ă', 'ᴀ'}


consonnes= "ztpqdfghjklmnbvcx"
voyelles="aeiouyäö"

same=[("à","a")]

NORMALFLAG=True

#Utiliser functools.partial pour avoir accès au mot d'avant et au mot d'après 
#prototype futur = matchfunction (mot, motprecedent,motsuivant,r1,r2) => partial =>matchfunc

standnorm=defaultdict(set)
normraw=defaultdict(set)
paradigm=defaultdict(set)

def matchfunc(stand,dial):
	r1=stand
	r2=dial

	if r1==r2:
		return 10
	elif (r1,r2) in same:
		return 10
	elif r1=="d" and r2=="r":
		return 7
	elif r1=="ö" and r2 in ["ä"]:
		return 7
	elif r1=="e" and r2=="i":
		return 4
	
	else:
		return -10

#matchfunc=functools.partial(matchfunction,"","","")

def gapfunction(x,y):
	return 0

if __name__=="__main__":
	#SKN=load_skn(sys.argv[1])
	
	b=0
	beginning=0
	end=860000
	tokennumber=0
	tokennumfile=0
	#ALPHABET=set()
	sentence_id=None
	fileoutput=""
	outfolder="sfeatures/"
	title=None
	sdict=None
	
	nassimil=0
	oassimil=0
	passit=True
	slength=0
	
	
	for rank in range(beginning,end,10000):
		if rank == 0:
			continue
		SKN=load_skn("skn_corpus_%s-%s.json" % (b,rank))
		tokennumfile=0
		b=rank
		for x,y in enumerate(SKN["kwic"]):
			word=y["tokens"][0]
			sentence=y["structs"]
			
			slength+=1
			
			try:
				nextword=SKN['kwic'][x+1]["tokens"][0]
			except IndexError as I:
				try:
					SKN2=load_skn("skn_corpus_%s-%s.json" % (b,rank+10000))
					nextword=SKN2['kwic'][0]["tokens"][0]
				except FileNotFoundError as e:
					passit=False
			
			if passit:
				if (word["normalized"][-1]==nextword["normalized"][0] ) and word["word"][-1]=="n" and word["normalized"] != "n":
					#print(word["normalized"],nextword["normalized"],word["word"])
					nassimil+=1
				elif (word["normalized"][-1]==nextword["normalized"][0] ) and word["word"][-1] != "n" and word["word"][-1]!=word["normalized"][-1]:
					print(word["normalized"],nextword["normalized"],word["word"])
					oassimil+=1
				
			
			
			if sentence_id != sentence['sentence_origid']:
				fileoutput+="\n"+"\t".join(map(str,[sentence_id,slength,nassimil,oassimil]))
				sentence_id=sentence['sentence_origid']
				oassimil=0
				slength=0
				nassimil=0

			
			if title != sentence['text_title']:
				if title:
					with open(outfolder+title+".txt","w") as out:
						out.write(fileoutput)
				fileoutput="id\tlen\tnass\toass"
				title=sentence['text_title']
				print(title)
			
			
			headid=int(word["dephead"])
			wordid=int(word["id"])
		

