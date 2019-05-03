#!/usr/bin/python

import json
#from Levenshtein import distance as lev
import sys
import logging
from Bio.pairwise2 import format_alignment, align
from functools import partial

from collections import defaultdict

from sknCONLL import load_skn

alphabet={'ᴱ', 'h', '’', 'ᴉ', 'ϑ', 'Ä', 'c', 'ᵏ', 'J', '˘', 'ˉ', 'E', 'B', 'a', 'ᵉ', '₍', '0', 'Ǹ', 'ŧ', '”', '“', '1', '»', 'ᵃ', 'ĕ', 'u', 'i', '%', 'q', '̬', 'd', 'ò', 'ê', 'P', '‿', 'R', 'ᶠ', 'ʸ', '"', 'ᴋ', 'Z', '–', 'ᵑ', 'Y', 'ʷ', '?', 'M', ' ', 'û', '͔', 'ẅ', 'ʜ', '̮', ']', 'ʲ', 'ŏ', 'Ŏ', '^', 'ᴹ', 'Ì', '\uf21d', ')', 'á', '̳', 'ẁ', '\u0ff1', 'w', '+', 'ᴍ', 'Ĺ', 'å', 'ᴡ', 'χ', 'A', ',', '.', 'é', '*', '§', 'ŕ', 'N', '_', 'I', 'ŋ', 'b', 'ᵐ', 'ỳ', 'e', 'ŷ', 'k', 'ɦ', 'p', 'Ŭ', 'ə', 'ʟ', '[', 'ɢ', 'V', 'H', 'r', '°', 'W', 'ä', '̇', 'ḙ', 'ᴰ', 'ɪ', '̰', 'ʏ', 'ð', '\u1cff', 'ᴅ', 'ᴇ', '̹', 'ᵗ', 'ȯ', 'ì', '̲', 'ʀ', 'ᴺ', 'C', 's', 'U', 'ḭ', 'ᴵ', 'ᴠ', 'ɔ', 'ʾ', 'ʼ', 'ᴬ', 'y', 'ᴜ', '͇', 'Ń', 'ᵎ', 'ʙ', 'ᵻ', 'ö', '!', '5', '́', 'è', 'f', 'î', 'И', 't', '̀', 'Ö', 'ù', 'ị', '̄', '=', 'ṋ', '∼', 'β', 'δ', '-', 'à', 'D', 'ᴼ', 'Ȯ', 'Î', 'ʳ', 'v', 'ĺ', '(', 'ⁱ', 'ᴛ', 'ė', '~', '̥', 'ᵒ', '̭', 'ᵛ', 'ᵊ', 'ᵁ', '̆', 'È', 'ⁿ', '̜', 'ˈ', 'ú', 'ʰ', 'Å', '̤', 'ᴊ', '\xa0', '̈', 'ʿ', 'ś', 'n', 'ń', 'ḱ', 'ɐ', 'Ò', 'F', 'ǹ', 'ˢ', 'í', 'š', 'ǰ', '̗', 'ᵘ', 'O', 'ŭ', 'ᴮ', 'T', 'ᴏ', 'o', 'ĭ', 'ḛ', 'ˡ', 'ó', 'ᴶ', 'À', '\u0fff', 'ᴎ', 'l', 'z', "'", 'g', 'ₒ', '#', 'j', 'ṵ', 'ô', '<', 'm', 'L', '´', 'K', '`', 'ɴ', 'G', '̂', '͕', 'ɛ', 'ˀ', 'S', 'ᵖ', 'â', '3', 'ă', 'ᴀ'}


consonnes= "zδtpqdfghjklmnbvcxðsrvŋw"
voyelles="aeiouyäöəɛᴏ"
punctuation=" &~#\"'({[-|`_\^@)]°=+}$^¨%µ*!:;,?./§1234567890"

same=[("à","a")]

NORMALFLAG=True

#Utiliser functools.partial pour avoir accès au mot d'avant et au mot d'après 
#prototype futur = matchfunction (mot, motprecedent,motsuivant,r1,r2) => partial =>matchfunc

standnorm=defaultdict(set)
normraw=defaultdict(set)
paradigm=defaultdict(set)

#def gapfunc(x,y):
	

def matchfunc(stand,dial):
	r1=stand.lower()
	r2=dial.lower()
	
	if "ᴏ" == r1:
		r1="o"
	elif "ᴏ"== r2:
		r2 = "o"
	elif "δ"==r1:
		r1="ð"
	elif r2=="δ":
		r2="ð"
	
	spaces=(" ","-")
	
	pairesproches=[("t","d"),("d","r"),("ä","e"),("ö","ä"),("ö","y"),("o","a"),("o","u"),("i","j"),("t","s"),("ð","d"),("i","e"),("ə","i"),("ə","e"),("d","r")]
	
	if r1==r2:
		return 20
	else:
		for p in pairesproches:
			if r1 in p and r2 in p:
				return 7
		else:
			if r1 == "n" and r2 in consonnes:
				return 2
			else:
				return -10
	
	return -100

#matchfunc=functools.partial(matchfunction,"","","")

def gap_functionA(word, normalized,x,y):
	#print("A",word, normalized, x,y)
	#print(word[x:],normalized[x:])
	
	if 3 < x < len(normalized)  and normalized[x-2] == normalized[x] and normalized[x-1] in ["h","l"] and normalized[x] in voyelles:
		return -10-(y-1) if y > 1 else -1
	
	if x < len(normalized) and normalized[x-1]=="i" and normalized[x]=="j":
		if y > 1:
			return -10-(y-1)
		else:
			return -1
	
	if x > 2:
		if x < len(normalized) and normalized[x] == normalized[x-1]:
			if y > 1:
				return -10-(y-1)
			else:
				return -1
	
	return (-10 + -1*y)

def gap_functionB(word, normalized,x,y):
	#print("B",word, normalized, x,y)
	if x > 2:
		if x < len(normalized) and normalized[x] == normalized[x-1]:
			if y > 1:
				return -7-(0.7*(y-1))
			else:
				return -1
	
	return (-7-(0.7*y))

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
	
	HEADERS=["id","longueur","t_ass","n_ass","mp_ass","n_n","e_ass","o_ass","word_t","word_n","word_e","VS","align","interrog"]
	headershelp="""id : id de la phrase
				  longueur : longueur en mot de la phrase
				  t_ass nombre de t en fin de mot assimilés totalement
				  n_ass nombre de n en fin de mot assimilés totalement
				  mp_ass nombre de n assmilés à m devant p 
				  n_n nombre de n en fin de mot suivi de n initial
				 e_ass nombre de mot finissant en -e dans le standard et ayant leur glottale assimilée
				 o_ass autres assimilations de sandhi
				 word_t mots finissants en t
				 word_n
				 word_e
				 
				 VS : trait binaire SV ou VS
				 align : alignement utilisé parmi VSO, OSV, OVS, etc"""
	eassimil=0
	nassimil=0
	mpassimil=0
	nn=0
	
	tassimil=0
	oassimil=0
	passit=True
	slength=0
	
	word_n=0
	word_t=0
	word_e=0
	
	alignment="_"
	VSsentence = False
	sinter=False

	
	for rank in range(beginning,end,10000):
		if rank == 0:
			continue
		SKN=load_skn("skn_corpus_%s-%s.json" % (b,rank))
		tokennumfile=0
		b=rank
		for x,y in enumerate(SKN["kwic"]):
			word=y["tokens"][0]
			
			if word["word"] and word["word"] != word["normalized"]:
				for elem in word["normalized"]:
					elem=elem.lower()
					if elem not in voyelles:
						if elem not in consonnes:
							if elem not in punctuation:
								logging.error("non classé : "+elem)
				W=list(word["word"])
				N=list(word["normalized"])
				gapfunctionA=partial(gap_functionA,W,N)
				gapfunctionB=partial(gap_functionB,W,N)
				
				alignement=align.globalcc(W,N,matchfunc,gapfunctionA,gapfunctionB,gap_char=['-'])
				print("=")
				for elem in alignement:
					align1, align2, score, begin, end=elem
					print("".join(align1)+"\n"+"".join(align2)+"\n"+str(score))
			
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
			
			if passit and word["word"]:
				if word["word"][-1]=="n":
					word_n+=1
					if (word["normalized"][-1]==nextword["normalized"][0] ) and word["normalized"][-1] != "n":
						#print("N",word["normalized"],nextword["normalized"],word["word"])
						nassimil+=1
					elif word["normalized"][-1]=="m" and nextword["normalized"][0] == "p":
						mpassimil+=1
					elif nextword["normalized"][-1] == "n":
						nn += 1
					
					
				if word["word"][-1]=="t":
					word_t+=1	
					if  (word["normalized"][-1]==nextword["normalized"][0] ) and word["normalized"][-1] != "t":
						#print("T",word["normalized"],nextword["normalized"],word["word"])
						tassimil+=1
				if word["word"][-1]=="t":
					word_e+=1	
					if  (word["normalized"][-1]==nextword["normalized"][0] ) and word["normalized"][-1] != "e":
						#print("E",word["normalized"],nextword["normalized"],word["word"])
						eassimil+=1
				
				if (word["normalized"][-1]==nextword["normalized"][0] ) and word["word"][-1]!=word["normalized"][-1]:
					#print("O",word["normalized"],nextword["normalized"],word["word"])
					oassimil+=1
			elif not word["word"]:
				logging.error("Problem at "+title+" "+str(word))
			
			if len(word["msd"]) > 1:
				#print(word["msd"])
				msd=dict([tuple(x.split("_",1)) for x in word["msd"].split("|")])
				
				if msd.get("SUBCAT",None) == "Interr":
					sinter=True
			
			if "subj" in word["deprel"]:
				VSsentence=(int(word["id"]) > int(word["dephead"] ))
					
				
			
			
			if sentence_id != sentence['sentence_origid']:
				fileoutput+="\n"+"\t".join(map(str,[sentence_id,slength,tassimil,nassimil,mpassimil,nn,eassimil,oassimil,word_t,word_n,word_e,str(VSsentence).upper(),alignment,str(sinter).upper()]))
				sentence_id=sentence['sentence_origid']
				tassimil=0
				oassimil=0
				nassimil=0
				eassimil=0
				mpassimil=0
				nn=0
				word_n=0
				word_t=0
				word_e=0
				slength=0
				
				alignment="_"
				VSsentence = False
				sinter=False

			
			if title != sentence['text_title']:
				if title:
					with open(outfolder+title+".csv","w") as out:
						out.write(fileoutput)
				fileoutput="\t".join(HEADERS)
				title=sentence['text_title']
				print(title)
			
			
			headid=int(word["dephead"])
			wordid=int(word["id"])
		

