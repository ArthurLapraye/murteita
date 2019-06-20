#!/usr/bin/python

import json
#from Levenshtein import distance as lev
import sys
import logging
import re
from Bio.pairwise2 import format_alignment, align
from functools import partial
from random import randint
from collections import defaultdict

from sknCONLL import load_skn

import pandas as pd

alphabet={'ᴱ', 'h', '’', 'ᴉ', 'ϑ', 'Ä', 'c', 'ᵏ', 'J', '˘', 'ˉ', 'E', 'B', 'a', 'ᵉ', '₍', '0', 'Ǹ', 'ŧ', '”', '“', '1', '»', 'ᵃ', 'ĕ', 'u', 'i', '%', 'q', '̬', 'd', 'ò', 'ê', 'P', '‿', 'R', 'ᶠ', 'ʸ', '"', 'ᴋ', 'Z', '–', 'ᵑ', 'Y', 'ʷ', '?', 'M', ' ', 'û', '͔', 'ẅ', 'ʜ', '̮', ']', 'ʲ', 'ŏ', 'Ŏ', '^', 'ᴹ', 'Ì', '\uf21d', ')', 'á', '̳', 'ẁ', '\u0ff1', 'w', '+', 'ᴍ', 'Ĺ', 'å', 'ᴡ', 'χ', 'A', ',', '.', 'é', '*', '§', 'ŕ', 'N', '_', 'I', 'ŋ', 'b', 'ᵐ', 'ỳ', 'e', 'ŷ', 'k', 'ɦ', 'p', 'Ŭ', 'ə', 'ʟ', '[', 'ɢ', 'V', 'H', 'r', '°', 'W', 'ä', '̇', 'ḙ', 'ᴰ', 'ɪ', '̰', 'ʏ', 'ð', '\u1cff', 'ᴅ', 'ᴇ', '̹', 'ᵗ', 'ȯ', 'ì', '̲', 'ʀ', 'ᴺ', 'C', 's', 'U', 'ḭ', 'ᴵ', 'ᴠ', 'ɔ', 'ʾ', 'ʼ', 'ᴬ', 'y', 'ᴜ', '͇', 'Ń', 'ᵎ', 'ʙ', 'ᵻ', 'ö', '!', '5', '́', 'è', 'f', 'î', 'И', 't', '̀', 'Ö', 'ù', 'ị', '̄', '=', 'ṋ', '∼', 'β', 'δ', '-', 'à', 'D', 'ᴼ', 'Ȯ', 'Î', 'ʳ', 'v', 'ĺ', '(', 'ⁱ', 'ᴛ', 'ė', '~', '̥', 'ᵒ', '̭', 'ᵛ', 'ᵊ', 'ᵁ', '̆', 'È', 'ⁿ', '̜', 'ˈ', 'ú', 'ʰ', 'Å', '̤', 'ᴊ', '\xa0', '̈', 'ʿ', 'ś', 'n', 'ń', 'ḱ', 'ɐ', 'Ò', 'F', 'ǹ', 'ˢ', 'í', 'š', 'ǰ', '̗', 'ᵘ', 'O', 'ŭ', 'ᴮ', 'T', 'ᴏ', 'o', 'ĭ', 'ḛ', 'ˡ', 'ó', 'ᴶ', 'À', '\u0fff', 'ᴎ', 'l', 'z', "'", 'g', 'ₒ', '#', 'j', 'ṵ', 'ô', '<', 'm', 'L', '´', 'K', '`', 'ɴ', 'G', '̂', '͕', 'ɛ', 'ˀ', 'S', 'ᵖ', 'â', '3', 'ă', 'ᴀ'}


consonnes= "zδtpqdfghjklmnbvcxðsrvŋw"
voyelles="aeiouyäöəɛᴏ"
punctuation=" &~#\"'({[-|`_\^@)]°=+}$^¨%µ*!:;,?./§1234567890"

same=[("à","a")]

NORMALFLAG=True
DEBUGOUTPUT=False
DEBUGSAMPLE=1

def makemsdic(word):
	morphofeats=word["msd"]
	
	return  dict([tuple(x.split("_",1)) for x in morphofeats.split("|")])
	
pairesproches=set((("t","d"),("d","r"),("ä","e"),("ö","ä"),("ö","y"),("o","a"),('ö','e'),("o","u"),("i","j"),("t","s"),("ð","d"),("i","e"),("ə","i"),("ə","e"),("d","r"),("d","w"),("d","j"),('y','w'),("ö","w"),("v","f"),("t","m"),("m","p"), ("u","w")))

re1=re.compile("ᴏ")
re2=re.compile("δ")

def betternorm(s):
	try:
		z=re2.sub("ð",re1.sub("o",s))
	except TypeError as e:
		logging.error(e,"at",s)
		return None
		
	return z
	

def matchfunc(stand,dial):
	r1=stand.lower()
	r2=dial.lower()
	
	spaces=(" ","-")
	
	if r1==r2:
		return 20
	else:
		for p in pairesproches:
			if r1 in p and r2 in p:
				return 7
		else:
			if r1 in "mn" and r2 in "mn":
				return 2
			elif r1 == "n" and r2 in consonnes:
				return 1.5
			
			elif r1 in voyelles and r2 in voyelles:		
				return 1
			else:
				return -20
		
	return -100

#matchfunc=functools.partial(matchfunction,"","","")

def gap_functionA(word, normalized,x,y):
	if 3 < x < len(normalized)  and normalized[x-2] == normalized[x] and normalized[x-1] in ["h","l"] and normalized[x] in voyelles:
		return -10-(y-1) if y > 1 else -1
	
	if x < len(normalized) and normalized[x-1]=="i" and normalized[x]=="j":
		if y > 1:
			return -2-(y-1)
		else:
			return -1
	
	if x > 2:
		if x < len(normalized) and normalized[x] == normalized[x-1]:
			if y > 1:
				return -2-(y-1)
			else:
				return -1
	
	return (-10 + -1*y)

def gap_functionB(word, normalized,x,y):
	if x > 2:
		if x < len(normalized) and normalized[x] == normalized[x-1]:
			if y > 1:
				return -1.5-(0.15*(y-1))
			else:
				return -1
	
	return (-1.5-(0.15*y))


def alignage(w,n):
	#print(w,n)
	Word=list(w.lower())
	Normalized=list(n.lower())
	gapfunctionA=partial(gap_functionA,Word,Normalized)
	gapfunctionB=partial(gap_functionB,Word,Normalized)
		
	phonalignement=align.globalcc(Word,Normalized,matchfunc,gapfunctionA,gapfunctionB,gap_char=['-'])
	
	
	if randint(0,1001) > 999:
		for elem in phonalignement:
			align1, align2, score, begin, end=elem
			print(str("".join(align1)+"\n"+"".join(align2)+"\t"+str(score)))
		
	return phonalignement

def traitsphonos(word, tabl,tableauindex):
	assimilation=False
	
	nextcons=nextword["normalized"][0]
	lastcons=word["word"][-1]
	asscons=word["normalized"][-1]
	
	if lastcons=="n":
		tabl.loc[tableauindex,"word_n"]+=1
		if (asscons==nextcons ) and asscons != "n":
			tabl.loc[tableauindex,"n_ass"]+=1
			assimilation=True
		elif asscons=="m" and nextcons == "p":
			tabl.loc[tableauindex,"mp_ass"]+=1
			assimilation=True
		elif nextcons == "n":
			tabl.loc[tableauindex,"n_n"] += 1
			assimilation=True
		
		
	elif lastcons=="t":
		tabl.loc[tableauindex,"word_t"]+=1
		if  (asscons==nextcons ) and asscons != "t":
			tabl.loc[tableauindex,"t_ass"]+=1
			assimilation=True
		
	elif lastcons=="e":
		tabl.loc[tableauindex,"word_glott"]+=1
		if  (asscons==nextcons ) and asscons != "e":
			tabl.loc[tableauindex,"glott_ass"]+=1
			assimilation=True
	
	elif word["pos"]== "V":
		vmsd=makemsdic(word)
		
		if vmsd.get("INF",None):
			if vmsd["INF"]=="Inf1":
				tabl.loc[tableauindex,"word_glott"]+=1
				if asscons == nextcons:
					tabl.loc[tableauindex, "glott_ass"] += 1
					assimilation=True
		
		if vmsd.get("MOOD",None):
			if vmsd["MOOD"]=="Imprt":
				tabl.loc[tableauindex,"word_glott"]+=1
				if	asscons == nextcons:
					tabl.loc[tableauindex, "glott_ass"] += 1
					assimilation=True
	
	elif (asscons==nextcons ) and lastcons !=asscons:
		tabl.loc[tableauindex,"otherass"]+=1
		assimilation=True
	
	return assimilation
	

def synstruc(syntacdic,tableau,tableauindex):
	tableau.loc[tableauindex, "props"] += 1
	
	if False and not (DEBUGSAMPLE % 101):
		print(' '.join([ syntacdic[x]["normalized"] for x in syntacdic]))
		print(' '.join([ syntacdic[x]["word"] for x in syntacdic]))
	
	proporders=defaultdict(dict)
	
	for elem in syntacdic:
		for e in syntacdic:
			if syntacdic[e]["dephead"] == elem:
				msd=makemsdic(syntacdic[e])
				
				if False and not (DEBUGSAMPLE % 101):
					print("\t",syntacdic[e]["normalized"],e,syntacdic[e]["deprel"],"child of",elem,syntacdic[elem]["normalized"])
				
				if "subj" in syntacdic[e]["deprel"]:
					proporders[elem]["subj"]=e
				
				elif "obj" in syntacdic[e]["deprel"]:
					proporders[elem]["obj"]=e
				
				if msd.get("SUBCAT",None) == "Interr":
					tableau.loc[tableauindex, "interrog"] += 1
				
				if msd.get("CLIT",None) and msd.get("CLIT",None).startswith("Qst"):
					tableau.loc[tableauindex, "interrog"] += 1
					
				if word["pos"] == "C":
					tableau.loc[tableauindex, "props"] += 1
					
	for elem in proporders:
		aligndict={"v":elem,
		"s":proporders[elem].get("subj",None),
		"o":proporders[elem].get("obj",None) }
		
		alignementstring=""
		
		for e in sorted((x for x in aligndict if aligndict[x]),key=lambda x : aligndict[x]):
			alignementstring+=e
		
		if not (DEBUGSAMPLE % 101) :
			print(alignementstring)
		
		tableau.loc[tableauindex, alignementstring] += 1
				
				
				
	if not (DEBUGSAMPLE % 101) :
		if False:
			input()
	
	return 0


if __name__=="__main__":
	#SKN=load_skn(sys.argv[1])
	
	
	
	b=0
	beginning=0
	end=860000
	passit=True
	tokennumber=0
	#ALPHABET=set()
	sentence_id=None
	fileoutput=""
	outfolder="sfeatures/"
	title=None
	assdict=defaultdict(str)
	
	HEADERS=["id","longueur","longueurpho","t_ass","n_ass","mp_ass","n_n","glott_ass","otherass","word_t","word_n","word_glott","interrog","props","vs","sv","os","so","vo","ov","svo","sov","ovs","osv","vso","vos","sv"]+list(map(lambda x : "var_"+ x, ['ö+w', 'y+o', 'n+w', 'o+ə', 'a+ə', 'y+a', 'ä+o', 'y+e', 'ö+e', 'u+y', 'e+ᴏ', 'u+e', 'e+y', 'i+ö', 'm+n', 'u+i', 'o+ö', 'o+e', 'a+ᴏ', 'y+i', 'ä+y', 'ö+y', 'ä+ə', 'u+ö', 'ö+o', 'ä+ö', 'y+u', 'i+ə', 'i+y', 'n+g', 'e+u', 'o+ᴏ', 'i+u', 'd+δ', 'd+w', 'o+i', 'n+p', 'y+ö', 'i+j', 'y+ä', 'a+ä', 'd+t', 'e+ə', 'i+o', 't+d', 'ö+ä', 'u+a', 'e+a', 'j+i', 'n+k', 'i+ä', 'o+u', 'a+e', 'n+s', 'n+r', 'd+j', 'a+u', 'n+t']))

	NONNUMHEADERS=["id"]
	headershelp="""id : id de la phrase
				  longueur : longueur en mot de la phrase
				  longueurpho : nombre de phonèmes de la phrase
				  t_ass nombre de t en fin de mot assimilés totalement
				  n_ass nombre de n en fin de mot assimilés totalement
				  mp_ass nombre de n assmilés à m devant p 
				  n_n nombre de n en fin de mot suivi de n initial
				 glott_ass nombre de mot finissant en occlusive glottale dans le standard et ayant leur glottale assimilée
				 o_ass autres assimilations de sandhi
				 word_t mots finissants en t
				 word_n
				 word_glott
				 
				 props : nombre de propositions
				 svo, osv, etc : nombre de propositions utilisant l'ordre svo, osv, etc
				 
				 """
	tableau=pd.DataFrame()
	
	syntacdict=dict()
	phonocorr=defaultdict(float)
	
	for rank in range(beginning,end,10000):
		if rank == 0:
			continue
		SKN=load_skn("skn_corpus_%s-%s.json" % (b,rank))
		b=rank
		for x,y in enumerate(SKN["kwic"]):
			sentence=y["structs"] 
			
			if title != sentence['text_title']:
				if title:
					
					#print(tableau)
					#input()
					#with open(,"w") as out:
					
					try:
						tableau.to_csv(outfolder+title+"-v2.csv", sep='\t', encoding='utf-8',index=False)
					except Exception as e:
						print(e)
						pass
				
				print(title)
				title=sentence['text_title']
				
				
				tableau=pd.DataFrame(dict( ( (x,[sentence_id]) if x == "id" else (x,0) for x in HEADERS  )) )
				tableauindex=0
				
			
			if sentence_id != sentence['sentence_origid']:
				if DEBUGOUTPUT:
					DEBUGSAMPLE+=1
				
				if sentence_id:
					tableau.loc[tableauindex, "longueurpho"] = len("".join([ syntacdict[x]["normalized"] for x in syntacdict]))
					synstruc(syntacdict,tableau,tableauindex)
					#print(tableau.iloc[-1:])
					tableauindex+=1
					sentence_id= sentence['sentence_origid']
					for elem in HEADERS:
						tableau.loc[tableauindex,elem] = sentence_id if elem == "id" else 0
					syntacdict=dict()
			
				else:
					sentence_id= sentence['sentence_origid']	
				
			
			word=y["tokens"][0]
			word["word"]=betternorm(word["word"])
			word["normalized"]=betternorm(word["normalized"])
			
			if word["word"] and word["word"] != word["normalized"]:
				alignement=alignage(word["word"],word["normalized"])
				if alignement:
					
					align1, align2,align3,align4,align5 = alignement[0]
					for orig,dialect in zip(align1, align2):
						if orig != dialect:
							featname="var_"+orig+"+"+dialect
							phonocorr[featname]+=1
							if featname in HEADERS:
								tableau.loc[tableauindex, featname] += 1
							
				else:
					logging.error("Problème d'alignement pour la paire "+word["word"]+" "+word["normalized"]) 
			
				
				#print(tableauindex,sentence_id)
			
						
			#print("caca")
			tableau.loc[tableauindex,"longueur"] += 1
			
			#Bloc qui recherche le mot suivant pour détecter les phénomènes de sandhi
			#Si le mot est en fin de phrase, la phrase suivante est chargée, voire le fichier suivant
			try:
				nextword=SKN['kwic'][x+1]["tokens"][0]
			except IndexError as I:
				try:
					SKN2=load_skn("skn_corpus_%s-%s.json" % (b,rank+10000))
					nextword=SKN2['kwic'][0]["tokens"][0]
				except FileNotFoundError as e:
					passit=False
			
			if passit and word["word"]:
				assimilation=traitsphonos(word,tableau,tableauindex)
				if assimilation:
					assdict[sentence_id]+=word["normalized"]+"_"+nextword["normalized"]
				
			elif not word["word"]:
				logging.error("Problem at "+title+" "+str(word))
			
			syntacdict[word["id"]]=word
			
	
	for elem in sorted(phonocorr,key=lambda x : phonocorr[x],reverse=False):
		print(elem,":",phonocorr[elem])
	
	
	sys.exit(0)		

