#!/usr/bin/python

import json
#from Levenshtein import distance as lev
import sys
import logging
from Bio.pairwise2 import format_alignment, align
from functools import partial

from collections import defaultdict

from sknCONLL import load_skn

import pandas as pd

alphabet={'ᴱ', 'h', '’', 'ᴉ', 'ϑ', 'Ä', 'c', 'ᵏ', 'J', '˘', 'ˉ', 'E', 'B', 'a', 'ᵉ', '₍', '0', 'Ǹ', 'ŧ', '”', '“', '1', '»', 'ᵃ', 'ĕ', 'u', 'i', '%', 'q', '̬', 'd', 'ò', 'ê', 'P', '‿', 'R', 'ᶠ', 'ʸ', '"', 'ᴋ', 'Z', '–', 'ᵑ', 'Y', 'ʷ', '?', 'M', ' ', 'û', '͔', 'ẅ', 'ʜ', '̮', ']', 'ʲ', 'ŏ', 'Ŏ', '^', 'ᴹ', 'Ì', '\uf21d', ')', 'á', '̳', 'ẁ', '\u0ff1', 'w', '+', 'ᴍ', 'Ĺ', 'å', 'ᴡ', 'χ', 'A', ',', '.', 'é', '*', '§', 'ŕ', 'N', '_', 'I', 'ŋ', 'b', 'ᵐ', 'ỳ', 'e', 'ŷ', 'k', 'ɦ', 'p', 'Ŭ', 'ə', 'ʟ', '[', 'ɢ', 'V', 'H', 'r', '°', 'W', 'ä', '̇', 'ḙ', 'ᴰ', 'ɪ', '̰', 'ʏ', 'ð', '\u1cff', 'ᴅ', 'ᴇ', '̹', 'ᵗ', 'ȯ', 'ì', '̲', 'ʀ', 'ᴺ', 'C', 's', 'U', 'ḭ', 'ᴵ', 'ᴠ', 'ɔ', 'ʾ', 'ʼ', 'ᴬ', 'y', 'ᴜ', '͇', 'Ń', 'ᵎ', 'ʙ', 'ᵻ', 'ö', '!', '5', '́', 'è', 'f', 'î', 'И', 't', '̀', 'Ö', 'ù', 'ị', '̄', '=', 'ṋ', '∼', 'β', 'δ', '-', 'à', 'D', 'ᴼ', 'Ȯ', 'Î', 'ʳ', 'v', 'ĺ', '(', 'ⁱ', 'ᴛ', 'ė', '~', '̥', 'ᵒ', '̭', 'ᵛ', 'ᵊ', 'ᵁ', '̆', 'È', 'ⁿ', '̜', 'ˈ', 'ú', 'ʰ', 'Å', '̤', 'ᴊ', '\xa0', '̈', 'ʿ', 'ś', 'n', 'ń', 'ḱ', 'ɐ', 'Ò', 'F', 'ǹ', 'ˢ', 'í', 'š', 'ǰ', '̗', 'ᵘ', 'O', 'ŭ', 'ᴮ', 'T', 'ᴏ', 'o', 'ĭ', 'ḛ', 'ˡ', 'ó', 'ᴶ', 'À', '\u0fff', 'ᴎ', 'l', 'z', "'", 'g', 'ₒ', '#', 'j', 'ṵ', 'ô', '<', 'm', 'L', '´', 'K', '`', 'ɴ', 'G', '̂', '͕', 'ɛ', 'ˀ', 'S', 'ᵖ', 'â', '3', 'ă', 'ᴀ'}


consonnes= "zδtpqdfghjklmnbvcxðsrvŋw"
voyelles="aeiouyäöəɛᴏ"
punctuation=" &~#\"'({[-|`_\^@)]°=+}$^¨%µ*!:;,?./§1234567890"

same=[("à","a")]

NORMALFLAG=True
DEBUGSAMPLE=1
#Utiliser functools.partial pour avoir accès au mot d'avant et au mot d'après 
#prototype futur = matchfunction (mot, motprecedent,motsuivant,r1,r2) => partial =>matchfunc

	

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
	if x > 2:
		if x < len(normalized) and normalized[x] == normalized[x-1]:
			if y > 1:
				return -7-(0.7*(y-1))
			else:
				return -1
	
	return (-7-(0.7*y))

def alignage(w,n):
	Word=list(w)
	Normalized=list(n)
	gapfunctionA=partial(gap_functionA,Word,Normalized)
	gapfunctionB=partial(gap_functionB,Word,Normalized)
					
	phonalignement=align.globalcc(Word,Normalized,matchfunc,gapfunctionA,gapfunctionB,gap_char=['-'],one_alignment_only=True)
	for elem in phonalignement:
		align1, align2, score, begin, end=elem
		#print("".join(align1)+"\n"+"".join(align2)+"\n"+str(score))
		
	return phonalignement

def traitsphonos(word, tabl,tableauindex):
	assimilation=False
	
	if word["word"][-1]=="n":
		tabl.loc[tableauindex,"word_n"]+=1
		if (word["normalized"][-1]==nextword["normalized"][0] ) and word["normalized"][-1] != "n":
			tabl.loc[tableauindex,"n_ass"]+=1
			assimilation=True
		elif word["normalized"][-1]=="m" and nextword["normalized"][0] == "p":
			tabl.loc[tableauindex,"mp_ass"]+=1
			assimilation=True
		elif nextword["normalized"][-1] == "n":
			tabl.loc[tableauindex,"n_n"] += 1
			assimilation=True
		
		
	elif word["word"][-1]=="t":
		tabl.loc[tableauindex,"word_t"]+=1
		if  (word["normalized"][-1]==nextword["normalized"][0] ) and word["normalized"][-1] != "t":
			tabl.loc[tableauindex,"t_ass"]+=1
			assimilation=True
		
	elif word["word"][-1]=="e":
		tabl.loc[tableauindex,"word_e"]+=1
		if  (word["normalized"][-1]==nextword["normalized"][0] ) and word["normalized"][-1] != "e":
			tabl.loc[tableauindex,"e_ass"]+=1
			assimilation=True
	
	elif (word["normalized"][-1]==nextword["normalized"][0] ) and word["word"][-1]!=word["normalized"][-1]:
		tabl.loc[tableauindex,"otherass"]+=1
		assimilation=True
	
	return assimilation
	

def synstruc(syntacdic,tableau,tableauindex):
	tableau.loc[tableauindex, "props"] += 1
	
	if not (DEBUGSAMPLE % 101):
		print(' '.join([ syntacdic[x]["normalized"] for x in syntacdic]))
		print(' '.join([ syntacdic[x]["word"] for x in syntacdic]))
	
	proporders=defaultdict(dict)
	
	for elem in syntacdic:
		for e in syntacdic:
			if syntacdic[e]["dephead"] == elem:
				msd=dict([tuple(x.split("_",1)) for x in syntacdic[e]["msd"].split("|")])
				
				if not (DEBUGSAMPLE % 101):
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
	
	HEADERS=["id","longueur","t_ass","n_ass","mp_ass","n_n","e_ass","otherass","word_t","word_n","word_e","interrog","props","vs","sv","os","so","vo","ov","svo","sov","ovs","osv","vso","vos","sv"]
	NONNUMHEADERS=["id"]
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
				 
				 props : nombre de propositions
				 svo, osv, etc : nombre de propositions utilisant l'ordre svo, osv, etc
				 
				 """
	tableau=pd.DataFrame()
	
	syntacdict=dict()
	
	for rank in range(beginning,end,10000):
		if rank == 0:
			continue
		SKN=load_skn("skn_corpus_%s-%s.json" % (b,rank))
		b=rank
		for x,y in enumerate(SKN["kwic"]):
			word=y["tokens"][0]
			
			if False and word["word"] and word["word"] != word["normalized"]:
				alignage(word["word"],word["normalized"])
			
			sentence=y["structs"]
			
			if sentence_id != sentence['sentence_origid']:
				#DEBUGSAMPLE+=1
				
				if sentence_id:
					synstruc(syntacdict,tableau,tableauindex)
					#print(tableau.iloc[-1:])
					tableauindex+=1
					sentence_id= sentence['sentence_origid']
					for elem in HEADERS:
						tableau.loc[tableauindex,elem] = sentence_id if elem == "id" else 0
					syntacdict=dict()
			
				else:
					sentence_id= sentence['sentence_origid']	
				
			if title != sentence['text_title']:
				if title:
					
					#print(tableau)
					#input()
					#with open(,"w") as out:
					try:
						tableau.to_csv(outfolder+title+"-v2.csv", sep='\t', encoding='utf-8')
					except Exception as e:
						print(e)
						pass
				
				print(title)
				title=sentence['text_title']
				
				
				tableau=pd.DataFrame(dict( ( (x,[sentence_id]) if x == "id" else (x,0) for x in HEADERS  )) )
				tableauindex=0
				sentence_id = None
				
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
			
			
			

