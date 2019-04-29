#!/usr/bin/python

import json
#from Levenshtein import distance as lev
import sys
import logging
from Bio.pairwise2 import format_alignment, align
import functools

from collections import defaultdict

NORMALFLAG=True

#Utiliser functools.partial pour avoir accès au mot d'avant et au mot d'après 
#prototype futur = matchfunction (mot, motprecedent,motsuivant,r1,r2) => partial =>matchfunc

def load_skn(sknpath):
	with open(sknpath) as sknfile:
		return json.loads(sknfile.read())

def prettify(SKN,n):
	word=SKN['kwic'][n]["tokens"][0]
	sentence=SKN['kwic'][n]["structs"]

	"""
	
    ADJ: adjective
    ADP: adposition
    ADV: adverb
    AUX: auxiliary
    CCONJ: coordinating conjunction
    DET: determiner
    INTJ: interjection
    NOUN: noun
    NUM: numeral
    PART: particle
    PRON: pronoun
    PROPN: proper noun
    PUNCT: punctuation
    SCONJ: subordinating conjunction
    SYM: symbol
    VERB: verb
    X: other

	"""
	
	uposdict={"A":"ADJ", "Adp":"ADP", "Adv":"ADV", "Foreign":"X", "Interj":"INTJ", "N":"NOUN", "Num":"NUM", "Pron":"PRON", "Punct":"PUNCT", "Symb":"SYMB", "V":"VERB", "CSUBCAT_CS":"SCONJ", "CSUBCAT_CC":"CCONJ" }
	
	ID=word["id"]
	FORM= word["normalized"] if NORMALFLAG else word["original"]
	LEMMA = word["lemma"]
	UPOS = uposdict[ word["pos"]+word["msd"].split("|")[0] ] if word["pos"] == "C" else uposdict[ word["pos"]]
	XPOS = word["pos"]
	FEATS = word["msd"]
	HEAD =word["dephead"] 
	DEPREL = word["deprel"]
	DEPS = HEAD+':'+DEPREL
	MISC = "_"+sentence['sentence_origid']";"+word["word"]
	
	
	
	
	return "\t".join(map(str,[ID,FORM,LEMMA, UPOS, XPOS, FEATS, HEAD, DEPREL, DEPS, MISC]))

if __name__=="__main__":
	
	b=0
	beginning=0
	end=860000
	tokennumber=0
	tokennumfile=0
	#ALPHABET=set()
	sentence_id=None
	fileoutput=""
	outfolder="conll/"
	title=None
	
	for rank in range(beginning,end,10000):
		if rank == 0:
			continue
		SKN=load_skn("skn_corpus_%s-%s.json" % (b,rank))
		tokennumfile=0
		b=rank
		for x,y in enumerate(SKN["kwic"]):
			word=y["tokens"][0]
			sentence=y["structs"]
			
			if title != sentence['text_title']:
				if title:
					with open(outfolder+title+".conll","w") as out:
						out.write(fileoutput)
				fileoutput=""
				title=sentence['text_title']
				print(title)
			
			if sentence_id != sentence['sentence_origid']:
				if fileoutput:
					fileoutput+="\n"
				#print(sentence)
				sentence_id=sentence['sentence_origid']
			
			fileoutput+=prettify(SKN,x)+"\n"
			
		
			

				
	
