#!/usr/bin/python

from urllib.request import urlopen
import json

beginning=0
end=860000



"""
0	"word"
1	"original"
2	"normalized"
3	"comment"
4	"id"
5	"ref"
6	"lemma"
7	"lemmacomp"
8	"pos"
9	"msd"
10	"dephead"
11	"deprel"
12	"nertag"
13	"lex"
14	"nerbio"


	
0	"text"
1	"text_name"
2	"text_title"
3	"text_editor"
4	"text_parish"
5	"text_dialect_group"
6	"text_dialect_region"
7	"text_date"
8	"text_datefrom"
9	"text_dateto"
10	"text_urlwav"
11	"text_urltextgrid"
12	"text_urleaf"
13	"text_timefrom"
14	"text_timeto"
15	"paragraph"
16	"paragraph_id"
17	"paragraph_speaker"
18	"paragraph_sex"
19	"paragraph_role"
20	"sentence"
21	"sentence_id"
22	"sentence_origid"
23	"sentence_beg"
24	"sentence_duration"
25	"sentence_urlview"
26	"ne"
27	"ne_name"
28	"ne_fulltype"
29	"ne_ex"
30	"ne_type"
31	"ne_subtype"
32	"ne_placename"
33	"ne_placename_source"
"""

b=beginning
for x in range(beginning,end,10000):
	SKN_URL = ("https://korp.csc.fi/cgi-bin/korp.cgi"+ 
		"?command=query&cqp=[]&corpus=skn"+
		"&show=word,original,normalized,comment,id,ref,lemma,lemmacomp,pos,msd,dephead,deprel,nertag,lex,nerbio&"+
		"show_struct=text,text_name,text_title,text_editor,text_parish,text_dialect_group,text_dialect_region,text_date,text_datefrom,text_dateto,text_urlwav,text_urltextgrid,text_urleaf,text_timefrom,"+
		"text_timeto,paragraph,paragraph_id,paragraph_speaker,paragraph_sex,paragraph_role,sentence,sentence_id,sentence_origid,sentence_beg,sentence_duration,sentence_urlview,ne,ne_name,ne_fulltype,"+
		"ne_ex,ne_type,ne_subtype,ne_placename,ne_placename_source&start=%s&end=%s&defaultcontext=0" % (b,x) )
	z=json.loads(urlopen(SKN_URL).read())
	
	with open("skn_corpus_%s-%s.json" % (b,x),"w") as out:
		out.write(json.dumps(z))
	
	print(x)	
		
	b=x


	#

#print(SKN_URL)

#"""text,text_name,text_title,text_editor,text_parish,text_dialect_group,text_dialect_region,text_date,text_datefrom,text_dateto,text_urlwav,text_urltextgrid,text_urleaf,text_timefrom,text_timeto,paragraph,paragraph_id,paragraph_speaker,paragraph_sex,paragraph_role,sentence,sentence_id,sentence_origid,sentence_beg,sentence_duration,sentence_urlview,ne,ne_name,ne_fulltype,ne_ex,ne_type,ne_subtype,ne_placename,ne_placename_source" """

#print(type(z))

#
