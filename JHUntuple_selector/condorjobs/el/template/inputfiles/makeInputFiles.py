#! /usr/bin/python

import glob
import os

all_txt = glob.glob('./*.txt')

for ifile in all_txt:
	# /eos/uscms/store/user/eminizer/TT_CT10_TuneZ2star_8TeV-powheg-tauola/TTBar_Powheg_v1/151116_191806/0000/*.root
	with open(ifile, 'r') as tmp:
		# find the eos dir of files
	    lines = tmp.readlines()
	    lines = lines[0].strip()
	    lines = lines.split('/eos/uscms')[1]
	    lines = lines.split('/*.root')[0]
	    cmd = 'eosls %s >> %s'%(lines,ifile)
	    print cmd



