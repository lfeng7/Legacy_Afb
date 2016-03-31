#! /usr/bin/python
import glob
import sys

argv = sys.argv[1:]
if len(argv)<1:
    print    """
    Generate a shell code that can be used to run makeTemp.py to get template root file from angles.root
    Usage: ./gen_tmps.py angles_files/signal/*.root angles_files/signal/MC.txt sideband
    """
    sys.exit(1)

file_path = argv.pop(0)
MC_txt = argv.pop(0)
if 'sideband' in argv:
    sideband = 'yes'
else:
    sideband = 'no'

allfiles= glob.glob(file_path)
fout = open('run_tmps.sh','w')
towrite = '' 
for ifile in allfiles:
    fname = ifile.split('_')
    postfix='.root'
    towrite += 'python makeTemp.py --inputfiles %s --evtsperjob -1 --MCinfo %s --sideband %s\n'%(ifile,MC_txt,sideband)
fout.write(towrite)
fout.close()
