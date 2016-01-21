import glob
allfiles= glob.glob('angles_files/signal/*.root')
fout = open('run_tmps.sh','w')
towrite = '' 
for ifile in allfiles:
    fname = ifile.split('_')
    postfix='.root'
    towrite += 'python makeTemp.py --inputfiles '+ifile+' --evtsperjob -1\n'
fout.write(towrite)
fout.close()
