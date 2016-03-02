import glob
#files = glob.glob('reco_files/signal_selection/*.root')
files = glob.glob('reco_files/signal_pscale_removed/*.root')
fout = open('run_makeAngles.sh','w')
for ifile in files:
    towrite = 'python makeAngles.py --inputfiles '+ifile+' --evtsperjob -1\n'
    fout.write(towrite)
fout.close()
