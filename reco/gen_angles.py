import glob
import sys

argv = sys.argv[1:]
if len(argv)==0:
    print 'Usage: python gen_angles.py reco_file_dir'
    sys.exit(1)
reco_path = argv.pop(0)
#files = glob.glob('reco_files/signal_selection/*.root')
files = glob.glob(reco_path)
fout = open('run_makeAngles.sh','w')
for ifile in files:
    towrite = 'python makeAngles.py --inputfiles '+ifile+' --evtsperjob -1\n'
    fout.write(towrite)
fout.close()
