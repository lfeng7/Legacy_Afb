fout = open('submitall.sh','w')
import glob
alldir = glob.glob('*')
for idir in alldir:
    idir = idir.strip()
    # ignore files in current dir, only loop over dirs
    if len(idir.split('.'))>1 : continue
    towrite = 'cd '+idir+'\n'
    towrite += 'rm -r '+idir+'\nsource grid_sub.csh \ncd ..\n\n'
    fout.write(towrite)
fout.close()
