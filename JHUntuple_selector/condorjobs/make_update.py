fout = open('update_all.sh','w')
import glob
alldir = glob.glob('*')
for idir in alldir:
    idir = idir.strip()
    # ignore files in current dir, only loop over dirs
    if len(idir.split('.'))>1 : continue
#    towrite = 'cp -r template/inputfiles '+idir+'\n'
    towrite = 'cp  template/grid_sub.csh '+idir+'\n'
    fout.write(towrite)
fout.close()
