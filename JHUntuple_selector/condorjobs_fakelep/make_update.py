fout = open('update_all.sh','w')
import glob
alldir = glob.glob('*')
towrite = ''
for idir in alldir:
    idir = idir.strip()
    # ignore files in current dir, only loop over dirs
    if len(idir.split('.'))>1 : continue
    #towrite += 'rm -r '+idir+'/inputfiles\n'
    #towrite += 'cp -r template/inputfiles '+idir+'\n'
    towrite += 'rm -r '+idir+'/grid_sub.csh\n'
    towrite += 'cp  template/grid_sub.csh '+idir+'\n'
    fout.write(towrite)
fout.close()
