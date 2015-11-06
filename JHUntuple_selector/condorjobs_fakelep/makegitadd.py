prepend = '/uscms_data/d3/lfeng7/Legacy_Afb/analysis/CMSSW_5_3_24/work/selection/JHUntuple_selector/output_rootfiles/v2_cutflow_modified/'
fout = open('gitaddall.sh','w')
import glob
alldir = glob.glob('*')
for idir in alldir:
    idir = idir.strip()
    # ignore files in current dir, only loop over dirs
    if len(idir.split('.'))>1 : continue
    towrite = 'git add '+idir+'/*.py '+idir+'/ana.* '+idir+'/*.*sh \n'
    fout.write(towrite)
fout.close()
