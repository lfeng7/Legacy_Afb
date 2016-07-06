import sys
argv = sys.argv[1:]
if len(argv)==0:
    print """ usage: python thisfile path/to/selected_files """
    sys.exit(1)

prepend = argv.pop(0) 
fout = open('cleanall.sh','w')
import glob
alldir = glob.glob('*')
for idir in alldir:
    idir = idir.strip()
    # ignore files in current dir, only loop over dirs
    if len(idir.split('.'))>1 : continue
    towrite = 'cd '+idir+'\nsource cleanup.sh\n'
    towrite += 'hadd ../%s/%s_selected.root %s/*.root \n'%(prepend,idir,idir)
    towrite += 'rm -r %s \ncd ..\n\n'%idir
    fout.write(towrite)
fout.close()
