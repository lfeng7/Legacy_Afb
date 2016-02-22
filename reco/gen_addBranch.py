import glob
from optparse import OptionParser

parser = OptionParser()

parser.add_option('--dir1', metavar='F', type='string', action='store',
                  default='',
                  dest='dir1',
                  help='the file to be cloned')
parser.add_option('--dir2', metavar='F', type='string', action='store',
                  default='',
                  dest='dir2',
                  help='the file with a new branch to add')

(options, args) = parser.parse_args()

dir1 = options.dir1
dir2 = options.dir2

# file1 should be angles files
allfiles1= glob.glob(dir1+'/*.root')
# file2 should be reco files
allfiles2= glob.glob(dir2+'/*.root')

# find pairs of file 1 and 2
pairs = [(f1,f2) for f1 in allfiles1 for f2 in allfiles2 if f1.split('/')[-1].split('reco')[0] == f2.split('/')[-1].split('reco')[0] ]
# write sh command
fout = open('run_addBranch.sh','w')
towrite = '' 
for item in pairs:
    towrite += 'python add_branch.py --file1 '+item[0]+' --file2 '+item[1]+'\n'
fout.write(towrite)
fout.close()
