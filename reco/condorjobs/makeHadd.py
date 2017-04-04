import sys 
from glob import glob

argv = sys.argv[1:]
if argv == []:
    print 'Usage: python thisfile.py input_files output_dir'
    sys.exit(1)
path = argv.pop(0)
output_dir = argv.pop(0)

input_dir = path.split('/')[-2]
print 'input_dir',input_dir
print 'output_dir',output_dir

all_files = glob(path)

unique_files = set()

for ifile in all_files:
    unique_files.add(ifile.split('reco')[0].split('/')[-1])

fout = open('hadd_all.sh','w')
to_write = []
for item in unique_files:
    tmp = 'hadd %s/%s_reco.root %s/%s*.root\n'%(output_dir,item,input_dir,item)
    to_write += [tmp]
to_write.sort()
to_write = ''.join(to_write)
fout.write(to_write)
fout.close() 
