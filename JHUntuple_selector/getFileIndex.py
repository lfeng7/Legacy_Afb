# This small macro will read in all edm files in a directory and count the total number of events 

from utility import *

from optparse import OptionParser

# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='If you want more information than you usually need.')

parser.add_option('--outname', metavar='F', type='string', action='store',
                  default = 'incomplete.txt',
                  dest='outname',
                  help='')

(options, args) = parser.parse_args()

argv = []

#debug
print 'Getting input files from '+ options.inputFiles

fout = open(options.outname,'w')

# Get the inputfiles.
if options.inputFiles:
    allfiles = glob.glob( options.inputFiles ) 

if options.verbose == 'yes' :   
    for ifile in files : print 'Getting these files',ifle

print 'Getting',len(files),'files'

all_index = [int(ifile.split('output')[1].split('.root')[0].strip('_')) for ifile in allfiles]
for i in all_index : fout.write(str(i)+',')    


