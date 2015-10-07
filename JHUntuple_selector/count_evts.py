# This small macro will read in all edm files in a directory and count the total number of events 

from fwlite_boilerplate import *

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

parser.add_option('--maxfiles', metavar='F', type='int', action='store',
                  default = -1,
                  dest='maxfiles',
                  help='max number of input ntuple files')

parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='If you want more information than you usually need.')
(options, args) = parser.parse_args()

argv = []

#debug
print 'Getting input files from '+ options.inputFiles

# Get the inputfiles.
if options.inputFiles:
    allfiles = glob.glob( options.inputFiles )
    # Only keep certain number of input files for fexibility
    if options.maxfiles > 0 :
        files = [allfiles[i] for i in range(options.maxFiles)]
    else : files = allfiles  

if options.verbose == 'yes' :   
    for ifile in files : print 'Getting these files',ifle

print 'Getting',len(files),'files'

# Read input files
events = Events(files)

print 'Getting',events.size(),'events'

# Wrtie output into a txt file
fout = open('evts_count.txt','a')
outputlines = 'Getting input files from '+ options.inputFiles+'\n'
outputlines += 'Getting '+str(len(files))+' files'+'\n'
outputlines += 'Getting '+str(events.size())+' events'+'\n\n'
fout.write(outputlines)
fout.close()

  
