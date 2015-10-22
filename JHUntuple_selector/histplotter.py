# Plot all histograms in input file and dump to my webpage

from utility import *

from optparse import OptionParser

# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputfiles',
                  help='Input files')

(options, args) = parser.parse_args()

argv = []

# Get files
if options.inputfiles != '':
    files = glob.glob(inputfiles)

for ifile in files:
    print 'Plotting for file',ifile
    # Get type of the file
    tmp = ifile.split('/')
    outputname = tmp[len(tmp)-1].split('.root')[0]
    # Loop over all histograms in the file
    tmpfile = ROOT.TFile(ifile)
    keys = tmpfile.GetListOfKeys()
    hlist = []
    for ikey in keys :
        if 'TH1' in ikey.GetClassName():
            hlist.append(ikey.ReadObj())
    # Plot and upload to web
    plotting(hlist,outputname,'dump','notlog',None,'','recreate')
    # Closure
    tmpfile.Close()
    