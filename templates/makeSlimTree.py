import ROOT
import glob
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
(options, args) = parser.parse_args()

argv = []

postfix = '_slim.root'

fin = ROOT.TFile(options.inputFiles)
oldtree =fin.Get('angles')
sample_name = options.inputFiles.split('/')
sample_name = sample_name[len(sample_name)-1].split('.root')[0]

fout = ROOT.TFile(sample_name+postfix,'recreate')
fout.SetCompressionLevel(9)
oldtree.SetBranchStatus('*w*',0)
oldtree.SetBranchStatus('*weight*',0)
newtree = oldtree.CloneTree()
fout.Write()
fout.Close()
fin.Close()
