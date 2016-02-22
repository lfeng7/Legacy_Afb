# add a branch to an existing tree and save as a new root file
import os
import glob
import math
import ROOT
from ROOT import *
import sys
from optparse import OptionParser
from array import array

parser = OptionParser()

parser.add_option('--file1', metavar='F', type='string', action='store',
                  default='',
                  dest='the file to be cloned',
                  help='')
parser.add_option('--file2', metavar='F', type='string', action='store',
                  default='',
                  dest='the file with a new branch to add',
                  help='')

(options, args) = parser.parse_args()

f1 = options.file1
f2 = options.file2
branch = 'final_chi2'

file1 = ROOT.TFile(f1)
file2 = ROOT.TFile(f2)

# Find the name of the ttree
keys = f1.GetListOfKeys()
for ikey in keys:
    if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
print 'Getting ttree',treename

# create a new root file
fname = f1.split('.root')[0]
fname += '_new.root'
f_new = ROOT.TFile(fname,'recreate')

# clone tree 1
tree1 = f1.Get(treename)
tree2 = f2.Get(treename)

newtree = tree1.CloneTree(0)


# add new branch from tree2 to newtree
br = array('f',[0.])
newtree.Branch(branch+'_new',br,'br/F')

# Loop over events
nevts = tree1.GetEntries()
if nevts != tree2.GetEntries():
    print 'num entries of two files are not the same. Cannot add branch to file1.'
    sys.exit()
for iev in range(nevts):
    tree2.GetEntry(iev)
    br[0] = tree2.final_chi2[0]
    newtree.Fill()

# write and save tree
f_new.Write()
f_new.Close()




