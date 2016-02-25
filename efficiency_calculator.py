# take many root files and make comparison plots using different colors
# if a MC.txt is provided, "stack" samples according to type

import os
import glob
import math
import ROOT
from ROOT import *
import sys
from optparse import OptionParser
from plot_tools import *

ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )

parser = OptionParser()

parser.add_option('--cut', metavar='F', type='string', action='store',
          default="",
                  dest='cut',
                  help='')
parser.add_option('--dir', metavar='F', type='string', action='store',
                  default='',
                  dest='dir',
                  help='')
parser.add_option('--output', metavar='F', type='string', action='store',
                  default='cut_results.txt',
                  dest='output',
                  help='output txt for cut results')
parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='')

(options, args) = parser.parse_args()


cut = options.cut
rundir = options.dir
verbose = options.verbose
output = options.output

# parse .txt to get the correct normalization 
data_lumi = 19700

# Open output txt file. If not exist, will create a new one
file_output = open(output,'a')

# Get input MC and data files according to txt file
txt_MC = open(rundir+'/MC_input_with_bkg.txt')
mc_samples = []
for line in txt_MC:
    items = line.split()
    #                TT_CT10_qq           qq              21560109              245.9   
    mc_samples += [{'file_name':items[1],'type':items[2],'nevts':int(items[3]),'xsec':float(items[4])}]

# Getting all file path
all_files = glob.glob(rundir+'/*.root')
for ifile in all_files:
    for item in mc_samples:
        if item['file_name'] in ifile:
            item['file_path'] = ifile

if verbose in ['yes','verbose']:
    print 'Plotting from these files!'
    for isample in mc_samples:
        print isample['file_path']   

# Get a list of types of samples 
# mc_samples is a list of dictionary like: {'file_name','type','nevts','xsec','file_path'}
all_types = [item['type'] for item in mc_samples]
all_types = list(set(all_types))

# Find tree name
tf = ROOT.TFile(mc_samples[0]['file_path'])
keys = tf.GetListOfKeys()
for ikey in keys:
    if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
print 'Getting ttree',treename+'\n'
tf.Close()

# add the weight of the tree to L*xsec/nevts
for isample in mc_samples:
    iweight = data_lumi*isample['xsec']/isample['nevts']
    isample['weight'] = iweight
    if verbose in ['yes','verbose']:
        print isample['file_name'],'weight=',iweight

# Loop over types to get nevts before and after cuts
cut_table = []
to_write = '\n----------------------------------------\n\n'
to_write += 'cut %s\n\n'%cut
for i,itype in enumerate(all_types):
    # Loop over samples, find samples of current type, and make a temp hist
    n_samples = 0
    nevts_no_cut = 0
    nevts_cut = 0
    for isample in mc_samples:
        if isample['type'] == itype :
            ifile = ROOT.TFile(isample['file_path'])
            ttree = ifile.Get(treename)
            iweight = isample['weight']
            nevts_no_cut += ttree.GetEntries()*iweight
            nevts_cut += ttree.GetEntries(cut)*iweight
    iefficiency = nevts_cut*1.0/nevts_no_cut
    ientry = {'type':itype,'nevts_no_cut':nevts_no_cut,'nevts_cut':nevts_cut,'efficiency':iefficiency}
    cut_table.append(ientry) 
    to_write += '%(type)s %(nevts_no_cut)i %(nevts_cut)i %(efficiency).3f\n'%ientry
 

# calculate total signal and bkg efficiency of cuts, and signal/bkg before and after cuts
n_total_signal, n_total_bkg , n_total_signal_cut, n_total_bkg_cut = 0,0,0,0
for itype in cut_table:
    if itype['type'] in ['signal','qq','gg'] :
        n_total_signal += itype['nevts_no_cut']
        n_total_signal_cut += itype['nevts_cut']
    else:
        n_total_bkg += itype['nevts_no_cut']
        n_total_bkg_cut += itype['nevts_cut']
to_write += '\nsig/bkg before cut: %.2f\n'%(n_total_signal*1.0/n_total_bkg)
to_write += 'sig/bkg after  cut: %.2f\n'%(n_total_signal_cut*1.0/n_total_bkg_cut)
to_write += 'sig efficiency: %.2f\n'%(n_total_signal_cut*1.0/n_total_signal)
to_write += 'bkg efficiency: %.2f\n'%(n_total_bkg_cut*1.0/n_total_bkg)
# Calcultate total number of events before and after cut
total_nevts = n_total_signal+n_total_bkg
total_nevts_cut = n_total_signal_cut+n_total_bkg_cut
to_write += '\ntotal nevts before cut: %i\n'%total_nevts
to_write += 'total nevts after  cut: %i\n'%total_nevts_cut
to_write += 'fraction of total events left: %.2f\n'%(total_nevts_cut*1.0/total_nevts)
print to_write

file_output.write(to_write)
file_output.close()




    
