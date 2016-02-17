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

parser.add_option('--plot', action='store_true', default=False,
                  dest='plot',
                  help='plot interactively')

parser.add_option('--cut', metavar='F', type='string', action='store',
          default="",
                  dest='cut',
                  help='')
parser.add_option('--var', metavar='F', type='string', action='store',
                  dest='var',
                  help='')
parser.add_option('--Min', metavar='F', type='float', action='store',
          default=0,
                  dest='Min',
                  help='')
parser.add_option('--Max', metavar='F', type='float', action='store',
                  dest='Max',
                  help='')
parser.add_option('--name', metavar='F', type='string', action='store',
              default = "blank",
                  dest='name',
                  help='')
parser.add_option('--log', action='store_true', default=False,
                  dest='log',
                  help='log sacle on y')
parser.add_option('--scale', action='store_true', default=False,
                  dest='scale',
                  help='scale to integral = 1')
parser.add_option('--bin', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin',
                  help='')
parser.add_option('--dir', metavar='F', type='string', action='store',
                  default='',
                  dest='dir',
                  help='')
parser.add_option('--save', action='store_true', default=True,
                  dest='save',
                  help='save plot')
parser.add_option('--title', metavar='F', type='string', action='store',
              default = "",
                  dest='title',
                  help='')
parser.add_option('--xaxis', metavar='F', type='string', action='store',
              default = "",
                  dest='xaxis',
                  help='')
parser.add_option('--yaxis', metavar='F', type='string', action='store',
              default = "",
                  dest='yaxis',
                  help='')
parser.add_option('--label', metavar='F', type='string', action='store',
              default = "",
                  dest='label',
                  help='')
parser.add_option('--config', action='store_true',
              default = False,
                  dest='config',
                  help='')
parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='')
parser.add_option('--weight', metavar='F', type='string', action='store',
                  default='',
                  dest='weight',
                  help='which event weight to use for MC')

(options, args) = parser.parse_args()


scale = options.scale
cut = options.cut
var = options.var
x = options.Min
y = options.Max
log = options.log
bin = options.bin
rundir = options.dir
name = options.name
title = options.title+' '+name
xaxis = options.xaxis
yaxis = options.yaxis
save = options.save
label = options.label
plot = options.plot
verbose = options.verbose
weight = options.weight

# Set root interactive or not
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
if not plot: ROOT.gROOT.SetBatch(True)

# some preset stuff
colors = [2,4,6,8,9]

# parse .txt to get the correct normalization 
data_lumi = 19700

# Some default setting for output names and title
if name == 'blank': name = var
if title == '': title=cut


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

if len(all_types)>len(colors):
    print 'Too many input files for now!'
    sys.exit(0)

# Find tree name
tf = ROOT.TFile(mc_samples[0]['file_path'])
keys = tf.GetListOfKeys()
for ikey in keys:
    if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
print 'Getting ttree',treename
tf.Close()

# add the weight of the tree to L*xsec/nevts
for isample in mc_samples:
    iweight = data_lumi*isample['xsec']/isample['nevts']
    isample['weight'] = iweight
    if verbose in ['yes','verbose']:
        print isample['file_name'],'weight=',iweight

# Loop over types to make plots
ymax = 0
hists = []

for i,itype in enumerate(all_types):
    # Define histograms
    newhist1 = ROOT.TH1F(name+'_'+itype, name+'_'+itype, bin, x, y) 
    newhist1.SetDirectory(0)
    # Loop over samples, find samples of current type, and make a temp hist
    n_samples = 0
    for isample in mc_samples:
        if isample['type'] == itype :
            ifile = ROOT.TFile(isample['file_path'])
            ttree = ifile.Get(treename)
            hist = ROOT.TH1F(name+isample['file_name'], name+isample['file_name'], bin, x, y) 
            ttree.Draw(var+">>"+name+isample['file_name'],""+ cut, "goff")
            newhist1.Add(newhist1,hist,1,isample['weight'])
            n_samples += 1
    newhist1 = overflow(newhist1)
    newhist1.SetDirectory(0)            
    # printing some information about current sample type
    print 'Processing type',itype,'with %i sample files'%n_samples
    # cosmetics
    newhist1.SetStats(0)
    newhist1.SetLineColor(colors[i])
    newhist1.SetLineWidth(1)
    newhist1.SetLineStyle(1)        
    newhist1.SetFillColor(0)

    if xaxis == "":
        newhist1.GetXaxis().SetTitle(var)
    else:
        newhist1.GetXaxis().SetTitle(xaxis)
    if yaxis == "":
        if not scale:
            newhist1.GetYaxis().SetTitle("Events")
        else:
            newhist1.GetYaxis().SetTitle("scaled to Integral = 1")
    else:
        newhist1.GetYaxis().SetTitle(yaxis)
    newhist1.GetYaxis().SetTitleSize(0.04)
    newhist1.GetYaxis().SetTitleOffset(1.6)
    newhist1.GetXaxis().SetTitleOffset(1.2)
    newhist1.GetXaxis().SetTitleSize(0.04)
    if scale:
        newhist1.Scale(1/newhist1.Integral())
        newhist1.SetMinimum(0)
    ymax = max(ymax,newhist1.GetMaximum())
    newhist1.SetTitle(title)
    hists.append(newhist1)

# actual plotting
c = TCanvas()
c.cd()
if log:
    c.SetLogy()
leg = ROOT.TLegend(0.6057047,0.6685824,0.8288591,0.8275862)
leg.SetFillColor(0)
leg.SetLineColor(0) 
leg.SetTextSize(0.03448276)
# Read labels
labels = []
labels = label.split(' ')
total_evts = 0
print '\n'
for i,ihist in enumerate(hists):
    ihist.SetMaximum(ymax*1.1)
    if i==0 : ihist.Draw('HIST')
    else : ihist.Draw("same HIST")
    #leg.SetTextSize(0.02)
    if len(labels) == len(hists):
        leg.AddEntry(ihist,labels[i], "l")
    else:
        leg.AddEntry(ihist,all_types[i], "l")
    print "Entries of type",all_types[i],':',ihist.GetEntries()
    total_evts += ihist.GetEntries()
print '\ntotal number of events:',total_evts

leg.Draw("same")

plotdir = 'plots/'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
    print 'Creating new dir '+plotdir
c.SaveAs(plotdir+name + ".png")
c.SaveAs(plotdir+name + ".root")

    
