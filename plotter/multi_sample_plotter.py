import os
import glob
import math
import ROOT
from ROOT import *
import sys
from optparse import OptionParser

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
parser.add_option('--files', metavar='F', type='string', action='store',
                  default='',
                  dest='file1',
                  help='')
parser.add_option('--save', metavar='F', type='string', action='store',default='no',
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

(options, args) = parser.parse_args()

if options.config:
    file1 = ""
    file2 = ""
    var = ""
    x = 0
    y = 100
    bin = 100
    cut = ""
    log = False
    scale = False
    title = ""
    xaxis = ""
    yaxis = ""
    label1 = ""
    label2 = ""
    save = False
    name = ""

else:
    scale = options.scale
    cut = options.cut
    var = options.var
    x = options.Min
    y = options.Max
    log = options.log
    bin = options.bin
    file1 = options.file1
    name = options.name
    title = options.title
    xaxis = options.xaxis
    yaxis = options.yaxis
    save = options.save
    label = options.label
    plot = options.plot

# Set root interactive or not
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
if not plot: ROOT.gROOT.SetBatch(True)

# some preset stuff
colors = [2,4,6,8,9,11]
# Set root interactive or not
if not plot: ROOT.gROOT.SetBatch(True)

# Getting all file names
all_files = file1.split()
if len(all_files)==1:
  all_files = glob.glob(file1)
files = []
print 'Plotting from these files!'
for ifile in all_files:
    files.append(ROOT.TFile(ifile))
    print ifile
# Find the name of the ttree
if len(files)>len(colors):
    print 'Too many input files for now! Only %i colors available'%len(colors)
    sys.exit(1)

tf = files[0]
keys = tf.GetListOfKeys()
for ikey in keys:
    if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
print 'Getting ttree',treename
# tf.Close()
# Some default setting for output names and title
if name == 'blank': name = var
if title == '': title=cut
# loop over files
ymax = 0
hists = []
for i,ifile in enumerate(files):
    ttree = ifile.Get(treename)
    newhist1 = ROOT.TH1F(name+str(i), name+str(i), bin, x, y) 
    ttree.Draw(var+">>"+name+str(i),""+ cut, "goff")

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
for i,ihist in enumerate(hists):
    ihist.SetMaximum(ymax*1.1)
    if i==0 : ihist.Draw('hist')
    else : ihist.Draw("same hist")
    #leg.SetTextSize(0.02)
    if len(labels) == len(hists):
        leg.AddEntry(ihist, labels[i], "l")
    else:
        leg.AddEntry(ihist, 'dist'+str(i), "l")
    print "entries file%i: "%i + str(ihist.GetEntries())

leg.Draw("same")


plotdir = 'plots/'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
    print 'Creating new dir '+plotdir

c.SaveAs(plotdir+name + ".png")

if save.lower() in ['true','yes']:
    rootdir = plotdir+'root/'
    if not os.path.exists(rootdir):
        os.mkdir(rootdir)
        print 'Creating new dir '+rootdir  
    c.SaveAs(rootdir+name+ ".root")

    
