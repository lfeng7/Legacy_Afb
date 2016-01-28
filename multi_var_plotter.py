import os
import glob
import math
import ROOT
from ROOT import *
import sys
from optparse import OptionParser

parser = OptionParser()

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
parser.add_option('--file1', metavar='F', type='string', action='store',
                  default='',
                  dest='file1',
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

# some preset stuff
colors = [2,4,6,8,9]

# Getting all file names
all_files = glob.glob(file1)
files = []
print 'Plotting from these files!'
for ifile in all_files:
    files.append(ROOT.TFile(ifile))
    print ifile
# Find the name of the ttree
if len(files)>len(colors):
    print 'Too many input files for now!'
    sys.exit(0)

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
    newhist1.GetYaxis().SetTitleOffset(1.05)
    newhist1.GetXaxis().SetTitleOffset(0.9)
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
leg = ROOT.TLegend(0.65, 0.75, 0.89, 0.89)
leg.SetFillColor(0)
leg.SetLineColor(0) 
# Read labels
labels = []
labels = label.split(' ')
for i,ihist in enumerate(hists):
    ihist.SetMaximum(ymax*1.1)
    if i==0 : ihist.Draw()
    else : ihist.Draw("same")
    #leg.SetTextSize(0.02)
    if len(labels) == len(hists):
        leg.AddEntry(ihist, labels[i], "l")
    else:
        leg.AddEntry(ihist, 'dist'+str(i), "l")
    print "entries file%i: "%i + str(ihist.GetEntries())

leg.Draw("same")

plotdir = 'plots/'
if save == True:
    c.SaveAs(plotdir+name + ".png")
if not save:
    print "Enter save/saveas, or other to close:"
    save = raw_input()
    if save == "save":
        c.SaveAs(name + ".png")
    if save == "saveas":
        print "enter file name:"
        savename = raw_input()
        c.SaveAs(savename + ".png")

