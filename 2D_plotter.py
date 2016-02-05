import os
import glob
import math
import ROOT
from ROOT import *
import sys

from optparse import OptionParser

parser = OptionParser()

parser.add_option('--cut', metavar='F', type='string', action='store',
                  dest='cut',
                  help='')
parser.add_option('--varx', metavar='F', type='string', action='store',
                  dest='var1',
                  help='')
parser.add_option('--vary', metavar='F', type='string', action='store',
                  dest='var2',
                  help='')
parser.add_option('--Minx', metavar='F', type='float', action='store',
                  dest='Min1',
                  help='')
parser.add_option('--Maxx', metavar='F', type='float', action='store',
                  dest='Max1',
                  help='')
parser.add_option('--Miny', metavar='F', type='float', action='store',
                  dest='Min2',
                  help='')
parser.add_option('--Maxy', metavar='F', type='float', action='store',
                  dest='Max2',
                  help='')
parser.add_option('--name', metavar='F', type='string', action='store',
	    	  default = "blahblahblah",
                  dest='name',
                  help='')
parser.add_option('--log', metavar='F', action='store_true',
                  default=False,
                  dest='log',
                  help='')
parser.add_option('--scale', metavar='F', action='store_true',
                  default=False,
                  dest='scale',
                  help='')
parser.add_option('--save', metavar='F', action='store_true',
                  default=False,
                  dest='save',
                  help='')
parser.add_option('--binx', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin1',
                  help='')
parser.add_option('--biny', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin2',
                  help='')
parser.add_option('--file', metavar='F', type='string', action='store',
                  default='no',
                  dest='file',
                  help='')
parser.add_option('--title', metavar='F', type='string', action='store',
                  default='blank',
                  dest='title',
                  help='')
parser.add_option('--xaxis', metavar='F', type='string', action='store',
                  default='blank',
                  dest='xaxis',
                  help='')
parser.add_option('--yaxis', metavar='F', type='string', action='store',
                  default='blank',
                  dest='yaxis',
                  help='')

(options, args) = parser.parse_args()

scale = options.scale
cut = options.cut
var1 = options.var1
var2 = options.var2
x = options.Min1
y = options.Max1
x2 = options.Min2
y2 = options.Max2
log = options.log
bin = options.bin1
bin2 = options.bin2
file = options.file
name = options.name

#f = ROOT.TFile( options.name + ".root", "recreate" )
#f.cd()
c = TCanvas()
c.cd()

# Find the name of the ttree
tf = ROOT.TFile(file)
keys = tf.GetListOfKeys()
for ikey in keys:
    if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
print 'Getting ttree',treename
tf.Close()

chain = ROOT.TChain(treename)
chain.Add(file)
newhist = ROOT.TH2F(name, name, bin, x, y, bin2, x2, y2)	
chain.Draw(var2+":"+var1+">>"+name,""+ cut, "Colz")
if scale:
  newhist.Scale(1/newhist.Integral())
newhist.SetMarkerColor(ROOT.kBlue)
newhist.SetFillColor(0)
newhist.SetLineWidth(2)
newhist.SetLineStyle(2)	
newhist.SetStats(0)
#f.Write()

newhist.SetTitle(options.title)
if options.xaxis == "blank":
  newhist.GetXaxis().SetTitle(var1)
else:
  newhist.GetXaxis().SetTitle(options.xaxis)
if options.yaxis == "blank":
  newhist.GetYaxis().SetTitle(var2)
else:
  newhist.GetYaxis().SetTitle(options.yaxis)

if log:
  c.SetLogz()

if options.save:
  c.SaveAs(name + ".png")

print str(newhist.GetEntries())

plotdir = 'plots/'
rootdir = 'plots/root/'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
    print 'Creating new dir '+plotdir
if not os.path.exists(rootdir):
    os.mkdir(rootdir)
    print 'Creating new dir '+rootdir

c.SaveAs(plotdir+name + ".png")
c.SaveAs(plotdir+'/root/'+name + ".root")

