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

parser.add_option('--cut', metavar='F', type='string', action='store',default='',
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

parser.add_option('--log', metavar='F', type='string', action='store',
                  default='no',
                  dest='log',
                  help='')

parser.add_option('--scale', metavar='F', type='float', action='store',
                  default='1.0',
                  dest='scale',
                  help='')

parser.add_option('--binx', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin1',
                  help='')

parser.add_option('--biny', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin2',
                  help='')

parser.add_option('--file1', metavar='F', type='string', action='store',
                  default='no',
                  dest='fi1',
                  help='')

parser.add_option('--file2', metavar='F', type='string', action='store',
                  default='no',
                  dest='fi2',
                  help='')

(options, args) = parser.parse_args()


plot = options.plot
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
fi1 = options.fi1
fi2 = options.fi2
name = options.name

#f = ROOT.TFile( options.name + ".root", "recreate" )
#f.cd()


# Some global root style 
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
if not plot: ROOT.gROOT.SetBatch(True)

# Find the name of the ttree
tf = ROOT.TFile(fi1)
keys = tf.GetListOfKeys()
for ikey in keys:
    if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
print 'Getting ttree',treename

chain = ROOT.TChain(treename)
chain.Add(fi1)
newhist = ROOT.TH2F(name, name, bin, x, y, bin2, x2, y2)	
chain.Draw(var2+":"+var1+">>"+name,""+ cut, "lego2")
#newhist.Scale(1/newhist.Integral())
newhist.SetMarkerColor(ROOT.kBlue)
newhist.SetFillColor(0)
newhist.SetLineWidth(2)
newhist.SetLineStyle(2)	
newhist.SetStats(0)
#f.Write()

chain2 = ROOT.TChain(treename)
chain2.Add(fi2)
newhist2 = ROOT.TH2F(name+"2", name, bin, x, y, bin2, x2, y2)	
chain2.Draw(var2+":"+var1+">>"+name+"2",""+ cut, "lego2")
newhist2.Scale(1/newhist.Integral())
newhist2.SetMarkerColor(ROOT.kRed)
newhist2.SetFillColor(0)
newhist2.SetLineWidth(2)
newhist2.SetLineStyle(2)	
newhist2.SetStats(0)

c = TCanvas()
c.cd()
newhist.SetTitle(name)
newhist.GetXaxis().SetTitle(var1)
newhist.GetYaxis().SetTitle(var2)

newhist.Draw()
newhist2.Draw("same")

print str(newhist.GetEntries())
print str(newhist2.GetEntries())

plotdir = 'plots/'
rootdir = 'plots/root/'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
    print 'Creating new dir '+plotdir
if not os.path.exists(rootdir):
    os.mkdir(rootdir)
    print 'Creating new dir '+rootdir

c.SaveAs(plotdir+name + ".png")


