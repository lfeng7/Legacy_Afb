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
                  dest='cut',
		  default = '',
                  help='')

parser.add_option('--var1', metavar='F', type='string', action='store',
                  dest='var1',
                  help='')

parser.add_option('--var2', metavar='F', type='string', action='store',
                  dest='var2',
                  help='')

parser.add_option('--Min', metavar='F', type='float', action='store',
                  dest='Min',
                  help='')

parser.add_option('--Max', metavar='F', type='float', action='store',
                  dest='Max',
                  help='')

parser.add_option('--name', metavar='F', type='string', action='store',
	    	  default = "blahblahblah",
                  dest='name',
                  help='')

parser.add_option('--log', metavar='F', type='string', action='store',
                  default='no',
                  dest='log',
                  help='')

parser.add_option('--scale', action='store_true', default=False,
                  dest='scale',
                  help='scale to integral = 1')

parser.add_option('--bin', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin',
                  help='')

parser.add_option('--file', metavar='F', type='string', action='store',
                  default='no',
                  dest='fi',
                  help='')

parser.add_option('--save', metavar='F', type='string', action='store',default='no',
                  dest='save',
                  help='save to root file')

parser.add_option('--title', metavar='F', type='string', action='store',
              default = "",
                  dest='title',
                  help='')

parser.add_option('--label', metavar='F', type='string', action='store',
              default = "",
                  dest='label1',
                  help='')

parser.add_option('--label2', metavar='F', type='string', action='store',
              default = "",
                  dest='label2',
                  help='')

parser.add_option('--xaxis', metavar='F', type='string', action='store',
              default = "",
                  dest='xaxis',
                  help='')

parser.add_option('--yaxis', metavar='F', type='string', action='store',
              default = "",
                  dest='yaxis',
                  help='')

(options, args) = parser.parse_args()

scale = options.scale
cut = options.cut
var1 = options.var1
var2 = options.var2
x = options.Min
y = options.Max
log = options.log
bin = options.bin
fi = options.fi
file_name = options.name
save = options.save
title = options.title
label1 = options.label1
label2 = options.label2
xaxis = options.xaxis
yaxis = options.yaxis
plot = options.plot

# Set root interactive or not
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
if not plot: ROOT.gROOT.SetBatch(True)

# Find the name of the ttree
tf = ROOT.TFile(fi)
keys = tf.GetListOfKeys()
for ikey in keys:
    if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
print 'Getting ttree',treename
tf.Close()

#f = ROOT.TFile( options.name + ".root", "recreate" )
#f.cd()

chain = ROOT.TChain(treename)
chain.Add(fi)

name = var1
newhist = ROOT.TH1F(name, name, bin, x, y)	
chain.Draw(var1+">>"+name,""+ cut, "goff")
newhist.SetLineColor(ROOT.kBlue)
newhist.SetFillColor(0)
newhist.SetLineWidth(2)
# newhist.SetLineStyle(2)	
# newhist.SetStats(0)
#f.Write()

name2 = var2
newhist2 = ROOT.TH1F(name2, name2, bin, x, y)	
chain.Draw(var2+">>"+name2,""+ cut, "goff")

newhist2.SetLineColor(ROOT.kRed)
newhist2.SetFillColor(0)
newhist2.SetLineWidth(2)

# newhist2.SetLineStyle(2)	
# newhist2.SetStats(0)
if scale:
      newhist.Scale(1/newhist.Integral())
      newhist2.Scale(1/newhist2.Integral())


if xaxis == "":
    newhist.GetXaxis().SetTitle(var)
else:
    newhist.GetXaxis().SetTitle(xaxis)
if yaxis == "":
    if not scale:
        newhist.GetYaxis().SetTitle("Events")
    else:
        newhist.GetYaxis().SetTitle("scaled to Integral = 1")
else:
    newhist.GetYaxis().SetTitle(yaxis)
newhist.GetYaxis().SetTitleSize(0.04)
# newhist1.GetYaxis().SetTitleOffset(1.05)
newhist.GetXaxis().SetTitleOffset(0.9)
newhist.GetXaxis().SetTitleSize(0.04)


c = TCanvas()
c.cd()
newhist.SetTitle(title)


if log == "yes":
	c.SetLogy()

newhist.SetMaximum(max(newhist.GetMaximum(), newhist2.GetMaximum()) * 1.05)
newhist.SetMinimum(1e-4)
newhist2.SetMinimum(1e-4)



newhist.Draw('hist')
gPad.Update()
statbox1 = newhist.FindObject("stats")
statbox1.SetTextColor(ROOT.kBlue)
statbox1.SetX1NDC(0.8)
statbox1.SetY1NDC(0.83)
statbox1.SetX2NDC(1)
statbox1.SetY2NDC(1)

newhist2.Draw("sames hist")
gPad.Update()
statbox2 = newhist2.FindObject("stats")  
statbox2.SetTextColor(ROOT.kRed)
statbox2.SetX1NDC(0.8)
statbox2.SetY1NDC(0.67)
statbox2.SetX2NDC(1)
statbox2.SetY2NDC(0.83)
statbox2.Draw('sames')

statbox1.Draw('sames')
# leg.Draw("same")

if label1 != "":
    leg = ROOT.TLegend(0.8, 0.450, 0.99, 0.66)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    #leg.SetTextSize(0.02)
    leg.AddEntry(newhist, label1, "l")
    leg.AddEntry(newhist2, label2, "l")
    leg.Draw("same")


print str(newhist.GetEntries())


plotdir = 'plots/'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
    print 'Creating new dir '+plotdir

c.SaveAs(plotdir+file_name + ".png")

if save.lower() in ['true','yes']:
    rootdir = plotdir+'root/'
    if not os.path.exists(rootdir):
        os.mkdir(rootdir)
        print 'Creating new dir '+rootdir  
    c.SaveAs(rootdir+file_name+ ".root")

