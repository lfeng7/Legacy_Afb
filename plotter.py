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
parser.add_option('--file', metavar='F', type='string', action='store',
                  default='',
                  dest='file',
                  help='')
parser.add_option('--save', action='store_true', default=False,
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
  file = ""
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
  label = ""
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
  file = options.file
  name = options.name
  title = options.title
  xaxis = options.xaxis
  yaxis = options.yaxis
  save = options.save
  label = options.label

chain = ROOT.TChain("tree")
chain.Add(file)
newhist = ROOT.TH1F(name, name, bin, x, y)	
chain.Draw(var+">>"+name,""+ cut, "goff")

newhist.SetStats(0)
newhist.SetLineColor(ROOT.kBlack)
newhist.SetLineWidth(1)
newhist.SetLineStyle(1)	
newhist.SetFillColor(0)

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
newhist.GetYaxis().SetTitleOffset(1.05)
newhist.GetXaxis().SetTitleOffset(0.9)
newhist.GetXaxis().SetTitleSize(0.04)

c = TCanvas()
c.cd()
newhist.SetTitle(title)

if log:
	c.SetLogy()
if scale:
	newhist.Scale(1/newhist.Integral())

newhist.Draw()
if label != "":
	leg = ROOT.TLegend(0.65, 0.75, 0.89, 0.89)
	leg.SetFillColor(0)
	leg.SetLineColor(0)
	#leg.SetTextSize(0.02)
	leg.AddEntry(newhist, label, "l")
	leg.Draw("same")
print "entries: " + str(newhist.GetEntries())

if save == True:
	c.SaveAs(name + ".png")
if not save:
	print "Enter save/saveas, or other to close:"
	save = raw_input()
	if save == "save":
		c.SaveAs(name + ".png")
	if save == "saveas":
		print "enter file name:"
		savename = raw_input()
		c.SaveAs(savename + ".png")

