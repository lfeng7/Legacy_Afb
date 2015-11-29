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
parser.add_option('--cut2', metavar='F', type='string', action='store',
          default="",
                  dest='cut2',
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
parser.add_option('--label2', metavar='F', type='string', action='store',
              default = "",
                  dest='label2',
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
  cut2 = options.cut2
  label2 = options.label2

# Find the name of the ttree
tf = ROOT.TFile(file)
keys = tf.GetListOfKeys()
for ikey in keys:
    if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
print 'Getting ttree',treename
tf.Close()
# Some default setting for output names and title
if name == 'blank': name = var
if title == '': title=cut

chain = ROOT.TChain(treename)
chain.Add(file)
newhist1 = ROOT.TH1F(name, name, bin, x, y)  
chain.Draw(var+">>"+name,""+ cut, "goff")
hists = [newhist1]
all_cuts = [cut]
if cut2!=cut and cut2!='':
    newhist2 = ROOT.TH1F(name+'2', name+'2', bin, x, y)  
    chain.Draw(var+">>"+name+'2',""+ cut2, "goff")
    hists.append(newhist2)
    all_cuts.append(cut2)

icolor = [ROOT.kBlue,ROOT.kRed,ROOT.kGreen]
ihist = 0
ymax = 0
for newhist in hists:
    if len(hists)>1:
        newhist.SetStats(0)
    newhist.SetLineColor(icolor[ihist])
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
    newhist.SetTitle(title)
    if scale:
        newhist.Scale(1/newhist.Integral())
    if newhist.GetMaximum()>ymax: ymax = newhist.GetMaximum()        
    ihist += 1
    
# set ymax for hists
for newhist in hists:
    newhist.SetMaximum(1.2*ymax)

c = TCanvas()
c.cd()

if log:
    c.SetLogy()

newhist1.Draw()
if len(hists)>1: hists[1].Draw('same')

labels=[label,label2]

if label != "" or len(all_cuts)>1:
    leg = ROOT.TLegend(0.65, 0.75, 0.89, 0.89)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    #leg.SetTextSize(0.02)
    ihist = 0
    for newhist in hists:
        if labels[ihist]=='': 
            label_ = all_cuts[ihist]
        else: 
            label_ = labels[ihist]
        leg.AddEntry(newhist, label_, "l")
        ihist +=1
    leg.Draw("same")
print "entries: " + str(newhist.GetEntries())

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

