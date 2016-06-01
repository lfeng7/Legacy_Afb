#! /usr/bin/python

# Calculate isolation extrapolation factor from a TH1
import ROOT
import sys
import glob
import numpy as np
import os

ROOT.gROOT.SetBatch(True)

argv = sys.argv[1:]
if len(argv) == 0:
    print """
    Usage:
    ./hist_plotter.py "*.root"
    """
    sys.exit(1)

argv = ' '.join(argv)
argv = argv.split('--')
# input parameter parsing
for item in argv:
    if 'files' in item:
        all_files = item.split('files')[-1].split()
    elif 'outname' in item:
        outname = item.split('outname')[-1].strip()   

do_norm = True
labels = []

colors = [2,4,6,8,9,11]


hists = []
ymax = 0
for i,item in enumerate(all_files):
    f1 = ROOT.TFile(item)
    keys = f1.GetListOfKeys()
    for ikey in keys:
        if 'TH1' in  ikey.GetClassName() : histname = ikey.GetName()
    ihist = f1.Get(histname)
    ihist.SetDirectory(0)
    htitle = ihist.GetTitle()
    # post process hist
    if do_norm:
        ihist.Scale(1.0/ihist.Integral())
    ihist.SetLineColor(colors[i])
    ymax = max(ihist.GetMaximum(),ymax)
    hists.append(ihist)
    f1.Close()

# plotting
# actual plotting
c = ROOT.TCanvas()
c.cd()

leg = ROOT.TLegend(0.6057047,0.6685824,0.8288591,0.8275862)
leg.SetFillColor(0)
leg.SetLineColor(0) 
leg.SetTextSize(0.03448276)

for i,ihist in enumerate(hists):
    ihist.SetMaximum(ymax*1.1)
    if i==0 : 
        ihist.Draw('hist')
    else : 
        ihist.Draw("same hist")
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

c.SaveAs(plotdir+outname + ".png")

