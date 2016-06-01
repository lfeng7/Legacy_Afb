#! /usr/bin/python

# Calculate isolation extrapolation factor from a TH1
import ROOT
import sys
import glob
import numpy as np

def hist_integral(hist,x1,x2):
    bin1 = hist.FindBin(x1)
    bin2 = hist.FindBin(x2)
    return hist.Integral(bin1,bin2)

argv = sys.argv[1:]
if len(argv) == 0:
    print """
    Usage:
    ./get_factor.py "*.root"
    """
    sys.exit(1)

inputfiles = argv.pop(0)
all_files = glob.glob(inputfiles)  

factors = []
for item in all_files:
    f1 = ROOT.TFile(item)
    keys = f1.GetListOfKeys()
    for ikey in keys:
        if 'TH1' in  ikey.GetClassName() : histname = ikey.GetName()
    ihist = f1.Get(histname)
    htitle = ihist.GetTitle()
    num1 = hist_integral(ihist,0.0,0.1)
    num2 = hist_integral(ihist,0.2,1.2)
    denom = ihist.Integral()
    ifactor = num1*1.0/num2
    factors.append(ifactor)
    print '%s , factor = %.4f'%(htitle,ifactor)
    f1.Close()

avg_fact = np.mean(factors)
stdev_fact = np.std(factors)
print 'avg = %.3f, stdev = %.3f'%(avg_fact,stdev_fact) 
