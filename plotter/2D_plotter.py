import os
import glob
import math
import ROOT
from ROOT import *
import sys

from optparse import OptionParser


common_text = 'CMS Simulation, 19.7 fb^{-1} at #sqrt{s} = 8 TeV'


parser = OptionParser()

parser.add_option('--plot', action='store_true', default=False,
                  dest='plot',
                  help='plot interactively')

parser.add_option('--cut', metavar='F', type='string', action='store',
                  default = '',
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
parser.add_option('--save', metavar='F', type='string', action='store',default='no',
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
files = options.file
name = options.name
save = options.save

# Some global root style 
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
if not plot: ROOT.gROOT.SetBatch(True)


#f = ROOT.TFile( options.name + ".root", "recreate" )
#f.cd()
c = TCanvas()
c.cd()

# Open all input files
files = files.split()
if len(files)>1:
  all_files = files
elif len(files)==1:
  all_files = glob.glob(files[0])
else:
  print 'zero files!'
  sys.exit(1)

# Find the name of the ttree
tf = ROOT.TFile(all_files[0])
keys = tf.GetListOfKeys()
for ikey in keys:
    if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
print 'Getting ttree',treename
tf.Close()

chain = ROOT.TChain(treename)
for ifile in all_files:
  chain.Add(ifile)

# Print out cut efficiency
total_data = chain.GetEntries()
selected_data = chain.GetEntries(cut)
data_cut_eff = selected_data*1.0/total_data
print 'Num Entries in data: %i, selected entries %i, cut efficiency %.3f'%(total_data,selected_data,data_cut_eff)

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

pt = ROOT.TLatex(.23,.80,common_text);
pt.SetNDC(ROOT.kTRUE);
pt.SetTextSize(0.03)
pt.Draw();

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

print str(newhist.GetEntries())

plotdir = 'plots/'
if not os.path.exists(plotdir):
    os.mkdir(plotdir)
    print 'Creating new dir '+plotdir

c.SaveAs(plotdir+name + ".pdf")

if save.lower() in ['true','yes']:
    rootdir = plotdir+'root/'
    if not os.path.exists(rootdir):
        os.mkdir(rootdir)
        print 'Creating new dir '+rootdir  
    c.SaveAs(rootdir+name+ ".root")

