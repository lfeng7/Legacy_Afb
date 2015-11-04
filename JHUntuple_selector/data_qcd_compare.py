# This small macro will read in all edm files in a directory and count the total number of events 
# v2. Will take a ttree, make histograms, and stack with data. And calculate correction weights for MC

from utility import *

from optparse import OptionParser

# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--maxfiles', metavar='F', type='int', action='store',
                  default = -1,
                  dest='maxfiles',
                  help='max number of input ntuple files')

parser.add_option('--upload', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='dumpplots',
                  help='if you want to dump plots to MY webpage! Unless you modify fwlite_boilerplate >_<')

parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='If you want more information than you usually need.')

parser.add_option('--MCplots', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='plotMConly',
                  help='If you want stack plots without data comparison')

parser.add_option('--makehists', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='makehists',
                  help='If you want to remake all histograms from selection files')

parser.add_option('--makeplots', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='makeplots',
                  help='If you want to make data MC comparison plots')

parser.add_option('--toppt', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='toppt',
                  help='If do top pt reweighting')

parser.add_option('--correction', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='correction',
                  help='If use correction')

parser.add_option('--tmptype', metavar='F', type='string', action='store',
                  default = 'scaled',
                  dest='tmptype',
                  help='If use correction')

parser.add_option('--applytrigger', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='applytrigger',
                  help='If apply trigger on MC')

parser.add_option('--fakelep', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='fakelep',
                  help='If run on selected events with fake lepton.')

parser.add_option('--lepisocut', metavar='F', type='int', action='store',
                  default = 0.15,
                  dest='lepisocut',
                  help='Lower bound for fake electron isolation.')

parser.add_option('--nbcut', metavar='F', type='int', action='store',
                  default = 1,
                  dest='nbcut',
                  help='Number of b-tagged jets cut')

parser.add_option('--nlepcut', metavar='F', type='int', action='store',
                  default = 1,
                  dest='nlepcut',
                  help='Number of selected leptons cut')

(options, args) = parser.parse_args()

argv = []

# Some preset constants
data_lumi = 19748 
csvm = 0.679

# Get input files
if options.fakelep == 'yes':
    prepend = './selected_files/v3_fakelep/all/'
else:
    prepend = './selected_files/v2_trigger_removed/all/'   # dir of output files to make histograms
postfix='_selected'
# Set up output histogram files
template_type = options.tmptype
tmptype_name = template_type
if options.fakelep == 'yes':
    tmptype_name += '_fakelep'
hist_prepend = './selected_hists/'+tmptype_name+'/' 
    
dir_name = 'data_qcd_'+tmptype_name # name of the dir for control plots, such as corrected, topPT, un_corrected
  
# initialization and declaration
hlist = ['MET','jets_pt','Njets','m3','lep_pt','jets_eta','lep_eta','Nbjets','lep_charge','npv','5jets_pt']
hlist+= ['lep_iso','Nleps']

#### Set up MC input files
# QCD
qcdfile = 'QCD_Pt-15to3000'
datafile = 'SingleEl_Run2012ABCD'

# Getting files
data_input = hist_prepend+datafile+'_control_plots.root'
print 'processing data file',data_input
fdata = ROOT.TFile(data_input)
qcd_input = hist_prepend+qcdfile+'_control_plots.root'
print 'processing qcd file',qcd_input
fqcd= ROOT.TFile(qcd_input)

# Make a legend
leg = ROOT.TLegend(0.85,0.65,1.0,1.0)

# Get histograms

data_hists = []
# Get histograms from data
icolor = ROOT.kYellow
data_scale = 1
for ihist in hlist :
    ih = fdata.Get(ihist)
    if ihist == 'MET' : data_scale = ih.Integral()
    ih.SetDirectory(0)
    ih.SetName(ih.GetName()+'_data')
    ih.SetFillColor(icolor)
    ih.SetLineColor(icolor)
    ih.SetMarkerStyle(21)  
    data_hists.append(ih)
# Add data entry to legend
leg.AddEntry(data_hists[0],'data')

qcd_hists = []
qcd_scale = 1
# Get histograms from qcd
for ihist in hlist :
    ih = fqcd.Get(ihist)
    if ihist == 'MET': qcd_scale = ih.Integral()
for ihist in hlist : 
    ih = fqcd.Get(ihist)   
    ih.SetDirectory(0)
    ih.Scale(1.0*data_scale/qcd_scale)
    ih.SetName(ih.GetName()+'_qcd')
    qcd_hists.append(ih)
# Add data entry to legend
leg.AddEntry(qcd_hists[0],'qcd')

# Make comparison plots
print 'Plotting comparison plots'

# Create the output root file
plotting([],dir_name,'not dump','not log',None,'','recreate')

# Make data MC comparison plots
data_mc = zip(data_hists,qcd_hists)

for item in data_mc :
    comparison_plot_v1(item[0],item[1],leg,dir_name)

# Dump plots to web
if options.dumpplots == 'yes':
    print 'Uploading all plots to web'
    plotting([],dir_name,'dump')

