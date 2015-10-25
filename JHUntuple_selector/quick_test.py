from utility import *

from optparse import OptionParser

import sys

# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "none",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--maxfiles', metavar='F', type='int', action='store',
                  default = -1,
                  dest='maxfiles',
                  help='max number of input ntuple files')

parser.add_option('--startfile', metavar='F', type='int', action='store',
                  default = 0,
                  dest='startfile',
                  help='starting file index of input ntuple files')

(options, args) = parser.parse_args()

argv = []

#debug
print options.inputFiles

# Get the inputfiles.
if options.inputFiles != 'none':
    allfiles = glob.glob( options.inputFiles )
    files = GetSomeFiles(allfiles,options.startfile,options.maxfiles)
    print 'Getting these files:'
    for ifile in files :    
        print ifile

files = ['ntuples/sample_jhudiffmo/TT_jhutester_numEvent1000_99.root']
#files = ['ntuples/sample_jhudiffmo/SingleEl_Run2012A_jhutester_numEvent1000_191.root']

# Read input files
events = Events(files)

# Control constants
nevt_cut = 10000
event_type  = 'test'

# Handles and labels
hndl1 = Handle('vector<double>')
label1 = ('jhuAk5','AK5JEC')
hndl2 = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >')
label2 =  ('jhuAk5','AK5')

jet_p4_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > ')
jet_p4_label = ("jhuAk5","AK5")

trig_hndl = Handle('edm::TriggerResults')
trig_label = ("TriggerResults","","HLT")

# MET
met_phi_hndl = Handle('double')
met_hndl = Handle('double')
met_phi_label = ("jhuGen","metphi")
met_label = ("jhuGen","metpt")

# JHU ntuple format
jet_csv_hndl = Handle('vector<double>')
jet_csv_label = ("jhuAk5","AK5csv")
jet_flavor_hndl = Handle('vector<int>')
jet_flavor_label =("jhuAk5","AK5PartonFlavour")

# Book histograms
h1 = ROOT.TH1D('jetspt','jetspt;pt GeV;events',50,0.,300.0)
# cutflows
h_cutflow = ROOT.TH1D('cutflow_TTbar',event_type+' cutflow;cuts;events',3,0.,3.)
h_cutflow.SetBit(ROOT.TH1.kCanRebin)

# h3 = ROOT.TH1D('csv_all_jets', type+' CSV of all jets;csv;events',100,0,1)


# Counter initiation 
n_evt = 0
n_evt_csv = 0

print 'Getting',events.size(),'events'
# Event loop
for evt in events:

    # counting and stuff
    if n_evt == nevt_cut: break
    #print 'loop over',n_evt,'events'
    n_evt += 1
    if n_evt%5000 == 1: print 'Loop over',n_evt,'event'

    h_cutflow.Fill('all',1)

    evt.getByLabel(trig_label,trig_hndl)
    trig_ = trig_hndl.product()
    iev = evt.object()
    triggerNames = iev.triggerNames(trig_)
    pathName='HLT_Ele27_WP80_v'
    trigName = ''
    for itrig in triggerNames.triggerNames():
        if pathName in itrig : trigName = itrig
    if pathName not in trigName :
        print 'No trigger',pathName,'found in evt',n_evt,'! Will skip this event.'
    passTrig=trig_.accept(triggerNames.triggerIndex(trigName))
    

    if not passTrig : continue
    h_cutflow.Fill('trigger',1)

# End of event loop

# Run summary
print 'break at event',n_evt


# Plotting and saving 
histlist = [h_cutflow]

plotting(histlist,event_type,"dump")
  
