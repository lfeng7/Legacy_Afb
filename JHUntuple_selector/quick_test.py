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

files = ['ntuples/sample_jhudiffmo/SingleEl_Run2012A.root']
# Read input files
events = Events(files)

# Control constants
nevt_cut = 10000
type = 'test'

# Handles and labels
hndl1 = Handle('vector<double>')
label1 = ('jhuAk5','AK5JEC')
hndl2 = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >')
label2 =  ('jhuAk5','AK5')
jet_p4_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > ')
jet_p4_label = ("jhuAk5","AK5")
trig_hndl = Handle('edm::TriggerResults')
trig_label = ("TriggerResults","","HLT")

# JHU ntuple format
jet_csv_hndl = Handle('vector<double>')
jet_csv_label = ("jhuAk5","AK5csv")
jet_flavor_hndl = Handle('vector<int>')
jet_flavor_label =("jhuAk5","AK5PartonFlavour")

# nick old ntuple format
#jet_csv_hndl = Handle('std::vector<float>')
#jet_csv_label = ("pfShyftTupleJets","csv")
#jet_flavor_hndl = Handle('std::vector<float>')
#jet_flavor_label =("pfShyftTupleJets","jetFlavor")

# Book histograms
h1 = ROOT.TH1D('jetspt','jetspt;pt GeV;events',50,0.,300.0)
h2 = ROOT.TH1D('njets',';njets;events',5,3,8)

h3 = ROOT.TH1D('csv_all_jets', type+' CSV of all jets;csv;events',100,0,1)
h4 = ROOT.TH1D('csv_b_jets', type+' CSV of b jets;csv;events',100,-20,40)
h5 = ROOT.TH1D('csv_light_jets', type+' CSV of light jets;csv;events',100,-20,40)
h6 = ROOT.TH1D('csv_gluon_jets', type+' CSV of gluon jets;csv;events',100,-20,40)

h_jets_eta = ROOT.TH1D('jets_eta',type+' jets eta;eta;events',100,-4.5,4.5)

# Make output file
fout = ROOT.TFile('test.root','recreate')
# Add and book ttree
testtree = ROOT.TTree('test','test')

jetsp4 = ROOT.vector('TLorentzVector')()
testtree.Branch('jetsp4',jetsp4)
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

    evt.getByLabel(trig_label,trig_hndl)
    trig_ = trig_hndl.product()

    # testing writing jets p4 into a vector of TLorentzVector  and then into the ttree

    jetsp4.Clear()
    
    evt.getByLabel(jet_p4_label, jet_p4_hndl)
    jets_p4 = jet_p4_hndl.product() 
    
    for ijet in jets_p4 :
        if ijet.pt()>30 :
            ijet_p4 = ROOT.TLorentzVector()
            ijet_p4.SetPtEtaPhiM(ijet.pt(),ijet.eta(),ijet.phi(),ijet.mass())
            jetsp4.push_back(ijet_p4)
    testtree.Fill()

# End of event loop

# Run summary
print 'break at event',n_evt

fout.cd()
testtree.Write()
fout.Close()

filein = ROOT.TFile('test.root')
ttree_ = filein.Get('test')

for i in range(ttree_.GetEntries()):
    ttree_.GetEntry(i)
    jetsp4 = ttree_.jetsp4
    #Fill numjets and pt of jets
    h2.Fill(jetsp4.size())
    for ijet in jetsp4:
        h1.Fill(ijet.Pt())

plotting([h1,h2],'test','dump')
filein.Close()


# Plotting and saving 
#histlist = [h3,h4,h5,h6]

#plotting(histlist,type,"dump")
  
