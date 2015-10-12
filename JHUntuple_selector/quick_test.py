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
                  default = 2,
                  dest='maxFiles',
                  help='max number of input ntuple files')

(options, args) = parser.parse_args()

argv = []

#debug
print options.inputFiles

# Get the inputfiles.
if options.inputFiles:
    allfiles = glob.glob( options.inputFiles )
    # Only keep certain number of input files for fexibility
    if len(allfiles)>= options.maxFiles and options.maxFiles > 0 :
        files = [allfiles[i] for i in range(options.maxFiles)]  
    else :
        files = allfiles
    
print 'getting files', files

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
h1 = ROOT.TH1D('met','met;GeV;events',100,0.,200.0)
h2 = ROOT.TH1D('nels',';nels;events',5,0.5,5.5)

h3 = ROOT.TH1D('csv_all_jets', type+' CSV of all jets;csv;events',100,0,1)
h4 = ROOT.TH1D('csv_b_jets', type+' CSV of b jets;csv;events',100,-20,40)
h5 = ROOT.TH1D('csv_light_jets', type+' CSV of light jets;csv;events',100,-20,40)
h6 = ROOT.TH1D('csv_gluon_jets', type+' CSV of gluon jets;csv;events',100,-20,40)

h_jets_eta = ROOT.TH1D('jets_eta',type+' jets eta;eta;events',100,-4.5,4.5)
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

    # get objects
    evt.getByLabel(jet_csv_label,jet_csv_hndl)
    evt.getByLabel(jet_flavor_label,jet_flavor_hndl)
    if not jet_csv_hndl.isValid() : continue 
    n_evt_csv += 1

    csv = jet_csv_hndl.product()
    flavor = jet_flavor_hndl.product()    
    jets = zip(csv,flavor)
    
    evt.getByLabel(jet_p4_label, jet_p4_hndl)
    jets_p4 = jet_p4_hndl.product() 

    # Fill histograms
    for icsv,iflavor in jets : 
        if icsv<0 or icsv>1 : continue
        h3.Fill(icsv)
        if abs(iflavor) == 5 : h4.Fill(icsv)
        if abs(iflavor) in [1,2,3,4] : h5.Fill(icsv)
        if abs(iflavor) == 21 : h6.Fill(icsv)
    
    for ijet in jets_p4 :
        if ijet.pt()>30 :
            h_jets_eta.Fill(ijet.eta())

# End of event loop

# Run summary
print 'break at event',n_evt
print 'number of events with valid jet csv',n_evt_csv

# Plotting and saving 
#histlist = [h3,h4,h5,h6]
histlist = [h_jets_eta]

plotting(histlist,type,"dump")
  
