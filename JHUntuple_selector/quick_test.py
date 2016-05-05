from Legacy_Afb.Tools.fwlite_boilerplate import *
# from fwlite_boilerplate import *
import ROOT

from optparse import OptionParser

import sys
from array import array

argv = sys.argv[1:]
if len(argv)<2:
    print"""
    Usage:
    python quick_test.py --inputfiles --maxfiles 1 --maxevts 500000 --startfile 0
    """
    sys.exit(1)

timer = ROOT.TStopwatch()
timer.Start()

iso_cut = 1.0
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

parser.add_option('--maxevts', metavar='F', type='int', action='store',
                  default = -1,
                  dest='maxevts',
                  help='max number of input ntuple events')

parser.add_option('--evts_passed', metavar='F', type='int', action='store',
                  default = -1,
                  dest='evts_passed',
                  help='max number of evts_passed to stop')

parser.add_option('--startfile', metavar='F', type='int', action='store',
                  default = 0,
                  dest='startfile',
                  help='starting file index of input ntuple files')

(options, args) = parser.parse_args()

argv = []

#debug
print options.inputFiles

# options
maxfiles = options.maxfiles
startfile = options.startfile

# Get the inputfiles.
if options.inputFiles != 'none':
    files = glob.glob( options.inputFiles )
    files = files[startfile:startfile+maxfiles]
    print 'Getting these files:'
    for ifile in files :    
        print ifile

# files = ['ntuples/sample_jhudiffmo/TT_jhutester_numEvent1000_99.root']
#files = ['ntuples/sample_jhudiffmo/SingleEl_Run2012A_jhutester_numEvent1000_191.root']

# Read input files
events = Events(files)

# Control constants
nevt_cut = options.maxevts
nevt_passed_cut = options.evts_passed
event_type  = 'test'

# Handles and labels
hndl1 = Handle('vector<double>')
label1 = ('jhuAk5','AK5JEC')
hndl2 = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >')
label2 =  ('jhuAk5','AK5')

jets_p4_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > ')
jets_p4_label = ("jhuAk5","AK5")

trig_hndl = Handle('edm::TriggerResults')
trig_label = ("TriggerResults","","HLT")

jets_csv_hndl = Handle('vector<double>')
jets_csv_label = ("jhuAk5"        ,       "AK5csv")

# electrons

el_p4_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >')
el_p4_label = ("jhuElePFlowLoose"    , "electronLoose" ,  "jhu" )

el_iso_hndl = Handle('vector<double>')
el_iso_label = ("jhuElePFlowLoose"   ,  "electronLooseiso" ,  "jhu" )

electronLooseispseudotight_hndl = Handle('vector<unsigned int>' )
electronLooseispseudotight_label = ("jhuElePFlowLoose"   ,  "electronLooseispseudotight" ,  "jhu")

electronLooseispseudoLoose_hndl = Handle('vector<unsigned int>' )
electronLooseispseudoLoose_label = ("jhuElePFlowLoose"   ,  "electronLooseispseudoloose" ,  "jhu")


electronLooseistight_hndl = Handle('vector<unsigned int>' )
electronLooseistight_label = ("jhuElePFlowLoose"  ,   "electronLooseistight" ,  "jhu")

electronLoosemodtight_hndl = Handle('vector<unsigned int>' )
electronLoosemodtight_label = ("jhuElePFlowLoose"  ,   "electronLoosemodtight" ,  "jhu")

# study electrons further
electronLooseMVA_hndl = Handle('vector<double>' )
electronLooseMVA_label = ("jhuElePFlowLoose"  ,   "electronLooseMVA" ,  "jhu")

electronLooseTransverseIP_hndl = Handle('vector<double>' )
electronLooseTransverseIP_label = ("jhuElePFlowLoose"  ,   "electronLooseTransverseIP" ,  "jhu")

electronLooseisEBEEGap_hndl = Handle('vector<unsigned int>' )
electronLooseisEBEEGap_label = ("jhuElePFlowLoose"  ,   "electronLooseisEBEEGap" ,  "jhu")

electronLoosenumberOfHits_hndl = Handle('vector<unsigned int>' )
electronLoosenumberOfHits_label = ("jhuElePFlowLoose"  ,   "electronLoosenumberOfHits" ,  "jhu")

electronLoosepassConversionVeto_hndl = Handle('vector<unsigned int>' )
electronLoosepassConversionVeto_label = ("jhuElePFlowLoose"  ,   "electronLoosepassConversionVeto" ,  "jhu")

#PDF
pdf_hndls = [Handle('vector<double>'),Handle('vector<double>'),Handle('vector<double>')]
pdf_label = [("pdfWeights"   ,        "cteq66"      ,      "jhu")]
pdf_label.append( ("pdfWeights"   ,        "CT10"     ,       "jhu") ) 
pdf_label.append( ("pdfWeights"   ,        "GJR08VFnloE"     ,       "jhu") ) 

pdf_types = ["cteq66","CT10","GJR08VFnloE"]

# Make output file with ttree
output_name = 'testtree_%i_%i.root'%(startfile,maxfiles+startfile)
fout = ROOT.TFile(output_name,'recreate')
outputtree = ROOT.TTree('selected','selected')

print 'output file is %s'%output_name

lep_pt = ROOT.vector('float')()
lep_eta = ROOT.vector('float')()

jets_csv_vec = ROOT.vector('float')()
lep_iso_vec = ROOT.vector('float')()
electronLooseispseudotight = ROOT.vector('int')()
electronLooseispseudoLoose = ROOT.vector('int')()
electronLooseistight = ROOT.vector('int')()

MVA = ROOT.vector('float')()
TransverseIP = ROOT.vector('float')()
isEBEEGap = ROOT.vector('int')()
numberOfHits = ROOT.vector('int')()
passConversionVeto = ROOT.vector('int')()



vecs = [jets_csv_vec,lep_iso_vec,electronLooseispseudotight,electronLooseispseudoLoose,electronLooseistight,lep_pt,lep_eta]
vecs += [MVA,TransverseIP,isEBEEGap,numberOfHits,passConversionVeto]
br_names = ['jets_csv','lep_iso','electronLooseispseudotight','electronLooseispseudoLoose','electronLooseistight','lep_pt','lep_eta']
br_names += ['MVA','TransverseIP','isEBEEGap','numberOfHits','passConversionVeto']

#pdf
pdf_vecs = [] 
for i,ipdf in enumerate(pdf_label): 
    pdf_vecs += [ROOT.vector('float')()]
    vecs += [pdf_vecs[i]]
    br_names += [ipdf[1]]

branches = zip(br_names,vecs)

for ibr in branches:
    outputtree.Branch(ibr[0],ibr[1])

# array formated stuff
br_defs = []

njets = array('i',[-1])
n_el = array('i',[-1])
br_defs += [('njets',njets,'njets/I')]
br_defs += [('n_el',n_el,'n_el/I')]

# Add branches contains array into ttree
for ibr in br_defs:
    outputtree.Branch(ibr[0],ibr[1],ibr[2])

# Counter initiation 
n_evt,n_evts_passed = 0,0
n_evt_csv = 0

# cutflows
h_cutflow = ROOT.TH1D('cutflow',event_type+' cutflow;cuts;events',7,0.,7.)
h_cutflow.SetBit(ROOT.TH1.kCanRebin)

print 'Getting',events.size(),'events'
# Event loop
for evt in events:
    for ivec in vecs: ivec.clear()

    # counting and stuff
    if n_evt == nevt_cut or n_evts_passed == nevt_passed_cut: 
        break

    #print 'loop over',n_evt,'events'
    n_evt += 1
    if n_evt%5000 == 0: 
        print 'Loop over %i evts, %i evts selected'%(n_evt,n_evts_passed)

    evt.getByLabel(jets_csv_label, jets_csv_hndl)
    evt.getByLabel(jets_p4_label, jets_p4_hndl)

    jets_csv = jets_csv_hndl.product()
    jets_p4 = jets_p4_hndl.product()


    evt.getByLabel(el_p4_label,el_p4_hndl)
    evt.getByLabel(el_iso_label,el_iso_hndl)
    evt.getByLabel(electronLooseispseudotight_label,electronLooseispseudotight_hndl)
    evt.getByLabel(electronLooseispseudoLoose_label,electronLooseispseudoLoose_hndl)

    evt.getByLabel(electronLooseistight_label,electronLooseistight_hndl)
    evt.getByLabel(electronLoosemodtight_label,electronLoosemodtight_hndl)

    # other stuff
    evt.getByLabel(electronLooseMVA_label,electronLooseMVA_hndl)
    evt.getByLabel(electronLooseTransverseIP_label,electronLooseTransverseIP_hndl)
    evt.getByLabel(electronLooseisEBEEGap_label,electronLooseisEBEEGap_hndl)
    evt.getByLabel(electronLoosenumberOfHits_label,electronLoosenumberOfHits_hndl)
    evt.getByLabel(electronLoosepassConversionVeto_label,electronLoosepassConversionVeto_hndl)

    # Get opbject from ntuple
    el_iso = el_iso_hndl.product()
    el_is_pseudotight = electronLooseispseudotight_hndl.product()
    el_istight = electronLooseistight_hndl.product() 
    el_is_pseudoLoose = electronLooseispseudoLoose_hndl.product()
    el_isModTight = electronLoosemodtight_hndl.product()
    el_p4 = el_p4_hndl.product()

    electronLooseMVA_prod = electronLooseMVA_hndl.product()
    electronLooseTransverseIP_prod = electronLooseTransverseIP_hndl.product()
    electronLooseisEBEEGap_prod = electronLooseisEBEEGap_hndl.product()
    electronLoosenumberOfHits_prod = electronLoosenumberOfHits_hndl.product()
    electronLoosepassConversionVeto_prod = electronLoosepassConversionVeto_hndl.product()



    # # cuts
    h_cutflow.Fill("no cut",1)

    if not el_iso.size()>0 : 
        continue
    h_cutflow.Fill("no el",1)

    selected_jets = []
    for i,ijet in enumerate(jets_p4):
        if ijet.pt()>30 and abs(ijet.eta())<2.5: 
            selected_jets.append(i)
    if len(selected_jets)!=2:
        continue
    h_cutflow.Fill("not 2 pT>30 jets",1)

    # #### PF electrons selection ####
    # el_cand = []
    # for i in range(len(el_p4)):
    #     if el_is_pseudoLoose[i]:
    #         el_cand.append(i)
    # if not len(el_cand)>0:
    #     continue
    # h_cutflow.Fill("PseudoLoose",1)

    # el_cand = []
    # for i in range(len(el_p4)):
    #     if el_is_pseudoLoose[i] and el_iso[i] < iso_cut :
    #         el_cand.append(i)
    # if not len(el_cand)>0:
    #     continue
    # h_cutflow.Fill("iso<%.1f"%iso_cut,1)

    # el_cand = []
    # for i in range(len(el_p4)):
    #     if el_is_pseudoLoose[i] and not el_isModTight[i] and el_iso[i] < iso_cut  :
    #         el_cand.append(i)
    # if not len(el_cand)>0:
    #     continue
    # h_cutflow.Fill("ECAL_gap",1)


    # el_cand = []
    # for i in range(len(el_p4)):
    #     el = el_p4[i]
    #     if el_is_pseudoLoose[i] and not el_isModTight[i] and el_iso[i] < iso_cut and el.pt()>30 :
    #         el_cand.append(i)
    # if not len(el_cand)>0:
    #     continue
    # h_cutflow.Fill("pt>30",1)

    # el_cand = []
    # for i in range(len(el_p4)):
    #     el = el_p4[i]
    #     if el_is_pseudoLoose[i] and not el_isModTight[i] and el_iso[i] < iso_cut and el.pt()>30 and abs(el.eta())<2.5 :
    #         el_cand.append(i)
    # if not len(el_cand)>0:
    #     continue
    # h_cutflow.Fill("|#eta|<2.5",1)

    # # extra loose leptons
    # el_loose = []
    # for i in range(len(el_p4)):
    #     el = el_p4[i]
    #     if el_is_pseudoLoose and el_iso[i]<0.15 and el.pt()>20 and math.fabs(el.eta())<2.5 : 
    #         el_loose.append(i)

    # el_extra = list( ipar for ipar in el_loose if ipar not in el_cand)
    # if len(el_extra)>0:
    #     continue
    # h_cutflow.Fill("dilep",1)


    # # pdf
    # for i in range(len(pdf_hndls)):
    #     continue
    #     evt.getByLabel(pdf_label[i],pdf_hndls[i])
    #     pdf_ws = pdf_hndls[i].product()
    #     pdf_w0 = pdf_ws[0]
    #     for item in pdf_ws:
    #         pdf_vecs[i].push_back(item/pdf_w0)


    # write some branches into ttree
    njets[0] = jets_p4.size()
    n_el[0] = el_iso.size()

    for i,item in enumerate(el_iso) : 
        el = el_p4[i]
        lep_iso_vec.push_back(el_iso[i])
        electronLooseispseudotight.push_back(el_is_pseudotight[i])
        electronLooseistight.push_back(el_istight[i])
        electronLooseispseudoLoose.push_back(el_is_pseudoLoose[i])
        lep_pt.push_back(el.pt())
        lep_eta.push_back(el.eta())
        # other electron stuff
        MVA.push_back(electronLooseMVA_prod[i])
        TransverseIP.push_back(electronLooseTransverseIP_prod[i])
        isEBEEGap.push_back(electronLooseisEBEEGap_prod[i])
        numberOfHits.push_back(electronLoosenumberOfHits_prod[i])
        passConversionVeto.push_back(electronLoosepassConversionVeto_prod[i])
    
    outputtree.Fill()
    n_evts_passed += 1

# End of event loop

# Run summary
print 'break at event %i, selected evts %i'%(n_evt,n_evts_passed)

fout.Write()
fout.Close()

# Stop our timer
timer.Stop()
# Print out our timing information
print '\n'
rtime = timer.RealTime(); # Real time (or "wall time")
ctime = timer.CpuTime(); # CPU time
print("RealTime={0:6.2f} seconds, CpuTime={1:6.2f} seconds").format(rtime,ctime)
print("{0:4.2f} events / RealTime second .").format( n_evt/rtime)
print("{0:4.2f} events / CpuTime second .").format( n_evt/ctime)
print("{0:4.2f} candidate events / RealTime second .").format( n_evts_passed/rtime)
print("{0:4.2f} candidate events / CpuTime second .").format( n_evts_passed/ctime)
# Run summary
print("Analyzed events: {0:6d}").format(n_evt)
print("Candidate events: {0:6d}").format(n_evts_passed)
print '\n'

