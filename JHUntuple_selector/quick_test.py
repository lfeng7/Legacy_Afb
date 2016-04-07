from Legacy_Afb.Tools.fwlite_boilerplate import *
# from fwlite_boilerplate import *
import ROOT

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
nevt_cut = -1 
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

jets_csv_hndl = Handle('vector<double>')
jets_csv_label = ("jhuAk5"        ,       "AK5csv")

el_iso_hndl = Handle('vector<double>')
el_iso_label = ("jhuElePFlowLoose"   ,  "electronLooseiso" ,  "jhu" )

electronLooseispseudotight_hndl = Handle('vector<unsigned int>' )
electronLooseispseudotight_label = ("jhuElePFlowLoose"   ,  "electronLooseispseudotight" ,  "jhu")

electronLooseispseudoLoose_hndl = Handle('vector<unsigned int>' )
electronLooseispseudoLoose_label = ("jhuElePFlowLoose"   ,  "electronLooseispseudoloose" ,  "jhu")


electronLooseistight_hndl = Handle('vector<unsigned int>' )
electronLooseistight_label = ("jhuElePFlowLoose"  ,   "electronLooseistight" ,  "jhu")

#PDF
pdf_hndls = [Handle('vector<double>'),Handle('vector<double>'),Handle('vector<double>')]
pdf_label = [("pdfWeights"   ,        "cteq66"      ,      "jhu")]
pdf_label.append( ("pdfWeights"   ,        "CT10"     ,       "jhu") ) 
pdf_label.append( ("pdfWeights"   ,        "GJR08VFnloE"     ,       "jhu") ) 

pdf_types = ["cteq66","CT10","GJR08VFnloE"]

# Make output file with ttree
fout = ROOT.TFile('testtree.root','recreate')
outputtree = ROOT.TTree('selected','selected')


jets_csv_vec = ROOT.vector('float')()
lep_iso_vec = ROOT.vector('float')()
electronLooseispseudotight = ROOT.vector('int')()
electronLooseispseudoLoose = ROOT.vector('int')()
electronLooseistight = ROOT.vector('int')()
vecs = [jets_csv_vec,lep_iso_vec,electronLooseispseudotight,electronLooseispseudoLoose,electronLooseistight]
br_names = ['jets_csv','lep_iso','electronLooseispseudotight','electronLooseispseudoLoose','electronLooseistight']
#pdf
pdf_vecs = [] 
for i,ipdf in enumerate(pdf_label): 
    pdf_vecs += [ROOT.vector('float')()]
    vecs += [pdf_vecs[i]]
    br_names += [ipdf[1]]

branches = zip(br_names,vecs)

for ibr in branches:
    outputtree.Branch(ibr[0],ibr[1])

# Counter initiation 
n_evt = 0
n_evt_csv = 0

print 'Getting',events.size(),'events'
# Event loop
for evt in events:
    for ivec in vecs: ivec.clear()

    # counting and stuff
    if n_evt == nevt_cut: break
    #print 'loop over',n_evt,'events'
    n_evt += 1
    if n_evt%5000 == 0: print 'Loop over',n_evt,'event'

    evt.getByLabel(jets_csv_label, jets_csv_hndl)
    evt.getByLabel(el_iso_label,el_iso_hndl)
    evt.getByLabel(electronLooseispseudotight_label,electronLooseispseudotight_hndl)
    evt.getByLabel(electronLooseistight_label,electronLooseistight_hndl)

    evt.getByLabel(electronLooseispseudoLoose_label,electronLooseispseudoLoose_hndl)


    jets_csv = jets_csv_hndl.product()
    el_iso = el_iso_hndl.product()
    el_is_pseudotight = electronLooseispseudotight_hndl.product()
    el_istight = electronLooseistight_hndl.product() 
    el_is_pseudoLoose = electronLooseispseudoLoose_hndl.product()


    # cuts
    if not el_iso.size()>0 and el_is_pseudotight.size()>0 and jets_csv.size()>0 : continue
    if len([item for item in el_iso if item < 0.15 ])==0: continue # we will see how many events has iso<0.2 electrons

    # pdf
    for i in range(len(pdf_hndls)):
        continue
        evt.getByLabel(pdf_label[i],pdf_hndls[i])
        pdf_ws = pdf_hndls[i].product()
        pdf_w0 = pdf_ws[0]
        for item in pdf_ws:
            pdf_vecs[i].push_back(item/pdf_w0)

#    for ijet in jets_csv: jets_csv_vec.push_back(ijet)
    for i,iel in enumerate(el_iso) : 
        if iel > 0.15: continue
        lep_iso_vec.push_back(iel)
        electronLooseispseudotight.push_back(el_is_pseudotight[i])
        electronLooseistight.push_back(el_istight[i])
        electronLooseispseudoLoose.push_back(el_is_pseudoLoose[i])


    outputtree.Fill()
# End of event loop

# Run summary
print 'break at event',n_evt

fout.Write()
fout.Close()

  
