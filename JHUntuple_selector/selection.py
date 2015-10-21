# Edit log
# This is version 1. Rewrite the structure of code to allow the main selection as a main function
# Enabling two different running modes. 
# (A) Plotting mode. Where all patuple files will be put in a file list and only one set of plots will be made.
# (B) Non-plotting mode. Where no plots will be made. Each input pattuple file will generate an output root file. All files 
# ... from same type of samples are put in the same directory.
# Edit log
# This is the grid version of v1 selection
# Add a grid option
# Mostly instead of outputting to local dir, directlly write outputfiles to current dir
# Maybe I can add some functionallity to rename to output dir as the --type?
# Edit log
# Start the v2 development.


from utility import *
import os
import glob
import math

# Some control var
evt_to_run = -1 
csv_cut = 0.679
#event_type = 'Powheg_TT_btag'
events_passed = -1

f_index = 0

from optparse import OptionParser

parser = OptionParser()


############################################
#            Job steering                  #
############################################

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 
parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--txtfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='txtfiles',
                  help='Input txt files')

parser.add_option('--mcordata', metavar='F', type='string', action='store',
                  default = "mc",
                  dest='mcordata',
                  help='If this run is on data or mc')

parser.add_option('--grid', metavar='F', type='string', action='store',
                  default = "no",
                  dest='grid',
                  help='If will run on grid using condor')

parser.add_option('--maxfiles', metavar='F', type='int', action='store',
                  default = -1,
                  dest='maxfiles',
                  help='max number of input ntuple files')

parser.add_option('--startfile', metavar='F', type='int', action='store',
                  default = 0,
                  dest='startfile',
                  help='starting file index of input ntuple files')

parser.add_option('--maxevts', metavar='F', type='int', action='store',
                  default = 100000,
                  dest='maxEvts',
                  help='max number of input ntuple files')

parser.add_option('--type', metavar='F', type='string', action='store',
                  default = 'test',
                  dest='evtType',
                  help='type of sample files')

# if want to make plots, use multiple input patfiles instead of looping over each files
parser.add_option('--makeplots', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='makeplots',
                  help='If we want to make plots directly. This is for testing purpose mostly, which enables plotting directly.')

(options, args) = parser.parse_args()

argv = []


def main():
    # Get the file list with all input files.
    if options.inputFiles:
        allfiles = glob.glob( options.inputFiles )
    elif options.txtfiles:
        allfiles = []
        with open(options.txtfiles, 'r') as input_:
            for line in input_:
                print 'Getting files from this dir '+line.strip()
                somefiles =  glob.glob(line.strip())
                allfiles.extend(somefiles)
    else:
        allfiles = []

    # Job splitting is done here
    # Only keep certain number of input files for fexibility
    files = GetSomeFiles(allfiles,options.startfile,options.maxfiles)

    # Print out information on the input files
    print 'Getting these files:'
    for ifile in files : print ifile
   
    # Run selection function to do selections
    # If we want to make plots, use many files input form
    # If not make plots, each PATtuple file will generate a ntuple files, with index go from 0 to maxFiles
    if options.makeplots == 'True' :
        selection(files)
    else :
        global f_index
        f_index = 0
        for ifile in files:
            # find the index of the input file
            index_ = ifile.split('jhutester_numEvent1000_')
            index_ = index_[1].split('.root')
            f_index = int( index_[0])
            print 'processing file  '+ifile
            #print 'current file index is',f_index
            selection(ifile)


# selection is the function to do selection. patfile should be EDM PATtuple files
def selection(patfile):

    # Get input root files
    files = patfile

    # Set sample type
    event_type = options.evtType

    # Get input
    events = Events(files)
    print 'Getting',events.size(),'events'

    ######## Define handles here ########

    # leptons
    el_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >')
    el_iso_hndl = Handle('vector<double>')
    el_charge_hndl = Handle('vector<int>')
    el_isLoose_hndl = Handle('vector<unsigned int>')
    el_isTight_hndl = Handle('vector<unsigned int>')
    el_isModTight_hndl = Handle('vector<unsigned int>')

    mu_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >')
    mu_label = ("jhuMuonPFlow", "muon")
    mu_iso_hndl = Handle('vector<double>')
    mu_iso_label = ("jhuMuonPFlow","muoniso" )
    mu_charge_hndl = Handle('vector<double>')
    mu_charge_label = ("jhuMuonPFlow","muoncharge")
    mu_isLoose_hndl = Handle('vector<unsigned int>')
    mu_isLoose_label = ("jhuMuonPFlow","muonisloose")
    mu_isTight_hndl = Handle('vector<unsigned int>')
    mu_isTight_label = ("jhuMuonPFlow","muonistight")

    # MET
    met_phi_hndl = Handle('double')
    met_hndl = Handle('double')
    met_phi_label = ("jhuGen","metphi")
    met_label = ("jhuGen","metpt")

    # AK5 Jets
    jet_p4_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > > ')
    jet_p4_label = ("jhuAk5","AK5")
    jet_csv_hndl = Handle('vector<double>')
    jet_csv_label = ("jhuAk5","AK5csv")
    jet_PartonFlavor_hndl = Handle('vector<int> ')
    jet_PartonFlavor_label = ( "jhuAk5","AK5PartonFlavour") 

    # gen info
    gen_hndl = Handle('vector<reco::GenParticle>  ') 
    gen_label = "prunedGenParticles"

    # define label module names here
    el_prefix = 'jhuElePFlow'
    el_loose_prefix = 'jhuElePFlowLoose'
    mu_prefix = 'jhuMuonPFlow'
    muloose_prefix = 'jhuMuonPFlowLoose'


    ## Initialization
    pdg_jets = [1,2,3,4,5,21]
    n_evt = 0
    n_evts_passed,n_PFel,n_PFel_isloose,n_PFel_istight = 0,0,0,0
    timer = ROOT.TStopwatch()
    timer.Start()

    # cutflows
    h_cutflow = ROOT.TH1D('cutflow',event_type+' cutflow;cuts;events',7,0.,7.)
    h_cutflow.SetBit(ROOT.TH1.kCanRebin)

    # generator level info
    h_num_gen_b = ROOT.TH1D('NGenbjets',event_type+' Num Gen bjets;Nbjets;events',5,1,6)
    h_num_gen_jets = ROOT.TH1D('NGenJets',event_type+' Num Gen jets;Njets;events',9,1,10) 

    # Making ttree
    outputtree = ROOT.TTree('selected','selected')
    # Setup container

    # branches for data and mc
    jets_pt = ROOT.vector('float')()
    jets_eta = ROOT.vector('float')()
    jets_phi = ROOT.vector('float')()
    jets_mass = ROOT.vector('float')()
    jets_csv_vec = ROOT.vector('float')()

    lep_pt = ROOT.vector('float')()
    lep_eta = ROOT.vector('float')()
    lep_phi = ROOT.vector('float')()
    lep_mass = ROOT.vector('float')()
    lep_charge = ROOT.vector('int')()

    met_pt_vec = ROOT.vector('float')()
    met_phi_vec = ROOT.vector('float')()

    # branches only for MC samples
    jets_flavor = ROOT.vector('int')()

    gen_whad = ROOT.vector('TLorentzVector')()
    gen_whad_pdgid = ROOT.vector('int')()

    gen_wlep = ROOT.vector('TLorentzVector')()
    gen_wlep_pdgid = ROOT.vector('int')()

    gen_thad = ROOT.vector('TLorentzVector')()
    gen_thad_pdgid = ROOT.vector('int')()

    gen_tlep = ROOT.vector('TLorentzVector')()
    gen_tlep_pdgid = ROOT.vector('int')()

    gen_q = ROOT.vector('TLorentzVector')()
    gen_q_pdgid = ROOT.vector('int')()

    gen_qbar = ROOT.vector('TLorentzVector')()
    gen_qbar_pdgid = ROOT.vector('int')()

    all_vecs = [jets_pt,jets_eta,jets_phi,jets_mass,jets_csv_vec,lep_pt,lep_eta,lep_phi,lep_mass,lep_charge,met_pt_vec,met_phi_vec]
    all_mc_vecs = [jets_flavor]

    # Set up branches
    branch_names = ['jets_pt','jets_eta','jets_phi','jets_mass','jets_csv','lep_pt','lep_eta','lep_phi','lep_mass','lep_charge','met_pt','met_phi']
    mc_branch_names = ['jets_flavor']

    all_branches = zip(branch_names,all_vecs)
    all_mc_branches = zip(mc_branch_names,all_mc_vecs)

    if options.mcordata == 'mc':
        all_vecs = all_vecs+all_mc_vecs
        all_branches = all_branches+all_mc_branches

    for ibranch in all_branches:
        outputtree.Branch(ibranch[0],ibranch[1])

    # Start main event loop
    for evt in events:
        # progrss reporting
        if n_evt%5000 == 1: print 'Loop over',n_evt,'event',', selected',n_evts_passed,'candidate events'
        if n_evts_passed == events_passed : 
            print 'reached',n_evts_passed,'candidate events'
            break
        if n_evt == options.maxEvts: 
            print 'reached',n_evt,'events looped'
            break
        n_evt += 1

        # Reset all vector containers
        for ivec in all_vecs: ivec.clear()

        # Read objects in nTuple
        evt.getByLabel(el_prefix,'electron',el_hndl)
        evt.getByLabel(el_prefix,'electroniso',el_iso_hndl)
        evt.getByLabel(el_prefix,'electronisloose',el_isLoose_hndl)
        evt.getByLabel(el_prefix,'electronistight',el_isTight_hndl)
        evt.getByLabel(el_prefix,'electronmodtight',el_isModTight_hndl)
        evt.getByLabel(el_prefix,'electroncharge',el_charge_hndl)

        evt.getByLabel('jhuMuonPFlow','muon',mu_hndl)
        evt.getByLabel('jhuMuonPFlow','muoniso',mu_iso_hndl)
        evt.getByLabel('jhuMuonPFlow','muonisloose',mu_isLoose_hndl)
      
        evt.getByLabel(met_label,met_hndl)
        evt.getByLabel(met_phi_label,met_phi_hndl)

        evt.getByLabel(jet_p4_label, jet_p4_hndl)
        evt.getByLabel(jet_csv_label, jet_csv_hndl)

        if options.mcordata == 'mc' :
            evt.getByLabel(jet_PartonFlavor_label, jet_PartonFlavor_hndl)  # not for data

        el_p4 = el_hndl.product()
        el_iso = el_iso_hndl.product()
        el_isLoose = el_isLoose_hndl.product()
        el_isTight = el_isTight_hndl.product()
        el_isModTight = el_isModTight_hndl.product()
        el_charge = el_charge_hndl.product()

        mu_p4 = mu_hndl.product()
        mu_is_loose = mu_isLoose_hndl.product()
        mu_iso = mu_iso_hndl.product()

        met_pt = met_hndl.product()
        met_phi = met_phi_hndl.product()

        jets_p4 = jet_p4_hndl.product()
        jets_csv = jet_csv_hndl.product()

        # Informations for MC only     
        if options.mcordata == 'mc' :

            jets_PartonFlavor = jet_PartonFlavor_hndl.product()

            # Get gen particles and find out the true identy of the PF electron collection
            evt.getByLabel(gen_label,gen_hndl)
            genpars = gen_hndl.product()
            # get all final state particles,status = 3,and particles with no daughters(final state partons) 
            # or particles from W's which may have daughters
            final_par = []
            for ipar in genpars:
                if ipar.status() == 3:
                    if ipar.numberOfDaughters()==0 or abs(ipar.pdgId())==5 : final_par.append(ipar)
                    elif ipar.numberOfMothers():
                        if abs(ipar.mother(0).pdgId())==24: final_par.append(ipar)
            # Find events with electron in final states
            gen_el = list(ipar for ipar in final_par if ipar.status() == 3 and abs(ipar.pdgId()) == 11)
            gen_mu = list(ipar for ipar in final_par if ipar.status() == 3 and abs(ipar.pdgId()) == 13)
            gen_tau = list(ipar for ipar in final_par if ipar.status() == 3 and abs(ipar.pdgId()) == 15)
            gen_b = list(ipar for ipar in final_par if abs(ipar.pdgId()) == 5)
            gen_jets = list(ipar for ipar in final_par if abs(ipar.pdgId()) in pdg_jets)
            # Determine if this event is e+jets event
            if len(gen_el) == 1 and len(gen_el)+len(gen_mu)+len(gen_tau) == 1 : 
                Is_ejets = True

        ################ Find all physics objects for reconstruction ###################

        #### PF electrons ####
        el_loose,el_cand = [],[]
        for i in range(len(el_p4)):
            el = el_p4[i]
            icharge = el_charge[i]
            # PFelectrons passed loose selection
            # https://twiki.cern.ch/twiki/bin/view/CMS/TopEGMRun1#Veto
            if el_isLoose[i] and el_iso[i]<0.15 and el.pt()>20 and math.fabs(el.eta())<2.5 : el_loose.append((el,icharge))
            # PFelectrons passed tight selection
            # https://twiki.cern.ch/twiki/bin/view/CMS/TopEGMRun1#Signal
            if el_isTight[i] and not el_isModTight[i] and el_iso[i]<0.1 and el.pt()>30 and abs(el.eta())<2.5: el_cand.append((el,icharge))
        el_extra = list( ipar for ipar in el_loose if ipar not in el_cand)

        #### PF muons ####
        mu_loose = []
        # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUORun1
        for i in range(len(mu_p4)):
            mu = mu_p4[i]
            if mu_is_loose[i] and mu_iso[i]< 0.2 and mu.pt()>10 and abs(mu.eta())<2.5: mu_loose.append(mu)

        ##### AK5 jets ####

        # jets Selection       https://twiki.cern.ch/twiki/bin/view/CMS/TopJMERun1#Jets
        jets_cand = []
        for i in range(len(jets_p4)):
            if jets_p4[i].pt()>30 and abs(jets_p4[i].eta())<2.5: 
                if options.mcordata == 'mc' : 
                    jets_cand.append((jets_p4[i].pt(),jets_p4[i],jets_csv[i],jets_PartonFlavor[i]))
                elif options.mcordata == 'data' :
                    jets_cand.append((jets_p4[i].pt(),jets_p4[i],jets_csv[i]))
                else :
                    print 'The sample is neither mc or data! Serious bug!'
                    break
        jets_cand_p4 = [ ijet[1] for ijet in jets_cand ]
 
        #### Signal events selection ####

        # Initialize cutflow histogram
        h_cutflow.Fill("step0",1)
 
        if not len(el_cand)==1 : continue # continue
        h_cutflow.Fill('el',1)
        if len(mu_loose) > 0 : continue
        h_cutflow.Fill('loose mu veto',1)
        if len(el_extra) > 0 : continue
        h_cutflow.Fill('dilep veto',1)
        if not len(jets_cand) >= 4 : continue
        h_cutflow.Fill('jets',1)

        # Do b-tagging. Work for both MC and data
        bjets = [ jet for jet in jets_cand if jet[2] > csv_cut ]

        # cut on btags
        if not len(bjets) >= 2: continue
        h_cutflow.Fill('b-tagging',1)

        n_evts_passed += 1

        ######## Finish event selection #########

        ######## Fill TTree ########

        # jets
        # First sort jets by pT
        jets_cand.sort()
        
        for ijet in jets_cand :
            icsv = ijet[2]
            ip4 = ijet[1]
            jets_csv_vec.push_back(icsv)
            jets_pt.push_back(ip4.pt()); jets_eta.push_back(ip4.eta()); jets_phi.push_back(ip4.phi()); jets_mass.push_back(ip4.mass())                        
            if options.mcordata == 'mc' :
                iflavor = ijet[3]
                jets_flavor.push_back(iflavor)
        # lep
        lepp4 = el_cand[0][0]
        lepcharge = el_cand[0][1]
        lep_pt.push_back(lepp4.pt());lep_eta.push_back(lepp4.eta());lep_phi.push_back(lepp4.phi());lep_mass.push_back(lepp4.mass())
        lep_charge.push_back(lepcharge)
        # MET
        met_pt_vec.push_back(met_pt[0])
        met_phi_vec.push_back(met_phi[0])

        # Fill all branches
        outputtree.Fill()

    ######## end main event loop ########

    h_cutflow_norm = norm(h_cutflow)
    h_cutflow_log = h_cutflow.Clone()
    h_cutflow_log.SetName(h_cutflow.GetName()+'_log')
    h_cutflow_norm_log = h_cutflow_norm.Clone()
    h_cutflow_norm_log.SetName(h_cutflow_norm.GetName()+'_log')
       
    ## Make and save plots

    # cutflows
    histlist = [h_cutflow,h_cutflow_norm]

    if options.makeplots == 'yes' and options.grid != 'yes':
        print '\nPlot and saven'
        plotting(histlist,event_type,'not dump')
        histlist1 = [h_cutflow_log,h_cutflow_norm_log]
        plotting(histlist1,event_type,"dump","setlogy")
        # Save to root files
        gridsaving(histlist,event_type,'hists')
    else :
        if options.grid == 'yes' :
            print '\nSaving output into root files for grid use\n'
            gridsaving(histlist,event_type,f_index)
        else :
            print '\nSaving output into root files to local dir \n'
            saving(histlist+[outputtree],event_type,f_index)

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



#########################  Actual process  ##############################
main()
