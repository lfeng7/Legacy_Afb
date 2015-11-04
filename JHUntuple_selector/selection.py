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
# Finish ttree making in selection
# Start to include all corrections


from utility import *
from fwlite_boilerplate import *
import os
import glob
import math

# Some predefined var
evt_to_run = -1 
csv_cut = 0.679
trigger_path='HLT_Ele27_WP80_v'
lepiso_cut = 0.2

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

parser.add_option('--applyHLT', metavar='F', type='string', action='store',
                  default = "no",
                  dest='applyHLT',
                  help='If apply HLT as first selection cut.')

parser.add_option('--mcordata', metavar='F', type='string', action='store',
                  default = "mc",
                  dest='mcordata',
                  help='If this run is on data or mc')

parser.add_option('--signal', metavar='F', type='string', action='store',
                  default = "no",
                  dest='isSignal',
                  help='If this run is on signal MC ( TTbar )')

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
                  default = 10000,
                  dest='maxEvts',
                  help='max number of input ntuple files')

parser.add_option('--type', metavar='F', type='string', action='store',
                  default = 'test',
                  dest='evtType',
                  help='type of sample files or the name of the sample, i.e., W4Jets, Data, TT_8TeV etc')

parser.add_option('--mctype', metavar='F', type='string', action='store',
                  default = 'test',
                  dest='sampletype',
                  help='type of sample files')

parser.add_option('--fakelep', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='fakelep',
                  help='If select fake leptons')

# if want to make plots, use multiple input patfiles instead of looping over each files
parser.add_option('--makeplots', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='makeplots',
                  help='If we want to make plots directly. This is for testing purpose mostly, which enables plotting directly.')

(options, args) = parser.parse_args()

argv = []

def main():
    # Get the file list with all input files.
    if options.inputFiles != '':
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

    # debug only
    #files = ['ntuples/sample_jhudiffmo/TT_jhutester_numEvent1000_99.root']
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
            if options.inputFiles == '':
                # find the index of the input file
                index_ = ifile.split('jhutester_numEvent1000_')
                index_ = index_[1].split('.root')
                f_index = int( index_[0])
            print 'processing file  '+ifile
            #print 'current file index is',f_index
            selection(ifile)


# selection is the function to do selection. patfile should be EDM PATtuple files
def selection(rootfiles):

    if options.fakelep == 'yes': 
        btag_cut = 1
        el_postfix = 'Loose'
    else : 
        btag_cut = 2
        el_postfix = ''

    # Get input files
    files = rootfiles
    events = Events(files)
    print 'Getting',events.size(),'events'    

    # Set sample type
    event_type = options.evtType

    # Make a output root file
    if options.grid == 'yes' :
        print '\nRunning in grid mode. Creating outputfile in current dir\n'
        outname = gridsaving([],event_type,f_index,'recreate')
    else :
        print '\nRunning in interactive mode. Creating outputfile to the output dir \n'
        outname = saving([],event_type,f_index,'recreate')
    fout = ROOT.TFile(outname,'update')

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

    # define label module names here
    el_prefix = 'jhuElePFlow'
    el_loose_prefix = 'jhuElePFlowLoose'
    #el_postfix = 'Loose'    
    mu_prefix = 'jhuMuonPFlow'
    muloose_prefix = 'jhuMuonPFlowLoose'

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

    # Trigger
    trig_hndl = Handle('edm::TriggerResults')
    trig_label = ("TriggerResults","","HLT")

    # PU
    dataPileupHandle = Handle('unsigned int')
    dataPileupLabel  = ('jhuGen','npv')
    # MC PU
    npvRealTrueHandle = Handle('unsigned int')
    npvRealTrueLabel  = ('jhuGen','npvTrue')   
      
    # gen info
    gen_hndl = Handle('vector<reco::GenParticle>  ') 
    gen_label = "prunedGenParticles"


    ## Initialization
    pdg_jets = [1,2,3,4,5,21]
    pdg_leps = [11,13,15]
    n_evt = 0
    n_evts_passed,n_PFel,n_PFel_isloose,n_PFel_istight = 0,0,0,0
    timer = ROOT.TStopwatch()
    timer.Start()

    # Book histograms
    # cutflows
    h_cutflow = ROOT.TH1D('cutflow',event_type+' cutflow;cuts;events',7,0.,7.)
    h_cutflow.SetBit(ROOT.TH1.kCanRebin)

    ################################################################ 
    #                   Making TTree                               #
    ################################################################

    outputtree = ROOT.TTree('selected','selected')

    # Data

    # set up vector containers
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
    lep_iso = ROOT.vector('float')()

    met_pt_vec = ROOT.vector('float')()
    met_phi_vec = ROOT.vector('float')()

    trigger_vec = ROOT.vector('bool')()   

    pileup_events = ROOT.vector('float')()

    data_vecs = [jets_pt,jets_eta,jets_phi,jets_mass,jets_csv_vec,lep_pt,lep_eta,lep_phi,lep_mass,lep_charge,met_pt_vec,met_phi_vec]
    data_vecs += [trigger_vec,pileup_events,lep_iso]

    # Set up branches
    branch_names = ['jets_pt','jets_eta','jets_phi','jets_mass','jets_csv','lep_pt','lep_eta','lep_phi','lep_mass','lep_charge','met_pt','met_phi']
    branch_names += ['trigger','pileup_events','lep_iso']

    all_branches = zip(branch_names,data_vecs)
    for ibranch in all_branches:
        outputtree.Branch(ibranch[0],ibranch[1])

    # Book all vectors for initiation in the event loop
    all_vecs = []
    all_vecs += data_vecs

    # MC samples

    if options.mcordata == 'mc' :

        # set up vector containers
        jets_flavor = ROOT.vector('int')()
        mc_pileup_events = ROOT.vector('float')()

        mc_vecs = [jets_flavor,mc_pileup_events]
        mc_branch_names = ['jets_flavor','mc_pileup_events']

        if options.isSignal == 'yes' or options.sampletype == 'ttbar':
            # The sequence of storage in each vector would be :
            #   0      1     2     3     4     5
            # init1, init2, top1, top2, w1, w2
            gen_pt = ROOT.vector('float')()
            gen_eta = ROOT.vector('float')()
            gen_phi = ROOT.vector('float')()
            gen_mass = ROOT.vector('float')()
            gen_pdgid = ROOT.vector('int')()
            # wether the W and/or top is on 'had' or 'lep' side,  or 'NA' 
            gen_side = ROOT.vector('string')() 
            # if this event is hadronic, dilep or e_jets,mu_jets,tau_jets        
            gen_type = ROOT.vector('string')() 
            # if this event is e+jets event
            gen_is_ejets = ROOT.vector('int')()
            # add gen info into MC branches
            mc_vecs += [gen_pt,gen_eta,gen_phi,gen_mass,gen_pdgid,gen_side,gen_type,gen_is_ejets]
            mc_branch_names += ['gen_pt','gen_eta','gen_phi','gen_mass','gen_pdgid','gen_side','gen_type','gen_is_ejets']
            
        ################################################################
        #                   All the corrections                        #
        ################################################################

        # Top pT weights
        weight_top_pT = ROOT.vector('float')()

        # Set vectors for corrections
        mc_vecs += [weight_top_pT]
        mc_branch_names += ['weight_top_pT']

        # Add MC branches to ttree
        all_mc_branches = zip(mc_branch_names,mc_vecs)
        for ibranch in all_mc_branches:
            outputtree.Branch(ibranch[0],ibranch[1])

        # Add MC vectors for intialization
        all_vecs += mc_vecs


    ################################################################
    #                       Start main event loop                  # 
    ################################################################

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
        evt.getByLabel(el_prefix+el_postfix,'electron'+el_postfix,el_hndl)
        evt.getByLabel(el_prefix+el_postfix,'electron'+el_postfix+'iso',el_iso_hndl)
        evt.getByLabel(el_prefix+el_postfix,'electron'+el_postfix+'isloose',el_isLoose_hndl)
        evt.getByLabel(el_prefix+el_postfix,'electron'+el_postfix+'istight',el_isTight_hndl)
        evt.getByLabel(el_prefix+el_postfix,'electron'+el_postfix+'modtight',el_isModTight_hndl)
        evt.getByLabel(el_prefix+el_postfix,'electron'+el_postfix+'charge',el_charge_hndl)

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

        # Get trigger bits
        evt.getByLabel(trig_label,trig_hndl)
        trig_ = trig_hndl.product()
        iev = evt.object()
        triggerNames = iev.triggerNames(trig_)        
        trigName = ''
        for itrig in triggerNames.triggerNames():
            if trigger_path in itrig : trigName = itrig
        if trigger_path not in trigName :
            print 'No trigger',trigger_path,'found in evt',n_evt,'! Will skip this event.'
        passTrig=trig_.accept(triggerNames.triggerIndex(trigName))   
        trigger_vec.push_back(passTrig)     

        # Initialize cutflow histogram
        h_cutflow.Fill("no cut",1)
        # Apply trigger 
        if options.applyHLT == 'yes' :
            if not passTrig : continue
        h_cutflow.Fill('HLT',1)

        # Informations for MC only     
        if options.mcordata == 'mc' :

            jets_PartonFlavor = jet_PartonFlavor_hndl.product()

            # Get gen particles and find out the true identy of the PF electron collection
            evt.getByLabel(gen_label,gen_hndl)
            if not gen_hndl.isValid() :
                print 'No Genparticles info available!'
                continue
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

        ################################################################ 
        #       Physics objects picking and event selections           # 
        ################################################################

        #### PF electrons ####
        el_loose,el_cand = [],[]
        for i in range(len(el_p4)):
            el = el_p4[i]
            icharge = el_charge[i]
            # PFelectrons passed loose selection
            # https://twiki.cern.ch/twiki/bin/view/CMS/TopEGMRun1#Veto
            if el_isLoose[i] and el_iso[i]<0.15 and el.pt()>20 and math.fabs(el.eta())<2.5 : 
                el_loose.append((el,icharge,el_iso[i]))
            # PFelectrons passed tight selection
            # https://twiki.cern.ch/twiki/bin/view/CMS/TopEGMRun1#Signal
            if el_isTight[i] and not el_isModTight[i] and el.pt()>30 and abs(el.eta())<2.5: 
                if options.fakelep == 'no' and el_iso[i]<0.1:
                    el_cand.append((el,icharge,el_iso[i]))
                if options.fakelep == 'yes' and lepiso_cut < el_iso[i] < 1.0 :
                    el_cand.append((el,icharge,el_iso[i]))
        el_extra = list( ipar for ipar in el_loose if ipar not in el_cand)

        #### PF muons ####
        mu_loose = []
        # https://twiki.cern.ch/twiki/bin/view/CMS/TopMUORun1
        for i in range(len(mu_p4)):
            mu = mu_p4[i]
            if mu_is_loose[i] and mu_iso[i]< 0.2 and mu.pt()>10 and abs(mu.eta())<2.5: mu_loose.append(mu)

        # Selection on leptons 
        if options.fakelep == 'yes' and not len(el_cand)>=1  : continue            
        elif not len(el_cand)==1 : continue # continue
        h_cutflow.Fill('el',1)
        if len(mu_loose) > 0 : continue
        h_cutflow.Fill('loose mu veto',1)
        if len(el_extra) > 0 : continue
        h_cutflow.Fill('dilep veto',1)            

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
 
        # Selection on jets
        if not len(jets_cand) >= 4 : continue
        h_cutflow.Fill('jets',1)

        # Do b-tagging. Work for both MC and data
        bjets = [ jet for jet in jets_cand if jet[2] > csv_cut ]

        # Selection on b-tagging
        if not len(bjets) >= btag_cut: continue
        h_cutflow.Fill('b-tagging',1)

        n_evts_passed += 1


        ################################################################
        #                    Fill TTree                                # 
        ################################################################

        # jets
        # First sort jets by pT
        jets_cand.sort(reverse = True)

        for ijet in jets_cand :
            icsv = ijet[2]
            ip4 = ijet[1]
            jets_csv_vec.push_back(icsv)
            jets_pt.push_back(ip4.pt()); jets_eta.push_back(ip4.eta()); jets_phi.push_back(ip4.phi()); jets_mass.push_back(ip4.mass())                        
            if options.mcordata == 'mc' :
                iflavor = ijet[3]
                jets_flavor.push_back(iflavor)
        # lep
        for iel in el_cand:
            lepp4 = iel[0]
            lepcharge = iel[1]
            lepiso = iel[2]
            lep_pt.push_back(lepp4.pt());lep_eta.push_back(lepp4.eta());lep_phi.push_back(lepp4.phi());lep_mass.push_back(lepp4.mass())
            lep_charge.push_back(lepcharge)
            lep_iso.push_back(lepiso)
        # MET
        met_pt_vec.push_back(met_pt[0])
        met_phi_vec.push_back(met_phi[0])

        # NPV
        evt.getByLabel(dataPileupLabel,dataPileupHandle)
        if not dataPileupHandle.isValid() :
            continue
        data_pileup_number = dataPileupHandle.product()
        pileup_events.push_back(1.0*data_pileup_number[0]) 

        if options.mcordata == 'mc' :
            #pileup reweighting
            evt.getByLabel(npvRealTrueLabel,npvRealTrueHandle)
            if not npvRealTrueHandle.isValid() :
                continue
            npvRealTrue = npvRealTrueHandle.product()
            mc_pileup_events.push_back(1.0*npvRealTrue[0])  


        ############################################################
        #               Get some correction weights for MC         #
        ############################################################

        # Initialize all weights
        w_top_pT = 1.0

        if options.mcordata == 'mc' :

            # GenParticles info for signal MC only

            # Look into genparticel info and get Gen Top, W's and initial particles (qqbar, gg etc)
            if options.isSignal == 'yes' or options.sampletype == 'ttbar':
                is_ejets = 0
                # Determine if this event is e+jets event
                if len(gen_el) == 1 and len(gen_el)+len(gen_mu)+len(gen_tau) == 1 : 
                    is_ejets = 1
                gen_is_ejets.push_back(is_ejets)

                # Another way to determine the status of this signal event  
                init_pars = []  
                gentops = []  
                genWs = []
                w_daughters = []       
                for ig in genpars:
                    # Get initial particles
                    if ig.pt()<0: continue
                    # Look through all the particles for protons; append their first daughters to the list
                    if ig.pdgId() == 2212: init_pars.append(ig.daughter(0))
                    # Look through particles for all ts
                    if math.fabs(ig.pdgId()) == 6 and ig.status() == 3 :
                        # By default t is in had side, unless the W from this t decays leptonically
                        whichside = 'had'
                        #look through all the daughters for top to find W.
                        for i in range(ig.numberOfDaughters()) :
                            dau = ig.daughter(i)
                            if math.fabs(dau.pdgId()) == 24 :
                                # Look into daughters of this W
                                # if the W doesn't have two daughters, I don't know what the hell happened.
                                if dau.numberOfDaughters() != 2 :
                                    info = 'W without two daughters. PARTICLE: ' + getId(ig.pdgId()) + ', DAUGHTERS : '
                                    for j in range(ig.numberOfDaughters()) :
                                        info = info + getId(ig.daughter(i)) + ' '
                                    print info
                                    continue  
                                # append all W daughters into the list and decide if it is leptonic or hadronic W
                                for j in range(dau.numberOfDaughters()):
                                    w_daughters.append(dau.daughter(j))
                                    if abs(dau.daughter(j).pdgId()) in pdg_leps: 
                                        whichside = 'lep'
                                # append W
                                genWs.append((dau,whichside))
                        # append to gentop 
                        gentops.append((ig,whichside))

                # Analyse the selected gen particles, write into ttree
                # The sequence of storage in each vector would be :
                #   0      1     2     3     4     5
                # init1, init2, thad, tlep, whad, wlep  

                # Initial particles
                if not len(init_pars) == 2:
                    print 'Events with ',len(init_pars),'initial partons..'
                    continue
                for i in range(2):
                    ig = init_pars[i]
                    gen_pt.push_back(ig.pt())
                    gen_eta.push_back(ig.eta())
                    gen_phi.push_back(ig.phi())
                    gen_mass.push_back(ig.mass())
                    gen_pdgid.push_back(ig.pdgId())
                    gen_side.push_back('NA')
                # GenTops
                if not len(gentops) == 2 :
                    print 'Events with ',len(gentops),'gen tops'
                    continue
                for i in range(2):
                    ig = gentops[i][0]
                    iside = gentops[i][1]
                    gen_pt.push_back(ig.pt())
                    gen_eta.push_back(ig.eta())
                    gen_phi.push_back(ig.phi())
                    gen_mass.push_back(ig.mass())
                    gen_pdgid.push_back(ig.pdgId())
                    gen_side.push_back(iside)
                # GenWs
                if not len(genWs) == 2 :
                    print 'Events with ',len(gentops),'gen Ws'
                    continue
                for i in range(2):
                    ig = genWs[i][0]
                    iside = genWs[i][1]
                    gen_pt.push_back(ig.pt())
                    gen_eta.push_back(ig.eta())
                    gen_phi.push_back(ig.phi())
                    gen_mass.push_back(ig.mass())
                    gen_pdgid.push_back(ig.pdgId())
                    gen_side.push_back(iside) 

                # Find the decay type of the signal ttbar events
                if not len(w_daughters) == 4:
                    print 'Events with ',len(w_daughters),'w daughters'
                    continue
                genleps = []
                for ig in w_daughters:
                    if abs(ig.pdgId()) in pdg_leps : genleps.append(ig.pdgId())
                if len(genleps) == 0 : gen_type.push_back('had')
                elif len(genleps) == 2  and 15 not in genleps and -15 not in genleps: gen_type.push_back('dilep') 
                elif len(genleps) == 2 and 15 in genleps or -15 in genleps : gen_type.push_back('tau_lep')
                elif len(genleps) == 1 :
                    if abs(genleps[0]) == 11 : gen_type.push_back('e_jets')
                    if abs(genleps[0]) == 13 : gen_type.push_back('mu_jets')
                    if abs(genleps[0]) == 15 : gen_type.push_back('tau_jets')
                else:
                    print 'This event has more than 2 leptons!'
                    continue
                    gen_type.push_back('Wrong gen info!')

                # Top pT reweighting
                if gen_type[0] in ['e_jets','mu_jets'] :
                    a = 0.159; b = -0.00141;
                elif gen_type[0] in ['dilep'] :
                    a = 0.148; b = -0.00129;
                else :  # for hadronic or weird situations...
                    a = 0.156; b = -0.00137;
                w_top_pT = GetTopPtWeights(a,b,gen_pt[2],gen_pt[3])

            weight_top_pT.push_back( w_top_pT )

        # Fill all branches
        outputtree.Fill()         

    ######## end main event loop ########
 
    ################################################################
    #                   Make and save plots                        # 
    ################################################################

    h_cutflow_norm = norm(h_cutflow)
       
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
    #        gridsaving(histlist+[outputtree],event_type,f_index,'update')
        else :
            print '\nSaving output into root files to local dir \n'
    #        saving(histlist+[outputtree],event_type,f_index)#,'update')

    for item in histlist+[outputtree]:
        item.SetDirectory(fout) 
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


################################################################
#                    Run main program                          #
################################################################
main()
