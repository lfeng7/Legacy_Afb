# Edit log
# making ttree with selected jets info only


from utility import *
import os
import glob
import math

# Some predefined var
evt_to_run = -1 
csv_cut = 0.679
trigger_path='HLT_Ele27_WP80_v'
lepiso_cut = 0.15

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

parser.add_option('--fakelep', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='fakelep',
                  help='If select fake leptons')

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
    # Each PATtuple file will generate a ntuple files, with index go from 0 to maxFiles
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
def selection(rootfiles):

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
    lep_iso = ROOT.vector('int')()

    met_pt_vec = ROOT.vector('float')()
    met_phi_vec = ROOT.vector('float')()

    trigger_vec = ROOT.vector('bool')()   

    data_vecs = [jets_pt,jets_eta,jets_phi,jets_mass,jets_csv_vec,lep_pt,lep_eta,lep_phi,lep_mass,lep_charge,met_pt_vec,met_phi_vec]
    data_vecs += [trigger_vec,lep_iso]

    # Set up branches
    branch_names = ['jets_pt','jets_eta','jets_phi','jets_mass','jets_csv','lep_pt','lep_eta','lep_phi','lep_mass','lep_charge','met_pt','met_phi']
    branch_names += ['trigger','lep_iso']

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

        mc_vecs = [jets_flavor]
        mc_branch_names = ['jets_flavor']
            
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
            if el_isLoose[i] and el_iso[i]<0.15 and el.pt()>20 and math.fabs(el.eta())<2.5 : el_loose.append((el,icharge,el_iso[i]))
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
        if not len(el_cand)==1 : continue # continue
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
        lepp4 = el_cand[0][0]
        lepcharge = el_cand[0][1]
        lepiso = el_cand[0][2]
        lep_pt.push_back(lepp4.pt());lep_eta.push_back(lepp4.eta());lep_phi.push_back(lepp4.phi());lep_mass.push_back(lepp4.mass())
        lep_charge.push_back(lepcharge)
        lep_iso.push_back(lepiso)
        # MET
        met_pt_vec.push_back(met_pt[0])
        met_phi_vec.push_back(met_phi[0])

        # Fill all branches
        outputtree.Fill()         

    ######## end main event loop ########
 
    ################################################################
    #                   Make and save plots                        # 
    ################################################################

    h_cutflow_norm = norm(h_cutflow)
       
    # cutflows
    histlist = [h_cutflow,h_cutflow_norm]

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
