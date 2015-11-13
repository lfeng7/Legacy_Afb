# Take in the selection output ttree, do reco reconstruction via kinfit

from Legacy_Afb.Tools.kinfit import *
from Legacy_Afb.Tools.ttbar_utility import *
import glob
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

parser.add_option('--evtsperjob', metavar='F', type='int', action='store',
                  default = 1000,
                  dest='evtsperjob',
                  help='number of events to run for each job')

parser.add_option('--evtstart', metavar='F', type='int', action='store',
                  default = 0,
                  dest='evtstart',
                  help='the evt to start')

parser.add_option('--upload', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='dumpplots',
                  help='if you want to dump plots to MY webpage! Unless you modify fwlite_boilerplate >_<')

parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='If you want more information than you usually need.')

parser.add_option('--fakelep', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='fakelep',
                  help='If run on selected events with fake lepton.')

parser.add_option('--lepisocut', metavar='F', type='float', action='store',
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

parser.add_option('--applytrigger', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='applytrigger',
                  help='If apply trigger on MC')

(options, args) = parser.parse_args()

argv = []

# Some preset constants
data_lumi = 19748 
csvm = 0.679

# Get input files
if options.fakelep == 'yes':
    prepend = './selected_files/v3_fakelep_updated/all/'
else:
    prepend = './selected_files/v2_trigger_removed/all/'   # dir of output files to make histograms
postfix='_selected'  

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

    timer = ROOT.TStopwatch()
    timer.Start()

    # Job splitting is done here
    for ifile in allfiles:
        sample_name = ifile.split('/')
        sample_name = sample_name[len(sample_name)-1].split('.root')[0]
        sample_type = GetTypeBtagging(sample_name)
        evt_start = options.evtstart
        evt_to_run = options.evtsperjob
        isFakeLep = options.fakelep
        tfile = ROOT.TFile(ifile)
        # Do reco
        print 'run options:',tfile.GetName(),sample_name,sample_type,evt_start,evt_to_run,isFakeLep
        reco_message = reconstruction(tfile,sample_name,sample_type,evt_start,evt_to_run,isFakeLep)   
        print reco_message

    print 'All done!'

    # Stop our timer
    timer.Stop()
    # Print out our timing information
    print '\n'
    rtime = timer.RealTime(); # Real time (or "wall time")
    ctime = timer.CpuTime(); # CPU time
    print("RealTime={0:6.2f} seconds, CpuTime={1:6.2f} seconds").format(rtime,ctime)



# Do reconstruction for the given file, process only a specific part of it, and output to a new file
# ex. tf1, 'Tbar_s','singletopbar' (for b-tagging correction purpose),2000,10000
def reconstruction(tfile,sample_name,sample_type,evt_start=0,evt_to_run=1000,isFakeLep='no'):
    # Get ttree
    tmptree = tfile.Get('selected')
    print 'Start to process file',tfile.GetName()
    if evt_to_run>0 :
        evt_end = evt_start+evt_to_run
    else : # (evt_to_run = -1 indicates run all events)
        evt_end = -1 
    print 'Prcess from event',evt_start,'to event',evt_end
    # Check if the starting event is more than total number of evts in current files
    if evt_start > tmptree.GetEntries():
        return 'evt_start > total_num_events in input ttree. Will stop here.'
        
    # Make output file
    fout_name = sample_name+'_reco_'+str(evt_start)+'_'+str(evt_end)+'.root'
    fout = ROOT.TFile(fout_name,'recreate')
    # Make output ttree
    tmptree.SetBranchStatus('jets*',0)
    tmptree.SetBranchStatus('lep*',0)
    tmptree.SetBranchStatus('lep_charge',1)
    if isFakeLep == 'yes':
        tmptree.SetBranchStatus('lep_iso',1)
    tmptree.SetBranchStatus('met*',0)
    tmptree.SetBranchStatus('pileup*',0)
    newtree = tmptree.CloneTree(0)
    tmptree.SetBranchStatus('*',1)

    # Add new branches to the output tree
    vecs = []
    br_names = []
    # Reconstruction results
    # reco momenta: t_lep, t_had, w_lep, w_had
    reco_pt = ROOT.vector('float')()
    reco_eta = ROOT.vector('float')()
    reco_phi = ROOT.vector('float')()
    reco_mass = ROOT.vector('float')()
    N_btag = ROOT.vector('int')()
    N_jets = ROOT.vector('int')()        
    vecs += [reco_pt,reco_eta,reco_phi,reco_mass,N_btag,N_jets]
    br_names += ['reco_pt','reco_eta','reco_phi','reco_mass','N_btag','N_jets']
    # kinfit results 
    kinfit_results = ROOT.vector('float')() #'pZv','scaleLep','scaleblep','scalebhad','scaleWsub1','scaleWsub2'
    final_chi2 = ROOT.vector('float')()
    final_nv_pz = ROOT.vector('float')()
    N_combos = ROOT.vector('int')() 
    N_combo_errs = ROOT.vector('int')()
    final_errflags = ROOT.vector('int')()
    vecs += [final_chi2,kinfit_results,final_nv_pz,N_combos,N_combo_errs,final_errflags]
    br_names += ['final_chi2','kinfit_results','final_nv_pz','N_combos','N_combo_errs','final_errflags']
    
    # Corrections
    if sample_type != 'data':
        # Load some SFs
        EleID_SFs = LoadEleSFs()
        pu_dists  = LoadPUfiles()
        btagEff_type = sample_type
        eff_hists = LoadBtagEfficiency(btagEff_type)
        # Book branches        
        w_PU_vec = ROOT.vector('float')()
        w_PU_up = ROOT.vector('float')()
        w_PU_down = ROOT.vector('float')()
        w_btag_vec = ROOT.vector('float')()
        w_btag_up = ROOT.vector('float')()
        w_btag_down = ROOT.vector('float')()
        w_eleID_vec = ROOT.vector('float')()
        w_eleID_up = ROOT.vector('float')()
        w_eleID_down = ROOT.vector('float')()
        w_trigger_vec = ROOT.vector('float')()
        w_trigger_up = ROOT.vector('float')()
        w_trigger_down = ROOT.vector('float')()
        vecs += [w_PU_vec,w_PU_up,w_PU_down,w_btag_vec,w_btag_up,w_btag_down,w_eleID_vec,w_eleID_up,w_eleID_down]
        vecs += [w_trigger_vec,w_trigger_up,w_trigger_down]
        br_names += ['w_PU','w_PU_up','w_PU_down','w_btag','w_btag_up','w_btag_down','w_eleID','w_eleID_up','w_eleID_down']
        br_names += ['w_trigger','w_trigger_up','w_trigger_down']
    # Add branches to the tree
    branches = zip(br_names,vecs)
    for ibr in branches:
        newtree.Branch(ibr[0],ibr[1])
    # Add cutflow diagrams
    h_cutflow = ROOT.TH1D('cutflow_extra',' cutflow extra;cuts;events',5,0.,5.)
    h_cutflow.SetBit(ROOT.TH1.kCanRebin)
    h_cutflow.SetDirectory(fout)


    # Loop over entries
    n_evt = 0
    for iev in range(evt_start,tmptree.GetEntries()):
        # Reset all vector containers
        for ivec in vecs: ivec.clear()    

        # Progress report
        if iev%500 == 0 : print 'processing event',iev 
        # Break at the given event number
        if iev == evt_end : print 'Finish processing. Quit at event',evt_end; break ;
        
        tmptree.GetEntry(iev)
        # Book the branches from input ttree
        jets_pt = tmptree.jets_pt
        jets_eta = tmptree.jets_eta
        jets_phi = tmptree.jets_phi
        jets_mass = tmptree.jets_mass
        jets_csv = tmptree.jets_csv 
        lep_pts = tmptree.lep_pt
        lep_etas = tmptree.lep_eta
        lep_phis = tmptree.lep_phi
        lep_mass = tmptree.lep_mass
        met_pt = tmptree.met_pt
        met_phi = tmptree.met_phi         

        #print 'iev',iev

        #######################################################
        #          Additional Selection Cuts                  #
        #######################################################
        h_cutflow.Fill('no cut',1)
        #trigger
        if options.applytrigger == 'yes':
            if tmptree.trigger[0] == 0 : continue  
        h_cutflow.Fill('trigger',1)      
        # jets
        bjets = []
        n_selected_jets = jets_pt.size()
        if not 4<= n_selected_jets <=5 : continue
        h_cutflow.Fill('n_jets',1)
        for ijet in jets_csv:
            if ijet>csvm : bjets.append(ijet)
        nbtags = len(bjets)
        if not nbtags >= options.nbcut : continue
        h_cutflow.Fill('n_btag',1)
        # leptons  
        if isFakeLep == 'yes':
            lep_isos = tmptree.lep_iso
            if not lep_isos.size() >= options.nlepcut : continue
            toskip = 0
            for ilep in lep_isos:
                if ilep < options.lepisocut : toskip = 1
            if toskip == 1 : continue 
        h_cutflow.Fill('lep_iso',1)

        #######################################################
        #          Calculate correction SFs for MC            #
        #######################################################
        if sample_type != 'data' :

            # Do electron corrections
            lep_pt_ = lep_pts[0]
            lep_eta_ = lep_etas[0]

            # Get Ele ID efficiency SF                
            sf_eleID = GetEleSFs(lep_pt_,lep_eta_,EleID_SFs)
            w_eleID = sf_eleID[0]
            err_eleID_up = sf_eleID[1]
            err_eleID_down = sf_eleID[2]
            w_eleID_vec.push_back(w_eleID)
            w_eleID_up.push_back(w_eleID + err_eleID_up)
            w_eleID_down.push_back(w_eleID - err_eleID_down)

            # Get HLT efficiency SF
            sf_trigger = GetTriggerSFs(lep_pt_,lep_eta_)
            w_trigger = sf_trigger[0]
            err_trigger_up = sf_trigger[1]
            err_trigger_down = sf_trigger[2]
            w_trigger_vec.push_back(w_trigger)
            w_trigger_up.push_back(w_trigger + err_trigger_up)
            w_trigger_down.push_back(w_trigger - err_trigger_down)

            # Pileup
            npvRealTrue = tmptree.mc_pileup_events[0]
            w_PU = GetPUWeights(npvRealTrue,pu_dists)
            err_PU_up = 0
            err_PU_down = 0
            w_PU_vec.push_back(w_PU)
            w_PU_up.push_back(w_PU+err_PU_up)
            w_PU_down.push_back(w_PU+err_PU_down)   

            # b-tagging efficiency
            jets_flavor = tmptree.jets_flavor
            selected_jets = []
            for i in range(jets_pt.size()):
                selected_jets.append((jets_pt[i],jets_eta[i],jets_flavor[i],jets_csv[i]))
            # Get b-tag weights
            w_result = get_weight_btag(selected_jets,eff_hists)
            w_btag   = w_result[0]
            err_btag = w_result[1]
            w_btag_vec.push_back(w_btag)
            w_btag_up.push_back(w_btag + err_btag)
            w_btag_down.push_back(w_btag - err_btag)

        #######################################################
        #          Do reconstruction                          #
        #######################################################
        # Set up inputs  
        lep_p4 = ROOT.TLorentzVector()
        lep_p4.SetPtEtaPhiM(lep_pts[0],lep_etas[0],lep_phis[0],lep_mass[0])
        jets_p4,jets_csv_list = [],[]    
        for i in range(jets_pt.size()):
            tmp_p4 = ROOT.TLorentzVector()
            tmp_p4.SetPtEtaPhiM(jets_pt[i],jets_eta[i],jets_phi[i],jets_mass[i])
            jets_p4.append(tmp_p4)
            jets_csv_list.append(jets_csv[i])
        metPt = met_pt[0]; metPhi = met_phi[0]
        # Do reco
        lep_type = 'el'
        if sample_type == 'data': 
            mcordata = 'data'
        else : 
            mcordata = 'mc'
        reco_result = DoReco(jets_p4,jets_csv_list,lep_p4,metPt,metPhi,lep_type,mcordata)
        if reco_result == 'none':
            print 'No valid reco done! Will skip this event #',iev
            continue
        # Add reco quantity to the tree
        N_btag.push_back(nbtags)
        N_jets.push_back(n_selected_jets)
        plot_final_chi = reco_result[0]
        reco_p4s = reco_result[1]
        bestFitParValues = reco_result[2]
        N_combos.push_back(reco_result[3])
        N_combo_errs.push_back(reco_result[4])
        final_errflags.push_back(reco_result[5])

        final_chi2.push_back(plot_final_chi)
        for ip4 in reco_p4s:
            reco_pt.push_back(ip4.Pt())
            reco_eta.push_back(ip4.Eta())
            reco_phi.push_back(ip4.Phi())
            reco_mass.push_back(ip4.M())
        for i in range(1,len(bestFitParValues)):
            kinfit_results.push_back(bestFitParValues[i])
        final_nv_pz.push_back(bestFitParValues[0])
        #Fill the newtree
        newtree.Fill()
    # End loop over entries

    h_cutflow.SetMinimum(0)
    # Write and close files
    fout.Write()
    fout.Close()
    tfile.Close()

    return 'Reconstruction finished!'

main()
