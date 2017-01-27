# Take in the selection output ttree, do reco reconstruction via kinfit

from Legacy_Afb.Tools.kinfit import *
from Legacy_Afb.Tools.ttbar_utility import *
import Legacy_Afb.Tools.lepHelper as lepHelper
import Legacy_Afb.Tools.jetHelper as jetHelper 
import glob
from optparse import OptionParser
from array import array
import time

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
                  default = 0.2,
                  dest='lepisocut',
                  help='Lower bound for fake electron isolation.')

parser.add_option('--nbcut', metavar='F', type='int', action='store',
                  default = 2,
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

parser.add_option('--lep_type', metavar='F', type='string', action='store',
                  default = 'ele',
                  dest='lep_type',
                  help='e+jets or mu+jets')

parser.add_option('--PDF_name', metavar='F', type='string', action='store',
                  default = 'ct10',
                  dest='PDF_name',
                  help='PDF used in the MC sample. ct10/cteq/gjr')

(options, args) = parser.parse_args()

argv = []

# Some preset constants
data_lumi = 19748 
csvm = 0.679

PDF_branch = {'ct10':'weight_pdf_ct10','cteq':'weight_pdf_cteq','gjr':'weight_pdf_gjr'}

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
    #tmptree.SetBranchStatus('jets*',0)
    tmptree.SetBranchStatus('lep*',0)
    tmptree.SetBranchStatus('lep_charge',1)
    if isFakeLep == 'yes':
        tmptree.SetBranchStatus('lep_iso',1)
    tmptree.SetBranchStatus('met*',0)
    # tmptree.SetBranchStatus('pileup*',0)
    tmptree.SetBranchStatus('weight_*',0)
    # keep weight_gen if exist
    for branch in tmptree.GetListOfBranches():
        if branch.GetName() == "weight_gen":
            tmptree.SetBranchStatus('weight_gen',1)
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
    # for chi2 study
    total_chi2 =  ROOT.vector('float')()
    chi2_mass = ROOT.vector('float')()
    chi2_scale = ROOT.vector('float')()
    chi2_csv = ROOT.vector('float')()
    chi2_mass_terms = ROOT.vector('float')()
    vecs+=[total_chi2,chi2_mass,chi2_scale,chi2_csv,chi2_mass_terms]
    br_names+=['total_chi2','chi2_mass','chi2_scale','chi2_csv','chi2_mass_terms']

    # Corrections
    if sample_type != 'data':
        # Load some SFs
        EleID_SFs = LoadEleSFs()
        pu_dists  = LoadPUfiles()
        btagEff_type = sample_type
        eff_hists = LoadBtagEfficiency(btagEff_type)
        if options.lep_type == 'mu':
            muon_helper = lepHelper.muon_helper()
        # Book branches        
        w_PU_vec = array('f',[1.])
        w_PU_up = array('f',[1.])
        w_PU_down = array('f',[1.])

        w_btag_vec = array('f',[1.])
        w_btag_up = array('f',[1.])
        w_btag_down = array('f',[1.])

        w_lepID_vec = array('f',[1.])
        w_lepID_up = array('f',[1.])
        w_lepID_down = array('f',[1.])

        w_trigger_vec = array('f',[1.])
        w_trigger_up = array('f',[1.])
        w_trigger_down = array('f',[1.])

        w_tracking = array('f',[1.])
        w_tracking_up = array('f',[1.])
        w_tracking_down = array('f',[1.])

        w_lepIso = array('f',[1.])
        w_lepIso_up = array('f',[1.])
        w_lepIso_down = array('f',[1.])

        w_PDF = array('f',[1.])
        w_PDF_up = array('f',[1.])
        w_PDF_down = array('f',[1.])

        br_defs = [('w_PU',w_PU_vec,'w_PU/F')]
        br_defs += [('w_PU_up',w_PU_up,'w_PU_up/F')]
        br_defs += [('w_PU_down',w_PU_down,'w_PU_down/F')]

        br_defs += [('w_PDF',w_PDF,'w_PDF/F')]
        br_defs += [('w_PDF_up',w_PDF_up,'w_PDF_up/F')]
        br_defs += [('w_PDF_down',w_PDF_down,'w_PDF_down/F')]

        br_defs += [('w_btag',w_btag_vec,'w_btag/F')]
        br_defs += [('w_btag_up',w_btag_up,'w_btag_up/F')]
        br_defs += [('w_btag_down',w_btag_down,'w_btag_down/F')]

        br_defs += [('w_lepID',w_lepID_vec,'w_lepID/F')]
        br_defs += [('w_lepID_up',w_lepID_up,'w_lepID_up/F')]
        br_defs += [('w_lepID_down',w_lepID_down,'w_lepID_down/F')]

        br_defs += [('w_trigger',w_trigger_vec,'w_trigger/F')]
        br_defs += [('w_trigger_up',w_trigger_up,'w_trigger_up/F')]
        br_defs += [('w_trigger_down',w_trigger_down,'w_trigger_down/F')]

        br_defs += [('w_tracking',w_tracking,'w_tracking/F')]
        br_defs += [('w_tracking_up',w_tracking_up,'w_tracking_up/F')]
        br_defs += [('w_tracking_down',w_tracking_down,'w_tracking_down/F')]

        br_defs += [('w_lepIso',w_lepIso,'w_lepIso/F')]
        br_defs += [('w_lepIso_up',w_lepIso_up,'w_lepIso_up/F')]
        br_defs += [('w_lepIso_down',w_lepIso_down,'w_lepIso_down/F')]
    
    # Add vec branches to the tree
    branches = zip(br_names,vecs)
    for ibr in branches:
        newtree.Branch(ibr[0],ibr[1])
    # Add float branch into ttree
    for ibr in br_defs:
        newtree.Branch(ibr[0],ibr[1],ibr[2])


    # Add cutflow diagrams
    h_cutflow = ROOT.TH1D('cutflow_extra',' cutflow extra;cuts;events',5,0.,5.)
    h_cutflow.SetBit(ROOT.TH1.kCanRebin)
    h_cutflow.SetDirectory(fout)

    # arguments for DoReco
    if options.lep_type in ['ele','el']:
        lep_type = 'el'
    else:
        lep_type = 'mu'
    print '(info) Run %s samples!'%lep_type
    if sample_type == 'data': 
        mcordata = 'data'
    else : 
        mcordata = 'mc'

    # check if PDF weight branch in input ttree. also decide which PDF to use
    has_PDF = False
    pdf_key = 'ct10'
    if PDF_branch.get(options.PDF_name,0):
        tmp_key = PDF_branch[options.PDF_name]
        if tmptree.FindBranch(tmp_key):
            print '(info) Found PDF_%s. Will get pdf weights.'%tmp_key
            has_PDF = True
            pdf_key = tmp_key

    # Loop over entries
    n_evt = 0
    for iev in range(evt_start,tmptree.GetEntries()):
        start_time = time.clock()

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
            if not tmptree.trigger[0] : continue  
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
            if not lep_isos.size() == options.nlepcut : continue
            toskip = 0
            for ilep in lep_isos:
                if ilep < options.lepisocut : toskip = 1
            if toskip == 1 : continue 
        h_cutflow.Fill('lep_iso',1)

        #######################################################
        #          Calculate correction SFs for MC            #
        #######################################################
        if sample_type != 'data' :

           # Pileup
            npvRealTrue = tmptree.mc_pileup_events[0]
            w_PU = GetPUWeights(npvRealTrue,pu_dists)
            err_PU_up = 0
            err_PU_down = 0
            w_PU_vec[0] = w_PU
            w_PU_up[0] =  w_PU+err_PU_up
            w_PU_down[0] = w_PU+err_PU_down

            # Do lepton corrections
            lep_pt_ = lep_pts[0]
            lep_eta_ = lep_etas[0]

            if lep_type == 'el':
                # Get Ele ID efficiency SF                
                sf_eleID = GetEleSFs(lep_pt_,lep_eta_,EleID_SFs)
                w_eleID = sf_eleID[0]
                err_eleID_up = sf_eleID[1]
                err_eleID_down = sf_eleID[2]

                w_lepID_vec[0] = w_eleID
                w_lepID_up[0] = w_eleID + err_eleID_up
                w_lepID_down[0] = w_eleID - err_eleID_down

                # Get HLT efficiency SF
                sf_trigger = GetTriggerSFs(lep_pt_,lep_eta_)
                w_trigger = sf_trigger[0]
                err_trigger_up = sf_trigger[1]
                err_trigger_down = sf_trigger[2]

                w_trigger_vec[0] = w_trigger
                w_trigger_up[0] = w_trigger + err_trigger_up
                w_trigger_down[0] = w_trigger - err_trigger_down
            else:
                muon_helper.getSF(lep_eta_,lep_pt_,npvRealTrue)

                w_trigger_vec[0] = muon_helper.weight_trig_eff
                w_trigger_up[0] = muon_helper.weight_trig_eff_hi
                w_trigger_down[0] = muon_helper.weight_trig_eff_low

                w_tracking[0] = muon_helper.weight_tracking
                w_tracking_up[0] = muon_helper.weight_tracking_hi
                w_tracking_down[0] = muon_helper.weight_tracking_low

                w_lepID_vec[0] = muon_helper.weight_lep_ID
                w_lepID_up[0] = muon_helper.weight_lep_ID_hi
                w_lepID_down[0] = muon_helper.weight_lep_ID_low

                w_lepIso[0] = muon_helper.weight_lep_iso
                w_lepIso_up[0] = muon_helper.weight_lep_iso_hi
                w_lepIso_down[0] = muon_helper.weight_lep_iso_low


            # b-tagging efficiency
            jets_flavor = tmptree.jets_flavor
            selected_jets = []
            for i in range(jets_pt.size()):
                selected_jets.append((jets_pt[i],jets_eta[i],jets_flavor[i],jets_csv[i]))
            # Get b-tag weights
            w_result = get_weight_btag(selected_jets,eff_hists)
            w_btag   = w_result[0]
            err_btag = w_result[1]
            w_btag_vec[0] =  w_btag 
            w_btag_up[0] =  w_btag + err_btag 
            w_btag_down[0] =  w_btag - err_btag

            # PDF weights, if exist the branch
            if has_PDF:
                pdf_vec = getattr(tmptree,pdf_key)
                #   return (pdf.max(),pdf.min())
                pdf_w = jetHelper.get_PDF_SF(pdf_vec)
                w_PDF_up[0] = pdf_w[0]
                w_PDF_down[0] = pdf_w[1]

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
        # study of chi2
        all_chis1 = reco_result[6]
        all_chis2 = reco_result[7]
        total_chi2.push_back(all_chis1[0])
        chi2_mass.push_back(all_chis1[1])
        chi2_scale.push_back(all_chis1[2])
        chi2_csv.push_back(all_chis1[3])
        if len(all_chis2)>0:
            total_chi2.push_back(all_chis2[0])
            chi2_mass.push_back(all_chis2[1])
            chi2_scale.push_back(all_chis2[2])
            chi2_csv.push_back(all_chis2[3])    
        # chi2 mass constraint terms: thad,tlep,whad,wlep
        for i in [4,5,6,7]:
            chi2_mass_terms.push_back(all_chis1[i])

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
