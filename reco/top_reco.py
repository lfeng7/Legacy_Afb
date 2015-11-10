# Take in the selection output ttree, do reco reconstruction via kinfit

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
                  default = -1,
                  dest='maxfiles',
                  help='max number of input ntuple files')

parser.add_option('--upload', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='dumpplots',
                  help='if you want to dump plots to MY webpage! Unless you modify fwlite_boilerplate >_<')

parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='If you want more information than you usually need.')

parser.add_option('--recopt', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='recopt',
                  help='If do reco pt reweighting')

parser.add_option('--tmptype', metavar='F', type='string', action='store',
                  default = 'scaled',
                  dest='tmptype',
                  help='If use correction')

parser.add_option('--applytrigger', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='applytrigger',
                  help='If apply trigger on MC')

parser.add_option('--fakelep', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='fakelep',
                  help='If run on selected events with fake lepton.')

parser.add_option('--lepisocut', metavar='F', type='int', action='store',
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

# # Set up output files
# template_type = options.tmptype
# tmptype_name = template_type
# if options.fakelep == 'yes':
#     tmptype_name += '_fakelep'
# hist_prepend = './selected_hists/'+tmptype_name+'/' 
    
# dir_name = 'controlplots_'+tmptype_name # name of the dir for control plots, such as corrected, recoPT, un_corrected
  

#### Set up MC input files

#    0,                 1         2          3         4                   5
# (MC_sample_name, sample_type, Nevts_gen, x-sec, nevts_total_ntuple, btag_type)
# Single reco
flist.append(['T_s','singlereco',259961,3.79,259176,'singlereco'] )
flist.append(['T_t','singlereco',3758227,56.4,3748155,'singlereco'] )
flist.append(['T_tW','singlereco',497658,11.1,495559,'singlereco'])
flist.append(['Tbar_s','singlereco',139974, 1.76,139604,'singlerecobar'])
flist.append(['Tbar_t','singlereco',1935072, 30.7,1930185,'singlerecobar'])
flist.append(['Tbar_tW','singlereco',493460,11.1,491463,'singlerecobar'])
# Wjets
flist.append(['W1JetsToLNu_TuneZ2Star_8TeV','wjets',23141598,6662.8,23038253,'wjets'])
flist.append(['W2JetsToLNu_TuneZ2Star_8TeV','wjets',34044921,2159.2,33993463,'wjets'])
flist.append(['W3JetsToLNu_TuneZ2Star_8TeV','wjets',15539503,640.4,15507852,'wjets'])
flist.append(['W4JetsToLNu_TuneZ2Star_8TeV','wjets',13382803,246.0,13326400,'wjets'])
# DYjets
flist.append(['DY1JetsToLL_M','zjets',24045248,660.6,23802736,'zjets'])
flist.append(['DY2JetsToLL_M','zjets',2352304,215.1,2345857,'zjets'])
flist.append(['DY3JetsToLL_M','zjets',11015445,65.79,10655325,'zjets'])
flist.append(['DY4JetsToLL_M','zjets',6402827,28.59,5843425,'zjets'])
# QCD
flist.append(['QCD_Pt-15to3000','qcd',9991674,6662.6,9940092,'qcd'])
# signal
flist.append(['TT_CT10_TuneZ2star_8TeV','ttbar',21675970,245.9,21560109,'ttbar'])

#### Set up data input files
#    0,         1                   2             
# (filepath,   type,   sample_integrated_lumi
#datafile = ['SingleEl_Run2012A_v1','data',888]
datafile = ['SingleEl_Run2012ABCD','data',19748]

def main():
    if options.makehists == 'yes':
        print 'Making histograms from the selection root files.'
        MakeHistograms()
    if options.makeplots == 'yes':
        print 'Making comparison plots of MC and data.'
        MakeComparisonPlots()
    print 'All done!'

# Do reconstruction for the given file, process only a specific part of it, and output to a new file
# ex. tf1, 'Tbar_s','singlerecobar' (for b-tagging correction purpose),2000,10000
def reco(tfile,sample_name,sample_type,evt_start=0,evt_to_run=1000,isFakeLep='no'):
    # Get ttree
    print 'Start to process file',tfile.GetName()
    evt_end = evt_start+evt_to_run
    print 'Prcess from event',evt_start,'to event',evt_end
    tmptree = tfile.Get('selected')
    # Make output file
    fout_name = sample_name+'_reco_'+evt_start+'_'+evt_end+'.root'
    fout = ROOT.TFile(fout_name,'recreate')
    # Make output ttree
    tmptree.SetBranchStatus('jets*',0)
    tmptree.SetBranchStatus('lep*',0)
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
    vecs += [final_chi2,kinfit_results]
    br_names += ['final_chi2','kinfit_results']
    
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
        vecs += [w_PU_vec,w_PU_up,w_PU_down,w_btag_vec,w_btag_up,w_btag_down,w_eleID,w_eleID_up,w_eleID_down]
        vecs += [w_trigger,w_trigger_up,w_trigger_down]
        br_names += ['w_PU','w_PU_up','w_PU_down','w_btag','w_btag_up','w_btag_down','w_eleID_vec','w_eleID_up','w_eleID_down']
        br_names += ['w_trigger_vec','w_trigger_up','w_trigger_down']
    # Add branches to the tree
    branches = zip(br_names,vecs)
    for ibr in branches:
        newtree.Branch(ibr[0],ibr[1])
    # Add cutflow diagrams
    h_cutflow = ROOT.TH1D('cutflow_extra',event_type+' cutflow extra;cuts;events',4,0.,4.)
    h_cutflow.SetBit(ROOT.TH1.kCanRebin)


    # Loop over entries
    n_evt = 0
    for iev in range(evt_start,oldtree.GetEntries()):
        # Progress report
        if iev%5000 == 0 : print 'processing event',i 
        # Break at the given event number
        if iev == evt_end : print 'Finish processing. Quit at event',evt_end; break ;

        # Book the branches from input ttree
        jets_pt = tmptree.jets_pt
        jets_eta = tmptree.jets_eta
        jets_phi = tmptree.jets_phi
        jets_mass = tmptree.jets_mass
        jets_csv = tmptree.jets_csv 
        lep_pts = tmptree.lep_pt
        lep_etas = tmptree.lep_eta
        lep_phis = tmptree.lep_phis
        lep_mass = tmptree.lep_mass
        met_pt = tmptree.met_pt
        met_phi = tmptree.met_phi         


        #######################################################
        #          Additional Selection Cuts                  #
        #######################################################
        h_cutflow.Fill('no cut',1)
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
        if isFakeLep == 'yes'
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
        lep_p4.SetPtEtaPhiM(lep_pt[0],lep_eta[0],lep_phi[0],lep_mass[0])
        jets_p4,jets_csv_list = [],[]    
        for i in range(jets_pt.size()):
            tmp_p4 = ROOT.TLorentzVector()
            tmp_p4.SetPtEtaPhiM(jets_pt[i],jets_eta[i],jets_phi[i],jets_phi[i],jets_mass[i])
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
        # Add reco quantity to the tree
        N_btag.push_back(nbtags)
        N_jets.push_back(n_selected_jets)
        plot_final_chi = reco_result[0]
        reco_p4s = reco_result[1]
        bestFitParValues = reco_result[2]
        final_chi2.push_back(plot_final_chi)
        for ip4 in reco_p4s:
            reco_pt.push_back(ip4.Pt())
            reco_eta.push_back(ip4.Eta())
            reco_phi.push_back(ip4.Phi())
            reco_mass.push_back(ip4.M())
        for ipar in bestFitParValues:
            kinfit_results.push_back(ipar)

        #Fill the newtree
        newtree.Fill()
    # End loop over entries

    # Write and close files
    fout.Write()
    fout.Close()
    tfile.Close()








main()