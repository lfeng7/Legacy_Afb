# Take in the angles ttree from reco, and make input ttree for original template making and fitting codes

from Legacy_Afb.Tools.ttbar_utility import *
from Legacy_Afb.Tools.angles_tools import *
import glob
from optparse import OptionParser
import ROOT
from array import array
#from angles_tools import *
import math
ROOT.gROOT.SetBatch(True)

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

parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='If you want more information than you usually need.')

parser.add_option('--slim', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='slim',
                  help='If you want slimmed ttree with angles and stuff')

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

parser.add_option('--ttbar_type', metavar='F', type='string', action='store',
                  default = 'none',
                  dest='ttbar_type',
                  help='type of ttbar templates. gg/qq/bkg')

(options, args) = parser.parse_args()

argv = []

# Some preset constants
data_lumi = 19748 
csvm = 0.679 
deltaR_matching = 0.4

def main():
    # Get the file list with all input files.
    if options.inputFiles != '':
        allfiles = glob.glob( options.inputFiles )
    else:
        allfiles = []

    timer = ROOT.TStopwatch()
    timer.Start()

    # Job splitting is done here
    for ifile in allfiles:
        sample_name = ifile.split('/')
        sample_name = sample_name[len(sample_name)-1].split('.root')[0]
        evt_start = options.evtstart
        evt_to_run = options.evtsperjob
        tfile = ROOT.TFile(ifile)
        # Do make angles
        print 'Do makeAngles now!'
        print 'run options:',tfile.GetName(),sample_name,evt_start,evt_to_run
        reco_message = makeTemps(tfile,sample_name,evt_start,evt_to_run)   
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
def makeTemps(tfile,sample_name,evt_start=0,evt_to_run=1000):
    if 'Single' in sample_name or 'Run' in sample_name:
        newtree_name ='angles_data'
        print 'Running on data!'
    else:
        newtree_name = 'angles'
        print 'Running on MC!'
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
    if options.ttbar_type == 'none':
        fout_name = sample_name+'_template_'+str(evt_start)+'.root'
    else :
        fout_name = sample_name+'_template_'+options.ttbar_type+'_'+str(evt_start)+'.root'
    fout = ROOT.TFile(fout_name,'recreate')
    # Make output ttree
    newtree = ROOT.TTree(newtree_name,newtree_name)

    # Add new branches to the output tree
    br_defs = []
    # vars
    ttbar_mass = array('f',[0.])
    Qt = array('f',[0.])
    cos_theta_cs = array('f',[0.])
    Feynman_x = array('f',[0.])
    Q_l = array('i',[0])
    cos_theta_mc = array('f',[-10])
    Feynman_x_mc = array('f',[-10])
    ttbar_mass_mc = array('f',[-10])
    lnL = array('f',[0.])
    n_valid_jets = array('i',[0])
    n_bTags = array('i',[0])
    fitParValues = array('f',6*[1.])
    w_a = array('f',[1.])
    w_a_opp = array('f',[1.])
    w_s_xi = array('f',[1.])
    w_s_xi_opp = array('f',[1.])
    w_a_xi = array('f',[1.])
    w_a_xi_opp = array('f',[1.])
    w_s_delta = array('f',[1.])
    w_s_delta_opp = array('f',[1.])
    w_a_delta = array('f',[1.])
    w_a_delta_opp = array('f',[1.])
    motherPIDs = array('i',2*[-100])
    Pdf_weights = array('d',100*[1.])
    pileup_reweight = array('f',[1.])
    top_pT_reweight = array('f',[1.])
    GJR_reweight = array('f',[1.])
    CT10_reweight = array('f',[1.])
    cteq_reweight = array('f',[1.])
    btag_eff_reweight = array('f',[1.])
    tracking_reweight = array('f',[1.])
    lepID_reweight = array('f',[1.])
    lepIso_reweight = array('f',[1.])
    trigger_reweight = array('f',[1.])
    btag_eff_reweight_low = array('f',[1.])
    tracking_reweight_low = array('f',[1.])
    lepID_reweight_low = array('f',[1.])
    lepIso_reweight_low = array('f',[1.])
    trigger_reweight_low = array('f',[1.])
    btag_eff_reweight_hi = array('f',[1.])
    tracking_reweight_hi = array('f',[1.])
    lepID_reweight_hi = array('f',[1.])
    lepIso_reweight_hi = array('f',[1.])
    trigger_reweight_hi = array('f',[1.])
    pileup = array('f',[0.])
    pileup_real = array('f',[0.])

    # extra stuff not used by template fitter
    correction_weight = array('f',[1.])
    reco_pt   = array('f',4*[0.])
    reco_eta  = array('f',4*[0.])
    reco_phi  = array('f',4*[0.])
    reco_mass = array('f',4*[0.])
    leadingJet_pt = array('f',[0.])
    leadingJet_mass = array('f',[0.])
    isMatched = array('i',4*[-1])


    # names 
    br_defs += [('ttbar_mass',ttbar_mass,'ttbar_mass/F')]
    br_defs += [('Qt',Qt,'Qt/F')]
    br_defs += [('cos_theta_cs',cos_theta_cs,'cos_theta_cs/F')]
    br_defs += [('Feynman_x',Feynman_x,'Feynman_x/F')]
    br_defs += [('Q_l',Q_l,'Q_l/I')]
    br_defs += [('cos_theta_mc',cos_theta_mc,'cos_theta_mc/F')]
    br_defs += [('Feynman_x_mc',Feynman_x_mc,'Feynman_x_mc/F')]
    br_defs += [('ttbar_mass_mc',ttbar_mass_mc,'ttbar_mass_mc/F')]
    br_defs += [('lnL',lnL,'lnL/F')]
    br_defs += [('n_valid_jets',n_valid_jets,'n_valid_jets/I')]
    br_defs += [('n_bTags',n_bTags,'n_bTags/I')]
    br_defs += [('fitParValues',fitParValues,'fitParValues/F')]
    br_defs += [('w_a',w_a,'w_a/F')]
    br_defs += [('w_a_opp',w_a_opp,'w_a_opp/F')]
    br_defs += [('w_s_xi',w_s_xi,'w_s_xi/F')]
    br_defs += [('w_s_xi_opp',w_s_xi_opp,'w_s_xi_opp/F')]
    br_defs += [('w_a_xi',w_a_xi,'w_a_xi/F')]
    br_defs += [('w_a_xi_opp',w_a_xi_opp,'w_a_xi_opp/F')]
    br_defs += [('w_s_delta',w_s_delta,'w_s_delta/F')]
    br_defs += [('w_s_delta_opp',w_s_delta_opp,'w_s_delta_opp/F')]
    br_defs += [('w_a_delta',w_a_delta,'w_a_delta/F')]
    br_defs += [('w_a_delta_opp',w_a_delta_opp,'w_a_delta_opp/F')]
    br_defs += [('motherPIDs',motherPIDs,'motherPIDs[2]/I')]
    br_defs += [("Pdf_weights",Pdf_weights,"Pdf_weights[100]/D")]
    br_defs += [('pileup_reweight',pileup_reweight,'pileup_reweight/F')]
    br_defs += [('top_pT_reweight',top_pT_reweight,'top_pT_reweight/F')]
    br_defs += [('GJR_reweight',GJR_reweight,'GJR_reweight/F')]
    br_defs += [('CT10_reweight',CT10_reweight,'CT10_reweight/F')]
    br_defs += [('cteq_reweight',cteq_reweight,'cteq_reweight/F')]
    br_defs += [('btag_eff_reweight',btag_eff_reweight,'btag_eff_reweight/F')]
    br_defs += [('tracking_reweight',tracking_reweight,'tracking_reweight/F')]
    br_defs += [('lepID_reweight',lepID_reweight,'lepID_reweight/F')]
    br_defs += [('lepIso_reweight',lepIso_reweight,'lepIso_reweight/F')]
    br_defs += [('trigger_reweight',trigger_reweight,'trigger_reweight/F')]
    br_defs += [('btag_eff_reweight_low',btag_eff_reweight_low,'btag_eff_reweight_low/F')]
    br_defs += [('tracking_reweight_low',tracking_reweight_low,'tracking_reweight_low/F')]
    br_defs += [('lepID_reweight_low',lepID_reweight_low,'lepID_reweight_low/F')]
    br_defs += [('lepIso_reweight_low',lepIso_reweight_low,'lepIso_reweight_low/F')]
    br_defs += [('trigger_reweight_low',trigger_reweight_low,'trigger_reweight_low/F')]
    br_defs += [('btag_eff_reweight_hi',btag_eff_reweight_hi,'btag_eff_reweight_hi/F')]
    br_defs += [('tracking_reweight_hi',tracking_reweight_hi,'tracking_reweight_hi/F')]
    br_defs += [('lepID_reweight_hi',lepID_reweight_hi,'lepID_reweight_hi/F')]
    br_defs += [('lepIso_reweight_hi',lepIso_reweight_hi,'lepIso_reweight_hi/F')]
    br_defs += [('trigger_reweight_hi',trigger_reweight_hi,'trigger_reweight_hi/F')]
    br_defs += [('pileup',pileup,'pileup/F')]
    br_defs += [('pileup_real',pileup_real,'pileup_real/F')]

    # extra stuff not used by template fitter
    br_defs += [('correction_weight',correction_weight,'correction_weight/F')]
    br_defs += [('reco_pt',reco_pt,'reco_pt[4]/F')]
    br_defs += [('reco_eta',reco_eta,'reco_eta[4]/F')]
    br_defs += [('reco_phi',reco_phi,'reco_phi[4]/F')]
    br_defs += [('reco_mass',reco_mass,'reco_mass[4]/F')]
    br_defs += [('isMatched',isMatched,'isMatched[4]/I')]
    br_defs += [('leadingJet_pt',leadingJet_pt,'leadingJet_pt/F')]
    br_defs += [('leadingJet_mass',leadingJet_mass,'leadingJet_mass/F')]


    # study beta
    beta_v0 = array('f',[-1.])
    br_defs += [('beta_v0',beta_v0,'beta_v0/F')]
    charge_ratio = array('i',[0])
    br_defs += [('charge_ratio',charge_ratio,'charge_ratio/I')]

    # add alternative chi2 if exist in tree
    has_chi2_new = 0
    if tmptree.FindBranch('final_chi2_new'):
        has_chi2_new = 1
        chi2_new = array('f',[0.])
        br_defs += [('chi2_new',chi2_new,'chi2_new/F')]

    # Add branches to the tree
    for ibr in br_defs:
        newtree.Branch(ibr[0],ibr[1],ibr[2])

    # Add cutflow diagrams
    h_cutflow = ROOT.TH1D('cutflow_extra',' cutflow extra;cuts;events',5,0.,5.)
    # h_cutflow.SetBit(ROOT.TH1.kCanRebin)
    h_cutflow.SetDirectory(fout)
    # Charge ratio hist
    h_charge_ratio = ROOT.TH1D('charge_ratio',' charge_ratio;;events',5,0.,5.)
    # h_charge_ratio.SetBit(ROOT.TH1.kCanRebin)
    h_charge_ratio.SetDirectory(fout)

    # Find out if the sample is TTbar MC
    ttbar_names = ['TT','signal']
    is_ttbar = 0
    for item in ttbar_names:
        if item in sample_name:
            is_ttbar = 1
    if is_ttbar == 1:
        print 'Is ttbar MC sample!'

    # a special normalization for top pT reweighting only
    if tmptree.FindBranch('weight_top_pT'):
        tmp_h0 = ROOT.TH1F('tmp_h0','tmp_h0',10,3,6)
        tmp_h1 = ROOT.TH1F('tmp_h1','tmp_h1',10,3,6)
        tmptree.Draw('N_jets>>tmp_h0')
        tmptree.Draw('N_jets>>tmp_h1','weight_top_pT')
        toppt_scale = tmp_h1.Integral()/tmp_h0.Integral()
        # new_w_top_pt = w_top_pt/toppt_scale
        print 'toppt_scale',toppt_scale
        tmp_h0.SetDirectory(0)
        tmp_h1.SetDirectory(0)

    # Check if it is a TTbar MC sample, an MC sample , or a data sample
    if tmptree.FindBranch('gen_pt'):
        is_TT_MC = 1
    else:
        is_TT_MC = 0
    if tmptree.FindBranch('w_PU') :
        is_MC = 1
    else:
        is_MC = 0

    # Loop over entries
    n_evt = 0
    for iev in range(evt_start,tmptree.GetEntries()):
        # Progress report
        if iev%5000 == 0 : print 'processing event',iev 
        # Break at the given event number
        if iev == evt_end : print 'Finish processing. Quit at event',evt_end; break ;
        
        tmptree.GetEntry(iev)

        h_cutflow.Fill('no cut',1)
        # Skip events that kinfit did not converge ( error flag == 4 )
        if not tmptree.final_errflags[0]==0 : continue
        h_cutflow.Fill('kinfit error',1)

        # separate TTbar into qq,gg,other templates
        # qq/gg are qq/gg->ttbar->ljets events
        # tt_other are anythong->ttbar->non ljets events        
        if options.ttbar_type=='bkg' and tmptree.gen_type[0] == 'e_jets' : continue
        if options.ttbar_type=='qq' and not (tmptree.gen_type[0]=='e_jets' and tmptree.init_type[0]=='qqbar'):continue
        if options.ttbar_type=='gg' and not (tmptree.gen_type[0]=='e_jets' and tmptree.init_type[0]!='qqbar'):continue

        ttbar_mass[0] = tmptree.mtt[0]
        # reco_p4 is a list of tlep_p4,thad_p4,wlep_p4,whad_p4
        tlep_pt = tmptree.reco_pt[0]
        tlep_eta = tmptree.reco_eta[0]
        tlep_phi = tmptree.reco_phi[0]
        tlep_mass = tmptree.reco_mass[0]
        thad_pt = tmptree.reco_pt[1]
        thad_eta = tmptree.reco_eta[1]
        thad_phi = tmptree.reco_phi[1]
        thad_mass = tmptree.reco_mass[1]        
        lep_charge = tmptree.lep_charge[0]
        # Make 4vec of tlep and thad
        tlep_p4 = ROOT.TLorentzVector()
        thad_p4 = ROOT.TLorentzVector()
        tlep_p4.SetPtEtaPhiM(tlep_pt,tlep_eta,tlep_phi,tlep_mass)
        thad_p4.SetPtEtaPhiM(thad_pt,thad_eta,thad_phi,thad_mass) 
        Q_Data = tlep_p4+thad_p4
        Qt[0] = math.sqrt(Q_Data.Px()*Q_Data.Px()+Q_Data.Py()*Q_Data.Py())
        cos_theta_cs[0] = tmptree.cos_theta[0]
        Feynman_x[0] = tmptree.xf[0]
        Q_l[0] = tmptree.lep_charge[0]
        cos_theta_mc[0] = -10
        if tmptree.FindBranch('cos_theta_mc'):
            if tmptree.cos_theta_mc.size()>0 : cos_theta_mc[0] = tmptree.cos_theta_mc[0]
        Feynman_x_mc[0] = -10
        if tmptree.FindBranch('xf_mc'):
            if tmptree.xf_mc.size()>0 : Feynman_x_mc[0] = tmptree.xf_mc[0]
        ttbar_mass_mc[0] = -10
        if tmptree.FindBranch('mtt_mc'):
            if tmptree.mtt_mc.size()>0 : ttbar_mass_mc[0] = tmptree.mtt_mc[0] 
        lnL[0] = tmptree.final_chi2[0]

        if has_chi2_new :
            chi2_new[0] = tmptree.final_chi2_new

        n_valid_jets[0] = tmptree.N_jets[0]
        n_bTags[0] = tmptree.N_btag[0]
        fitParValues[0] = tmptree.final_nv_pz[0]
        for i in range(tmptree.kinfit_results.size()):
            fitParValues[i+1] = tmptree.kinfit_results[i] 

        # Get template reweighting factors for qq->ttbar MC and must be ejets events
        weight_is_valid = 0
        # set default weight to be 1
        w_a[0],w_a_opp[0],w_s_xi[0],w_s_xi_opp[0],w_a_xi[0],w_a_xi_opp[0] = 1,1,1,1,1,1
        w_s_delta[0],w_s_delta_opp[0],w_a_delta[0],w_a_delta_opp[0] = 1,1,1,1
        # decide if want to use temp reweight
        if is_TT_MC: 
            if tmptree.gen_type[0] == 'e_jets' and tmptree.init_type[0] == 'qqbar':
                weight_is_valid = 1
        beta_v0[0] = -1            
        if weight_is_valid == 1 :
            # print 'getting special ws!'
            # Make 4vecs of gen tlep_p4,thad_p4,wlep_p4,whad_p4
            gen_pt = tmptree.gen_pt
            gen_eta = tmptree.gen_eta
            gen_phi = tmptree.gen_phi
            gen_mass = tmptree.gen_mass
            gen_side = tmptree.gen_side
            gen_pdgid = tmptree.gen_pdgid       
            top_MC,Atop_MC = ROOT.TLorentzVector(),ROOT.TLorentzVector()            
            for i in range(gen_pdgid.size()):
                if gen_pdgid[i] == 6 :
                    top_MC.SetPtEtaPhiM(gen_pt[i],gen_eta[i],gen_phi[i],gen_mass[i])
                if gen_pdgid[i] == -6 :
                    Atop_MC.SetPtEtaPhiM(gen_pt[i],gen_eta[i],gen_phi[i],gen_mass[i])

            # Get weights
            tmp_w = GetAnglesWeights(top_MC,Atop_MC,tmptree.cos_theta_mc[0])
            w_a[0],w_a_opp[0],w_s_xi[0],w_s_xi_opp[0] = tmp_w[0],tmp_w[1],tmp_w[2],tmp_w[3]
            w_a_xi[0],w_a_xi_opp[0],w_s_delta[0],w_s_delta_opp[0] = tmp_w[4],tmp_w[5],tmp_w[6],tmp_w[7]
            w_a_delta[0],w_a_delta_opp[0] = tmp_w[8],tmp_w[9]
            # beta study
            beta_v0[0]=tmp_w[10]

        # gen branches only exist for ttbar MC, which makes PERFECT sense
        motherPIDs[0],motherPIDs[1] = -100,-100
        if is_TT_MC:
            motherPIDs[0],motherPIDs[1] = tmptree.gen_pdgid[0],tmptree.gen_pdgid[1]

        # a bunch of MC correction weights. Applies only for MC, literaly
        # So far, only weight_top_pT,w_btag,w_eleID,w_trigger,PU are there.
        top_pT_reweight[0],btag_eff_reweight[0],btag_eff_reweight_hi[0],btag_eff_reweight_low[0] = 1,1,1,1
        lepID_reweight[0],lepID_reweight_hi[0],lepID_reweight_low[0],trigger_reweight[0] = 1,1,1,1
        trigger_reweight_hi[0],trigger_reweight_low[0] = 1,1
        pileup_real[0] = -10
        if is_MC:
            pileup_reweight[0] = tmptree.w_PU[0]
            if tmptree.FindBranch('weight_top_pT'):
                top_pT_reweight[0]  = tmptree.weight_top_pT[0]/toppt_scale
            btag_eff_reweight[0]    = tmptree.w_btag[0]
            btag_eff_reweight_hi[0] = tmptree.w_btag_up[0]
            btag_eff_reweight_low[0]= tmptree.w_btag_down[0]
            lepID_reweight[0]       = tmptree.w_eleID[0]
            lepID_reweight_hi[0]    = tmptree.w_eleID_up[0]
            lepID_reweight_low[0]   = tmptree.w_eleID_down[0]
            trigger_reweight[0]     = tmptree.w_trigger[0]
            trigger_reweight_hi[0]  = tmptree.w_trigger_up[0]
            trigger_reweight_low[0] = tmptree.w_trigger_down[0]
            # PU
            pileup_real[0] = tmptree.mc_pileup_events[0]
            # total weight
            correction_weight[0] = top_pT_reweight[0]*btag_eff_reweight[0]*lepID_reweight[0]*trigger_reweight[0]*pileup_reweight[0]
            correction_weight[0]*= lepIso_reweight[0]*tracking_reweight[0]

        # Added branch
        if n_valid_jets[0] == 4 and Q_l[0]==  1 : 
            charge_ratio[0]=1
            h_charge_ratio.Fill('4jets,l+',1);h_charge_ratio.Fill('4jets,l-',0);h_charge_ratio.Fill('5jets,l+',0);h_charge_ratio.Fill('5jets,l-',0)
        if n_valid_jets[0] == 4 and Q_l[0]== -1 : 
            charge_ratio[0]=2
            h_charge_ratio.Fill('4jets,l+',0);h_charge_ratio.Fill('4jets,l-',1);h_charge_ratio.Fill('5jets,l+',0);h_charge_ratio.Fill('5jets,l-',0)
        if n_valid_jets[0] == 5 and Q_l[0]==  1 : 
            charge_ratio[0]=3
            h_charge_ratio.Fill('4jets,l+',0);h_charge_ratio.Fill('4jets,l-',0);h_charge_ratio.Fill('5jets,l+',1);h_charge_ratio.Fill('5jets,l-',0)
        if n_valid_jets[0] == 5 and Q_l[0]== -1 : 
            charge_ratio[0]=4
            h_charge_ratio.Fill('4jets,l+',0);h_charge_ratio.Fill('4jets,l-',0);h_charge_ratio.Fill('5jets,l+',0);h_charge_ratio.Fill('5jets,l-',1)

        ################################################################
        #               reconstruction study                           #
        ################################################################

        # Find and store reco p-4, leading jets pt and mass, for all MC and data samples
        for i in range(4):
            reco_pt[i]   = tmptree.reco_pt[i]
            reco_eta[i]  = tmptree.reco_eta[i]
            reco_phi[i]  = tmptree.reco_phi[i]
            reco_mass[i] = tmptree.reco_mass[i]

        # sort jets by pt
        max_jetpt,max_jet_index = 0,0
        for index,item in enumerate(tmptree.jets_pt):
            if item>max_jetpt:
                max_jet_index = index
        # Keep leading jet pt and mass
        leadingJet_pt[0] = tmptree.jets_pt[max_jet_index]
        leadingJet_mass[0] = tmptree.jets_mass[max_jet_index]

        # Match reco t's and w's with gen objects using deltaR for TT semilep events only
        do_Matching = 0
        if is_TT_MC:
            if tmptree.gen_type[0] == 'e_jets':
                do_Matching = 1

        # initialization
        for i in range(len(isMatched)):
            isMatched[i] = -1

        if do_Matching:
            # Make 4-vec for gen t's and w's
            # Make 4vec for reco t and W's
            reco_p4 = []
            for i in range(4):
                ip4 = ROOT.TLorentzVector()
                ip4.SetPtEtaPhiM(reco_pt[i],reco_eta[i],reco_phi[i],reco_mass[i])
                reco_p4.append(ip4)
            # Make 4vecs of gen tlep_p4,thad_p4,wlep_p4,whad_p4
            gen_pt = tmptree.gen_pt
            gen_eta = tmptree.gen_eta
            gen_phi = tmptree.gen_phi
            gen_mass = tmptree.gen_phi
            gen_side = tmptree.gen_side
            gen_pdgid = tmptree.gen_pdgid
            gen_index = [0,0,0,0]
            for i in range(len(gen_side)):
                if gen_side[i] == 'lep' and abs(gen_pdgid[i])==6 : gen_index[0]=i 
                if gen_side[i] == 'had' and abs(gen_pdgid[i])==6 : gen_index[1]=i 
                if gen_side[i] == 'lep' and abs(gen_pdgid[i])==24 : gen_index[2]=i 
                if gen_side[i] == 'had' and abs(gen_pdgid[i])==24 : gen_index[3]=i 
            gen_p4 = []
            for j in range(4):
                ip4 = ROOT.TLorentzVector()
                i = gen_index[j]
                ip4.SetPtEtaPhiM(gen_pt[i],gen_eta[i],gen_phi[i],gen_mass[i])
                gen_p4.append(ip4)  
            # Do matching
            for i in range(4):
                delR_ = reco_p4[i].DeltaR(gen_p4[i])
                if delR_< deltaR_matching: isMatched[i] = 1
                else : isMatched[i] = 0                          

        # Do matching for all signal events

        # Fill events that pass additional cuts
        newtree.Fill()

    # Write and close files
    fout.Write()
    fout.Close()
    tfile.Close()

    return 'Make templates finished!'

main()
