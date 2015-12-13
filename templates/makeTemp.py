# Take in the angles ttree from reco, and make input ttree for original template making and fitting codes
# This is the debug version
# 12-12-15

# from Legacy_Afb.Tools.ttbar_utility import *
# from Legacy_Afb.Tools.angles_tools import *
import glob
from optparse import OptionParser
import ROOT
from array import array
from angles_tools import *
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

    correction_weight = array('f',[1.])
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

    br_defs += [('correction_weight',correction_weight,'correction_weight/F')]

    # study beta
    beta_v0 = array('f',[0.])
    beta_v1 = array('f',[0.])
    beta_v2 = array('f',[0.])
    mtt_mc = array('f',[0.])
    mt_mc = array('f',[0.])
    m_tbar_mc = array('f',[0.])

    br_defs += [('beta_v0',beta_v0,'beta_v0/F')]
    br_defs += [('beta_v1',beta_v1,'beta_v1/F')]
    br_defs += [('beta_v2',beta_v2,'beta_v2/F')]
    br_defs += [('mtt_mc',mtt_mc,'mtt_mc/F')]
    br_defs += [('mt_mc',mt_mc,'mt_mc/F')]
    br_defs += [('m_tbar_mc',m_tbar_mc,'m_tbar_mc/F')]

    # study cos_theta
    positive_z = array('i',[0])
    cs_bisector = array('f',[0.])
    beta_mc = array('f',3*[0.])
    beta1 = array('f',3*[0.])
    beta2 = array('f',3*[0.])
    delta_mag = array('f',[0.])
    br_defs += [('positive_z',positive_z,'positive_z/I')]
    br_defs += [('cs_bisector',cs_bisector,'cs_bisector/F')]
    br_defs += [('beta_mc',beta_mc,'beta_mc[3]/F')]
    br_defs += [('beta1',beta1,'beta1[3]/F')]
    br_defs += [('beta2',beta2,'beta2[3]/F')]
    br_defs += [('delta_mag',delta_mag,'delta_mag/F')]

      
    # Add branches to the tree
    for ibr in br_defs:
        newtree.Branch(ibr[0],ibr[1],ibr[2])
    # Some added branches
    all_vecs = []
    charge_ratio = ROOT.vector('int')()
    newtree.Branch('charge_ratio',charge_ratio)
    all_vecs += [charge_ratio]

    # Add cutflow diagrams
    h_cutflow = ROOT.TH1D('cutflow_extra',' cutflow extra;cuts;events',5,0.,5.)
    # h_cutflow.SetBit(ROOT.TH1.kCanRebin)
    h_cutflow.SetDirectory(fout)

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


    # Loop over entries
    n_evt = 0
    for iev in range(evt_start,tmptree.GetEntries()):
        # Progress report
        if iev%50 == 0 : print 'processing event',iev 
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

        # Fill branches
        for ivec in all_vecs: ivec.clear()

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
        if tmptree.FindBranch('gen_type'): 
            if tmptree.gen_type[0] == 'e_jets' and tmptree.init_type[0] == 'qqbar':
                weight_is_valid = 1

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
            # beta v0, the original version in Nick's code
            # beta_v0[0]=tmp_w[10]
            # Alternative beta definitions

        # Some study of cos_theta_mc for all TTbar MC samples
        # init
        beta_v0[0]=-1
        beta_v1[0]=-1
        beta_v2[0]=-1
        mtt_mc[0] = -1
        mt_mc[0],m_tbar_mc[0] = -1,-1        
        positive_z[0],cs_bisector[0],delta_mag[0] = -10,-10,-1000
        for i in range(3):
            beta_mc[i],beta1[i],beta2[i] = 0,0,0

        if tmptree.FindBranch('gen_type'): 
            # Make 4vecs of gen tlep_p4,thad_p4,wlep_p4,whad_p4
            gen_pt = tmptree.gen_pt
            gen_eta = tmptree.gen_eta
            gen_phi = tmptree.gen_phi
            gen_mass = tmptree.gen_mass
            gen_side = tmptree.gen_side
            gen_pdgid = tmptree.gen_pdgid       
            top_MC,Atop_MC = ROOT.TLorentzVector(),ROOT.TLorentzVector()
            init1,init2 = ROOT.TLorentzVector(),ROOT.TLorentzVector()            
            for i in range(gen_pdgid.size()):
                if gen_pdgid[i] == 6 :
                    top_MC.SetPtEtaPhiM(gen_pt[i],gen_eta[i],gen_phi[i],gen_mass[i])
                if gen_pdgid[i] == -6 :
                    Atop_MC.SetPtEtaPhiM(gen_pt[i],gen_eta[i],gen_phi[i],gen_mass[i]) 
            init1.SetPtEtaPhiM(gen_pt[0],gen_eta[0],gen_phi[0],gen_mass[0])
            init2.SetPtEtaPhiM(gen_pt[1],gen_eta[1],gen_phi[1],gen_mass[1])
            # Study beta
            mtt_mc[0]= (top_MC+Atop_MC).M()
            mt_mc[0] = top_MC.M() 
            m_tbar_mc[0] = Atop_MC.M()

            beta_v1[0] = Get_beta_v0(mt_mc[0],m_tbar_mc[0],mtt_mc[0])
            # beta study v1. Use simple definition assuming mt1=mt2  
            beta_v1[0] = Get_beta_v1(mt_mc[0],m_tbar_mc[0],mtt_mc[0])      
            # beta study v2. Use beta = sqrt(1-2*(m1*m2+m2*m2)/Mtt*Mtt+(m1*m1-m2*m2)^2/Mtt^4)
            beta_v2[0] = Get_beta_v2(mt_mc[0],m_tbar_mc[0],mtt_mc[0])            
            # Study reco cos_theta_cs first
            ttbar_cm = tlep_p4+thad_p4
            if ttbar_cm.Pz()>0 : 
                positive_z[0] = 1
            else : 
                positive_z[0] = -1

            if lep_charge >0 :
                t_p4 = tlep_p4.Clone()
                tbar_p4 = thad_p4.Clone()
            elif lep_charge <0 :
                t_p4 = thad_p4.Clone()
                tbar_p4 = tlep_p4.Clone() 

            cs_results = get_angles_v1(t_p4,tbar_p4)
            # return [x_f,ttbar_mass_data,cos_theta_cs_data,cs_bisector,beta_mc,beta1,beta2,delta_mag]            
            cs_bisector[0] = cs_results[3]
            delta_mag[0] = cs_results[7]
            for i in range(3):
                beta_mc[i] = cs_results[4][i]
                beta1[i] = cs_results[5][i]
                beta2[i] = cs_results[6][i]

        # gen branches only exist for ttbar MC, which makes PERFECT sense
        motherPIDs[0],motherPIDs[1] = -100,-100
        if tmptree.FindBranch('gen_pdgid'):
            motherPIDs[0],motherPIDs[1] = tmptree.gen_pdgid[0],tmptree.gen_pdgid[1]
        # a bunch of MC correction weights. Applies only for MC, literaly
        # So far, only weight_top_pT,w_btag,w_eleID,w_trigger,PU are there.
        top_pT_reweight[0],btag_eff_reweight[0],btag_eff_reweight_hi[0],btag_eff_reweight_low[0] = 1,1,1,1
        lepID_reweight[0],lepID_reweight_hi[0],lepID_reweight_low[0],trigger_reweight[0] = 1,1,1,1
        trigger_reweight_hi[0],trigger_reweight_low[0] = 1,1
        pileup_real[0] = -10
        if tmptree.FindBranch('w_PU'):
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
        if n_valid_jets[0] == 4 and Q_l[0]==  1 : charge_ratio.push_back(1)
        if n_valid_jets[0] == 4 and Q_l[0]== -1 : charge_ratio.push_back(2)
        if n_valid_jets[0] == 5 and Q_l[0]==  1 : charge_ratio.push_back(3)
        if n_valid_jets[0] == 5 and Q_l[0]== -1 : charge_ratio.push_back(4)

        # Fill events that pass additional cuts
        newtree.Fill()

    # Write and close files
    fout.Write()
    fout.Close()
    tfile.Close()

    return 'Make templates finished!'
# m1 = m_t, m2 = m_tbar
def Get_beta_v2(m1,m2,mtt):
    return math.sqrt(1-2*(m1*m1+m2*m2)/pow(mtt,2)+pow((m1*m1-m2*m2),2)/pow(mtt,4))

def Get_beta_v1(m1,m2,mtt):
    if 1-pow(2*m1/mtt,2) < 0 : 
        return 0
    else: 
        return math.sqrt(1-pow(2*m1/mtt,2))

def Get_beta_v0(m1,m2,mtt):
    M2_1_mc   = m1*m2
    M2_2_mc   = m2*m2
    ttbar_mass_mc = mtt
    num_mc    = 1. - 2.*(M2_1_mc+M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc) + (M2_1_mc-M2_2_mc)*(M2_1_mc-M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc*ttbar_mass_mc*ttbar_mass_mc);
    denom_1_mc   = (1. + (M2_1_mc-M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc))*(1. + (M2_1_mc-M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc));
    denom_2_mc   = (1. + (M2_2_mc-M2_1_mc)/(ttbar_mass_mc*ttbar_mass_mc))*(1. + (M2_2_mc-M2_1_mc)/(ttbar_mass_mc*ttbar_mass_mc));
    beta_mc   = sqrt(sqrt((num_mc*num_mc)/(denom_1_mc*denom_1_mc) * (num_mc*num_mc)/(denom_2_mc*denom_2_mc)));    
    return beta_mc

main()
