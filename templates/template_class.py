# class def of template

# Here, a definition of MCinfo class should be imported 
# The idea is , all the information stored in MC.txt would be read in, processed, and 
# stored in an object of the class MCinfo. Later down the line, if we need info e.g.
# N_evts_gen we can simply do sth like self.mc_file.get_nevtsGen()

import ROOT
from array import array
#from angles_tools import *
import math
from Legacy_Afb.Tools.ttbar_utility import *
from Legacy_Afb.Tools.angles_tools import *


# import class that handles MC_info txt file
from MCinfo_class import MC_info

alpha = -0.26
# a helper function to handle vector and float at the same time
def take(item):
    if type(item) in [int,float]: return item
    elif len(item)!=0: return item[0]
    else: return None 

class template:
    """
    class to make nicely formatted template files , given a input root angles root file

    """
    csvm = 0.679 
    deltaR_matching = 0.4

    def __init__(self):
        # attributes 
        self.template_type = 'MC'   # data or MC
        self.evt_start = 0
        self.evt_to_run = 1000
        self.sample_name = None     # 'Tbar_s'
        self.input_file = None
        self.input_ttree_name = 'selected'
        self.ttbar_type = 'none'
        self.lep_type = 'e_jets'

        # for calculating the normalization
        self.totalLumi = 19700
        self.isSideband = False
        self.norm_w = 1.0
        # all attributes but totalLumi are from class MC_entry        
        self.xsec = None
        self.nevts_gen = None
        self.sample_key = None
        self.title = None
        self.type = None

        # new_MCinfo is a object of class MC_info
        # finding the entry can be done by new_MCinfo.get_entry('T_s')
        # getting any info in that entry can be done as new_MCinfo.get_entry('T_s').nevts_gen ( sample_key, title,type,xsec )
        self.MCinfo_txt = None
        self.new_MCinfo = None
      
    def set_MCinfo(self,MCinfo_txt):
        # create a new object of class MC_info as a container for MC.txt entries
        if MCinfo_txt is not None and MCinfo_txt is not '':
            self.new_MCinfo = MC_info(MCinfo_txt)
            self.MCinfo_txt = MCinfo_txt
            print 'Importing MC_info!'
        else:
            print 'No MC_info.txt provided. Will not import MC_info!'


    def set_templateType(self):
        # determine the type of input template to be either 'data' or 'MC'
        if 'Single' in self.sample_name or 'Run' in self.sample_name:
            self.template_type ='data'
            print 'Input file is data!'
        else:
            self.template_type = 'mc'
            print 'Input file is MC!'

    def set_newTreeName(self):
        # if input file is QCD sideband, the output templates should store 'angles' so as all other regular MC files
        # if input is data file, and not QCD sideband, should store ttree as 'angles_data'
        # the ttree naming convention is to be compatible with fitting code

        if self.isSideband: 
            self.newtree_name = 'angles'
        else:
            if self.template_type == 'data':
                self.newtree_name = 'angles_data'
            else:
                self.newtree_name = 'angles'
        print 'set_newTreeName: %s'%self.newtree_name


    def store_MC_entry(self,entry):
        # given an entry object of type MC_entry, set self.atributes from MC_entry.attributes

        self.xsec = float(entry.xsec)
        self.nevts_gen = float(entry.nevts_gen)
        self.type = entry.type
        self.sample_key = entry.sample_key  
        self.title = entry.title


    def set_normalizationWeight(self):
        # return normalization weight if MC_info is setup and template_type is MC

        if self.template_type == 'data':
            self.norm_w = 1.0
        elif self.new_MCinfo is not None :
            # set infos for this MC
            entry = self.new_MCinfo.get_entry(self.sample_name)
            if entry is None:
                print 'No entry found!'
                self.norm_w = 1.0
            else:
                self.store_MC_entry(entry)
                # calculate the normalization weight using xsec and IntegratedLumi
                self.norm_w = self.totalLumi*self.xsec/self.nevts_gen
                # if it is a MC in sideband region, it is a contamination of data-driven QCD templates,
                # In this case, we set the event weight to be negative
                if self.isSideband:
                    self.norm_w = -1.0*self.norm_w
        else:
            self.norm_w = 1.0
        print '(info) sample name %s, normalization weight is %.3f'%(self.sample_key,self.norm_w)    


    # Main method for template class
    def makeTemps(self):
    # Do reconstruction for the given file, process only a specific part of it, and output to a new file
    # self.MC_info has been added if MC.txt is provided, otherwise self.MC_info is None

        # some processing of input files

        # Determine if the input is data or MC
        self.set_templateType()
        # Determine the name of the output new tree, to be either 'angles' or 'angles_data'
        self.set_newTreeName()
        # Calculate normalziation weight
        self.set_normalizationWeight()

        # set the values of some class global var
        evt_to_run = self.evt_to_run
        evt_start = self.evt_start
        # Get ttree
        tfile = self.input_file
        tmptree = tfile.Get(self.input_ttree_name)
        print 'Start to process file',tfile.GetName()
        if evt_to_run>0 :
            evt_end = evt_start+evt_to_run
        else : # (evt_to_run = -1 indicates run all events)
            evt_end = -1 
        print '\nPrcess from event',evt_start,'to event',evt_end
        # Check if the starting event is more than total number of evts in current files
        if evt_start > tmptree.GetEntries():
            return 'evt_start > total_num_events in input ttree. Will stop here.'
            
        # Make output file
        if self.ttbar_type == 'none':
            fout_name = '%s__template_%i.root'%(self.sample_name,self.evt_start)
        else :
            fout_name = '%s_%s__%s__template_%i.root'%(self.sample_name,self.ttbar_type,self.lep_type,self.evt_start)
        fout = ROOT.TFile(fout_name,'recreate')
        # Make output ttree
        newtree = ROOT.TTree(self.newtree_name,self.newtree_name)

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

        Pdf_weights = array('f',[1.])
        Pdf_weights_low = array('f',[1.])
        Pdf_weights_hi = array('f',[1.])

        pileup_reweight = array('f',[1.])
        pileup_reweight_low = array('f',[1.])
        pileup_reweight_hi = array('f',[1.])


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

        # extra stuff might not used by template fitter
        correction_weight = array('f',[1.])
        reco_pt   = array('f',4*[0.])
        reco_eta  = array('f',4*[0.])
        reco_phi  = array('f',4*[0.])
        reco_mass = array('f',4*[0.])
        leadingJet_pt = array('f',[0.])
        leadingJet_mass = array('f',[0.])
        isMatched = array('i',4*[-1])
        # for qcd_sideband normalizarion
        normalization_weight = array('f',[1.])
	# for weighted MC
        gen_weight = array('f',[1.])


        # names 
        br_defs += [('cos_theta_cs',cos_theta_cs,'cos_theta_cs/F')]
        br_defs += [('Feynman_x',Feynman_x,'Feynman_x/F')]
        br_defs += [('ttbar_mass',ttbar_mass,'ttbar_mass/F')]
        br_defs += [('cos_theta_mc',cos_theta_mc,'cos_theta_mc/F')]
        br_defs += [('Feynman_x_mc',Feynman_x_mc,'Feynman_x_mc/F')]
        br_defs += [('ttbar_mass_mc',ttbar_mass_mc,'ttbar_mass_mc/F')]
        br_defs += [('Q_l',Q_l,'Q_l/I')]
        br_defs += [('n_bTags',n_bTags,'n_bTags/I')]

        br_defs += [('motherPIDs',motherPIDs,'motherPIDs[2]/I')]
        br_defs += [('Qt',Qt,'Qt/F')]
        br_defs += [('lnL',lnL,'lnL/F')]
        br_defs += [('n_valid_jets',n_valid_jets,'n_valid_jets/I')]
        br_defs += [('fit_MET',fitParValues,'fit_MET/F')]
        br_defs += [('pileup',pileup,'pileup/F')]
        br_defs += [('pileup_real',pileup_real,'pileup_real/F')]

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

        br_defs += [('top_pT_reweight',top_pT_reweight,'top_pT_reweight/F')]
        br_defs += [('GJR_reweight',GJR_reweight,'GJR_reweight/F')]
        br_defs += [('CT10_reweight',CT10_reweight,'CT10_reweight/F')]
        br_defs += [('cteq_reweight',cteq_reweight,'cteq_reweight/F')]

        br_defs += [("Pdf_weights",Pdf_weights,"Pdf_weights/F")]
        br_defs += [("Pdf_weights_low",Pdf_weights_low,"Pdf_weights_low/F")]
        br_defs += [("Pdf_weights_hi",Pdf_weights_hi,"Pdf_weights_hi/F")]

        br_defs += [('pileup_reweight',pileup_reweight,'pileup_reweight/F')]
        br_defs += [('pileup_reweight_low',pileup_reweight_low,'pileup_reweight_low/F')]
        br_defs += [('pileup_reweight_hi',pileup_reweight_hi,'pileup_reweight_hi/F')]


        br_defs += [('btag_eff_reweight',btag_eff_reweight,'btag_eff_reweight/F')]
        br_defs += [('btag_eff_reweight_low',btag_eff_reweight_low,'btag_eff_reweight_low/F')]
        br_defs += [('btag_eff_reweight_hi',btag_eff_reweight_hi,'btag_eff_reweight_hi/F')]

        br_defs += [('tracking_reweight',tracking_reweight,'tracking_reweight/F')]
        br_defs += [('tracking_reweight_low',tracking_reweight_low,'tracking_reweight_low/F')]
        br_defs += [('tracking_reweight_hi',tracking_reweight_hi,'tracking_reweight_hi/F')]

        br_defs += [('lepID_reweight',lepID_reweight,'lepID_reweight/F')]
        br_defs += [('lepID_reweight_low',lepID_reweight_low,'lepID_reweight_low/F')]
        br_defs += [('lepID_reweight_hi',lepID_reweight_hi,'lepID_reweight_hi/F')]

        br_defs += [('lepIso_reweight',lepIso_reweight,'lepIso_reweight/F')]
        br_defs += [('lepIso_reweight_low',lepIso_reweight_low,'lepIso_reweight_low/F')]
        br_defs += [('lepIso_reweight_hi',lepIso_reweight_hi,'lepIso_reweight_hi/F')]

        br_defs += [('trigger_reweight',trigger_reweight,'trigger_reweight/F')]        
        br_defs += [('trigger_reweight_low',trigger_reweight_low,'trigger_reweight_low/F')]
        br_defs += [('trigger_reweight_hi',trigger_reweight_hi,'trigger_reweight_hi/F')]

        # extra stuff not used by template fitter
        br_defs += [('correction_weight',correction_weight,'correction_weight/F')]
        br_defs += [('reco_pt',reco_pt,'reco_pt[4]/F')]
        br_defs += [('reco_eta',reco_eta,'reco_eta[4]/F')]
        br_defs += [('reco_phi',reco_phi,'reco_phi[4]/F')]
        br_defs += [('reco_mass',reco_mass,'reco_mass[4]/F')]
        br_defs += [('isMatched',isMatched,'isMatched[4]/I')]
        br_defs += [('leadingJet_pt',leadingJet_pt,'leadingJet_pt/F')]
        br_defs += [('leadingJet_mass',leadingJet_mass,'leadingJet_mass/F')]

        br_defs += [('normalization_weight',normalization_weight,'normalization_weight/F')]
        br_defs += [('gen_weight',gen_weight,'gen_weight/F')]

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
            if item in self.sample_name:
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
        # Check if has_genW
        if tmptree.FindBranch('weight_gen'):
            has_genW = True
        else:
            has_genW = False

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
            if not take(tmptree.final_errflags)==0 : continue
            h_cutflow.Fill('kinfit error',1)

            # separate TTbar into qq,gg,other templates
            # qq/gg are qq/gg->ttbar->ljets events
            # tt_other are anythong->ttbar->non ljets events        
            if self.ttbar_type=='bkg' and tmptree.gen_type[0] == self.lep_type : continue
            if self.ttbar_type=='qq' and not (tmptree.gen_type[0]==self.lep_type and tmptree.init_type[0]=='qqbar'):continue
            if self.ttbar_type=='gg' and not (tmptree.gen_type[0]==self.lep_type and tmptree.init_type[0]!='qqbar'):continue

            # if gen_w is in angles files, copied it over
            if has_genW:
                gen_weight[0]=tmptree.weight_gen
            ttbar_mass[0] = take(tmptree.mtt)
            # reco_p4 is a list of tlep_p4,thad_p4,wlep_p4,whad_p4
            tlep_pt = tmptree.reco_pt[0]
            tlep_eta = tmptree.reco_eta[0]
            tlep_phi = tmptree.reco_phi[0]
            tlep_mass = tmptree.reco_mass[0]
            thad_pt = tmptree.reco_pt[1]
            thad_eta = tmptree.reco_eta[1]
            thad_phi = tmptree.reco_phi[1]
            thad_mass = tmptree.reco_mass[1]        
            lep_charge = take(tmptree.lep_charge)
            # Make 4vec of tlep and thad
            tlep_p4 = ROOT.TLorentzVector()
            thad_p4 = ROOT.TLorentzVector()
            tlep_p4.SetPtEtaPhiM(tlep_pt,tlep_eta,tlep_phi,tlep_mass)
            thad_p4.SetPtEtaPhiM(thad_pt,thad_eta,thad_phi,thad_mass) 
            Q_Data = tlep_p4+thad_p4
            Qt[0] = math.sqrt(Q_Data.Px()*Q_Data.Px()+Q_Data.Py()*Q_Data.Py())
            cos_theta_cs[0] = take(tmptree.cos_theta)
            Feynman_x[0] = take(tmptree.xf)
            Q_l[0] = take(tmptree.lep_charge)

            cos_theta_mc[0] = -10
            if tmptree.FindBranch('cos_theta_mc'):
                if take(tmptree.cos_theta_mc) is not None: cos_theta_mc[0] = take(tmptree.cos_theta_mc) 
            Feynman_x_mc[0] = -10
            if tmptree.FindBranch('xf_mc'):
                if take(tmptree.xf_mc) is not None : Feynman_x_mc[0] = take(tmptree.xf_mc) 
            ttbar_mass_mc[0] = -10
            if tmptree.FindBranch('mtt_mc'):
                if take(tmptree.mtt_mc) is not None : ttbar_mass_mc[0] = take(tmptree.mtt_mc)
            lnL[0] = take(tmptree.final_chi2)

            if has_chi2_new :
                chi2_new[0] = tmptree.final_chi2_new

            n_valid_jets[0] = take(tmptree.N_jets)
            n_bTags[0] = take(tmptree.N_btag)
            fitParValues[0] = take(tmptree.final_nv_pz)
            for i in range(tmptree.kinfit_results.size()):
                fitParValues[i+1] = tmptree.kinfit_results[i] 

            # Get template reweighting factors for qq->ttbar MC and must be ejets events
            weight_is_valid = 0
            # set default weight to be 1
            w_a[0],w_a_opp[0],w_s_xi[0],w_s_xi_opp[0],w_a_xi[0],w_a_xi_opp[0] = 1,1,1,1,1,1
            w_s_delta[0],w_s_delta_opp[0],w_a_delta[0],w_a_delta_opp[0] = 1,1,1,1
            # decide if want to use temp reweight
            if is_TT_MC: 
                if tmptree.gen_type[0] == self.lep_type and tmptree.init_type[0] == 'qqbar':
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
                tmp_w = GetAnglesWeights(top_MC,Atop_MC,tmptree.cos_theta_mc,alpha_input=alpha,beta=True)
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
            # total correction weight and normalization weight
            correction_weight[0],normalization_weight[0] = 1,1
            if is_MC:
                pileup_reweight[0] = take(tmptree.w_PU)
                pileup_reweight_low[0] = take(tmptree.w_PU_down)
                pileup_reweight_hi[0] = take(tmptree.w_PU_up)

                if tmptree.FindBranch('weight_top_pT'):
                    top_pT_reweight[0]  = take(tmptree.weight_top_pT)/toppt_scale

                if tmptree.FindBranch('w_PDF_up'):
                    Pdf_weights_hi[0] = tmptree.w_PDF_up    
                    Pdf_weights_low[0] = tmptree.w_PDF_down  

                btag_eff_reweight[0]    = take(tmptree.w_btag)
                btag_eff_reweight_hi[0] = take(tmptree.w_btag_up)
                btag_eff_reweight_low[0]= take(tmptree.w_btag_down)

                if tmptree.FindBranch('w_eleID'):
                    lepID_reweight[0]       = take(tmptree.w_eleID)
                    lepID_reweight_hi[0]    = take(tmptree.w_eleID_up)
                    lepID_reweight_low[0]   = take(tmptree.w_eleID_down)
                elif tmptree.FindBranch('w_lepID'):
                    lepID_reweight[0]       = take(tmptree.w_lepID)
                    lepID_reweight_hi[0]    = take(tmptree.w_lepID_up)
                    lepID_reweight_low[0]   = take(tmptree.w_lepID_down)

                trigger_reweight[0]     = take(tmptree.w_trigger)
                trigger_reweight_hi[0]  = take(tmptree.w_trigger_up)
                trigger_reweight_low[0] = take(tmptree.w_trigger_down)

                if tmptree.FindBranch('w_tracking_up'):
                    tracking_reweight[0] = tmptree.w_tracking
                    tracking_reweight_hi[0] = tmptree.w_tracking_up
                    tracking_reweight_low[0] = tmptree.w_tracking_down

                if tmptree.FindBranch('w_lepIso_up'):
                    lepIso_reweight[0] = tmptree.w_lepIso
                    lepIso_reweight_hi[0] = tmptree.w_lepIso_up
                    lepIso_reweight_low[0] = tmptree.w_lepIso_down


                # PU
                pileup_real[0] = take(tmptree.mc_pileup_events)
                # total weight
                correction_weight[0] = top_pT_reweight[0]*btag_eff_reweight[0]*lepID_reweight[0]*trigger_reweight[0]*pileup_reweight[0]
                correction_weight[0]*= lepIso_reweight[0]*tracking_reweight[0]
                # normalization weight
                normalization_weight[0] = self.norm_w


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
                if tmptree.gen_type[0] == take(self.lep_type):
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
                    if delR_< self.deltaR_matching: isMatched[i] = 1
                    else : isMatched[i] = 0                          

            # Do matching for all signal events

            # Fill events that pass additional cuts
            newtree.Fill()

        # Write and close files
        fout.Write()
        fout.Close()
        tfile.Close()

        return '\nMake templates finished!'
