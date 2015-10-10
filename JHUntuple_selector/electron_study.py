from utility import *

# Some control var
evt_to_run = -1 
do_plot = True
# Definition of global parameters
dR = 0.15
lepID = (11,12,13,14,15,16)
jetID = (1,2,3,4,5,21)
# Get input
events = Events("test_input.root")

# define handles here
el_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >')
el_Loose_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >')
el_iso_hndl = Handle('vector<double>')
el_Loose_iso_hndl = Handle('vector<double>')

el_isLoose_hndl = Handle('vector<unsigned int>')
el_isTight_hndl = Handle('vector<unsigned int>')
el_isModTight_hndl = Handle('vector<unsigned int>')

el_Loose_isLoose_hndl = Handle('vector<unsigned int>')
el_Loose_isTight_hndl = Handle('vector<unsigned int>')
el_Loose_isModTight_hndl = Handle('vector<unsigned int>')

gen_hndl = Handle('vector<reco::GenParticle>  ')
gen_label = 'prunedGenParticles'

mu_hndl = Handle('vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > >')
mu_label = ('jhuMuonPFlow','muon')
mu_iso_hndl =  Handle('vector<double>')
mu_iso_label = ('jhuMuonPFlow','muoniso')
mu_tight_hndl = Handle('vector<unsigned int>')
mu_tight_label = ('jhuMuonPFlow','muonistight')
# define label module names here
el_prefix = 'jhuElePFlow'
elLoose_prefix = 'jhuElePFlowLoose'
# define and initiate counters here
n_evt = 0
n_PFel,n_PFel_isloose,n_PFel_istight = 0,0,0
n_PFel_all,n_PFcommon = 0,0
n_gen_el, n_ejets = 0,0
# book histograms here
event_type = 'ttbar_powheg'
h_el_istight_pt = ROOT.TH1D('el_istight_pt',event_type+';pT;events',50,0.,200.)
h_el_isloose_pt = ROOT.TH1D('el_isloose_pt',event_type+';pT;events',50,0.,200.)
h_el_loose_pt = ROOT.TH1D('el_Loose_pt',event_type+';pT;events',50,0.,200.)
h_el_status = ROOT.TH1D('el_status',event_type+';type;events',5,0.,5.)
h_el_status.SetBit(ROOT.TH1.kCanRebin)
h_el_loose_status = ROOT.TH1D('el_Loose_status',event_type+';type;events',5,0.,5.)
h_el_loose_status.SetBit(ROOT.TH1.kCanRebin)
h_el_iso = ROOT.TH1D('el_iso',event_type+';RelIso;events',50,0.,0.3)
h_el_tight_iso = ROOT.TH1D('el_tight_iso',event_type+';RelIso;events',50,0.,0.3)
 
h_gen_ttbar_channel = ROOT.TH1D('ttbar_channel',event_type+';channel;events',7,0.,7.)
h_el_status.SetBit(ROOT.TH1.kCanRebin)
h_el_selected_pdgid =  ROOT.TH1D('h_el_selected_pdgid',event_type+';PDG ID;events',50,0.,25.)
h_el_unselected_pdgid =  ROOT.TH1D('h_el_unselected_pdgid',event_type+';PDG ID;events',50,0.,25.)
h_el_eff = ROOT.TH1D('el_gen_efficiency',event_type+';status;events',4,0.,4.)
h_el_fakerate = ROOT.TH1D('el_selection_fake_rate',event_type+';status;events',4,0.,4.)
h_el_eff.SetBit(ROOT.TH1.kCanRebin)
h_el_fakerate.SetBit(ROOT.TH1.kCanRebin)
h_el_match_multiplicity = ROOT.TH1D('h_el_match_multiplicity',event_type+';n_gen;events',50,0.,4.)

h_mu_eff = ROOT.TH1D('mu_gen_efficiency',event_type+';status;events',4,0.,4.)
h_mu_fakerate = ROOT.TH1D('mu_selection_fake_rate',event_type+';status;events',4,0.,4.)
h_mu_eff.SetBit(ROOT.TH1.kCanRebin)
h_mu_fakerate.SetBit(ROOT.TH1.kCanRebin)
h_mu_match_multiplicity = ROOT.TH1D('h_mu_match_multiplicity',event_type+';n_gen;events',50,0.,4.)

hist_gen_pdgID = ROOT.TH1D('hist_gen_pdgID',event_type+';PDG ID;events',50,0.,25.)
h_gen_multiplicity = ROOT.TH1D('jet_multiplicity',event_type+';n_jets;events',50,0.,8.)

# Set the y-axis minimum to 0 for counting histograms
hlist_Yzero = [h_gen_ttbar_channel,h_el_eff,h_el_fakerate,h_mu_eff,h_mu_fakerate]
for ihist in hlist_Yzero:
     ihist.SetMinimum(0)


# start main event loop
for evt in events:
    if n_evt%1000 == 1: print 'Loop over',n_evt,'event'
    if n_evt ==  evt_to_run : break
    n_evt += 1

    evt.getByLabel(el_prefix,'electron',el_hndl)
    evt.getByLabel(elLoose_prefix,'electronLoose',el_Loose_hndl)
    evt.getByLabel(el_prefix,'electroniso',el_iso_hndl)
    evt.getByLabel(el_prefix,'electronisloose',el_isLoose_hndl)
    evt.getByLabel(el_prefix,'electronistight',el_isTight_hndl)
    evt.getByLabel(el_prefix,'electronmodtight',el_isModTight_hndl)
    evt.getByLabel(elLoose_prefix,'electronLooseisloose',el_Loose_isLoose_hndl)
    evt.getByLabel(elLoose_prefix,'electronLooseistight',el_Loose_isTight_hndl)
    evt.getByLabel(elLoose_prefix,'electronLoosemodtight',el_Loose_isModTight_hndl)

    el_p4 = el_hndl.product()
    el_iso = el_iso_hndl.product()
    el_isLoose = el_isLoose_hndl.product()
    el_isTight = el_isTight_hndl.product()
    el_isModTight = el_isModTight_hndl.product()

    el_loose_p4 = el_Loose_hndl.product()
    el_loose_isLoose = el_Loose_isLoose_hndl.product()
    el_loose_isTight = el_Loose_isTight_hndl.product()
    el_loose_isModTight = el_Loose_isModTight_hndl.product()

    # Loop over PF electrons
    el_selected = []
    for i in range(len(el_p4)):
        # plot the number of tight and loose electrons
        h_el_status.Fill("PFelectrons",1)
        h_el_iso.Fill(el_iso[i])
        n_PFel += 1
        # PFelectrons passed loose selection
        if el_isLoose[i]:
            h_el_status.Fill("loose",1)   
            h_el_isloose_pt.Fill(el_p4[i].pt())
            n_PFel_isloose += 1
        # PFelectrons passed tight selection
        if el_isTight[i] and not el_isModTight[i] and el_iso[i]<0.1: 
            h_el_status.Fill("tight",1)
            h_el_tight_iso.Fill(el_iso[i])
            h_el_istight_pt.Fill(el_p4[i].pt())
            n_PFel_istight += 1
            el_selected.append(el_p4[i])
    # Find all PFlooseelectrons that are not selected based on tight condition in PFelectrons
    el_unselected = list( el for el in el_loose_p4 if el not in el_selected)
    # Loop over PFelectronLoose
    PFloose_el_isloose = []
    for i in range(len(el_loose_p4)):
        # plot the number of tight and loose electrons
        h_el_loose_status.Fill("PFlooseElectrons",1)
        if el_loose_isLoose[i]:h_el_loose_status.Fill("loose",1)
        if el_loose_isTight[i] and not el_loose_isModTight[i] : h_el_loose_status.Fill("tight",1)
        h_el_loose_pt.Fill(el_loose_p4[i].pt())
        if not (el_loose_isTight[i] and not el_loose_isModTight[i]) : PFloose_el_isloose.append(el_loose_p4[i])

    # Find the electrons that are in common in both collections 
    common_el = list(i for i in el_loose_p4 if i in el_p4)
    # counting
    n_PFel_all += len(el_p4)
    n_PFcommon += len(common_el)

    # Get gen particles and find out the true identy of the PF electron collection
    evt.getByLabel(gen_label,gen_hndl)
    genpars = gen_hndl.product()
    # get all final state particles,status = 3,and particles with no daughters(final state partons) or leptons or b's (which will further decay)
    final_par = [] 
    for ipar in genpars: 
        if ipar.status() == 3:
            if ipar.numberOfDaughters()==0 or abs(ipar.pdgId())==5 : final_par.append(ipar)
            elif ipar.numberOfMothers(): 
                if abs(ipar.mother(0).pdgId())==24: final_par.append(ipar) 
    
    for igen in final_par:
        hist_gen_pdgID.Fill(igen.pdgId())
    # Find events with electron in final states
    gen_el = list(ipar for ipar in final_par if ipar.status() == 3 and abs(ipar.pdgId()) == 11)
    gen_mu = list(ipar for ipar in final_par if ipar.status() == 3 and abs(ipar.pdgId()) == 13)
    gen_tau = list(ipar for ipar in final_par if ipar.status() == 3 and abs(ipar.pdgId()) == 15)
    # find the matched electron and genparticle based on deltaR
    el_selected_gen = list( (ipar,j) for ipar in final_par for j in el_selected if DeltaR2(ipar.p4(),j) < dR ) 
    el_unselected_gen = list( (ipar,j) for ipar in final_par for j in el_unselected if DeltaR2(ipar.p4(),j) < dR )      
    # print out the information of matched particles
    # Check the matching multiplicity 
    for icand in el_selected:
        n_matched_gen = 0
        for igen in final_par:
            if DeltaR2(icand,igen.p4())< dR : n_matched_gen += 1
        h_el_match_multiplicity.Fill(n_matched_gen)
    # pdgID of selected and unselected PFelectrons
    for ipar,j in el_selected_gen:
        h_el_selected_pdgid.Fill(abs(ipar.pdgId()))   
    for ipar,j in el_unselected_gen:
        h_el_unselected_pdgid.Fill(abs(ipar.pdgId()))
    # count how many events with electrons and the total number of electrons
    if len(gen_el)+len(gen_mu)+len(gen_tau) == 0 : h_gen_ttbar_channel.Fill('hadronic',1)
    if len(gen_el) == 1 and len(gen_el)+len(gen_mu)+len(gen_tau) == 1 : h_gen_ttbar_channel.Fill('e+jets',1)
    if len(gen_mu) == 1 and len(gen_el)+len(gen_mu)+len(gen_tau) == 1 : h_gen_ttbar_channel.Fill('mu+jets',1)
    if len(gen_tau) == 1 and len(gen_el)+len(gen_mu)+len(gen_tau) == 1 : h_gen_ttbar_channel.Fill('tau+jets',1)
    if len(gen_el)+len(gen_mu)+len(gen_tau) == 2 : h_gen_ttbar_channel.Fill('dilep',1)
    if len(gen_el)+len(gen_mu)+len(gen_tau) > 2  : h_gen_ttbar_channel.Fill('wtf',1)
    # check jet multiplicity for semi-lep ttbar events
    gen_jets = list(ipar for ipar in final_par if abs(ipar.pdgId()) in jetID)
    if len(gen_el)+len(gen_mu) == 1 and len(gen_tau) == 0:
        h_gen_multiplicity.Fill(len(gen_jets))
        # debugging
        if len(gen_jets) < 4:
            final_parid = list(ipar.pdgId() for ipar in final_par)
            gen_par_look = list((ipar.pdgId(),ipar.status()) for ipar in genpars)
            print final_parid
            print gen_par_look
            break
    # Check the selection efficiency and fake rate
    el_selected_matched = list(j for ipar,j in el_selected_gen)
    el_selected_true = list(j for ipar,j in el_selected_gen if abs(ipar.pdgId())==11)
    gen_el_selected_true = list(ipar for ipar,j in el_selected_gen if abs(ipar.pdgId())==11)
    el_selected_fake = list(j for j in el_selected if j not in el_selected_true)
    # gen-electron selection efficiency 
    for ipar in gen_el : 
        h_el_eff.Fill('gen_el',1) 
        # gen electrons passed the same kin cuts as PFelectrons
        if ipar.pt()>10 and ipar.eta()<2.5 and not (ipar.eta()>1.44 and ipar.eta()<1.56) :
             h_el_eff.Fill('gen_el_kin_cuts',1)
        if ipar in gen_el_selected_true: h_el_eff.Fill('selected',1) 
    # PFelectron selection fake rate
    for ipar in el_selected: 
        h_el_fakerate.Fill('selected PFel',1)
        if ipar in el_selected_matched: h_el_fakerate.Fill('gen_matched',1)
        if ipar in el_selected_true: h_el_fakerate.Fill('gen_el_matched',1)

    # Check muon selection efficiency 
    evt.getByLabel(mu_label,mu_hndl)
    evt.getByLabel(mu_tight_label,mu_tight_hndl)
    evt.getByLabel(mu_iso_label,mu_iso_hndl)
    mu_p4 = mu_hndl.product()
    mu_is_tight = mu_tight_hndl.product()
    mu_iso = mu_iso_hndl.product()
    # select tight muon following the twiki recommendations
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopMUORun1
    mu_cand = []
    for i in range(len(mu_p4)):
        if mu_is_tight[i] and mu_iso[i]<0.12 : mu_cand.append(mu_p4[i])
    mu_cand_gen =  list( (ipar,j) for ipar in final_par for j in mu_cand if DeltaR2(ipar.p4(),j) < dR )
    mu_cand_matched = list( j for ipar,j in mu_cand_gen)
    mu_cand_true = list( j for ipar,j in mu_cand_gen if abs(ipar.pdgId())==13)
    gen_mu_cand_true = list( ipar for ipar,j in mu_cand_gen if abs(ipar.pdgId())==13)
    # Check the matching multiplicity 
    for icand in mu_cand:
        n_matched_gen = 0
        for igen in final_par:
            if DeltaR2(icand,igen.p4())< dR : n_matched_gen += 1 
        h_mu_match_multiplicity.Fill(n_matched_gen)
        
    # selection efficiency
    for ipar in gen_mu : 
        h_mu_eff.Fill('gen_mu',1) 
        # gen muons passed the same kin cuts as PFelectrons
        if ipar.pt()>10 and ipar.eta()<2.5 :
             h_mu_eff.Fill('gen_mu_kin_cuts',1)
        # gen muons that are matched with PFmuon candidates
        if ipar in gen_mu_cand_true: h_mu_eff.Fill('selected',1) 
    # PFmuon selection fake rate
    for ipar in mu_cand: 
        h_mu_fakerate.Fill('selected PFmuons',1)
        if ipar in mu_cand_matched : h_mu_fakerate.Fill('gen matched',1)
        if ipar in mu_cand_true : h_mu_fakerate.Fill('gen_mu_matched',1)
 
   # end event loop
print 'last event looped',n_evt
# Some outputs
#print 'num of PFel, PFel_loose, PFel_tight :',n_PFel,n_PFel_isloose,n_PFel_istight
#print 'num of PFel,PFcommon :',n_PFel_all,n_PFcommon

# Make and save plots
#histlist = [h_el_status,h_el_iso,h_el_isloose_pt,h_el_tight_iso,h_el_istight_pt,h_el_loose_status,h_el_loose_pt]
histlist = [h_el_selected_pdgid,h_el_unselected_pdgid,h_gen_ttbar_channel,h_el_eff,h_el_fakerate,h_mu_eff,h_mu_fakerate,h_mu_match_multiplicity,h_el_match_multiplicity,hist_gen_pdgID,h_gen_multiplicity]
if do_plot : plotting(histlist,event_type,True)

