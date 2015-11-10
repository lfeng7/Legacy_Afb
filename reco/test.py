from Legacy_Afb.Tools.kinfit import *
fin = ROOT.TFile('selected_files/v2_trigger_removed/all/TT_CT10_TuneZ2star_8TeV_selected.root')
tmptree = fin.Get('selected')

for iev in range(tmptree.GetEntries()):
    tmptree.GetEntry(iev)

    n_selected_jets = jets_pt.size()
    if not 4<= n_selected_jets <=5 : continue

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

    # Set up inputs  
    lep_p4 = ROOT.TLorentzVector()
    lep_p4.SetPtEtaPhiM(lep_pts[0],lep_etas[0],lep_phis[0],lep_mass[0])
    jets_p4,jets_csv_list = [],[]    
    for i in range(jets_pt.size()):
        tmp_p4 = ROOT.TLorentzVector()
        tmp_p4.SetPtEtaPhiM(jets_pt[i],jets_eta[i],jets_phi[i],jets_phi[i],jets_mass[i])
        jets_p4.append(tmp_p4)
        jets_csv_list.append(jets_csv[i])
    metPt = met_pt[0]; metPhi = met_phi[0]
    # Do reco
    lep_type = 'el'
    mcordata = 'mc'

    reco_result = DoReco(jets_p4,jets_csv_list,lep_p4,metPt,metPhi,lep_type,mcordata)
    break
