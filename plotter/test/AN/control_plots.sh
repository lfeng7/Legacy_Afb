python ../masterplot.py dir ../mu_nominal/, var1 lep_pt, bin 50, min 20, max 250, xaxis "p_{T}^{lep} (GeV)", title "mu lep_pt", plotname lep_pt_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 met_pt, bin 50, min 0, max 300, xaxis "MET (GeV)", plotname MET_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 jets_pt, bin 50, min 30, max 300, xaxis "p_{T}^{Jets} (GeV)", plotname jets_pt_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 N_btag, bin 5, min 0, max 5, xaxis "Num (b-jets)", plotname nbtag_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"
