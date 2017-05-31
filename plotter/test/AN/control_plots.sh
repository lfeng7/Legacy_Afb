# control plots for selection
python ../masterplot.py dir ../mu_nominal/, var1 lep_pt, bin 50, min 20, max 250, xaxis "p_{T}^{lep} (GeV)", title "mu lep_pt", plotname lep_pt_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 met_pt, bin 50, min 0, max 300, xaxis "MET (GeV)", plotname MET_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 jets_pt, bin 50, min 30, max 300, xaxis "p_{T}^{Jets} (GeV)", plotname jets_pt_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 N_btag, bin 5, min 0, max 5, xaxis "Num (b-jets)", plotname nbtag_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 pileup_events, bin 40, min 0, max 40, xaxis "Numeber of Vertices (Rewighted)", plotname PU_corr_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT", scale yes

# PU plots

python ../masterplot.py dir ../mu_nominal/, var1 pileup_events, bin 40, min 0, max 40, xaxis "Numeber of Vertices (Unweighted)", plotname PU_uncorr_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT", scale yes, weight_ignore w_PU

../masterplot.py file ../mu_nominal/TT_8TeV-mcatnlo_selected__reco.root , var mc_pileup_events, bin 50, min 0, max 50, xaxis "Numeber of PU Events (Reweighted)", plotname PU_init_cor_mu, scale yes, weight w_PU, title "t#bar{t}#rightarrow#mu+jets"

../masterplot.py file ../mu_nominal/TT_8TeV-mcatnlo_selected__reco.root , var mc_pileup_events, bin 50, min 0, max 50, xaxis "Numeber of PU Events (Unweighted)", plotname PU_init_uncor_mu, scale yes, title "t#bar{t}#rightarrow#mu+jets"

# kin reco vs gen 2D plots
../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --varx mtt_mc --vary mtt --binx 30 --biny 30 --Minx 350 --Maxx 1000 --Miny 350 --Maxy 1000 --xaxis "M_{t#bar{t}} (generated)[GeV]" --yaxis "M_{t#bar{t}} (reconstructed)[GeV]" --title "t#bar{t}#rightarrow#mu+jets" --plotname Mtt_2D_mu

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --varx cos_theta_mc --vary cos_theta --binx 30 --biny 30 --Minx -1 --Maxx 1 --Miny -1 --Maxy 1 --xaxis "c* (generated)" --yaxis "c* (reconstructed)" --plotname cstar_2D_mu --cut 'gen_type[0]==\"mu_jets\"&&init_type==\"qqbar\" ' --title "q#bar{q}#rightarrow t#bar{t}#rightarrow#mu+jets"

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --varx 'fabs(xf_mc)' --vary 'fabs(xf)' --binx 30 --biny 30 --Minx 0 --Maxx 0.2 --Miny 0 --Maxy 0.2 --xaxis "|x_{F}| (generated)" --yaxis "|x_{F}| (reconstructed)" --plotname xf_2D_mu --cut --title "t#bar{t}#rightarrow#mu+jets"

# kin reco residue plots
../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --var 'mtt-mtt_mc' --bin 50 --min -500 --max 600 --plotname Mtt_res_mu --xaxis "M_{t#bar{t}}(Reco) - M_{t#bar{t}}(Gen) [GeV]" --title "t#bar{t}#rightarrow#mu+jets" --weight "weight_gen*w_PU*w_btag*w_lepID*w_trigger*weight_top_pT" --scale True

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --var 'cos_theta-cos_theta_mc' --bin 50 --min -2 --max 2 --plotname cstar_res_mu --xaxis "c*(Reco) - c*(Gen) " --title "q#bar{q}#rightarrow t#bar{t}#rightarrow#mu+jets" --weight "weight_gen*w_PU*w_btag*w_lepID*w_trigger*weight_top_pT" --scale True --cut '(gen_type[0]==\"mu_jets\"&&init_type==\"qqbar\") '

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --var 'fabs(xf)-fabs(xf_mc)' --bin 50 --min -1 --max 1 --plotname xf_res_mu --xaxis "|x_{F}|(Reco) - |x_{F}|(Gen) " --title "t#bar{t}#rightarrow#mu+jets" --weight "weight_gen*w_PU*w_btag*w_lepID*w_trigger*weight_top_pT" --scale True

# top_pT SF vs kin 2D plots
../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --varx mtt --vary weight_top_pT --binx 50 --biny 50 --Minx 350 --Maxx 1000 --Miny 0.6 --Maxy 1.25 --xaxis "M_{t#bar{t}} [GeV]" --yaxis "top p_{T} SF" --title "t#bar{t}#rightarrow#mu+jets" --plotname Mtt_SF_mu

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --varx cos_theta --vary weight_top_pT --binx 50 --biny 50 --Minx -1 --Maxx 1 --Miny 0.6 --Maxy 1.25 --xaxis "c*" --yaxis "top p_{T} SF" --plotname cstar_SF_mu --title "t#bar{t}#rightarrow#mu+jets"

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --varx 'fabs(xf_mc)' --vary 'weight_top_pT' --binx 50 --biny 50 --Minx 0 --Maxx 0.2 --Miny 0.6 --Maxy 1.25 --xaxis "|x_{F}|" --yaxis "top p_{T} SF" --plotname xf_SF_mu  --title "t#bar{t}#rightarrow#mu+jets"


