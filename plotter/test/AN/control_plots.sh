# control plots for selection
python ../masterplot.py dir ../mu_nominal/, var1 lep_pt, bin 50, min 20, max 250, xaxis "p_{T}^{lep} (GeV)", title "mu lep_pt", plotname lep_pt_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 met_pt, bin 50, min 0, max 300, xaxis "MET (GeV)", plotname MET_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 jets_pt, bin 50, min 30, max 300, xaxis "p_{T}^{Jets} (GeV)", plotname jets_pt_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

python ../masterplot.py dir ../mu_nominal/, var1 N_btag, bin 5, min 0, max 5, xaxis "Num (b-jets)", plotname nbtag_mu, stack yes, mcinfo ../payloads/MC_input_muon_mcNLO.txt, weight "weight_gen*weight_top_pT"

# kin reco vs gen 2D plots
../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --varx mtt_mc --vary mtt --binx 30 --biny 30 --Minx 350 --Maxx 1000 --Miny 350 --Maxy 1000 --xaxis "M_{t#bar{t}} (generated)[GeV]" --yaxis "M_{t#bar{t}} (reconstructed)[GeV]" --title "t#bar{t}#rightarrow#mu+jets" --plotname Mtt_2D

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --varx cos_theta_mc --vary cos_theta --binx 30 --biny 30 --Minx -1 --Maxx 1 --Miny -1 --Maxy 1 --xaxis "c* (generated)" --yaxis "c* (reconstructed)" --plotname cstar_2D --cut 'gen_type[0]==\"mu_jets\"&&init_type==\"qqbar\" ' --title "q#bar{q}#rightarrow t#bar{t}#rightarrow#mu+jets"

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --varx 'fabs(xf_mc)' --vary 'fabs(xf)' --binx 30 --biny 30 --Minx 0 --Maxx 0.2 --Miny 0 --Maxy 0.2 --xaxis "|x_{F}| (generated)" --yaxis "|x_{F}| (reconstructed)" --plotname xf_2D_mu --cut --title "t#bar{t}#rightarrow#mu+jets"

# kin reco residue plots
../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --var 'mtt-mtt_mc' --bin 50 --min -500 --max 600 --plotname Mtt_res --xaxis "M_{t#bar{t}}(Reco) - M_{t#bar{t}}(Gen) [GeV]" --title "t#bar{t}#rightarrow#mu+jets" --weight "weight_gen*w_PU*w_btag*w_lepID*w_trigger*weight_top_pT" --scale True

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --var 'cos_theta-cos_theta_mc' --bin 50 --min -2 --max 2 --plotname cstar_res --xaxis "c*(Reco) - c*(Gen) " --title "q#bar{q}#rightarrow t#bar{t}#rightarrow#mu+jets" --weight "weight_gen*w_PU*w_btag*w_lepID*w_trigger*weight_top_pT" --scale True --cut '(gen_type[0]==\"mu_jets\"&&init_type==\"qqbar\") '

../masterplot.py --file ../angles_files/mu/signal/nominal/TT_8TeV-mcatnlo_reco_angles_0.root --var 'fabs(xf)-fabs(xf_mc)' --bin 50 --min -1 --max 1 --plotname xf_res_mu --xaxis "|x_{F}|(Reco) - |x_{F}|(Gen) " --title "t#bar{t}#rightarrow#mu+jets" --weight "weight_gen*w_PU*w_btag*w_lepID*w_trigger*weight_top_pT" --scale True
