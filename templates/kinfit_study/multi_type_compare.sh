#! /bin/bash
#python multi_sample_plotter_v2.py --var ttbar_mass --weight correction_weight --Min 300 --Max 1600 --name "Mtt_$2" --bin 100 --dir "$1" --xaxis "M_{tt}(GeV)" --scale True --title "Mtt $2"

# tlep Mass
python multi_sample_plotter_v2.py --dir "$1"  --var "reco_mass[0]" --Min 50 --Max 250 --name "mtlep_compare_$2" --scale True --bin 100 --title "m_tlep $2" --xaxis "M_{tlep}(GeV)"  --log ${3:-'no'} --weight correction_weight
# thad mass
python multi_sample_plotter_v2.py --dir "$1" --var "reco_mass[1]" --Min 50 --Max 500 --name "mthad_compare_$2" --scale True --bin 100 --title "m_thad $2" --xaxis "M_{thad}(GeV)"  --log ${3:-'no'} --weight correction_weight
# wlep mass
python multi_sample_plotter_v2.py --dir "$1"   --var "reco_mass[2]" --Min 70 --Max 200 --name "mwlep_compare_$2" --scale True --bin 100 --title "m_wlep $2" --xaxis "M_{wlep}(GeV)"  --log ${3:-'no'} --weight correction_weight
# whad mass
python multi_sample_plotter_v2.py --dir "$1"   --var "reco_mass[3]" --Min 20 --Max 300 --name "mwhad_compare_$2" --scale True --bin 100 --title "m_whad $2" --xaxis "M_{whad}(GeV)"   --log ${3:-'no'} --weight correction_weight
# chi2
python multi_sample_plotter_v2.py --dir "$1"   --var "lnL" --Min -20 --Max 100 --name "chi2_compare_$2" --scale True --bin 100 --title "chi2 $2" --xaxis "#Chi2"  --weight correction_weight
# N_jets
python multi_sample_plotter_v2.py --dir "$1"   --var "n_valid_jets" --Min 3 --Max 6 --name "njets_compare_$2" --bin 10 --scale True --title "N_jets $2" --xaxis "N_{j}"  --weight correction_weight
# leading jet pT
python multi_sample_plotter_v2.py --dir "$1"   --var "leadingJet_pt" --Min 20 --Max 100 --name "leadingJetPt_compare_$2" --scale True --bin 100 --title "leadingJet_pT $2" --xaxis "leading_jet pT (GeV)" --weight correction_weight  
# leading jet mass
python multi_sample_plotter_v2.py --dir "$1"   --var "leadingJet_mass" --Min 0 --Max 25 --name "leadingJetMass_compare_$2" --scale True --bin 100 --title "leadingJet_mass $2" --xaxis "leading_jet mass (GeV)"  --weight correction_weight
# ttbar mass
python multi_sample_plotter_v2.py --dir "$1"   --var "ttbar_mass" --Min 300 --Max 1500 --name "mtt_compare_$2" --scale True --bin 100 --title "ttbar_mass $2" --xaxis "ttbar mass (GeV)" --weight correction_weight
