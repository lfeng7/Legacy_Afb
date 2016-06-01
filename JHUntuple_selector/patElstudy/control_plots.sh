#! /bin/bash
cut1="el_j_delR>2.0&&(el_j_mass[0]<60||el_j_mass[0]>100)&&(el_j_mass[1]<60||el_j_mass[1]>100)&&met_pt_vec<20&&trigger_vec&&lep_isTight"
cut12="el_j_delR>2.0&&(el_j_mass[0]<40||el_j_mass[0]>120)&&(el_j_mass[1]<40||el_j_mass[1]>120)&&met_pt_vec<20&&trigger_vec&&lep_isTight"

cut21="lep_isLoose==0&&cos(jets_phi-lep_phi)< -0.95&&trigger_vec"
cut22="lep_isLoose==0&&el_j_delR>2.0&&trigger_vec"
cut23="lep_isLoose==0&&trigger_vec"
cut24="lep_isTight==0&&trigger_vec"

qcd_file="ABCD_v2/QCD/QCD_Pt-15to3000_selected.root"

# baseline
./masterplot.py --var lep_iso --Min 0 --Max 1.2 --xaxis RelIso --yaxis Events  --title "Baseline w/ 2 jets > 30GeV and 1 PF-el > 30GeV" --dir ABCD_v2  --stack yes --plotname eliso_baseline --cut "trigger_vec"
# Data QCD w/ substraction
./masterplot.py --var lep_iso --Min 0 --Max 1.2 --xaxis RelIso --yaxis Events  --title "Data QCD Baseline w/ 2 jets > 30GeV and 1 PF-el > 30GeV" --files "ABCD_v2/*.root"  --plotname eliso_baseline_dataQCD --mode plotter --saveroot yes --cut "weight_tmp"
# MC_QCD
./masterplot.py --var lep_iso --Min 0 --Max 1.2 --xaxis RelIso --yaxis Events  --title "MC QCD Baseline w/ 2 jets > 30GeV and 1 PF-el > 30GeV" --files "$qcd_file"  --plotname eliso_baseline_mcQCD --mode plotter --saveroot yes --cut "weight_norm"

# QCD enriched region 1, with good PATel selection
./masterplot.py --var lep_iso --Min 0 --Max 1.2 --xaxis RelIso --yaxis Events  --title "tight elID, MET<20, e-jet delR>2, e-jet mass<60 |>100" --dir ABCD_v2  --stack yes --plotname eliso_cut1 --cut "$cut1"
# Data QCD w/ substraction
./masterplot.py --var lep_iso --Min 0 --Max 1.2 --xaxis RelIso --yaxis Events  --title "Data QCD w/tight elID, MET<20, e-jet delR>2, e-jet mass<60 |>100" --files "ABCD_v2/*.root"  --plotname eliso_cut1_dataQCD --cut "($cut1)*weight_tmp" --mode plotter --saveroot yes 
#<< comment
# v2
./masterplot.py --var lep_iso --Min 0 --Max 1.2 --xaxis RelIso --yaxis Events  --title "tight elID, MET<20, e-jet delR>2, e-jet mass<40 |>
120" --dir ABCD_v2  --stack yes --plotname eliso_cut12 --cut "$cut12"
# Data QCD w/ substraction
./masterplot.py --var lep_iso --Min 0 --Max 1.2 --xaxis RelIso --yaxis Events  --title "Data QCD w/tight elID, MET<20, e-jet delR>2, e-jet
 mass<40 |>120" --files "ABCD_v2/*.root"  --plotname eliso_cut12_dataQCD --cut "($cut12)*weight_tmp" --mode plotter --saveroot yes 

# QCD enriched region 2, with bad (anti Loose ID) PATel
./masterplot.py --var lep_iso --Min 0  --Max 1.2  --xaxis RelIso --yaxis Events --title "2 jets, 1 anti-Loose El, cos(jet_phi-el_phi)< -0.95"  --dir ABCD_v2  --stack yes --plotname eliso_cut21 --cut "$cut21"
./masterplot.py --var lep_iso --Min 0  --Max 1.2  --xaxis RelIso --yaxis Events --title "Data QCD w/ 2 jets, 1 anti-Loose El, cos(jet_phi-el_phi)< -0.95"  --files "ABCD_v2/*.root"  --plotname eliso_cut21_dataQCD --cut "($cut21)*weight_tmp" --mode plotter --saveroot yes
# version2
./masterplot.py --var lep_iso --Min 0  --Max 1.2  --xaxis RelIso --yaxis Events --title "2 jets, 1 anti-Loose El, el_j_delR>2.0"  --dir ABCD_v2  --stack yes --plotname eliso_cut22 --cut "$cut22"
./masterplot.py --var lep_iso --Min 0  --Max 1.2  --xaxis RelIso --yaxis Events --title "Data QCD w/ 2 jets, 1 anti-Loose El, el_j_delR>2.0"  --files "ABCD_v2/*.root"  --plotname eliso_cut22_dataQCD --cut "($cut22)*weight_tmp" --mode plotter --saveroot yes
# version3 
./masterplot.py --var lep_iso --Min 0  --Max 1.2  --xaxis RelIso --yaxis Events --title "2 jets, 1 anti-Loose El"  --dir ABCD_v2  --stack yes --plotname eliso_cut23 --cut "$cut23"
./masterplot.py --var lep_iso --Min 0  --Max 1.2  --xaxis RelIso --yaxis Events --title "Data QCD w/ 2 jets, 1 anti-Loose El"  --files "ABCD_v2/*.root"  --plotname eliso_cut23_dataQCD --cut "($cut23)*weight_tmp" --mode plotter --saveroot yes
# version4
./masterplot.py --var lep_iso --Min 0  --Max 1.2  --xaxis RelIso --yaxis Events --title "2 jets, 1 anti-Tight El"  --dir ABCD_v2  --stack yes --plotname eliso_cut24 --cut "$cut24"
./masterplot.py --var lep_iso --Min 0  --Max 1.2  --xaxis RelIso --yaxis Events --title "Data QCD w/ 2 jets, 1 anti-Tight El"  --files "ABCD_v2/*.root"  --plotname eliso_cut24_dataQCD --cut "($cut24)*weight_tmp" --mode plotter --saveroot yes
# MC-QCD version
./masterplot.py --var lep_iso --Min 0  --Max 1.2  --xaxis RelIso --yaxis Events --title "MC QCD w/ 2 jets, 1 anti-Tight El"  --files "$qcd_file"  --plotname eliso_cut24_mcQCD --cut "($cut24)*weight_norm" --mode plotter --saveroot yes
#comment
