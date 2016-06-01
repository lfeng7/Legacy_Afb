# compare the shape of data and QCD MC
#! /bin/bash
cut0="trigger_vec"

cut11="el_j_delR>2.0&&(el_j_mass[0]<60||el_j_mass[0]>100)&&(el_j_mass[1]<60||el_j_mass[1]>100)&&met_pt_vec<20&&trigger_vec&&lep_isTight"
cut12="el_j_delR>2.0&&(el_j_mass[0]<40||el_j_mass[0]>120)&&(el_j_mass[1]<40||el_j_mass[1]>120)&&met_pt_vec<20&&trigger_vec&&lep_isTight"

cut21="lep_isLoose==0&&cos(jets_phi-lep_phi)< -0.95&&trigger_vec"
cut22="lep_isLoose==0&&el_j_delR>2.0"
cut23="lep_isLoose==0"
cut24="lep_isTight==0"

data_file="ABCD_v2/SingleEl_Run2012A_selected.root"
qcd_file="ABCD_v2/QCD/QCD_Pt-15to3000_selected.root"

# iso
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "Baseline w/ 2jets" --plotname lep_iso_baseline_shape --var lep_iso --scale yes --label1 "Data" --label2 "MC QCD" --min 0 --max 1.2  --bin 30
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "Baseline w/ 2jets and pass trigger " --plotname lep_iso_trigger_shape --var lep_iso --scale yes --label1 "Data" --label2 "MC QCD" --min 0 --max 1.2 --cut "$cut0" --bin 30
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "$cut21" --plotname lep_iso_cut21_shape --var lep_iso --scale yes --label1 "Data" --label2 "MC QCD" --min 0 --max 1.2  --cut "$cut21" --bin 30

# lep_pt
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "Baseline w/ 2jets" --plotname lep_pt_baseline_shape --var lep_pt --scale yes --label1 "Data" --label2 "MC QCD" --min 25 --max 200 --bin 30
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "Baseline w/ 2jets and pass trigger" --plotname lep_pt_trigger_shape --var lep_pt --scale yes --label1 "Data" --label2 "MC QCD" --min 25 --max 200 --cut "$cut0" --bin 30
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "$cut21" --plotname lep_pt_cut21_shape --var lep_pt --scale yes --label1 "Data" --label2 "MC QCD" --min 25 --max 200  --cut "$cut21" --bin 30

# jets_pt
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "Baseline w/ 2jets" --plotname jets_pt_baseline_shape --var jets_pt --scale yes --label1 "Data" --label2 "MC QCD" --min 25 --max 300 --bin 30
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "Baseline w/ 2jets and pass trigger" --plotname jets_pt_trigger_shape --var jets_pt --scale yes --label1 "Data" --label2 "MC QCD" --min 25 --max 300 --cut "$cut0" --bin 30
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "$cut21" --plotname jets_pt_cut21_shape --var jets_pt --scale yes --label1 "Data" --label2 "MC QCD" --min 25 --max 300  --cut "$cut21" --bin 30

# met_pt_vec
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "Baseline w/ 2jets" --plotname met_pt_vec_baseline_shape --var met_pt_vec --scale yes --label1 "Data" --label2 "MC QCD" --min 0 --max 200 --bin 30
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "Baseline w/ 2jets and pass trigger" --plotname met_pt_vec_trigger_shape --var met_pt_vec --scale yes --label1 "Data" --label2 "MC QCD" --min 0 --max 200 --cut "$cut0" --bin 30
./masterplot.py --file1 "$data_file"  --file2 "$qcd_file" --title "$cut21" --plotname met_pt_vec_cut21_shape --var met_pt_vec --scale yes --label1 "Data" --label2 "MC QCD" --min 0 --max 200  --cut "$cut21" --bin 30



