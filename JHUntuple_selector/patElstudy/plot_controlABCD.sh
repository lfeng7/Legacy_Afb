./masterplot.py --file ABCD_v4/DY1JetsToLL_M_selected.root ABCD_v4/DY2JetsToLL_M_selected.root ABCD_v4/DY3JetsToLL_M_selected.root ABCD_v4/DY4JetsToLL_M_selected.root --title "DYJets njets>=4&&n_btags==2&&lep_isTight" --plotname iso_met_signal_cutDYJets --cut  "(njets>=4&&n_btags==2&&lep_isTight)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/DY1JetsToLL_M_selected.root ABCD_v4/DY2JetsToLL_M_selected.root ABCD_v4/DY3JetsToLL_M_selected.root ABCD_v4/DY4JetsToLL_M_selected.root --title "DYJets lep_isLoose" --plotname iso_met_loose_cutDYJets --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/DY1JetsToLL_M_selected.root ABCD_v4/DY2JetsToLL_M_selected.root ABCD_v4/DY3JetsToLL_M_selected.root ABCD_v4/DY4JetsToLL_M_selected.root --title "DYJets lep_isLoose" --plotname iso_transIP_loose_cutDYJets --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx TransverseIP --binx 20 --biny 20 --Minx 0 --Maxx 0.1 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "TransverseIP" 

./masterplot.py --file ABCD_v4/W1JetsToLNu_TuneZ2Star_8TeV_selected.root ABCD_v4/W2JetsToLNu_TuneZ2Star_8TeV_selected.root ABCD_v4/W3JetsToLNu_TuneZ2Star_8TeV_selected.root ABCD_v4/W4JetsToLNu_TuneZ2Star_8TeV_selected.root --title "WJets njets>=4&&n_btags==2&&lep_isTight" --plotname iso_met_signal_cutWJets --cut  "(njets>=4&&n_btags==2&&lep_isTight)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/W1JetsToLNu_TuneZ2Star_8TeV_selected.root ABCD_v4/W2JetsToLNu_TuneZ2Star_8TeV_selected.root ABCD_v4/W3JetsToLNu_TuneZ2Star_8TeV_selected.root ABCD_v4/W4JetsToLNu_TuneZ2Star_8TeV_selected.root --title "WJets lep_isLoose" --plotname iso_met_loose_cutWJets --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/W1JetsToLNu_TuneZ2Star_8TeV_selected.root ABCD_v4/W2JetsToLNu_TuneZ2Star_8TeV_selected.root ABCD_v4/W3JetsToLNu_TuneZ2Star_8TeV_selected.root ABCD_v4/W4JetsToLNu_TuneZ2Star_8TeV_selected.root --title "WJets lep_isLoose" --plotname iso_transIP_loose_cutWJets --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx TransverseIP --binx 20 --biny 20 --Minx 0 --Maxx 0.1 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "TransverseIP" 

./masterplot.py --file ABCD_v4/TT_CT10_TuneZ2star_8TeV_selected.root --title "TTbar njets>=4&&n_btags==2&&lep_isTight" --plotname iso_met_signal_cutTTbar --cut  "(njets>=4&&n_btags==2&&lep_isTight)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/TT_CT10_TuneZ2star_8TeV_selected.root --title "TTbar lep_isLoose" --plotname iso_met_loose_cutTTbar --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/TT_CT10_TuneZ2star_8TeV_selected.root --title "TTbar lep_isLoose" --plotname iso_transIP_loose_cutTTbar --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx TransverseIP --binx 20 --biny 20 --Minx 0 --Maxx 0.1 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "TransverseIP" 

./masterplot.py --file ABCD_v4/T_s_selected.root ABCD_v4/T_t_selected.root ABCD_v4/T_tW_selected.root ABCD_v4/Tbar_s_selected.root ABCD_v4/Tbar_t_selected.root ABCD_v4/Tbar_tW_selected.root --title "SingleTop njets>=4&&n_btags==2&&lep_isTight" --plotname iso_met_signal_cutSingleTop --cut  "(njets>=4&&n_btags==2&&lep_isTight)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/T_s_selected.root ABCD_v4/T_t_selected.root ABCD_v4/T_tW_selected.root ABCD_v4/Tbar_s_selected.root ABCD_v4/Tbar_t_selected.root ABCD_v4/Tbar_tW_selected.root --title "SingleTop lep_isLoose" --plotname iso_met_loose_cutSingleTop --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/T_s_selected.root ABCD_v4/T_t_selected.root ABCD_v4/T_tW_selected.root ABCD_v4/Tbar_s_selected.root ABCD_v4/Tbar_t_selected.root ABCD_v4/Tbar_tW_selected.root --title "SingleTop lep_isLoose" --plotname iso_transIP_loose_cutSingleTop --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx TransverseIP --binx 20 --biny 20 --Minx 0 --Maxx 0.1 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "TransverseIP" 

./masterplot.py --file ABCD_v4/QCD_Pt-15to3000_selected.root --title "QCD njets>=4&&n_btags==2&&lep_isTight" --plotname iso_met_signal_cutQCD --cut  "(njets>=4&&n_btags==2&&lep_isTight)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/QCD_Pt-15to3000_selected.root --title "QCD lep_isLoose" --plotname iso_met_loose_cutQCD --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/QCD_Pt-15to3000_selected.root --title "QCD lep_isLoose" --plotname iso_transIP_loose_cutQCD --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx TransverseIP --binx 20 --biny 20 --Minx 0 --Maxx 0.1 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "TransverseIP" 

./masterplot.py --file ABCD_v4/SingleEl_Run2012B_selected.root ABCD_v4/SingleEl_Run2012C_part1_selected.root ABCD_v4/SingleEl_Run2012C_part2_selected.root ABCD_v4/SingleEl_Run2012D_selected.root --title "Data njets>=4&&n_btags==2&&lep_isTight" --plotname iso_met_signal_cutData --cut  "(njets>=4&&n_btags==2&&lep_isTight)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/SingleEl_Run2012B_selected.root ABCD_v4/SingleEl_Run2012C_part1_selected.root ABCD_v4/SingleEl_Run2012C_part2_selected.root ABCD_v4/SingleEl_Run2012D_selected.root --title "Data lep_isLoose" --plotname iso_met_loose_cutData --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" 

./masterplot.py --file ABCD_v4/SingleEl_Run2012B_selected.root ABCD_v4/SingleEl_Run2012C_part1_selected.root ABCD_v4/SingleEl_Run2012C_part2_selected.root ABCD_v4/SingleEl_Run2012D_selected.root --title "Data lep_isLoose" --plotname iso_transIP_loose_cutData --cut  "(lep_isLoose)*(weight_norm)"   --vary lep_iso --varx TransverseIP --binx 20 --biny 20 --Minx 0 --Maxx 0.1 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "TransverseIP" 

