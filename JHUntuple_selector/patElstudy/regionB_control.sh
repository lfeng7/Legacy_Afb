# ele-jet deltaR for 2jets, good ele region
./masterplot.py --var el_j_delR --Min 0 --Max 5 --xaxis "DelR(e-jets)" --yaxis Events --title "delR 2jets good el" --dir ABCD_v4 --stack yes --plotname study_delR_goodEl_2jets --cut "njets==2&&(el_j_mass[0]<60||el_j_mass[0]>100)&&(el_j_mass[1]<60||el_j_mass[1]>100)&&(el_j_delR[0]>2||el_j_delR[1]>2) && lep_iso<0.1 && lep_isTight"
./masterplot.py --var "el_j_delR[0]" --Min 0 --Max 5 --xaxis "DelR(e-jets0)" --yaxis Events --title "delR 2jets good el delR1>2" --dir ABCD_v4 --stack yes --plotname study_delR0_goodEl_2jets --cut "njets==2&&(el_j_mass[0]<60||el_j_mass[0]>100)&&(el_j_mass[1]<60||el_j_mass[1]>100)&&(el_j_delR[1]>2) && lep_iso<0.1 && lep_isTight"
./masterplot.py --var "el_j_delR[1]" --Min 0 --Max 5 --xaxis "DelR(e-jets1)" --yaxis Events --title "delR 2jets good el delR1>2" --dir ABCD_v4 --stack yes --plotname study_delR1_goodEl_2jets --cut "njets==2&&(el_j_mass[0]<60||el_j_mass[0]>100)&&(el_j_mass[1]<60||el_j_mass[1]>100)&&(el_j_delR[0]>2) && lep_iso<0.1 && lep_isTight"

# cos(delPhi) between ele and MET
./masterplot.py --var "cos(lep_phi-met_phi_vec)" --Min -1 --Max 1 --xaxis "cos ( DelPhi_e_MET ) " --yaxis Events --title "2jets good el" --dir ABCD_v4 --stack yes --plotname study_cosDelPhi_e_MET_goodEl_2jets --cut "njets==2 && ( el_j_mass[0] < 60 || el_j_mass[0] > 100 ) && ( el_j_mass[1] < 60 || el_j_mass[1] > 100 ) && ( el_j_delR[0] > 2 || el_j_delR[1] > 2 ) && lep_iso < 0.1 && lep_isTight"
