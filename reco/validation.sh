# Validate the resolution improvement by cutting on chi2
python -i ../plotter.py --file angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0_kinfit.root  --save True --name mtt_res_chi_cut35 --var "mtt-mtt_mc" --xaxis "Mtt(GeV) reco-mc"  --cut "" --cut2 "final_chi2<35" --title "qqbar->ttbar MC"  --label "no cut" --log True
python -i ../plotter.py --file angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0_kinfit.root  --save True --name mthad_res_chi_cut --var "reco_mass[1]-gen_mass_v2[1]" --xaxis "thad mass reco-mc"  --cut "" --cut2 "final_chi2<35" --title "TTbar MC"   --log True --Min -200 --Max 1500 --label "no cut"
# Study best chi2 cut
python -i ../plotter.py --file angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0_kinfit.root  --save True --name mthad_res_chi25 --var "reco_mass[1]-gen_mass_v2[1]" --xaxis "thad mass reco-mc"  --cut "final_chi2<25" --cut2 "final_chi2>25" --title "TTbar MC"   --log True --Min -200 --Max 1500
# validate chi2 cut on data
python -i ../plotter.py --file angles_files/SingleEl_Run2012ABCD_reco_angles_0.root  --save True --name mtt_data_chi35 --var "mtt" --xaxis "Mtt(GeV)"  --cut "" --cut2 "final_chi2<35" --title "SingleEl Run2012A-D"  --label "no cut" --bin 40 --Min 300 --Max 1750 --scale True
python -i ../plotter.py --file angles_files/SingleEl_Run2012ABCD_reco_angles_0.root  --save True --name cs_data_chi35 --var "cos_theta" --xaxis "cos#theta*"  --cut "" --cut2 "final_chi2<35" --title "SingleEl Run2012A-D"  --label "no cut" --scale True --Min -1 --Max 1 --bin 20
 python -i ../plotter.py --file angles_files/SingleEl_Run2012ABCD_reco_angles_0.root  --save True --name xf_data_chi35 --var "abs(xf)" --xaxis "abs(xf)"  --cut "" --cut2 "final_chi2<35" --title "SingleEl Run2012A-D"  --label "no cut" --scale True --Min 0 --Max 0.6 --bin 30
python -i ../plotter.py --file angles_files/SingleEl_Run2012ABCD_reco_angles_0.root  --save True --name mthad_data_chi35 --var "reco_mass[1]" --xaxis "m_thad (GeV)"  --cut "" --cut2 "final_chi2<35" --title "SingleEl Run2012A-D"  --label "no cut" --bin 100 --log True --Min 0 --Max 1500
