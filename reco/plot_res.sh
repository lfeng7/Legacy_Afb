#!/bin/bash
# source plot_res.sh "kinfit_study/*angles*" "no_pscale pscale_constrained small_res nominal"

# resolution plots comparing different kinfit set up

python  ../multi_sample_plotter.py --file "$1" --label "$2" --save True --name mtt_res --var "mtt-mtt_mc" --xaxis "Mtt-Mtt_{MC}(GeV)" --Min -300  --Max 600  --title "TTbar MC" --scale True
python  ../multi_sample_plotter.py --file "$1" --label "$2"  --save True --name xf_res --var "xf-xf_mc" --xaxis "xf-xf_mc" --Min -0.2 --Max 0.2  --title "TTbar MC" --scale True
python  ../multi_sample_plotter.py --file "$1" --label "$2"  --save True --name cos_res --var "cos_theta-cos_theta_mc" --xaxis "cos#theta* - (cos#theta*)_{MC}" --Min -2 --Max 2 --title "TTbar MC" --scale True


# resolution plots comparing different chi2 cut

#python -i ../plotter.py --file fname  --save True --name mtt_res_new --var "mtt/mtt_mc" --xaxis "Mtt/Mtt_{MC}" --Min 0.6  --Max 2.2 --cut "final_chi2<25" --cut2 "final_chi2>25" --title "TTbar MC" --scale True  
#python -i plotter.py --file reco/angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root  --save True --name xf_res_new --var "xf/xf_mc" --xaxis "xf/xf_{MC}" --Min -1 --Max 3 --cut "final_chi2<25" --cut2 "final_chi2>25" --title TTbar --scale True
#python -i plotter.py --file reco/angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root  --save True --name cos_res --var "cos_theta/cos_theta_mc" --xaxis "cos#theta*/(cos#theta*)_{MC}" --Min -3 --Max 3 --cut "final_chi2<25" --cut2 "final_chi2>25" --title "TTbar MC" --scale True 


# mass of tlep,thad,wlep,whad
python  ../multi_sample_plotter.py --file "$1" --label "$2"  --save True --name m_tlep --var "reco_mass[0]" --xaxis "M_{tlep}(GeV)" --Min 100 --Max 250 --title "TTbar MC" --scale True --bin 30 
python  ../multi_sample_plotter.py --file "$1" --label "$2"  --save True --name m_thad --var "reco_mass[1]" --xaxis "M_{thad}(GeV)" --Min 100 --Max 300 --title "TTbar MC" --scale True --bin 30
python  ../multi_sample_plotter.py --file "$1" --label "$2"  --save True --name m_wlep --var "reco_mass[2]" --xaxis "M_{Wlep}(GeV)" --Min 50 --Max 150 --title "TTbar MC" --scale True  --bin 30
python  ../multi_sample_plotter.py --file "$1" --label "$2"  --save True --name m_whad --var "reco_mass[3]" --xaxis "M_{Whad}(GeV)" --Min 20 --Max 200 --title "TTbar MC" --scale True  --bin 30


# Study cos_theta,Mtt,xf by cutting on chi2
#python -i plotter.py --file reco/angles_files/SingleEl_Run2012ABCD_reco_angles_0.root  --save True --name mtt_chi2_data --var "mtt" --xaxis "Mtt" --Min 300  --Max 1500 --cut "final_chi2<25" --cut2 "final_chi2>25" --title "Data" --scale True
#python -i plotter.py --file reco/angles_files/SingleEl_Run2012ABCD_reco_angles_0.root  --save True --name cos_chi2_data --var "cos_theta" --xaxis "cos#theta*" --Min -1  --Max 1 --cut "final_chi2<25" --cut2 "final_chi2>25" --title "Data" --scale True

