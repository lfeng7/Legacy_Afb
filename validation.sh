# resolution plots 
python -i plotter.py --file reco/angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root  --save True --name mtt_res_new --var "mtt/mtt_mc" --xaxis "Mtt/Mtt_{MC}" --Min 0.6  --Max 2.2 --cut "final_chi2<25" --cut2 "final_chi2>25" --title TTbar MC --scale True
python -i plotter.py --file reco/angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root  --save True --name mtt_res --var "mtt-mtt_mc" --xaxis "Mtt-Mtt_{MC}(GeV)" --Min -300  --Max 600 --cut "final_chi2<25" --cut2 "final_chi2>25" --title TTbar MC --scale True 

python -i plotter.py --file reco/angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root  --save True --name xf_res --var "xf-xf_mc" --xaxis "xf-xf_mc" --Min -0.2 --Max 0.2 --cut "final_chi2<25" --cut2 "final_chi2>25" --title TTbar MC --scale True
python -i plotter.py --file reco/angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root  --save True --name xf_res_new --var "xf/xf_mc" --xaxis "xf/xf_{MC}" --Min -1 --Max 3 --cut "final_chi2<25" --cut2 "final_chi2>25" --title TTbar --scale True

python -i plotter.py --file reco/angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root  --save True --name cos_res_new --var "cos_theta-cos_theta_mc" --xaxis "cos#theta* - (cos#theta*)_{MC}" --Min -2 --Max 2 --cut "final_chi2<25" --cut2 "final_chi2>25" --title TTbar MC --scale True
python -i plotter.py --file reco/angles_files/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root  --save True --name cos_res --var "cos_theta/cos_theta_mc" --xaxis "cos#theta*/(cos#theta*)_{MC}" --Min -3 --Max 3 --cut "final_chi2<25" --cut2 "final_chi2>25" --title TTbar MC --scale True 


# two vars plots
python -i double_var_plotter.py --file templates/TT_CT10_TuneZ2star_8TeV_reco_angles_0_template_qq_0.root --var1 beta_v0 --var2 beta_v1 --cut "beta_v1>-1" --name beta_qqbar --Min 0 --Max 1 
