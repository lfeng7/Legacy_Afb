#! /bin/bash 
python plotter.py --name "cos_res_cut_${2:-}" --var "abs(cos_theta_cs)-abs(cos_theta_mc)" --xaxis "|cos#theta*| - |cos#theta*_{MC}|" --Min -1 --Max 1 --cut '' --cut2 "$1"  --label 'no cut' --label2 "cut" --file template_files/pscale_removed_v2/TT_CT10_signal_TuneZ2star_8TeV_reco_no_p_scale_angles_0_new_template.root --scale True --title "TTbar signal ${2:-}"	
python plotter.py --name "mtt_res_cut_${2:-}" --var "ttbar_mass-ttbar_mass_mc" --xaxis "Mtt-Mtt_{MC}(GeV)" --Min -300  --Max 600 --cut '' --cut2 "$1"  --label 'no cut' --label2 "cut" --file template_files/pscale_removed_v2/TT_CT10_signal_TuneZ2star_8TeV_reco_no_p_scale_angles_0_new_template.root --scale True --title "TTbar signal ${2:-}"
python plotter.py --name "xf_res_cut_${2:-}" --var "Feynman_x-Feynman_x_mc" --xaxis "xf-xf_mc" --Min -0.2 --Max 0.2 --cut '' --cut2 "$1"  --label 'no cut' --label2 "cut" --file template_files/pscale_removed_v2/TT_CT10_signal_TuneZ2star_8TeV_reco_no_p_scale_angles_0_new_template.root --scale True --title "TTbar signal ${2:-}"


