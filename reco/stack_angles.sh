python  ../stack_plotter.py  --var cos_theta   --Min -1 --Max 1 --name sideband --bin 20 --dir angles_files/sideband/ --xaxis "cos#theta*"
python   ../stack_plotter.py --var ttbar_mass   --Min 300 --Max 1750 --name sideband --bin 40 --dir angles_files/sideband/ --xaxis "M_{tt}(GeV)"
python   ../stack_plotter.py --var "abs(x_f)"   --Min 0 --Max 0.6 --name sideband --bin 30 --dir angles_files/sideband/ --xaxis "|x_{F}|"
python   ../stack_plotter.py --var "final_chi2"   --Min -30 --Max 100 --name sideband --bin 40 --dir angles_files/sideband/ --xaxis "kinfit -2lnL"
python ../stack_plotter.py --var "reco_mass[0]" --Min 100 --Max 250 --name m_tlep_$2 --dir $1 --xaxis "M_{tlep}(GeV)"
python ../stack_plotter.py --var "reco_mass[1]" --Min 50 --Max 500 --name m_thad_$2 --dir $1 --xaxis "M_{thad}(GeV)"
python ../stack_plotter.py --var "reco_mass[2]" --Min 50 --Max 120 --name m_wlep_$2 --dir $1 --xaxis "M_{wlep}(GeV)"
python ../stack_plotter.py --var "reco_mass[3]" --Min 20 --Max 300 --name m_whad_$2 --dir $1 --xaxis "M_{whad}(GeV)"

