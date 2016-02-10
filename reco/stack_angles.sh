python ../stack_plotter.py --var "reco_mass[0]" --Min 100 --Max 250 --name m_tlep_$2 --dir $1 --xaxis "M_{tlep}(GeV)"
python ../stack_plotter.py --var "reco_mass[1]" --Min 50 --Max 500 --name m_thad_$2 --dir $1 --xaxis "M_{thad}(GeV)"
python ../stack_plotter.py --var "reco_mass[2]" --Min 50 --Max 120 --name m_wlep_$2 --dir $1 --xaxis "M_{wlep}(GeV)"
python ../stack_plotter.py --var "reco_mass[3]" --Min 20 --Max 300 --name m_whad_$2 --dir $1 --xaxis "M_{whad}(GeV)"

