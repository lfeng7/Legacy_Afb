python -i ../stack_plotter.py --yields yes --var cos_theta_cs --weight correction_weight  --Min -1 --Max 1 --name corrected --bin 30 --dir temp_angles --xaxis "cos#theta*"
python -i  ../stack_plotter.py --var ttbar_mass --weight correction_weight  --Min 270 --Max 1700 --name corrected --bin 50 --dir temp_angles --xaxis "M_{tt}(GeV)"
python -i  ../stack_plotter.py --var "abs(Feynman_x)" --weight correction_weight  --Min 0 --Max 0.6 --name corrected --bin 40 --dir temp_angles --xaxis "|x_{F}|"
python -i ../stack_plotter.py --var charge_ratio --weight correction_weight  --Min 0 --Max 5  --name corrected --bin 10 --dir temp_angles --xaxis "charge_ratio" 
