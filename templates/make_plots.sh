#! /bin/bash
if ["$1" eq ""];then
    echo make template fitting three vars plots
    echo usage:  ./make_plots.sh template_files/sideband/ sideband
    exit 1
fi
python  ../stack_plotter.py --yields yes --var cos_theta_cs --weight correction_weight  --Min -1 --Max 1 --name $2_cos_theta_cs --bin 20 --dir $1 --xaxis "cos#theta*"
python   ../stack_plotter.py --var ttbar_mass --weight correction_weight  --Min 350 --Max 1750 --name $2_ttbar_mass --bin 40 --dir $1 --xaxis "M_{tt}(GeV)"
python   ../stack_plotter.py --var "abs(Feynman_x)" --weight correction_weight  --Min 0 --Max 0.6 --name $2_xf --bin 30 --dir $1 --xaxis "|x_{F}|"
python  ../stack_plotter.py --var charge_ratio --weight correction_weight  --Min 0 --Max 6  --name $2_charge_ratio --bin 6 --dir $1 --xaxis "charge_ratio"
python   ../stack_plotter.py --var "lnL" --weight correction_weight  --Min -30 --Max 100 --name $2_lnL --bin 40 --dir $1 --xaxis "kinfit -2lnL"
