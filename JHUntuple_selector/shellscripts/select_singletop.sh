prepend=ntuples/jhu_diffmo/v1_complete
echo processing T_s-channel_TuneZ2star_8TeV-powheg-tauola_all
python selection.py --inputfiles $prepend/T_s-channel_TuneZ2star_8TeV-powheg-tauola_all/\*.root  --maxfiles -1 --maxevts -1 --type T_s_v2 --makeplots no 
echo processing T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_all
python selection.py --inputfiles $prepend/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_all/\*.root  --maxfiles -1 --maxevts -1 --type T_tW_v2  --makeplots no # >> logs/singleTop.log
echo processing Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_all
python selection.py --inputfiles $prepend/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_all/\*.root  --maxfiles -1 --maxevts -1 --type Tbar_s_v2  --makeplots no # >> logs/singleTop.log
echo processing Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_all
python selection.py --inputfiles $prepend/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_all/\*.root  --maxfiles -1 --maxevts -1 --type Tbar_tW_v2 --makeplots no  # >> logs/singleTop.log

# python selection.py --inputfiles ntuples/jhu_diffmo/v1_complete/T_t-channel_TuneZ2star_8TeV-powheg-tauola_all/\*.root  --maxfiles 10 --maxevts -1 --type T_t_v2 --makeplots no --startfile 0 
# python selection.py --inputfiles ntuples/jhu_diffmo/v1_complete/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_all/\*.root  --maxfiles 10 --maxevts -1 --type Tbar_t_v2 --makeplots no --startfile 0

