prepend=ntuples/jhu_diffmo/v1_complete
echo processing W4jets 

#python selection.py --inputfiles $prepend/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_all/jhutester_numEvent1000_10.root  --maxfiles 1 --maxevts -1 --type W4Jets_v1 --makeplots no >> logs/missedfiles.log
#python selection.py --inputfiles $prepend/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_all/jhutester_numEvent1000_23.root  --maxfiles 1 --maxevts -1 --type W4Jets_v1 --makeplots no >> logs/missedfiles.log
#python selection.py --inputfiles $prepend/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_all/jhutester_numEvent1000_6.root  --maxfiles 1 --maxevts -1 --type W4Jets_v1 --makeplots no >> logs/missedfiles.log

echo processing data
python selection.py --inputfiles $prepend/SingleEl_Run2012A_v1_all/jhutester_numEvent1000_133.root  --maxfiles 1 --maxevts -1 --type  SingleEl_Run2012A_all --makeplots no --mcordata data >> logs/missedfiles.log
python selection.py --inputfiles $prepend/SingleEl_Run2012A_v1_all/jhutester_numEvent1000_17.root  --maxfiles 1 --maxevts -1 --type  SingleEl_Run2012A_all --makeplots no  --mcordata data >> logs/missedfiles.log
python selection.py --inputfiles $prepend/SingleEl_Run2012A_v1_all/jhutester_numEvent1000_186.root  --maxfiles 1 --maxevts -1 --type  SingleEl_Run2012A_all --makeplots no   --mcordata data >> logs/missedfiles.log
python selection.py --inputfiles $prepend/SingleEl_Run2012A_v1_all/jhutester_numEvent1000_98.root  --maxfiles 1 --maxevts -1 --type  SingleEl_Run2012A_all --makeplots no  --mcordata data  >> logs/missedfiles.log
python selection.py --inputfiles $prepend/SingleEl_Run2012A_v1_all/jhutester_numEvent1000_99.root  --maxfiles 1 --maxevts -1 --type  SingleEl_Run2012A_all --makeplots no  --mcordata data  >> logs/missedfiles.log
