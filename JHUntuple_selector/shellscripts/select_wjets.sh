prepend=ntuples/jhu_diffmo/v1_complete
echo processing W4JetsToLNu_TuneZ2Star_8TeV 
python selection.py --inputfiles $prepend/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_all/\*.root  --maxfiles -1 --maxevts -1 --type W4Jets_v1 --makeplots no >> logs/wjets.log
echo processing W3JetsToLNu_TuneZ2Star_8TeV
python selection.py --inputfiles $prepend/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_all/\*.root  --maxfiles -1 --maxevts -1 --type W3Jets_v1 --makeplots no >> logs/wjets.log
