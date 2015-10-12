prepend=ntuples/jhu_diffmo/v1_complete
#echo processing W4JetsToLNu_TuneZ2Star_8TeV 
#python selection.py --inputfiles $prepend/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_all/\*.root  --maxfiles -1 --maxevts -1 --type W4Jets_v1 --makeplots no >> logs/wjets.log
echo processing W3JetsToLNu_TuneZ2Star_8TeV
python selection.py --inputfiles $prepend/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_all/\*.root  --maxfiles -1 --maxevts -1 --type W3Jets_v2 --makeplots no 
#python selection.py --inputfiles ./ntuples/jhu_diffmo/v1_complete/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_all/\*.root --type W4Jets_v2 --maxevts -1 --makeplots no --startfile 30  --maxfiles 30
