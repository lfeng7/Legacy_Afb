prepend=ntuples/jhu_diffmo/v1_complete

python selection.py --inputfiles ntuples/jhu_diffmo/v1_complete/W3JetsToLNu_TuneZ2Star_8TeV-madgraph_all/\*.root  --maxfiles 30 --maxevts -1 --type W3Jets_v2 --makeplots no --startfile 0 
# python selection.py --inputfiles ./ntuples/jhu_diffmo/v1_complete/W4JetsToLNu_TuneZ2Star_8TeV-madgraph_all/\*.root --type W4Jets_v2 --maxevts -1 --makeplots no  --maxfiles 20 --startfile 0
