prepend=ntuples/jhu_diffmo/
echo processing DY4JetsToLNu_TuneZ2Star_8TeV
python selection.py --inputfiles $prepend/DY4JetsToLL_M/\*.root  --maxfiles -1 --maxevts -1 --type DY4Jets_v2 --makeplots no 
echo processing DY3JetsToLNu_TuneZ2Star_8TeV
python selection.py --inputfiles $prepend/DY3JetsToLL_M/\*.root  --maxfiles -1 --maxevts -1 --type DY3Jets_v2 --makeplots no
