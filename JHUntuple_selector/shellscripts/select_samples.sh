prepend=ntuples/jhu_diffmo/
echo processing DY1JetsToLL_M 
python selection.py --inputfiles $prepend/DY1JetsToLL_M/\*.root  --maxfiles -1 --maxevts -1 --type DY1Jets --makeplots no >> logs/othersamples.log
echo processing DY2JetsToLL_M 
python selection.py --inputfiles $prepend/DY2JetsToLL_M/\*.root  --maxfiles -1 --maxevts -1 --type DY2Jets --makeplots no >> logs/othersa
mples.log
echo processing W1JetsToLNu_TuneZ2Star_8TeV 
python selection.py --inputfiles $prepend/W1JetsToLNu_TuneZ2Star_8TeV/\*.root  --maxfiles -1 --maxevts -1 --type W1Jets --makeplots no >> logs/othersa
mples.log
echo processing W2JetsToLNu_TuneZ2Star_8TeV 
python selection.py --inputfiles $prepend/W2JetsToLNu_TuneZ2Star_8TeV/\*.root  --maxfiles -1 --maxevts -1 --type W2Jets --makeplots no >> logs/othersa
mples.log

