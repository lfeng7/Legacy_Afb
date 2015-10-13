# -*- sh -*- # for font lock mode
# variable definitions
- env = cd /uscms_data/d3/lfeng7/Legacy_Afb/analysis/CMSSW_5_3_24; eval `scramv1 runtime -sh`; cd -
- tag = 
- output = outputFile=
- tagmode = none
- tarfile = /uscms_data/d3/lfeng7/Legacy_Afb/analysis/CMSSW_5_3_24/work/selection/JHUntuple_selector/condorjobs/tarball.tgz
- untardir = tardir
- copycommand = cp

# Sections listed
output_$(JID)        python ./tardir/selection.py --inputfiles ntuples/jhu_diffmo/v1_complete/SingleEl_Run2012A_v1_all/\*.root --makeplots no --mcordata data --maxevts -1 --type condor_test --maxfiles 2 --startfile 0
output_$(JID)        python ./tardir/selection.py --inputfiles ntuples/jhu_diffmo/v1_complete/SingleEl_Run2012A_v1_all/\*.root --makeplots no --mcordata data --maxevts -1 --type condor_test --maxfiles 2 --startfile 2

