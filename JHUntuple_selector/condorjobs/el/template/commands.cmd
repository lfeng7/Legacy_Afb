# -*- sh -*- # for font lock mode
# variable definitions
- env = cd /uscms_data/d3/lfeng7/B2G_FW/CMSSW_7_2_0; eval `scramv1 runtime -sh`; cd -
- tag = 
- output = outputFile=
- tagmode = none
- tarfile = /uscms_data/d3/lfeng7/B2G_FW/CMSSW_7_2_0/work/Legacy_Afb/JHUntuple_selector/condorjobs/template/tarball.tgz
- untardir = tardir
- copycommand = cp

# Sections listed
output_$(JID)        python ./tardir/selection.py --txtfiles tardir/inputfiles/SingleEl_Run2012A_v1.txt --makeplots no --mcordata data --maxevts -1 --type condor_test --grid yes --maxfiles 2 --startfile 0
output_$(JID)        python ./tardir/selection.py --txtfiles tardir/inputfiles/SingleEl_Run2012A_v1.txt --makeplots no --mcordata data --maxevts -1 --type condor_test --grid yes --maxfiles 2 --startfile 2

