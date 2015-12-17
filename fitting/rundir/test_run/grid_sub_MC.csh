#! /bin/sh
tar czvf tarball.tgz ../../../*.*  *.txt
/uscms_data/d3/eminizer/runManySections/runManySections.py --createCommandFile --cmssw --addLog --setTarball=tarball.tgz \ana.listOfJobs commands.cmd
/uscms_data/d3/eminizer/runManySections/runManySections.py --submitCondor commands.cmd
condor_q lfeng7 
