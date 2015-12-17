#! /bin/sh
tar czvf tarball.tgz ../../../../*.* histo_files/*.root *.txt
/uscms_data/d3/eminizer/runManySections/runManySections.py --createCommandFile --cmssw --addLog --setTarball=tarball.tgz \ana.listOfJobs_combined commands.cmd_combined
/uscms_data/d3/eminizer/runManySections/runManySections.py --submitCondor commands.cmd_combined
condor_q lfeng7 
