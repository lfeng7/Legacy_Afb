rm -r *.cmd  *.tgz notneeded condor-* templates condorlog
mkdir templates
mv *.root templates
mkdir condorlog
mv *.log condorlog
