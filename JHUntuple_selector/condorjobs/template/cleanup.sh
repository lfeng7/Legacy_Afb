rm -r *.cmd  *.tgz notneeded condor-* condorlog template
mkdir template
mv *.root template
mkdir condorlog
mv *.log condorlog
