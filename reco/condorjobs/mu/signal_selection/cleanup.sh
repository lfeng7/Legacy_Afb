mv *.log logs 
rm *.tgz  *.cmd condor-*
: << 'END'
hadd all_reco/T_s_reco.root T_s*.root
hadd all_reco/T_t_reco.root T_t*.root
hadd all_reco/T_tW_reco.root T_tW*.root
hadd all_reco/Tbar_s_reco.root Tbar_s*.root
hadd all_reco/Tbar_t_reco.root Tbar_t*.root
hadd all_reco/Tbar_tW_reco.root Tbar_tW*.root
hadd all_reco/W1JetsToLNu_TuneZ2Star_8TeV_reco.root W1JetsToLNu_TuneZ2Star_8TeV*.root
hadd all_reco/W2JetsToLNu_TuneZ2Star_8TeV_reco.root W2JetsToLNu_TuneZ2Star_8TeV*.root
hadd all_reco/W3JetsToLNu_TuneZ2Star_8TeV_reco.root W3JetsToLNu_TuneZ2Star_8TeV*.root
hadd all_reco/W4JetsToLNu_TuneZ2Star_8TeV_reco.root W4JetsToLNu_TuneZ2Star_8TeV*.root
hadd all_reco/DY1JetsToLL_M_reco.root DY1JetsToLL_M*.root
hadd all_reco/DY2JetsToLL_M_reco.root DY2JetsToLL_M*.root
hadd all_reco/DY3JetsToLL_M_reco.root DY3JetsToLL_M*.root
hadd all_reco/DY4JetsToLL_M_reco.root DY4JetsToLL_M*.root
hadd all_reco/TT_CT10_TuneZ2star_8TeV_reco.root TT_CT10_TuneZ2star_8TeV*.root
hadd all_reco/SingleEl_Run2012ABCD_reco.root SingleEl_Run2012ABCD*.root
hadd all_reco/SingleMu_Run2012ABCD_reco.root SingleMu_Run2012ABCD*.root
END
