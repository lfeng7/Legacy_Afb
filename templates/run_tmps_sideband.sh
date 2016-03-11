#! /bin/bash
python makeTemp.py --inputfiles angles_files/sideband_triggered/DY1JetsToLL_M_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/DY2JetsToLL_M_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/DY3JetsToLL_M_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/DY4JetsToLL_M_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/SingleEl_Run2012ABCD_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/T_s_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/T_tW_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/T_t_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/Tbar_s_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/Tbar_tW_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/Tbar_t_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/W1JetsToLNu_TuneZ2Star_8TeV_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/W2JetsToLNu_TuneZ2Star_8TeV_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/W3JetsToLNu_TuneZ2Star_8TeV_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/W4JetsToLNu_TuneZ2Star_8TeV_reco_angles_0.root --evtsperjob -1
python makeTemp.py --inputfiles angles_files/sideband_triggered/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root --evtsperjob -1 --ttbar_type qq
python makeTemp.py --inputfiles angles_files/sideband_triggered/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root --evtsperjob -1 --ttbar_type gg
python makeTemp.py --inputfiles angles_files/sideband_triggered/TT_CT10_TuneZ2star_8TeV_reco_angles_0.root --evtsperjob -1 --ttbar_type bkg

