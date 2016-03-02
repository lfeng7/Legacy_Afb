#! /bin/bash
# tlep Mass
python plotter.py --file "$1" --cut "isMatched[0]==1&&isMatched[1]==1" --cut2 "!(isMatched[0]==1&&isMatched[1]==1)" --var "reco_mass[0]" --Min 50 --Max 250 --name "mtlep_$2" --scale True --bin 100 --title "m_tlep $2" --xaxis "M_{tlep}(GeV)" --label "top matched" --label2 "not matched" --log ${3:-'no'} 
# thad mass
python plotter.py --file "$1" --cut "isMatched[0]==1&&isMatched[1]==1" --cut2 "!(isMatched[0]==1&&isMatched[1]==1)" --var "reco_mass[1]" --Min 50 --Max 500 --name "mthad_$2" --scale True --bin 100 --title "m_thad $2" --xaxis "M_{thad}(GeV)" --label "top matched" --label2 "not matched" --log ${3:-'no'}
# wlep mass
python plotter.py --file "$1" --cut "isMatched[0]==1&&isMatched[1]==1" --cut2 "!(isMatched[0]==1&&isMatched[1]==1)" --var "reco_mass[2]" --Min 70 --Max 200 --name "mwlep_$2" --scale True --bin 100 --title "m_wlep $2" --xaxis "M_{wlep}(GeV)" --label "top matched" --label2 "not matched" --log ${3:-'no'}
# whad mass
python plotter.py --file "$1" --cut "isMatched[0]==1&&isMatched[1]==1" --cut2 "!(isMatched[0]==1&&isMatched[1]==1)" --var "reco_mass[3]" --Min 20 --Max 300 --name "mwhad_$2" --scale True --bin 100 --title "m_whad $2" --xaxis "M_{whad}(GeV)" --label "top matched" --label2 "not matched" --log ${3:-'no'}
# chi2
python plotter.py --file "$1" --cut "isMatched[0]==1&&isMatched[1]==1" --cut2 "!(isMatched[0]==1&&isMatched[1]==1)" --var "lnL" --Min -20 --Max 100 --name "chi2_$2" --scale True --bin 100 --title "chi2 $2" --xaxis "#Chi2" --label "top matched" --label2 "not matched"
# N_jets
python plotter.py --file "$1" --cut "isMatched[0]==1&&isMatched[1]==1" --cut2 "!(isMatched[0]==1&&isMatched[1]==1)" --var "n_valid_jets" --Min 3 --Max 6 --name "njets_$2" --bin 10 --scale True --title "N_jets $2" --xaxis "N_{j}" --label "top matched" --label2 "not matched"
# leading jet pT
python plotter.py --file "$1" --cut "isMatched[0]==1&&isMatched[1]==1" --cut2 "!(isMatched[0]==1&&isMatched[1]==1)" --var "leadingJet_pt" --Min 20 --Max 100 --name "leadingJetPt_$2" --scale True --bin 100 --title "leadingJet_pT $2" --xaxis "leading_jet pT (GeV)" --label "top matched" --label2 "not matched"
# leading jet mass
python plotter.py --file "$1" --cut "isMatched[0]==1&&isMatched[1]==1" --cut2 "!(isMatched[0]==1&&isMatched[1]==1)" --var "leadingJet_mass" --Min 0 --Max 25 --name "leadingJetMass_$2" --scale True --bin 100 --title "leadingJet_mass $2" --xaxis "leading_jet mass (GeV)" --label "top matched" --label2 "not matched"
# ttbar mass
python plotter.py --file "$1" --cut "isMatched[0]==1&&isMatched[1]==1" --cut2 "!(isMatched[0]==1&&isMatched[1]==1)" --var "ttbar_mass" --Min 300 --Max 1500 --name "mtt_$2" --scale True --bin 100 --title "ttbar_mass $2" --xaxis "ttbar mass (GeV)" --label "top matched" --label2 "not matched"

#

