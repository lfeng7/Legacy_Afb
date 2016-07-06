cd SingleEl_Run2012A
source cleanup.sh
hadd ../../selected_files/qcd_sideband/SingleEl_Run2012A_selected.root SingleEl_Run2012A/*.root 
rm -r SingleEl_Run2012A 
cd ..

cd SingleEl_Run2012B
source cleanup.sh
hadd ../../selected_files/qcd_sideband/SingleEl_Run2012B_selected.root SingleEl_Run2012B/*.root 
rm -r SingleEl_Run2012B 
cd ..

cd SingleEl_Run2012C_part1
source cleanup.sh
hadd ../../selected_files/qcd_sideband/SingleEl_Run2012C_part1_selected.root SingleEl_Run2012C_part1/*.root 
rm -r SingleEl_Run2012C_part1 
cd ..

cd SingleEl_Run2012C_part2
source cleanup.sh
hadd ../../selected_files/qcd_sideband/SingleEl_Run2012C_part2_selected.root SingleEl_Run2012C_part2/*.root 
rm -r SingleEl_Run2012C_part2 
cd ..

cd SingleEl_Run2012D
source cleanup.sh
hadd ../../selected_files/qcd_sideband/SingleEl_Run2012D_selected.root SingleEl_Run2012D/*.root 
rm -r SingleEl_Run2012D 
cd ..

cd T_t
source cleanup.sh
hadd ../../selected_files/qcd_sideband/T_t_selected.root T_t/*.root 
rm -r T_t 
cd ..

cd T_tW
source cleanup.sh
hadd ../../selected_files/qcd_sideband/T_tW_selected.root T_tW/*.root 
rm -r T_tW 
cd ..

cd Tbar_s
source cleanup.sh
hadd ../../selected_files/qcd_sideband/Tbar_s_selected.root Tbar_s/*.root 
rm -r Tbar_s 
cd ..

cd Tbar_t
source cleanup.sh
hadd ../../selected_files/qcd_sideband/Tbar_t_selected.root Tbar_t/*.root 
rm -r Tbar_t 
cd ..

cd Tbar_tW
source cleanup.sh
hadd ../../selected_files/qcd_sideband/Tbar_tW_selected.root Tbar_tW/*.root 
rm -r Tbar_tW 
cd ..
