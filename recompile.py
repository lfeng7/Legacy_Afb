import os
print("CLEARING OLD CODE")
os.system('rm -rf *.so *.d')
print("RECOMPILING")
os.system("""echo 'gROOT->ProcessLine(".L main_MC.C++"); gROOT->ProcessLine(".L main_data.C++"); ' | root -b -l""")
print("DONE")
