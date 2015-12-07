import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--n', metavar='F', type='string', action='store', dest='name', help='') ## Sets the run name
(options, args) = parser.parse_args()
mc_run = """root -l -b -q '../../main_MC.C(0,0,\""""
mc_run = mc_run + options.name
mc_run = mc_run + """\","MC_input_with_bkg.txt",0)'"""
print mc_run
print("RUNNING MC")
os.system(mc_run)
data_run = """root -l -b -q '../../main_data.C(0,0,1,\""""
data_run = data_run + options.name
data_run = data_run + """\","MC_input_with_bkg.txt","data_input.txt")'"""
print data_run
print("RUNNING DATA")
os.system(data_run)
os.system('rm -rf final_fit_stuff')
os.system('rm -rf data_histo_files')
print("ALLDONE")
os.system("cat summary_separated.txt")
