import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--n', metavar='F', type='string', action='store', dest='name', help='') ## Sets the run name
(options, args) = parser.parse_args()
mc_run = """root -l -b -q '../../main_MC.C(0,1,\""""
mc_run = mc_run + options.name
mc_run = mc_run + """\","MC_input_with_bkg.txt")'"""
print mc_run
print("RUNNING MC")
os.system(mc_run)
