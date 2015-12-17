import os
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--n', metavar='F', type='string', action='store', dest='name', help='') ## Sets the run name
parser.add_option('--pdf', metavar='F', type='int', action='store', dest='pdf_index', help='') ## Sets the pdf member
(options, args) = parser.parse_args()

mc_run = """root -l -b -q '../../main_MC.C(0,0,\""""
mc_run = mc_run + options.name
mc_run = mc_run + """\","MC_input_with_bkg.txt","""
mc_run = mc_run + str(options.pdf_index)
mc_run = mc_run +""")'"""
print mc_run
print("RUNNING MC")
os.system(mc_run)
