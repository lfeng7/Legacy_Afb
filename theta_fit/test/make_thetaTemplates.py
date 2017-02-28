"""
main code that produce templates for theta
"""
import Legacy_Afb.theta_fit.thetaTempMaker as thetaTempMaker 
import sys
import glob
import argparse

parser = argparse.ArgumentParser(description="""
make theta-templates
Usage: python make_thetaTemplates.py -input "temp_angles/*.root" -output ele_mcNLO -nevts -1 -MC_info temp_angles/MC_input_with_bkg_mcNLO.txt
or   : python make_thetaTemplates.py -input "temp_angles/*.root" -output ele_mcNLO_mc_data -nevts -1 -MC_info temp_angles/MC_input_with_bkg_mcNLO.txt -mc_data
""")
parser.add_argument('-input',type=str, help='input is smoothed histograms')
parser.add_argument('-output',type=str, default='out', help='input is smoothed histograms')
parser.add_argument('-lep_type', type=str, help='el or mu, deciding what correction SF to use for templates')
parser.add_argument('-MC_info',type=str, default='MC_input_with_bkg.txt', help='path of MC_info.txt')
parser.add_argument('-smoothed', help='input is smoothed histograms', action='store_true')
parser.add_argument('-nevts', type=int, default=1000, help='nevts to run on each sample. default=1000')
parser.add_argument('-fixed', help='use fixed binning', action='store_true')
parser.add_argument('-verbose', help='if verbose', action='store_true')
parser.add_argument('-mc_data', help='use mc data sample for theoretical prediction.', action='store_true')
parser.add_argument('-top_reweight', help='use top_pT_reweighting', action='store_true')
parser.add_argument('-use_sys', help='create sys templates', action='store_true')


args = parser.parse_args()

print 'Run on %d evts.' % args.nevts

if args.fixed:
    bin_type = 'fxied'
else:
    bin_type = 'variable'

if not args.smoothed:
    """ Processing ttree files"""
    inputfile = glob.glob(args.input)
else:
    inputfile = args.input
output_name = '%s'%args.output

def bool_to_str(var):
    if var:
        return 'true'
    else:
        return 'false'

print 'options: thetaTemp(outputName=%s,lep_type=%s, use_smoothed_temp=%s, bin_type=%s, use_MC_Data=%s,top_weight=%s,use_sys = %s, MC_info=%s)'%(output_name,args.lep_type,bool_to_str(args.smoothed),bin_type,bool_to_str(args.mc_data),bool_to_str(args.top_reweight),bool_to_str(args.use_sys) , args.MC_info)

self = thetaTempMaker.thetaTemp(outputName=output_name, inputFile=inputfile,lep_type = args.lep_type, isTTree=(not args.smoothed), txtfile=args.MC_info, bin_type=bin_type, use_MC_DATA=args.mc_data, top_reweight=args.top_reweight, use_sys=args.use_sys, verbose = args.verbose, nevts=args.nevts)
self.main()
print 'All finished.'


