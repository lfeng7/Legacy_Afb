"""
main code that produce templates for theta
"""
import Legacy_Afb.theta_fit.thetaTempMaker as thetaTempMaker 
import sys
import glob
import argparse

parser = argparse.ArgumentParser(description='make theta-templates')
parser.add_argument('-input',type=str, help='input is smoothed histograms')
parser.add_argument('-output',type=str, default='out', help='input is smoothed histograms')
parser.add_argument('-MC_info',type=str, default='MC_input_with_bkg.txt', help='path of MC_info.txt')
parser.add_argument('-smoothed', help='input is smoothed histograms', action='store_true')
parser.add_argument('-nevts', type=int, default=1000, help='nevts to run on each sample. default=1000')
parser.add_argument('-fixed', help='use fixed binning', action='store_true')
parser.add_argument('-verbose', help='if verbose', action='store_true')
parser.add_argument('-mc_data', help='use mc data sample for theoretical prediction.', action='store_true')

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

print 'options: thetaTemp(outputName=%s, use_smoothed_temp=%s, bin_type=%s, use_MC_Data=%s, MC_info=%s)'%(output_name,bool_to_str(args.smoothed),bin_type,bool_to_str(args.mc_data),args.MC_info)

self = thetaTempMaker.thetaTemp(outputName=output_name, inputFile=inputfile, isTTree=(not args.smoothed), txtfile=args.MC_info, bin_type=bin_type, use_MC_DATA=args.mc_data, verbose = args.verbose, nevts=args.nevts)
self.main()
print 'All finished.'


