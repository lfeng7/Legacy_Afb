# Take in the angles ttree from reco, and make input ttree for original template making and fitting codes
import ROOT
import glob
from optparse import OptionParser
from template_class import *

ROOT.gROOT.SetBatch(True)

# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--evtsperjob', metavar='F', type='int', action='store',
                  default = 1000,
                  dest='evtsperjob',
                  help='number of events to run for each job')

parser.add_option('--evtstart', metavar='F', type='int', action='store',
                  default = 0,
                  dest='evtstart',
                  help='the evt to start')

parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='If you want more information than you usually need.')

parser.add_option('--slim', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='slim',
                  help='If you want slimmed ttree with angles and stuff')

parser.add_option('--fakelep', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='fakelep',
                  help='If run on selected events with fake lepton.')

parser.add_option('--lepisocut', metavar='F', type='float', action='store',
                  default = 0.15,
                  dest='lepisocut',
                  help='Lower bound for fake electron isolation.')

parser.add_option('--nbcut', metavar='F', type='int', action='store',
                  default = 1,
                  dest='nbcut',
                  help='Number of b-tagged jets cut')

parser.add_option('--nlepcut', metavar='F', type='int', action='store',
                  default = 1,
                  dest='nlepcut',
                  help='Number of selected leptons cut')

parser.add_option('--applytrigger', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='applytrigger',
                  help='If apply trigger on MC')

parser.add_option('--ttbar_type', metavar='F', type='string', action='store',
                  default = 'none',
                  dest='ttbar_type',
                  help='type of ttbar templates. gg/qq/bkg')

parser.add_option('--lep_type', metavar='F', type='string', action='store',
                  default = 'ele',
                  dest='lep_type',
                  help='e+jets or mu+jets')

# need special treatment for QCD bkg tempates
parser.add_option('--sideband', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='sideband',
                  help='if we are making templates for sideband')

parser.add_option('--MCinfo', metavar='F', type='string', action='store',
                  default = '',
                  dest='MCinfo',
                  help='dir path to MCinfo.txt')

(options, args) = parser.parse_args()

argv = []


def str_to_bool(_str):
    # convert 'yes','no' etc to Bool True or False
    if _str in ['yes','Yes','true','True','T']:
        return True
    elif _str in ['no','No','false','False','F']:
        return False

# Some preset constants
csvm = 0.679 
deltaR_matching = 0.4

# Get the file list with all input files.
if options.inputFiles != '':
    allfiles = glob.glob( options.inputFiles )
else:
    allfiles = []

timer = ROOT.TStopwatch()
timer.Start()

# create a template class object
new_template = template()
if options.lep_type in ['mu','mu_jets']:
  new_template.lep_type = 'mu_jets'
else:
  new_template.lep_type = 'e_jets'
# load MC.txt info
if options.MCinfo=='':
    MCinfo_txt='/uscms_data/d3/lfeng7/Payloads/run1/MC_input_with_bkg.txt'
else:
    MCinfo_txt= options.MCinfo

# Job splitting is done here
for ifile in allfiles:
    sample_name = ifile.split('/')
    sample_name = sample_name[len(sample_name)-1].split('.root')[0]
    evt_start = options.evtstart
    evt_to_run = options.evtsperjob
    tfile = ROOT.TFile(ifile)
    # Do make templates
    print 'run options:',tfile.GetName(),sample_name,evt_start,evt_to_run

    # Set template class attributes
    new_template.sample_name = sample_name
    new_template.evt_start = evt_start
    new_template.evt_to_run = evt_to_run
    new_template.input_file = tfile
    new_template.ttbar_type = options.ttbar_type
    new_template.isSideband = str_to_bool(options.sideband)
    new_template.set_MCinfo(MCinfo_txt)
    # Make template
    reco_message = new_template.makeTemps()   
    print reco_message

print 'All done!'

# Stop our timer
timer.Stop()
# Print out our timing information
print '\n'
rtime = timer.RealTime(); # Real time (or "wall time")
ctime = timer.CpuTime(); # CPU time
print("RealTime={0:6.2f} seconds, CpuTime={1:6.2f} seconds").format(rtime,ctime)



