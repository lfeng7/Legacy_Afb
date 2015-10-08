# This small macro will read in all edm files in a directory and count the total number of events 

from fwlite_boilerplate import *
from ttbar_utitilty import *

from optparse import OptionParser

# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--maxfiles', metavar='F', type='int', action='store',
                  default = -1,
                  dest='maxfiles',
                  help='max number of input ntuple files')

parser.add_option('--upload', metavar='F', type='string', action='store',
                  default = 'dump',
                  dest='dumpplots',
                  help='if you want to dump plots to MY webpage! Unless you modify fwlite_boilerplate >_<')

parser.add_option('--verbose', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='verbose',
                  help='If you want more information than you usually need.')
(options, args) = parser.parse_args()

argv = []

# Some preset constants
data_lumi = 19.7*1000 

# Get input files
prepend = './output_rootfiles/all_channels/'   # dir of output files to stack 
  
# initialization and declaration
sample_types = []
flist = []
hlist_ = []

# stucture of a list of files
#    0,         1         2        3         4                 5
# (filepath, nevts_gen, xsec_NLO, type, nevts_total_ntuple, nevts_used_ntuple)
flist.append(['T_s_selection_output_all.root',259176,3.79,'singletop',259176,259176] )
flist.append(['TT_CT10_selection_output_all.root',21560109,245.9,'ttbar',21560109,1150618])

# list of histogram to make stack plots
hlist = ['cutflow','jets_pt']

# Booking stack histograms 
for hist in hlist:
    hlist_.append(ROOT.THStack(hist,'undefined'))
# put the name of histogram in the stack and the stack histograms in a list
stacklist = zip(hlist,hlist_)

# Make a legend
leg = ROOT.TLegend(0.3,0.6,0.6,0.8)

# Loop over files
for ifile in flist :
    # Getting files
    print 'processing file',ifile[0]
    f_ = ROOT.TFile(prepend+ifile[0])
    # Find the color of the sample type
    sample_type = ifile[3]
    icolor = sample_color(sample_type) # need to write a sub routine that return color given sample type

    # Calculate weight for this channel
    fraction = ifile[5]*1.0/ifile[4]
    cross_section_NLO = ifile[2]
    nevts_gen = ifile[1]*fraction
    weight = data_lumi*cross_section_NLO/nevts_gen  

    # Determine if we want to add an entry to legend
    if sample_type not in sample_types: 
        add_legend = 1
        sample_types.append(sample_type)
        # debugg
        if options.verbose == 'yes' : print 'Adding new sample type',sample_type 
    else : add_legend = 0

    # Loop over histograms and make stacks 
    for ihist in stacklist:
        ihist_name = ihist[0]
        istack = ihist[1]
        ih = f_.Get(ihist_name)     
        ih.Scale(weight) 
        ih.SetFillColor(icolor)
        ih.SetLineColor(icolor)
        ih.SetMarkerStyle(21)
        istack.Add(ih) 
        # Add entries to legend
        if add_legend : leg.AddEntry(ih,sample_type,"F")

    # for debugging
    if options.verbose == 'yes' :
        print 'weight is:',weight   
    
    # close temp file
    f_.Close()

##### end loop over files ######

# Plot and save
plotlist = [istack for name,istack in stacklist]  
plotting(plotlist,'stackplots',options.dumpplots,'not log',leg) 
  
