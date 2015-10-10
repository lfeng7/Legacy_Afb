# This small macro will read in all edm files in a directory and count the total number of events 

from utility import *

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

parser.add_option('--MCplots', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='plotMConly',
                  help='If you want stack plots without data comparison')

(options, args) = parser.parse_args()

argv = []

# Some preset constants
data_lumi = 19.7*1000 
dir_name = 'stackplots'

# Get input files
prepend = './output_rootfiles/all_channels/'   # dir of output files to stack 

# Make an txt files for some information output
f_info = open('info_stackplots.txt','w')
  
# initialization and declaration
sample_types = []
flist = []
hlist_ = []
all_hists = []
data_hists = []

# stucture of a list of files
#    0,         1         2        3         4                 5
# (filepath, nevts_gen, xsec_NLO, type, nevts_total_ntuple, nevts_used_ntuple)
# backgrounds
flist.append(['T_s_selection_output_all.root',259176,3.79,'singletop',259176,259176] )
flist.append(['T_t_selection_output_all.root',3758227,56.4,'singletop',3748155,1536211] )
flist.append(['T_tW_selection_output_all.root',497658,11.1,'singletop',495559,495559])
flist.append(['Tbar_s_selection_output_all.root',139974,1.76,'singletop',139604,139604])
flist.append(['Tbar_t_selection_output_all.root',1935072,30.7,'singletop',1930185,1530347])
flist.append(['Tbar_tW_selection_output_all.root',493460,11.1,'singletop',491463,491463])
flist.append(['W3Jets_selection_output_all.root',15522558,640.4,'wjets',15507852,2263637])
flist.append(['W4Jets_selection_output_all.root',13326400,246.0,'wjets',13326400,1534294])
flist.append(['DY4_selection_output_all.root',6387032,28.59,'zjets',5843425,1508497])
flist.append(['DY3_selection_output_all.root',10997628,65.79,'zjets',10655325,1533162])
# signal
flist.append(['TT_CT10_selection_output_all.root',21560109,245.9,'ttbar',21560109,1150618])

######## data
#    0,         1                        2           3                  4                5    
# (filepath, sample_integrated_lumi, total_data_L, type, nevts_total_ntuple, nevts_used_ntuple)
datafile = ['SingleEl_Run2012A_all_selection_output_all.root',888,19748,'data',205,200]

# list of histogram to make stack plots
hlist = ['cutflow','jets_pt','Njets','m3','csv_all_jets','el_cand_pt','MET']

# Booking stack histograms 
for hist in hlist:
    hlist_.append(ROOT.THStack(hist,hist))
# put the name of histogram in the stack and the stack histograms in a list
stacklist = zip(hlist,hlist_)

# Make a legend
leg = ROOT.TLegend(0.3,0.6,0.6,0.8)

############ Get histograms from data files

# Getting files
print 'processing data file',datafile[0]
fdata = ROOT.TFile(prepend+datafile[0])
# Calculate weight for data 
fraction_ = datafile[5]*1.0/datafile[4]
weight_ = data_lumi/(datafile[1]*fraction_)
# Get histograms from data
for ihist in hlist :
    ih = fdata.Get(ihist)
    ih.SetDirectory(0)
    ih.Scale(weight_)
    ih.SetName(ih.GetName()+'_data')
    data_hists.append(ih)
# Add data entry to legend
leg.AddEntry(data_hists[0])

# write some informations about current sample
info_ = 'Sample type : data'+'\n'+'Events weight : '+str(weight_)+'\n'
info_ += 'Fraction of sample used : '+str(fraction_)+'\n\n'
f_info.write(info_)

############ Make stack histograms for MC samples

# Loop over files
for ifile in flist :
    # Getting files
    print 'processing MC file',ifile[0]
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
        # debug
        if options.verbose == 'yes' : print 'Adding new sample type',sample_type 
    else : add_legend = 0

    many_hists = []

    # Loop over histograms and make stacks 
    for ihist in stacklist:
        ihist_name = ihist[0]
        istack = ihist[1]
        ih = f_.Get(ihist_name)     
        # Delink the hist from the tmp file f_
        ih.SetDirectory(0)
        ih.Scale(weight) 
        ih.SetFillColor(icolor)
        ih.SetLineColor(icolor)
        ih.SetMarkerStyle(21)
        istack.Add(ih) 
        # Keep scaled histgrams in a list for later use
        many_hists.append(ih)

    all_hists.append(many_hists)
        
    # Add entries to legend
    if add_legend : 
        leg.AddEntry(many_hists[0],sample_type,"F")
        if options.verbose == 'yes' : print 'Adding a new legend entry for type:',sample_type

    # write some informations about current sample
    info_ = 'Sample type : '+sample_type+'\n'+'Events weight : '+str(weight)+'\n'
    info_ += 'Nevts : '+str(ifile[5])+'\n'+'Fraction of sample used : '+str(fraction)+'\n\n'
    f_info.write(info_)

    # for debugging
    if options.verbose == 'yes' :
        print 'weight is:',weight   
    

##### end loop over files ######

# Plot and save
mc_stacks = [istack for name,istack in stacklist]   # This is a list of stackplots

stack_cutflow = mc_stacks[0]
data_cutflow = data_hists[0]

stack_cutflow_norm = normstack(stack_cutflow)
data_cutflow_norm = norm(data_cutflow)

mc_stacks.append(stack_cutflow_norm)
data_hists.append(data_cutflow_norm)

# Plot MC stacks
if options.plotMConly == 'yes' :
    print 'Plotting MC stackplots without data comparison'
    plotting(mc_stacks,dir_name,'not dump','not log',leg) 
    plotting([stack_cutflow,stack_cutflow_norm],dir_name,options.dumpplots,'log',leg)

print 'Plotting comparison plots'

# Make data MC comparison plots
data_mc = zip(mc_stacks,data_hists)

for item in data_mc :
    comparison_plot(item[0],item[1],leg,dir_name,'not dump')

data_mc_log = ([stack_cutflow,data_cutflow],[stack_cutflow_norm,data_cutflow_norm])

for item in data_mc_log :
    comparison_plot(item[0],item[1],leg,dir_name,options.dumpplots,'log')

# Save MC stackplots and data histograms into an root files    
savelist = mc_stacks+data_hists+[leg]
saving(savelist,dir_name)
# file closure
f_info.close()
  
