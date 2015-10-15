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
data_lumi = 19748 
dir_name = 'stackplots'

Nbin = 2

# Get input files
prepend = './output_rootfiles/all_channels/'   # dir of output files to stack 

# Make an txt files for some information output
f_info = open('info_stackplots.txt','w')
f_yields = open('yields.txt','w')
  
# initialization and declaration
sample_types = []
flist = []
hlist_ = []
all_hists = []
data_hists = []

# stucture of a list of files
#    0,         1         2        3         4           
# (filepath, nevts_gen, xsec_NLO, type, nevts_total_ntuple)
# Single Top
flist.append(['T_s_v2_selection_output_all.root',259961,3.79,'singletop',259176] )
flist.append(['T_t_v2_selection_output_all.root',3758227,56.4,'singletop',3748155] )
flist.append(['T_tW_v2_selection_output_all.root',497658,11.1,'singletop',495559])
flist.append(['Tbar_s_v2_selection_output_all.root',139974,1.76,'singletop',139604])
flist.append(['Tbar_t_v2_selection_output_all.root',1935072,30.7,'singletop',1930185])
flist.append(['Tbar_tW_v2_selection_output_all.root',493460,11.1,'singletop',491463])
# Wjets
flist.append(['W1JetsToLNu_v2.root',23141598,6662.8,'wjets',23038253])
flist.append(['W2JetsToLNu_v2.root',34044921,2159.2,'wjets',33993463])
flist.append(['W3Jets_v2_selection_output_all.root',15539503,640.4,'wjets',15507852])
flist.append(['W4Jets_v2_selection_output_all.root',13382803,246.0,'wjets',13326400])
# DYjets
flist.append(['DY1JetsToLL_v2.root',24045248,660.6,'zjets',23802736])
flist.append(['DY2JetsToLL_v2.root',2352304,215.1,'zjets',2345857])
flist.append(['DY3Jets_v2_selection_output_all.root',11015445,65.79,'zjets',10655325])
flist.append(['DY4Jets_v2_selection_output_all.root',6402827,28.59,'zjets',5843425])
# signal
flist.append(['TT_CT10_v2_selection_output_all.root',21675970,245.9,'ttbar',21560109])

######## data
#    0,         1                        2           3                  4                5    
# (filepath, sample_integrated_lumi, total_data_L, type, nevts_total_ntuple, nevts_used_ntuple)
#datafile = ['SingleEl_Run2012A_v2_selection_output_all.root',888,19748,'data',11212832]
datafile = ['SingleEl_Run2012ABCD_v2.root',19748,19748,'data',11212832]

# list of histogram to make stack plots
hlist = ['cutflow','jets_pt','Njets','m3','csv_all_jets','el_cand_pt','MET','jets_eta','el_cand_eta']
rebinlist = ['jets_pt','m3','el_cand_pt','MET','jets_eta','el_cand_eta']

# Booking stack histograms 
for hist in hlist:
    hlist_.append(ROOT.THStack(hist,hist))
# put the name of histogram in the stack and the stack histograms in a list
stacklist = zip(hlist,hlist_)

# Make a legend
leg = ROOT.TLegend(0.7,0.65,1.0,0.85)

############ Get histograms from data files

# Getting files
print 'processing data file',datafile[0]
fdata = ROOT.TFile(prepend+datafile[0])
# Calculate weight for data 
nevts_data = fdata.Get('cutflow').GetBinContent(1)
fraction_ = 1.0 # nevts_data*1.0/datafile[4]
weight_ = data_lumi/(datafile[1]*fraction_)
# Get histograms from data
for ihist in hlist :
    ih = fdata.Get(ihist)
    ih.SetDirectory(0)
    ih.Scale(weight_)
    ih.SetName(ih.GetName()+'_data')
    # rebin some histograms
    if ihist in rebinlist: ih.Rebin(Nbin)

    data_hists.append(ih)
# Add data entry to legend
leg.AddEntry(data_hists[0],'data')

# write some informations about current sample
info_ = 'Sample type : data'+'\n'+'Events weight %.2f : '%weight_ +'\n'
info_ += 'Fraction of sample used : %.2f'%fraction_ +'\n\n'
f_info.write(info_)

############ Make stack histograms for MC samples

# Loop over files
for ifile in flist :
    # Getting files
    print 'processing MC file',ifile[0]
    f_ = ROOT.TFile(prepend+ifile[0])
    # Find the color of the sample type
    sample_type = ifile[3]
    icolor = GetSampleColor(sample_type) # need to write a sub routine that return color given sample type

    # Calculate weight for this channel
    nevts_total = int(f_.Get('cutflow').GetBinContent(1))
    fraction = nevts_total*1.0/ifile[4]
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

        # rebin some histograms
        if ihist in rebinlist: ih.Rebin(Nbin)

        istack.Add(ih) 
        # Keep scaled histgrams in a list for later use
        many_hists.append(ih)

    all_hists.append(many_hists)
        
    # Add entries to legend
    if add_legend : 
        leg.AddEntry(many_hists[0],sample_type,"F")
        if options.verbose == 'yes' : print 'Adding a new legend entry for type:',sample_type

    # write some informations about current sample
    info_ = 'Sample name : '+ifile[0]+'\n'+'Type: '+sample_type+'\n'+'Events weight : '+str(weight)+'\n'
    info_ += 'Nevts : '+str(nevts_total)+'\n'+'Fraction of sample used : %.2f'%fraction+'\n\n'
    f_info.write(info_)

    # for debugging
    if options.verbose == 'yes' :
        print 'weight is:',weight   
    

##### end loop over files ######

# Plot and save
mc_stacks = [istack for name,istack in stacklist]   # This is a list of stackplots

stack_cutflow = mc_stacks[0]
data_cutflow = data_hists[0]
data_cutflow.SetMinimum(2000)

# this is an ugly fix in order to get correct yields....
stack_cutflow_0 = stack_cutflow.Clone()
data_cutflow_0 = data_cutflow.Clone()

stack_cutflow_norm = normstack(stack_cutflow)
data_cutflow_norm = norm(data_cutflow)
data_cutflow_norm.SetMinimum(0.00002)

#mc_stacks.append(stack_cutflow_norm)
#data_hists.append(data_cutflow_norm)

# Create the output root file
plotting([],dir_name,'not dump','not log',None,'','recreate')
# Plot MC stacks
if options.plotMConly == 'yes' :
    print 'Plotting MC stackplots without data comparison'
    plotting(mc_stacks,dir_name,'not dump','not log',leg) 
    plotting([stack_cutflow,stack_cutflow_norm],dir_name,options.dumpplots,'log',leg)

print 'Plotting comparison plots'

# Make data MC comparison plots
data_mc = zip(mc_stacks,data_hists)
data_mc_norm = [[stack_cutflow_norm,data_cutflow_norm]]

for item in data_mc :
    comparison_plot(item[0],item[1],leg,dir_name)
for item in data_mc_norm :
    comparison_plot(item[0],item[1],leg,dir_name,'not dump','notlog','p')

data_mc_log = ([stack_cutflow,data_cutflow],[stack_cutflow_norm,data_cutflow_norm])

for item in data_mc_log :
    comparison_plot(item[0],item[1],leg,dir_name,'not dump','log','p')

# Dump plots to web
print 'Uploading all plots to web'
plotting([],dir_name,'dump')

############ Save MC stackplots and data histograms into an root files    
savelist = mc_stacks+data_hists+[leg]
saving(savelist,dir_name)

# Getting yields for MC and data
sample_yields = []
for ihist in stack_cutflow_0.GetHists():
    type_ = GetSampleType(ihist.GetFillColor())
    sample_yields.append([type_]+GetBinEntry(ihist))

MC_yields = []

all_types = GetListChoices([ sample[0] for sample in sample_yields])

for itype in all_types :
    ilist = [ sample for sample in sample_yields if sample[0] == itype]
    MC_yields.append(SumColumn(ilist))

mc_total_yields = SumColumn(MC_yields)
mc_total_yields.pop(0)
MC_yields.append(['MC']+mc_total_yields)

# Yields for data
data_yields = ['data']+GetBinEntry(data_cutflow_0)
# Write into yields output files
f_yields.write('sample,nocut, el selection, loose mu veto, dilep veto, jets selection, b-tagging, MET \n')
for row in MC_yields :
    for item in row :
        f_yields.write(str(item)+',')
    f_yields.write('\n')
for item in data_yields :
    f_yields.write(str(item)+',')

# file closure
f_info.close()
f_yields.close()
fdata.Close()  
