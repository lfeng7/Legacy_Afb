# This small macro will read in all edm files in a directory and count the total number of events 
# v2. Will take a ttree, make histograms, and stack with data.

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
                  default = 'yes',
                  dest='verbose',
                  help='If you want more information than you usually need.')

parser.add_option('--MCplots', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='plotMConly',
                  help='If you want stack plots without data comparison')

parser.add_option('--makehists', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='makehists',
                  help='If you want to remake all histograms from selection files')

parser.add_option('--makeplots', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='makeplots',
                  help='If you want to make data MC comparison plots')

(options, args) = parser.parse_args()

argv = []

# Some preset constants
data_lumi = 19748 
dir_name = 'controlplots'

Nbin = 2

# Get input files
prepend = './output_rootfiles/all_selected/'   # dir of output files to make histograms
hist_prepend = './selected_hists/test/' 
  
# initialization and declaration
flist = []

#### Set up MC input files
# stucture of a list of files
#    0,         1         2        3         4           
# (filepath, nevts_gen, xsec_NLO, type, nevts_total_ntuple)
# Single Top
# flist.append(['T_s_v2_selection_output_all','singletop',259961,3.79,259176] )
# flist.append(['T_t_v2_selection_output_all.root','singletop',3758227,56.4,3748155] )
# flist.append(['T_tW_v2_selection_output_all.root','singletop',497658,11.1,495559])
# flist.append(['Tbar_s_v2_selection_output_all.root','singletop',139974, 1.76,139604])
# flist.append(['Tbar_t_v2_selection_output_all.root','singletop',1935072, 30.7,1930185])
# flist.append(['Tbar_tW_v2_selection_output_all.root','singletop',493460,11.1,491463])
# # Wjets
# flist.append(['W1JetsToLNu_v2.root','wjets',23141598,6662.8,23038253])
# flist.append(['W2JetsToLNu_v2.root','wjets',34044921,2159.2,33993463])
# flist.append(['W3Jets_v2_selection_output_all.root','wjets',15539503,640.4,15507852])
# flist.append(['W4Jets_v2_selection_output_all.root','wjets',13382803,246.0,13326400])
# # DYjets
# flist.append(['DY1JetsToLL_v2.root','zjets',24045248,660.6,23802736])
# flist.append(['DY2JetsToLL_v2.root','zjets',2352304,215.1,2345857])
# flist.append(['DY3Jets_v2_selection_output_all.root','zjets',11015445,65.79,10655325])
# flist.append(['DY4Jets_v2_selection_output_all.root','zjets',6402827,28.59,5843425])
# signal
flist.append(['TT_CT10_TuneZ2star_8TeV','ttbar',21675970,245.9,21560109])

#### Set up data input files
#    0,         1                   2             
# (filepath,   type,   sample_integrated_lumi
datafile = ['SingleEl_Run2012A_v1',888,19748,'data',11212832]
#datafile = ['SingleEl_Run2012ABCD_v2.root','data',19748]

def main():
    if options.makehists == 'yes':
        print 'Making histograms from the selection root files.'
        MakeHistograms()
    if options.makeplots == 'yes':
        print 'Making comparison plots of MC and data.'
        MakeComparisonPlots()
    print 'All done!'

def MakeHistograms():
    #### Making histograms for MC samples
    nbins = 50
    # Get a list of MC and data input file names
    selected_samples = [ifile[0] for ifile in flist]+[datafile[0]]
    # debug
    if options.verbose =='yes':
        print 'Making histograms for these samples:'
        for ifile in selected_samples : print ifile
    # Loop over all selection root files
    for ifile in flist:
        event_type = ifile[0]
        ifilename = prepend+event_type+'.root'
        print 'Reading root file:',ifilename
        tmpfile = ROOT.TFile(ifilename)
        tmptree = tmpfile.Get('selected')
        h_cutflow = tmpfile.Get('cutflow')
        h_cutflow_norm = tmpfile.Get('cutflow_norm')
        # Book Histograms
        h_lep_pt = ROOT.TH1D('lep_pt',event_type+' selected lepton pT;pT (GeV);events',nbins,0.,200.)
        h_lep_eta = ROOT.TH1D('lep_eta',event_type+' selected lepton eta;eta;events',nbins,-2.7,2.7)
        h_lep_charge = ROOT.TH1D('lep_charge',event_type+' charge of the selected lepton ;charge;events',10,-2,2)
        h_m3 = ROOT.TH1D('m3',event_type+' M3;m3;events',nbins,0.,500.)
        h_Njets = ROOT.TH1D('Njets',event_type+' Num selected jets;Njets;events',5,3,8)
        h_jets_pt = ROOT.TH1D('jets_pt',event_type+' selected jets pT;pT (GeV);events',nbins,0.,300.)
        h_jets_eta = ROOT.TH1D('jets_eta',event_type+' selected jets eta;eta;events',nbins,-2.7,2.7)
        h_MET = ROOT.TH1D('MET',event_type+' MET;MET;events',nbins,0.,200.)
        h_Nbjets = ROOT.TH1D('Nbjets',event_type+' Num tagged bjets;Nbjets;events',5,1,6)
        # Make a list of histograms for write
        tmplist = [h_cutflow,h_cutflow_norm,h_lep_pt,h_lep_eta,h_lep_charge, h_m3,h_Njets,h_jets_pt,h_jets_eta,h_MET,h_Nbjets]
        # Fill Histograms
        for i in range(tmptree.GetEntries()):
            tmptree.GetEntry(i)
            # Jets
            num_jets = tmptree.jets_pt.size()
            h_Njets.Fill(num_jets)
            for ijet in tmptree.jets_pt: h_jets_pt.Fill(ijet)
            for ijet in tmptree.jets_eta: h_jets_eta.Fill(ijet)
            bjets = []
            for ijet in tmptree.jets_csv:
                if ijet>csvm : bjets.append(ijet)
            h_Nbjets.Fill(len(bjets))
            h_m3.Fill(M3(tmptree.jets_pt,tmptree.jets_eta,tmptree.jets_phi,tmptree.mass))
            # leptons
            h_lep_pt.Fill(tmptree.lep_pt[0])
            h_lep_eta.Fill(tmptree.lep_eta[0])
            h_lep_charge.Fill(tmptree.lep_charge[0])
            #MET
            h_MET.Fill(tmptree.met_pt[0])

        # Save histograms into root files
        print 'Saving',isample_name,' histograms into root file..'
        savetoroot(tmplist,'selected_hists','test','_validation_plots')

def MakeComparisonPlots():

    ########################################################
    #               Making comparison plots                #
    ########################################################

    # list of histogram to make stack plots
    hlist = ['cutflow','jets_pt','Njets','m3','lep_pt','MET','jets_eta','lep_eta','Nbjets','lep_charge']
    rebinlist = ['jets_pt','m3','el_cand_pt','MET','jets_eta','el_cand_eta']

    # Booking stack histograms 
    for hist in hlist:
        hlist_.append(ROOT.THStack(hist,hist))
    # put the name of histogram in the stack and the stack histograms in a list
    stacklist = zip(hlist,hlist_)

    # Make a legend
    leg = ROOT.TLegend(0.7,0.65,1.0,0.85)

    # Process data files

    # Getting files
    data_input = hist_prepend+datafile[0]+'_validation_plots.root'
    print 'processing data file',data_input
    fdata = ROOT.TFile(data_input)
    # Calculate weight for data 
    data_fraction = 20./205
    nevts_data = fdata.Get('cutflow').GetBinContent(1)
    data_weight = data_lumi*1.0/(datafile[2]*data_fraction)

    data_hists = []
    # Get histograms from data
    for ihist in hlist :
        ih = fdata.Get(ihist)
        ih.SetDirectory(0)
        ih.Scale(data_weight)
        ih.SetName(ih.GetName()+'_data')
        # rebin some histograms
        # if ihist in rebinlist: ih.Rebin(Nbin)
        data_hists.append(ih)
    # Add data entry to legend
    leg.AddEntry(data_hists[0],'data')

    # write some informations about current sample
    info_ = 'Sample type : data'+'\n'+'Events weight %.2f : '%data_weight +'\n\n'
    f_info.write(info_)

    ############ Make stack histograms for MC samples

    sample_types = []
    # Loop over files
    for ifile in flist :
        # Getting files
        tmp_fname = hist_prepend+ifile[0]+'_validation_plots.root'
        print 'processing MC file',tmp_fname
        tmp_file = ROOT.TFile(tmp_fname)
        # Find the color of the sample type
        sample_type = ifile[1]
        icolor = GetSampleColor(sample_type) # GetSampleColor is defined in ttbar_utility

        # Calculate weight for this channel
        nevts_total = int(tmp_file.Get('cutflow').GetBinContent(1))
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
            ih = tmp_file.Get(ihist_name)    
            # rebin some histograms
            # if ihist in rebinlist: ih.Rebin(Nbin)
            # Delink the hist from the tmp_file
            ih.SetDirectory(0)
            ih.Scale(weight) 
            ih.SetFillColor(icolor)
            ih.SetLineColor(icolor)
            ih.SetMarkerStyle(21)
            istack.Add(ih) 
            
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

    for item in data_mc :
        comparison_plot(item[0],item[1],leg,dir_name)

    data_mc_log = ([stack_cutflow,data_cutflow])

    for item in data_mc_log :
        comparison_plot(item[0],item[1],leg,dir_name,'not dump','log','p')

    # Dump plots to web
    print 'Uploading all plots to web'
    plotting([],dir_name,'dump')

    ############ Save MC stackplots and data histograms into an root files    
    savelist = mc_stacks+data_hists+[leg]
    saving(savelist,dir_name)

    #################################################################
    #               Making yields table                             #
    #################################################################

    # Make an txt files for some information output
    f_info = open('info_stackplots.txt','w')
    f_yields = open('yields.txt','w')

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


main()
