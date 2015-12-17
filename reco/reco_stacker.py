# This small macro will read in all edm files in a directory and count the total number of events 
# v2. Will take a ttree, make histograms, and stack with data. And calculate correction weights for MC

from Legacy_Afb.Tools.ttbar_utility import * 
from Legacy_Afb.Tools.root_utility import *

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
                  default = 'yes',
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

parser.add_option('--makehists', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='makehists',
                  help='If you want to remake all histograms from selection files')

parser.add_option('--makeplots', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='makeplots',
                  help='If you want to make data MC comparison plots')

parser.add_option('--toppt', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='toppt',
                  help='If do top pt reweighting')

parser.add_option('--correction', metavar='F', type='string', action='store',
                  default = 'yes',
                  dest='correction',
                  help='If use correction')

parser.add_option('--tmptype', metavar='F', type='string', action='store',
                  default = 'scaled',
                  dest='tmptype',
                  help='If use correction')

parser.add_option('--applytrigger', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='applytrigger',
                  help='If apply trigger on MC')

parser.add_option('--fakelep', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='fakelep',
                  help='If run on selected events with fake lepton.')

parser.add_option('--lepisocut', metavar='F', type='int', action='store',
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

(options, args) = parser.parse_args()

argv = []

# Some preset constants
data_lumi = 19748 
csvm = 0.679

# Get input files
if options.fakelep == 'no':
    prepend = './angles_files/signal/'
else:
    prepend = './angles_files/sideband/'   # dir of output files to make histograms
postfix='_reco_angles_0.root'
# Set up output histogram files
template_type = options.tmptype
tmptype_name = template_type
if options.fakelep == 'yes':
    tmptype_name += '_fakelep'
hist_prepend = './selected_hists/'+tmptype_name+'/' 
    
dir_name = 'reco_plots_'+tmptype_name # name of the dir for control plots, such as corrected, topPT, un_corrected
  
# initialization and declaration
flist = []

#### Set up MC input files

#    0,                 1         2          3         4                   5
# (MC_sample_name, sample_type, Nevts_gen, x-sec, nevts_total_ntuple, btag_type)
# Single Top
flist.append(['T_s','singletop',259961,3.79,259176,'singletop'] )
flist.append(['T_t','singletop',3758227,56.4,3748155,'singletop'] )
flist.append(['T_tW','singletop',497658,11.1,495559,'singletop'])
flist.append(['Tbar_s','singletop',139974, 1.76,139604,'singletopbar'])
flist.append(['Tbar_t','singletop',1935072, 30.7,1930185,'singletopbar'])
flist.append(['Tbar_tW','singletop',493460,11.1,491463,'singletopbar'])
# Wjets
flist.append(['W1JetsToLNu_TuneZ2Star_8TeV','wjets',23141598,6662.8,23038253,'wjets'])
flist.append(['W2JetsToLNu_TuneZ2Star_8TeV','wjets',34044921,2159.2,33993463,'wjets'])
flist.append(['W3JetsToLNu_TuneZ2Star_8TeV','wjets',15539503,640.4,15507852,'wjets'])
flist.append(['W4JetsToLNu_TuneZ2Star_8TeV','wjets',13382803,246.0,13326400,'wjets'])
# DYjets
flist.append(['DY1JetsToLL_M','zjets',24045248,660.6,23802736,'zjets'])
flist.append(['DY2JetsToLL_M','zjets',2352304,215.1,2345857,'zjets'])
flist.append(['DY3JetsToLL_M','zjets',11015445,65.79,10655325,'zjets'])
flist.append(['DY4JetsToLL_M','zjets',6402827,28.59,5843425,'zjets'])
# QCD
#flist.append(['QCD_Pt-15to3000','qcd',9991674,6662.6,9940092,'qcd'])
# signal
flist.append(['TT_CT10_TuneZ2star_8TeV','ttbar',21675970,245.9,21560109,'ttbar'])

#### Set up data input files
#    0,         1                   2             
# (filepath,   type,   sample_integrated_lumi
#datafile = ['SingleEl_Run2012A_v1','data',888]
datafile = ['SingleEl_Run2012ABCD','data',19748]

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
    htitle = '19.7 fb^{-1} @ 8 TeV'
    nbins = 50
    # Get a list of MC and data input file names
    selected_samples = [(ifile[0],ifile[1],ifile[5]) for ifile in flist]+[(datafile[0],datafile[1])]
    # debug
    if options.verbose =='yes':
        print 'Making histograms for these samples:'
        for ifile in selected_samples : print ifile

    # Loop over all selection root files
    for ifile in selected_samples:
        event_type = ifile[0] # Name of the sample
        sample_type = ifile[1] # type of the sample, data, ttbar ,etc
        ifilename = prepend+event_type+postfix
        print '\nReading root file:',ifilename
        tmpfile = ROOT.TFile(ifilename)
        tmptree = tmpfile.Get('selected')
        # Make an output file for histograms
        savetoroot([],'selected_hists',tmptype_name,event_type+'_control_plots')

        # Physics objects
        h_cos = ROOT.TH1D('cos',htitle+';cos(theta);events',nbins,-1.,1.)
        h_xf = ROOT.TH1D('xf',htitle+';x_{f};events',nbins,-0.6,0.6)
        h_mtt = ROOT.TH1D('mtt',htitle+';M_{tt} (GeV);events',nbins,345,1600)
        h_chi2 = ROOT.TH1D('chi2',htitle+';chi^2;events',nbins,-45,100.)
        h_nz = ROOT.TH1D('nz',htitle+';nv P_{z} (GeV);events',nbins,-700,700.)

        # Make a list of histograms for write
        tmplist = [h_cos,h_xf,h_mtt,h_chi2,h_nz]

        # Remove the attachement of histograms from input root file, debug only
        for ihist in tmplist : ihist.SetDirectory(0)

        if sample_type == 'data':
            # Fill Histograms
            for i in range(tmptree.GetEntries()):
                # Progress report
                if i%10000 == 0 : print 'processing event',i
                tmptree.GetEntry(i)

                #######################################################
                #          Fill the control hists                     #
                #######################################################
                h_cos.Fill(tmptree.cos_theta[0])
                h_xf.Fill(tmptree.xf[0])
                h_mtt.Fill(tmptree.mtt[0])
                h_chi2.Fill(tmptree.final_chi2[0])
                h_nz.Fill(tmptree.final_nv_pz[0])
        else:
            ########################################################
            #               Make histograms for MC                 #
            ########################################################       

            #######################################################
            #     Calculate the normalization for top pT SFs      #
            #######################################################
            # To avoid the top pt weight change the overall normalization 
            topPtScale = 1.0
            if sample_type == 'ttbar' :
                # Get the correct normalization of topPt reweights
                h_cos_0 = h_cos.Clone()
                h_cos_1 = h_cos.Clone()
                h_cos_0.SetDirectory(0)
                h_cos_1.SetDirectory(0)

                for i in range(tmptree.GetEntries()):
                    tmptree.GetEntry(i)

                    #######################################################
                    #          Calculate Normalzation                     #
                    #######################################################                                    
                    w_top_pt = tmptree.weight_top_pT[0]
                    h_cos_0.Fill(tmptree.cos_theta[0])
                    h_cos_1.Fill(tmptree.cos_theta[0],w_top_pt)

                topPtScale = h_cos_1.Integral()*1.0/h_cos_0.Integral()
                print 'topPtScale = ',topPtScale

            # Loop over entries for MC
            for i in range(tmptree.GetEntries()):
                # Progress report
                if i%10000 == 0 : print 'processing event',i
                tmptree.GetEntry(i)

                #######################################################
                #          Get       correction SFs                   #
                #######################################################

                # top pT
                w_top_pt = tmptree.weight_top_pT[0]/topPtScale
                w_PU = tmptree.w_PU[0]
                w_btag = tmptree.w_btag[0]
                w_eleID = tmptree.w_eleID[0]
                w_trigger = tmptree.w_trigger[0]
                # book all weights
                all_weights = []
                if options.tmptype == 'scaled':
                    all_weights += [w_PU,w_btag,w_eleID,w_top_pt]
                if options.tmptype == 'no_top_pt':
                    all_weights += [w_PU,w_btag,w_eleID]
                if options.applytrigger == 'yes':
                    all_weights += [w_trigger]
                # Get final weight
                event_weight = 1.0
                for iw in all_weights : event_weight *= iw

                #######################################################
                #          Fill control plots                         #
                #######################################################             
                h_cos.Fill(tmptree.cos_theta[0],event_weight)
                h_xf.Fill(tmptree.xf[0],event_weight)
                h_mtt.Fill(tmptree.mtt[0],event_weight)
                h_chi2.Fill(tmptree.final_chi2[0],event_weight)
                h_nz.Fill(tmptree.final_nv_pz[0],event_weight)

        # Save histograms into root files
        print 'Saving',event_type,' histograms into root file..'
        #        will save to selected_hists/test/event_type_control_plots.root
        savetoroot(tmplist,'selected_hists',tmptype_name,event_type+'_control_plots','update')
        # Close selected root file
        tmpfile.Close()

def MakeComparisonPlots():

    ########################################################
    #               Making comparison plots                #
    ########################################################

    # list of histogram to make stack plots
    hlist = ['cos','xf','mtt','chi2','nz']

    stacks = []
    # Booking stack histograms 
    for hist in hlist:
        stacks.append(ROOT.THStack(hist,hist))
    # put the name of histogram in the stack and the stack histograms in a list
    stacklist = zip(hlist,stacks)

    # Make a legend
    leg = ROOT.TLegend(0.85,0.65,1.0,1.0)

    # Process data files

    # Getting files
    data_input = hist_prepend+datafile[0]+'_control_plots.root'
    print 'processing data file',data_input
    fdata = ROOT.TFile(data_input)
    # Calculate weight for data 
    data_fraction = 1.0
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
    f_info = open('info_stackplots.txt','w')

    ############ Make stack histograms for MC samples

    sample_types = []
    # Loop over files
    for ifile in flist :
        # Skip QCD MC samples for good reasons..
        if ifile[1] == 'qcd': continue
        # Getting files
        tmp_fname = hist_prepend+ifile[0]+'_control_plots.root'
        print 'processing MC file',tmp_fname
        tmp_file = ROOT.TFile(tmp_fname)
        # Find the color of the sample type
        sample_type = ifile[1]
        icolor = GetSampleColor(sample_type) # GetSampleColor is defined in ttbar_utility
        #print sample_type,icolor

        # Calculate weight for this channel
        fraction = 1.0
        cross_section_NLO = ifile[3]
        nevts_gen = ifile[2]*fraction
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
            many_hists.append(ih)
            
        # Add entries to legend
        if add_legend : 
            leg.AddEntry(many_hists[0],sample_type,"F")
            if options.verbose == 'yes' : print 'Adding a new legend entry for type:',sample_type

        # for debugging
        if options.verbose == 'yes' :
            print 'weight is:',weight           

    ##### end loop over files ######

    # Plot and save
    mc_stacks = [istack for name,istack in stacklist]   # This is a list of stackplots

    # Create the output root file
    plotting([],dir_name,'not dump','not log',None,'','recreate')
    # Plot MC stacks
    if options.plotMConly == 'yes' :
        print 'Plotting MC stackplots without data comparison'
        plotting(mc_stacks,dir_name,'not dump','not log',leg) 

    print 'Plotting comparison plots'

    # Make data MC comparison plots
    data_mc = zip(mc_stacks,data_hists)

    for item in data_mc :
        comparison_plot_v1(item[0],item[1],leg,dir_name)

    data_mc_log = ([])

    for item in data_mc_log :
        comparison_plot(item[1],item[0],leg,dir_name,'not dump','log','p')

    # Dump plots to web
    if options.dumpplots == 'yes':
        print 'Uploading all plots to web'
        plotting([],dir_name,'dump')

    ############ Save MC stackplots and data histograms into an root files    
    savelist = mc_stacks+data_hists+[leg]
    fstacks = ROOT.TFile('fstacks.root','recreate')
    for item in savelist:
    #    item.SetDirectory(fstacks)
        item.Write()
    fstacks.Close()

    # file closure
    f_info.close()
    fdata.Close()  


main()
