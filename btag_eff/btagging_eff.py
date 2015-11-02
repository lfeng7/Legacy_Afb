##############################################################################################################
##                                      Imports and Other Header Stuff                                      ##
##############################################################################################################
# This code will read in selected events root file and generate b-tagging efficiency root files. The efficiency files are
# the input of tagger.
# Hey new version! 
# 11-1-15
from utility import *
from optparse import OptionParser

parser = OptionParser()


############################################
#            Job steering                  #
############################################

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 
parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='inputFiles',
                  help='Input files')

parser.add_option('--txtfiles', metavar='F', type='string', action='store',
                  default = "",
                  dest='txtfiles',
                  help='Input txt files')

parser.add_option('--applyHLT', metavar='F', type='string', action='store',
                  default = "no",
                  dest='applyHLT',
                  help='If apply HLT as first selection cut.')

parser.add_option('--maxevts', metavar='F', type='int', action='store',
                  default = -1,
                  dest='maxevts',
                  help='max number of input ntuple files')

parser.add_option('--type', metavar='F', type='string', action='store',
                  default = 'regular',
                  dest='run_type',
                  help='the type of the run')

parser.add_option('--fakelep', metavar='F', type='string', action='store',
                  default = 'no',
                  dest='fakelep',
                  help='If select fake leptons')

(options, args) = parser.parse_args()

argv = []

# Get input files
prepend = 'selected_files/regular/all/'   # dir of output files to make histograms
postfix='_selected.root'
input_tree_name = 'selected'
# Set up output root file dir
template_type = options.run_type
output_dir = 'btagging_efficiency_files/'
 
# initialization and declaration

#### Set up MC input files
# stucture of a list of files
#    0,         1         2        3         4           
# (filepath, nevts_gen, xsec_NLO, type, nevts_total_ntuple)
flist = []
# Single Top
flist.append(['T_s','singletop',259961,3.79,259176] )
flist.append(['T_t','singletop',3758227,56.4,3748155] )
flist.append(['T_tW','singletop',497658,11.1,495559])
flist.append(['Tbar_s','singletopbar',139974, 1.76,139604])
flist.append(['Tbar_t','singletopbar',1935072, 30.7,1930185])
flist.append(['Tbar_tW','singletopbar',493460,11.1,491463])
# Wjets
flist.append(['W1JetsToLNu_TuneZ2Star_8TeV','wjets',23141598,6662.8,23038253])
flist.append(['W2JetsToLNu_TuneZ2Star_8TeV','wjets',34044921,2159.2,33993463])
flist.append(['W3JetsToLNu_TuneZ2Star_8TeV','wjets',15539503,640.4,15507852])
flist.append(['W4JetsToLNu_TuneZ2Star_8TeV','wjets',13382803,246.0,13326400])
# DYjets
flist.append(['DY1JetsToLL_M','zjets',24045248,660.6,23802736])
flist.append(['DY2JetsToLL_M','zjets',2352304,215.1,2345857])
flist.append(['DY3JetsToLL_M','zjets',11015445,65.79,10655325])
flist.append(['DY4JetsToLL_M','zjets',6402827,28.59,5843425])
# signal
flist.append(['TT_CT10_TuneZ2star_8TeV','ttbar',21675970,245.9,21560109])

def main():
    print 'Making btagging efficiency files.'
    MakeBtaggingEfficiency()

def MakeBtaggingEfficiency():
    fout_postfix = '_CSVM_bTaggingEfficiencyMap.root'
    # Get list of selected files for each type of events
    all_types = ['singletop','singletopbar','wjets','zjets','ttbar']

    all_type_names = ['T_star-channel_TuneZ2star_8TeV-powheg-tauola','Tbar_star-channel_TuneZ2star_8TeV-powheg-tauola']  
    all_type_names+= ['WnJetsToLNu_TuneZ2Star_8TeV-madgraph','DYnJetsToLL_M-50_TuneZ2Star_8TeV-madgraph']
    all_type_names+= ['TT_CT10_TuneZ2star_8TeV-powheg-tauola']

    selected_files = []  
    for itype in all_types:
        selected_files.append([ifile[0] for ifile in flist if ifile[1] == itype ])

    # Loop over type of samples
    for itype in range(len(all_types)):
        # Make output file
        # will save to btagging_efficiency_files/template_type/event_type_control_plots.root
        typename = all_type_names[itype]        
        MakeDirectory(output_dir+template_type)
        foutname = output_dir+template_type+'/'+typename+fout_postfix
        fout = ROOT.TFile(foutname,'recreate')  
        print 'Making:',foutname
        ################################################################
        #           Set up 2D histograms for efficiency files          #
        ################################################################
        # constants
        pt_min = 10
        pt_max = 500
        eta_min = 0
        eta_max = 2.4
        CSVM = 0.679
        # Binning
        # bins_pt = array('d',[0.5,30.0,40.0,50.0,60.0,70.0,80.0,100.0,120.0,160.0,210.0,260.0,320.0,400.0,500.0,600.0,4000.0]) # For signal
        if all_types[itype] == 'ttbar':
            bins_pt_b = array('d',[0., 40., 60., 80., 100., 150., 200., 300., 400., 500., 1000.]) # For ttbar # the same as Sal's code
            bins_eta_b = array('d',[0., 0.6, 1.2, 2.4])
            bins_pt_c = array('d',[0., 40., 60., 80., 100., 150., 200., 300., 400., 1000.]) # For ttbar # the same as Sal's code
            bins_eta_c = array('d',[0., 0.6, 1.2, 2.4])
            bins_pt_udsg = array('d',[0., 40., 60., 80., 100., 150., 200., 300., 400., 1000.]) # For ttbar # the same as Sal's code
            bins_eta_udsg = array('d',[0., 0.6, 1.2, 2.4])
        elif all_types[itype] in ['singletop','singletopbar'] :
            bins_pt_b = array('d',[0., 40., 60., 80., 100., 150., 200., 300., 1000.]) # For ttbar # the same as Sal's code
            bins_eta_b = array('d',[0., 0.6, 1.2, 2.4])
            bins_pt_c = array('d',[0., 40., 60., 80., 100., 150., 1000.]) # For ttbar # the same as Sal's code
            bins_eta_c = array('d',[0., 0.6, 1.2, 2.4])
            bins_pt_udsg = array('d',[0., 40., 60., 80., 100., 150., 1000.]) # For ttbar # the same as Sal's code
            bins_eta_udsg = array('d',[0., 0.6, 1.2, 2.4]) 
        else :
            bins_pt_b = array('d',[0., 40., 60., 80., 100., 150., 1000.]) # For ttbar # the same as Sal's code
            bins_eta_b = array('d',[0., 0.6, 1.2, 2.4])
            bins_pt_c = array('d',[0., 40., 60., 80., 100., 150., 1000.]) # For ttbar # the same as Sal's code
            bins_eta_c = array('d',[0., 0.6, 1.2, 2.4])
            bins_pt_udsg = array('d',[0., 40., 60., 80., 100., 150., 1000.]) # For ttbar # the same as Sal's code
            bins_eta_udsg = array('d',[0., 0.6, 1.2, 2.4])                       

        pt_binning = [bins_pt_b,bins_pt_c,bins_pt_udsg]
        eta_binning = [bins_eta_b,bins_eta_c,bins_eta_udsg]

        pt_bin = []
        # len(bins_pt)-1 # this is what should be
        eta_bin = []
        # len(bins_eta)-1
        for i in range(len(pt_binning)) :
            pt_bin.append(len(pt_binning[i])-1)
            eta_bin.append(len(eta_binning[i])-1)
        # Create 2D histograms
        h_name_num = ['b_num','c_num','udsg_num']
        h_name_denom = ['b_denom','c_denom','udsg_denom']
        h_name_eff = ['efficiency_b','efficiency_c','efficiency_udsg']
        eff_num_list = []
        eff_denom_list = []
        eff_list = []
        for i in range(len(h_name_num)):
            eff_num_list.append(ROOT.TH2D(h_name_num[i],h_name_num[i],pt_bin[i],pt_binning[i],eta_bin[i],eta_binning[i]))
            eff_denom_list.append(ROOT.TH2D(h_name_denom[i],h_name_denom[i],pt_bin[i],pt_binning[i],eta_bin[i],eta_binning[i]))
        # Set the address of the histogram as the output root file
        for i in range(len(eff_num_list)):
            eff_num_list[i].SetDirectory(fout)
            eff_denom_list[i].SetDirectory(fout)

        ################################################################            
        #           Loop over input files in current type              #
        ################################################################
        for ifile in selected_files[itype] :
            tmp_name =  prepend+ifile+postfix
            print 'Processing',tmp_name            
            tmp_f = ROOT.TFile(tmp_name)
            tmp_tree = tmp_f.Get(input_tree_name)   
            # Loop over entries
            nev = tmp_tree.GetEntries()
            # print 'num entries is',nev
            jets_count,evt_count,num_all_jets = 0 , 0 , 0
            for iev in range(nev):
                evt_count += 1
                if iev == options.maxevts:
                    break
                # Report progress
                if iev%20000 == 1 :
                    # print 'finishing event # ',iev
                    print 'Progress ',int(100.*iev/nev),'%'
                tmp_tree.GetEntry(iev)

                njets_denom = 0
                njets_num = 0

                num_all_jets+=len(tmp_tree.jets_pt)

                for j in range(len(tmp_tree.jets_pt)):
                    j_csv = tmp_tree.jets_csv[j]
                    j_PID = abs(tmp_tree.jets_flavor[j])
                    j_pt  = tmp_tree.jets_pt[j]
                    j_eta = abs(tmp_tree.jets_eta[j])
                    # Skip empty events
                    if not ( j_pt > 0 and j_PID > 0 ):
                        continue
                    # debug
                    if j_pt == 0 :
                        print 'pt and PID '+ str(j_pt)+','+str(j_PID)
                    jets_count += 1
                    # Fill numerator
                    if j_PID==5 and j_csv>CSVM:
                        eff_num_list[0].Fill(j_pt,j_eta)
                    if j_PID==4 and j_csv>CSVM:
                        eff_num_list[1].Fill(j_pt,j_eta)
                    if j_PID!=4 and j_PID!=5 and j_csv>CSVM:
                        eff_num_list[2].Fill(j_pt,j_eta)
                    # Fill denominator
                    if j_PID==5:
                        eff_denom_list[0].Fill(j_pt,j_eta)
                    if j_PID==4:
                        eff_denom_list[1].Fill(j_pt,j_eta)
                    if j_PID!=5 and j_PID!=4:
                        eff_denom_list[2].Fill(j_pt,j_eta)
                # print 'Number of valid jets is : '+str(valid_jets) # For debug
            print 'Total number of events : '+str(evt_count)
            print 'Total number of jets :',num_all_jets
            print 'Total number of good Jets   : '+str(jets_count)
            print 'Average number of jets per event is : %.2f'%(float(jets_count)/float(evt_count))+'\n'
            tmp_f.Close()
        # Finish Loop over files in current sample type

        # Get efficiency
        for i in range(len(eff_num_list)):
            h_tmp = eff_num_list[i].Clone()
            h_tmp.SetDirectory(fout)
            h_tmp.SetName(h_name_eff[i])
            h_tmp.Divide(eff_denom_list[i])
            eff_list.append(h_tmp)
        # Save rootfiles
        fout.Write()
        fout.Close()
        print 'Finished making efficiency file for sample',all_types[itype],'\n'

    # Finish loop over sample types
    print 'All done!'

main()
