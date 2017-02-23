# This small macro will do stacker plot and compare with data 
# just give var and working directory, with a txt file containing MC channel name, type, xsec and nevts_gen
# will make plots/ in current dir and put all .png and .root file in plots/
# 12-7-2015

import ROOT
from optparse import OptionParser
import os
import glob
from helper.plot_tools import overflow
from helper import helper


# pavetex.SetShadowColor(0)
# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

parser.add_option('--plot', action='store_true', default=False,
                  dest='plot',
                  help='plot interactively')

parser.add_option('--var', metavar='F', type='string', action='store',
                  default = "",
                  dest='var',
                  help='var to plot')

parser.add_option('--cut', metavar='F', type='string', action='store',
                  default="",
                  dest='cut',
                  help='')

parser.add_option('--Min', metavar='F', type='float', action='store',
                  default=0,
                  dest='Min',
                  help='')
parser.add_option('--Max', metavar='F', type='float', action='store',
                  default=0,
                  dest='Max',
                  help='')

parser.add_option('--name', metavar='F', type='string', action='store',
              default = "blank",
                  dest='name',
                  help='')
parser.add_option('--log', action='store_true', default=False,
                  dest='log',
                  help='log sacle on y')

parser.add_option('--bin', metavar='F', type='int', action='store',
                  default=100,
                  dest='bin',
                  help='')

parser.add_option('--weight', metavar='F', type='string', action='store',
                  default='1',
                  dest='weight',
                  help='which event weight to use for MC')

parser.add_option('--xaxis', metavar='F', type='string', action='store',
              default = "",
                  dest='xaxis',
                  help='')
parser.add_option('--yaxis', metavar='F', type='string', action='store',
              default = "",
                  dest='yaxis',
                  help='')

parser.add_option('--title', metavar='F', type='string', action='store',
              default = "",
                  dest='title',
                  help='')

parser.add_option('--dir', metavar='F', type='string', action='store',
              default = "",
                  dest='dir',
                  help='')

parser.add_option('--yields', metavar='F', type='string', action='store',
              default = "no",
                  dest='yields',
                  help='If you want to make a yields table')

parser.add_option('--overflow', metavar='F', type='string', action='store',
              default = "yes",
                  dest='overflow',
                  help='If you want to plot overflow bin')

parser.add_option('--MCinfo', metavar='F', type='string', action='store',
              default = 'MC_input.txt',
                  dest='MCinfo',
                  help='MCinfo.txt file')

parser.add_option('--verbose', action='store_true',
              default = False,
                  dest='verbose',
                  help='addional output')


(options, args) = parser.parse_args()

argv = []

#canvas_title = 'CMS Private Work, 19.7 fb^{-1} at #sqrt{s} = 8 TeV'

h_qcd = []

QCD_SF = 0.13

def main():
    global xaxis_name , fout , canvas_title, h_qcd
    plot = options.plot
    cut = options.cut
    var = options.var
    xmin = options.Min
    xmax = options.Max
    dolog = options.log
    bin = options.bin
    hname = options.name
    htitle = options.title
    xaxis_name = options.xaxis
    yaxis_name = options.yaxis
    rundir = options.dir
    weight = options.weight
    makeyields = options.yields
    plot_overflow = options.overflow
    verbose = options.verbose 

    if htitle != '':
        canvas_title = htitle
    else :
        canvas_title = hname

    # Some global root style 
    ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
    if not plot: ROOT.gROOT.SetBatch(True)

    if xaxis_name == '': xaxis_name = var

    dir_hists = '' # 'histo_files/'
    data_lumi = 19700

    # Get input MC and data files according to txt file
    txt_MC = open(options.MCinfo)
    mc_samples = []
    for line in txt_MC:
        if line.startswith('#'): continue
        items = line.split()
#        if '#' in items[0] or len(items)<3 : continue
        # TT_CT10_qq                      qq      21560109        245.9   
        mc_samples += [(items[1],items[2],int(items[3]),float(items[4]))]

    # Get data histogram
    dir_hists =rundir+'/'+dir_hists
    all_inputfiles = glob.glob(dir_hists+'*.root')
    # find data root file
    data_path = []
    for ifile in all_inputfiles:
        if 'data' in ifile or 'Run' in ifile : 
            data_path.append(ifile)
    if not len(data_path)>0 :
        print 'Cannot find any data root file! Will exit.'
        sys.exit(1)   
    fdata = ROOT.TFile(data_path[0])
    keys = fdata.GetListOfKeys()
    for ikey in keys:
        if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
    print 'Getting ttree',treename
    fdata.Close()
    # put all data root files in a TChain
    ttree_data = ROOT.TChain(treename)
    for item in data_path:
        ttree_data.Add(item)

    # Make histogram
    hname_data = hname+'_data'    
    if xmin != xmax:
        h_data = ROOT.TH1F(hname_data, hname_data, bin, xmin, xmax)        
    ttree_data.Draw(var+">>"+hname_data,""+ cut, "goff")
    # Print out cut efficiency
    total_data = ttree_data.GetEntries()
    selected_data = ttree_data.GetEntries(cut)
    data_cut_eff = selected_data*1.0/total_data
    print 'Num Entries in data: %i, selected entries %i, cut efficiency %.3f'%(total_data,selected_data,data_cut_eff)
    # if we want overflow bin
    if plot_overflow == 'yes':
        h_data = overflow(h_data) 

    h_data.SetDirectory(0)
    data_integral = h_data.Integral()
    fdata.Close()

    # Loop over all MC files to make MC stack and legend
    mc_stack = ROOT.THStack(hname+'_stack',var+' comparison')
    sample_types = []
    leg = ROOT.TLegend(0.7550143,0.4991304,0.9054441,0.8469565)
    leg.SetName('legend')
    hlist_mc = []
    for i,isample in enumerate(mc_samples):
        sample_type = isample[1]
        # find input root file
        mc_path = ''
        for ifile in all_inputfiles:
            if isample[0] in ifile  : mc_path = ifile   
        if mc_path == '':
            print isample[0],'cannot be found! Will skip this sample!'
            continue   
        tmpf = ROOT.TFile(mc_path)
        treename = helper.GetTTreeName(tmpf)
        ttree_mc = tmpf.Get(treename)
        # Get nominal correction weights
        corr_w = helper.GetWeightNames(ttree_mc)
        # Get valid weights in arg
        val_w = helper.GetValidWeight(ttree_mc,weight)
        if val_w != '':
            tmp_weight = '%s*%s'%(weight,corr_w)
        else:
            tmp_weight = corr_w
        # hard coded top_pT reweight for now ...
        if 'CT10' in mc_path and 'TT' in mc_path: 
            if ttree_mc.FindBranch('top_pT_reweight'):
                tmp_weight = '%s*top_pT_reweight'%tmp_weight
            if ttree_mc.FindBranch('weight_top_pT'):
                tmp_weight = '%s*weight_top_pT'%tmp_weight            
        if verbose: print '(info) weight: %s'%tmp_weight
        # Make histogram
        hname_mc = hname+'_'+isample[0]
        if xmin!=xmax:
            h_mc = ROOT.TH1F(hname_mc, hname_mc, bin, xmin, xmax) 
        if cut == '': 
            ttree_mc.Draw(var+">>"+hname_mc,tmp_weight, "goff")
        else :
            ttree_mc.Draw(var+">>"+hname_mc,"(%s)*(%s)"%(cut,tmp_weight), "goff")
        # decide if we want a last bin for overflow
        if plot_overflow == 'yes':
            h_mc = overflow(h_mc) 

        h_mc.SetDirectory(0)    
        cross_section_NLO = isample[3]
        nevts_gen = isample[2]
        w_scale = data_lumi*cross_section_NLO/nevts_gen 
        # here w_scale uncertainty can be implemented

        # special normalization for QCD MC. 
        if sample_type == 'qcd':
            w_scale *=  QCD_SF
        h_mc.Scale(w_scale)
        print '%s, norm_w = %.3f'%(isample[0],w_scale)

        # Add into stack
        icolor = helper.getColors(sample_type)
        h_mc.SetFillColor(icolor)
        if sample_type == 'qcd':
            h_mc.SetFillColor(0)
            h_qcd.append(h_mc)
        h_mc.SetLineColor(icolor)
        h_mc.SetMarkerStyle(0)
        h_mc.SetYTitle('events')
        # Special setting for charge ratio plot
        if var == 'charge_ratio':
            h_mc.SetBarWidth(0.70);
            h_mc.SetBarOffset(0.15);
            h_mc.SetMarkerStyle();
            h_mc.SetMarkerColor(icolor)
        hlist_mc.append(h_mc)
        mc_samples[i] += (h_mc,)
        # Determine if we want to add an entry to legend
        if sample_type not in sample_types: 
            print '(info)Adding %s\n'%sample_type 
            leg.AddEntry(h_mc,GetSampleName(sample_type),"F")
            sample_types.append(sample_type)
        tmpf.Close()
    # Add hists into stack in certain order
    alltypes = ['bck','zjets','WJets','singletop','tt_bkg','other_bkg','gg','qq','signal','ttbar']
    for itype in alltypes :
        isamples = [item for item in mc_samples if item[1]==itype]
        if len(isamples)>0:
            print 'Include MC sample type %s'%itype
        else:
            continue
        for i,item in enumerate(isamples):
            if len(item)<5: continue
            tmp_h = item[4]
            if i == len(isamples)-1 : tmp_h.SetLineColor(1)
            mc_stack.Add(tmp_h)
    # Write out yields
    yield_types = ['ttbar','qq','gg','WJets','other_bkg','qcd']
    if makeyields == 'yes':
        table_yields = []
        f_yds = open('yields.csv','w')
        nevts_data = int(h_data.Integral())
        nevts_mc_total = 0
	#print mc_samples
        for item in mc_samples:
            if len(item)<5: continue
            tmp_h = item[4]
            tmp_type = item[1]
            tmp_name = item[0]
            nevts = tmp_h.Integral()
            nevts_mc_total += int(nevts)
            table_yields +=[(tmp_name,tmp_type,int(nevts))]
        for itype in yield_types:
            tmp_list = [item[2] for item in table_yields if item[1]==itype]
            if len(tmp_list)==0 : continue
            nevts_mc = sum(tmp_list)
            per_mc = nevts_mc*100./nevts_mc_total
            f_yds.write('%15s,%15i,%15.1f%%\n'%(itype,nevts_mc,per_mc))
        f_yds.write('%15s,%15i,%15.1f%%\n'%('mc',nevts_mc_total,100.0))
        f_yds.write('%15s,%15i,%15.1f%%\n'%('data',nevts_data,nevts_data*100./nevts_mc_total))
        f_yds.close()

        # for i,isample in enumerate(mc_samples):
        #     sample_type = isample[1]
        #     if sample_type == itype :
        #         mc_stack.Add(hlist_mc[i])     

    # Make data/MC comparison plot
    leg.AddEntry(h_data,'data')
    if var != 'charge_ratio':
        #def comparison_plot_v1(mc_,data_,legend,event_type='plots',bin_type = 'fixed',logy=False,draw_option = 'hist'):
        c_plot,final_hist = helper.comparison_plot_v1(mc_stack,h_data,leg,hname)
    else :
        fout = ROOT.TFile('plots/'+hname+'_'+var+'_plots.root','recreate')        
        c_plot = ROOT.TCanvas()
        mc_stack.Draw('bar1')
        mc_stack.GetXaxis().SetBinLabel(mc_stack.GetXaxis().FindFixBin(1.5),"4jets, l+");
        mc_stack.GetXaxis().SetBinLabel(mc_stack.GetXaxis().FindFixBin(2.5),"4jets, l-");
        mc_stack.GetXaxis().SetBinLabel(mc_stack.GetXaxis().FindFixBin(3.5),"5jets, l+");
        mc_stack.GetXaxis().SetBinLabel(mc_stack.GetXaxis().FindFixBin(4.5),"5jets, l-");
        h_data.Draw('same PE1X0')
        leg.Draw() 
        c_plot.Update()
        c_plot.SaveAs('plots/'+hname+'_compare.png') 

    c_plot.SetName('cplot')
    # Write out
    plotdir = 'plots/'
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)
        os.system('cp ~/index.php '+plotdir)
        print 'Creating new dir '+plotdir
    fout = ROOT.TFile('%s/%s.root'%(plotdir,hname),'recreate')
    fout.cd()
    mc_stack.Write()
    final_hist.Write()
    h_data.Write()
    leg.Write()
    fout.Write()
    # fout.Close()
    c_plot.SaveAs('%s/%s.png'%(plotdir,hname))

def GetSampleName(itype):
    if itype=='gg'                  : return 'gg/qg->t#bar{t}'
    if itype=='qq'                  : return 'q#bar{q}->t#bar{t}'
    if itype in ['ttbar','signal']     : return 't#bar{t}'
    if itype in ['bck','bkg','other_bkg']       : return 'other bkg'
    if itype=='WJets'               : return 'W+jets'
    if itype=='tt_bkg'              : return 'dilep/had t#bar{t}'
    if itype=='qcd'                 : return 'QCD MJ'
    if itype in ['singletop','T+x'] : return 'single top'
    if itype in ['zjets','DY']      : return 'z+jets'
    return 0



main()

