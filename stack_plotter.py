# This small macro will do stacker plot and compare with data 
# just give var and working directory, with a txt file containing MC channel name, type, xsec and nevts_gen
# will make plots/ in current dir and put all .png and .root file in plots/
# 12-7-2015

import ROOT
from optparse import OptionParser
import os
import glob
from plot_tools import *

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
                  default='',
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


(options, args) = parser.parse_args()

argv = []

#canvas_title = 'CMS Private Work, 19.7 fb^{-1} at #sqrt{s} = 8 TeV'


def main():
    global xaxis_name , fout , canvas_title
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
    txt_MC = open(rundir+'/MC_input_with_bkg.txt')
    mc_samples = []
    for line in txt_MC:
        items = line.split()
        if '#' in items[0] or len(items)<3 : continue
        # TT_CT10_qq                      qq      21560109        245.9   
        mc_samples += [(items[1],items[2],int(items[3]),float(items[4]))]

    # Get data histogram
    dir_hists =rundir+'/'+dir_hists
    all_inputfiles = glob.glob(dir_hists+'*.root')
    # find data root file
    for ifile in all_inputfiles:
        if 'data' in ifile or 'Run' in ifile : data_path = ifile    
    fdata = ROOT.TFile(data_path)
    keys = fdata.GetListOfKeys()
    for ikey in keys:
        if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
    print 'Getting ttree',treename
    ttree_data = fdata.Get(treename)
    # Make histogram
    hname_data = hname+'_data'    
    if xmin != xmax:
        h_data = ROOT.TH1F(hname_data, hname_data, bin, xmin, xmax)        
    ttree_data.Draw(var+">>"+hname_data,""+ cut, "goff")
    # if we want overflow bin
    if plot_overflow == 'yes':
        h_data = overflow(h_data) 

    h_data.SetDirectory(0)
    fdata.Close()

    # Loop over all MC files to make MC stack and legend
    mc_stack = ROOT.THStack(hname+'_stack',var+' comparison')
    sample_types = []
    leg = ROOT.TLegend(0.7550143,0.4991304,0.9054441,0.8469565)
    leg.SetName('legend')
    hlist_mc = []
    for i,isample in enumerate(mc_samples):
        # find input root file
        mc_path = ''
        for ifile in all_inputfiles:
            if isample[0] in ifile  : mc_path = ifile   
        if mc_path == '':
            print isample[0],'cannot be found! Will skip this sample!'
            continue
        tmpf = ROOT.TFile(mc_path)
        keys = tmpf.GetListOfKeys()
        for ikey in keys:
            if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()   
        ttree_mc = tmpf.Get(treename)
        # Make histogram
        hname_mc = hname+'_'+isample[0]
        if xmin!=xmax:
            h_mc = ROOT.TH1F(hname_mc, hname_mc, bin, xmin, xmax) 
        if cut == '': 
            ttree_mc.Draw(var+">>"+hname_mc,weight, "goff")
        else :
            ttree_mc.Draw(var+">>"+hname_mc,"(%s)*(%s)"%(cut,weight), "goff")
        # decide if we want a last bin for overflow
        if plot_overflow == 'yes':
            h_mc = overflow(h_mc) 

        h_mc.SetDirectory(0)    
        cross_section_NLO = isample[3]
        nevts_gen = isample[2]
        w_scale = data_lumi*cross_section_NLO/nevts_gen 
        h_mc.Scale(w_scale)
        print 'sample name %s norm_w = %.3f'%(hname_mc,w_scale)

        # Add into stack
        sample_type = isample[1]
        icolor = GetSampleColor(sample_type)
        h_mc.SetFillColor(icolor)
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
            leg.AddEntry(h_mc,GetSampleName(sample_type),"F")
            sample_types.append(sample_type)
        tmpf.Close()
    # Add hists into stack in certain order
    alltypes = ['qcd','bck','zjets','WJets','singletop','tt_bkg','gg','qq','signal']
    for itype in alltypes :
        isamples = [item for item in mc_samples if item[1]==itype]
        if len(isamples)>0:
            print 'Include MC sample type %s'%itype
        for i,item in enumerate(isamples):
            tmp_h = item[4]
            if i == len(isamples)-1 : tmp_h.SetLineColor(1)
            mc_stack.Add(tmp_h)
    # Write out yields
    if makeyields == 'yes':
        table_yields = []
        f_yds = open('plots/yields.csv','w')
        nevts_data = int(h_data.Integral())
        nevts_mc_total = 0
        for item in mc_samples:
            tmp_h = item[4]
            tmp_type = item[1]
            tmp_name = item[0]
            nevts = tmp_h.Integral()
            nevts_mc_total += int(nevts)
            table_yields +=[(tmp_name,tmp_type,int(nevts))]
        for itype in alltypes:
            tmp_list = [item[2] for item in table_yields if item[1]==itype]
            if len(tmp_list)==0 : continue
            nevts_mc = sum(tmp_list)
            per_mc = nevts_mc*1.0/nevts_mc_total
            f_yds.write(itype+'    '+str(nevts_mc)+'    %.3f\n'%per_mc)
        f_yds.write('mc      '+str(int(nevts_mc_total))+'    1.0\n')
        f_yds.write('data    '+str(nevts_data)+'     1.0')
        f_yds.close()

        # for i,isample in enumerate(mc_samples):
        #     sample_type = isample[1]
        #     if sample_type == itype :
        #         mc_stack.Add(hlist_mc[i])     

    # Make data/MC comparison plot
    leg.AddEntry(h_data,'data')
    if var != 'charge_ratio':
        c_plot = comparison_plot_v1(mc_stack,h_data,leg,hname)
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
    fout.cd()
    mc_stack.Write()
    h_data.Write()
    leg.Write()
    fout.Write()
    # fout.Close()

def GetSampleColor(itype):
    if itype in ['gg','w1jet','z1jet']                          : return 38
    if itype in ['qq','signal','w2jet','z2jet']             : return 46
    if itype in ['bck','bkg','w3jet','z3jet']               : return 41
    if itype in ['singletop','T+x','w4jet','z4jet']         : return 41
    if itype=='WJets'                       : return ROOT.kGreen-3
    if itype=='tt_bkg'                      : return ROOT.kMagenta
    if itype=='qcd'                         : return ROOT.kYellow
    if itype in ['zjets','DY']              : return ROOT.kBlue
    return 0

def GetSampleName(itype):
    if itype=='gg'                  : return 'gg/qg->t#bar{t}'
    if itype=='qq'                  : return 'q#bar{q}->t#bar{t}'
    if itype in ['tt','signal']     : return 't#bar{t}'
    if itype in ['bck','bkg']       : return 'other bkg'
    if itype=='WJets'               : return 'W+jets'
    if itype=='tt_bkg'              : return 'dilep/had t#bar{t}'
    if itype=='qcd'                 : return 'QCD MJ'
    if itype in ['singletop','T+x'] : return 'single top'
    if itype in ['zjets','DY']      : return 'z+jets'
    return 0

# This is specifically for comparing the stacked MC plots with data adding residule plots too
import math
def comparison_plot_v1(mc_,data_,legend,event_type='plots',draw_option = 'h',logy=False):
    global fout,var
    prefix = 'plots'
    # check if plotting dir is made. If not , make it now
    if not os.path.exists(prefix):
        os.mkdir(prefix)
        print 'Making '+prefix
    # Set the dir to put all plots
    plotdir = prefix+'/'
    fout = ROOT.TFile(plotdir+event_type+'_plots.root','recreate')

    # plotting
    c1 = ROOT.TCanvas(data_.GetName()+'_compare')
    if logy == "log" : 
        c1.SetLogy()
        name = plotdir+event_type+'_compare_log.png'
    else :
        name = plotdir+event_type+'_compare.png'
    # Find the max of both histograms
    max_mc = mc_.GetMaximum()
    max_data = data_.GetMaximum()
    max_ = max(max_mc,max_data)*1.1
    data_.SetMaximum(max_)
    mc_.SetMaximum(max_)
    # Calculating the residual of each bin
    h_data = data_
    if mc_.ClassName() == 'THStack':
        h_stack = mc_.GetStack().Last() # This is the combined histogram in stack
    else : 
        h_stack = mc_
    # Make residual histogram
    h_res = ROOT.TH1D(event_type+'_residuals',';; Data/MC',h_data.GetNbinsX(),h_data.GetXaxis().GetXmin(),h_data.GetXaxis().GetXmax())
    # h_res.GetXaxis().SetName(h_data.GetXaxis().GetName())
    # h_res.SetXTitle(h_data.GetXaxis().GetName())

    maxxdeviations = 0.0
    for ibin in range(h_data.GetNbinsX()):
        databin = h_data.GetBinContent(ibin)
        mcbin = h_stack.GetBinContent(ibin)
        # Calculate residual
        if mcbin != 0 and databin != 0:
            res = databin*1.0/mcbin
            # Calculate error of residual, delta(res) = residual*sqrt(1/data+1/mc)
            res_err = res*math.sqrt(1.0/databin+1.0/mcbin)
            # Find maximum residual
            maxxdeviations = max(maxxdeviations,max(abs(res+res_err-1.0),abs(res-res_err-1.0)))
        else :
            res = 0 ; res_err = 0
        # Set residual histograms
        h_res.SetBinContent(ibin,res)
        h_res.SetBinError(ibin,res_err)
    # print 'maxxdeviations',maxxdeviations
    # Setup residual histograms
    h_res.SetStats(0)
    h_res.GetXaxis().SetLabelSize((0.05*0.72)/0.28); h_res.GetXaxis().SetTitleOffset(0.8)
    h_res.GetYaxis().SetLabelSize((0.05*0.72)/0.28); h_res.GetYaxis().SetTitleOffset(0.4)
    h_res.GetXaxis().SetTitleSize((0.72/0.28)*h_res.GetXaxis().GetTitleSize())
    h_res.GetYaxis().SetTitleSize((0.72/0.28)*h_res.GetYaxis().GetTitleSize())
    maxx = 1.0+1.1*maxxdeviations
    minx = 1.0-1.1*maxxdeviations
    h_res.GetYaxis().SetRangeUser(minx,maxx)
    h_res.GetYaxis().SetNdivisions(503)
    h_res.GetXaxis().SetTitle(xaxis_name)
    # cosmetics
    h_res.SetLineStyle(0)
    h_res.SetMarkerStyle(20)
    h_res.SetMarkerSize(0.5) 

    # Some cosmetics for data
    data_.GetYaxis().SetTitle('Events')
    data_.GetYaxis().SetTitleOffset(1.2)
    data_.SetMarkerStyle(20)
    data_.SetMarkerSize(0.5)
    data_.SetLineStyle(0)

    #Build the lines that go at 1 on the residuals plots
    xline = ROOT.TLine(h_res.GetXaxis().GetXmin(),1.0,h_res.GetXaxis().GetXmax(),1.0); xline.SetLineWidth(2); xline.SetLineStyle(2)
    #plot stacks with data overlaid and residuals. Totally stole from Nick :)
    c1.cd()
    channame = event_type
    # Make and adjust pads
    x_histo_pad=ROOT.TPad(channame+'_x_histo_pad',channame+'_x_histo_pad',0,0.25,1,1)
    if logy == 'log': x_histo_pad.SetLogy()
    x_resid_pad=ROOT.TPad(channame+'_x_residuals_pad',channame+'_x_residuals_pad',0,0,1.,0.25)
    x_histo_pad.SetCanvas(c1); x_resid_pad.SetCanvas(c1)
    x_histo_pad.SetLeftMargin(0.16); x_histo_pad.SetRightMargin(0.05) 
    x_histo_pad.SetTopMargin(0.11);  x_histo_pad.SetBottomMargin(0.02)
    x_histo_pad.SetBorderMode(0)
    x_resid_pad.SetLeftMargin(0.16); x_resid_pad.SetRightMargin(0.05)
    x_resid_pad.SetTopMargin(0.0);   x_resid_pad.SetBottomMargin(0.3)
    x_resid_pad.SetBorderMode(0)
    x_resid_pad.Draw(); x_histo_pad.Draw()
    x_histo_pad.cd(); 
    mc_.Draw(draw_option);
    mc_.GetYaxis().SetTitle("events")
    mc_.GetYaxis().SetTitleOffset(1.0)
    mc_.GetXaxis().SetLabelOffset(999)   
    mc_.SetTitle(canvas_title)
    data_.Draw('SAME PE1X0'); 
    # Draw data stat box
    ROOT.gStyle.SetOptStat("e");
    data_.SetStats(1)
    data_.Draw('SAMEs PE1X0'); 
    x_histo_pad.Update()
    statbox1 = data_.FindObject("stats")
    statbox1.SetLineColor(0)
    statbox1.Draw('sames')

    legend.SetShadowColor(0)
    legend.SetLineColor(0)
    legend.Draw()
    x_resid_pad.cd(); 
    h_res.Draw('PE1X0 '); xline.Draw()
    h_res.GetXaxis().SetTitle(xaxis_name)

    obj_title = c1.FindObject("title")
    obj_title.SetShadowColor(0)
    obj_title.SetLineColor(0)
    c1.Update()    
    # Saving
    c1.SaveAs(name)
    c1.Write()
    # fout.Close()
    return (c1)

main()

