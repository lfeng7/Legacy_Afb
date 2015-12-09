# This small macro will read in all edm files in a directory and count the total number of events 
# v2. Will take a ttree, make histograms, and stack with data. And calculate correction weights for MC

import ROOT
from optparse import OptionParser
import os

# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

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
(options, args) = parser.parse_args()

argv = []

def main():
    global xaxis_name 

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

    if xaxis_name == '': xaxis_name = var

    dir_hists = 'histo_files/'
    data_lumi = 19700

    # Get input MC and data files according to txt file
    txt_MC = open(rundir+'/MC_input_with_bkg.txt')
    mc_samples = []
    for line in txt_MC:
        items = line.split()
        # TT_CT10_qq                      qq      21560109        245.9   
        mc_samples += [(items[1],items[2],int(items[3]),float(items[4]))]

    # Get data histogram
    dir_hists =rundir+'/'+dir_hists
    fdata = ROOT.TFile(dir_hists+'all_data_histos.root')
    keys = fdata.GetListOfKeys()
    for ikey in keys:
        if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
    print 'Getting ttree',treename
    ttree_data = fdata.Get(treename)
    # Make histogram
    hname_data = hname+'_data'
    h_data = ROOT.TH1F(hname_data, hname_data, bin, xmin, xmax)  
    ttree_data.Draw(var+">>"+hname_data,""+ cut, "goff")
    h_data.SetDirectory(0)
    fdata.Close()

    # Loop over all MC files to make MC stack and legend
    mc_stack = ROOT.THStack(hname+'_stack',var+' comparison')
    sample_types = []
    leg = ROOT.TLegend(0.85,0.65,1.0,1.0)
    hlist_mc = []
    for i,isample in enumerate(mc_samples):
        tmpf = ROOT.TFile(dir_hists+isample[0]+'_histos.root')
        keys = tmpf.GetListOfKeys()
        for ikey in keys:
            if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()   
        ttree_mc = tmpf.Get(treename)
        # Make histogram
        hname_mc = hname+'_'+isample[0]
        h_mc = ROOT.TH1F(hname_mc, hname_mc, bin, xmin, xmax) 
        if cut == '': 
            ttree_mc.Draw(var+">>"+hname_mc,weight, "goff")
        else :
            ttree_mc.Draw(var+">>"+hname_mc,'('+cut+')*('+weight+')', "goff")
        h_mc.SetDirectory(0)    
        cross_section_NLO = isample[3]
        nevts_gen = isample[2]
        w_scale = data_lumi*cross_section_NLO/nevts_gen 
        h_mc.Scale(w_scale)
        # Add into stack
        sample_type = isample[1]
        icolor = GetSampleColor(sample_type)
        h_mc.SetFillColor(icolor)
        h_mc.SetLineColor(icolor)
        h_mc.SetMarkerStyle(0)
        h_mc.SetYTitle('events')
        hlist_mc.append(h_mc)
        mc_samples[i] += (h_mc,)
        # Determine if we want to add an entry to legend
        if sample_type not in sample_types: 
            leg.AddEntry(h_mc,sample_type,"F")
            sample_types.append(sample_type)
        tmpf.Close()
    # Add hists into stack in certain order
    alltypes = ['tt_bkg','bck','WJets','gg','qq']
    for itype in alltypes :
        isamples = [item for item in mc_samples if item[1]==itype]
        for i,item in enumerate(isamples):
            tmp_h = item[4]
            if i == len(isamples)-1 : tmp_h.SetLineColor(1)
            mc_stack.Add(tmp_h)


        # for i,isample in enumerate(mc_samples):
        #     sample_type = isample[1]
        #     if sample_type == itype :
        #         mc_stack.Add(hlist_mc[i])     


    # Make data/MC comparison plot
    leg.AddEntry(h_data,'data')
    c_plot = comparison_plot_v1(mc_stack,h_data,leg,hname)
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
    fout.Write()
    fout.Close()

def GetSampleColor(itype):
    if itype=='gg' : return 38
    if itype=='qq' : return 46
    if itype=='bck' : return 41
    if itype=='WJets': return ROOT.kGreen-3
    if itype=='tt_bkg': return ROOT.kMagenta
    return 0

# This is specifically for comparing the stacked MC plots with data adding residule plots too
import math
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )

def comparison_plot_v1(mc_,data_,legend,event_type='plots',logy=False):
    global fout
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
        name = plotdir+event_type+'_'+mc_.GetName()+'_compare_log.png'
    else :
        name = plotdir+event_type+'_'+mc_.GetName()+'_compare.png'
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
    data_.GetYaxis().SetTitleOffset(1.0)
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
    mc_.Draw('h');
    mc_.GetYaxis().SetTitle("events")
    mc_.GetYaxis().SetTitleOffset(1.0)
    mc_.GetXaxis().SetLabelOffset(999)   
    data_.Draw('SAME PE1X0'); 
    legend.Draw()
    x_resid_pad.cd(); 
    h_res.Draw('PE1X0 '); xline.Draw()
    # h_res.GetXaxis().SetTitle(xaxis_name)

    c1.Update()    
    # Saving
    c1.SaveAs(name)
    c1.Write()
    # fout.Close()
    return (c1)

main()

