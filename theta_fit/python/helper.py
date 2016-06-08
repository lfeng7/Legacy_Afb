import ROOT
"""
common helper functions
"""

def GetTTreeName(tfile):
    # Find the name of the ttree
    keys = tfile.GetListOfKeys()
    for ikey in keys:
        if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
    print 'Getting ttree',treename
    return treename

def GetListTH1D(tfile):
    # Find the names of all TH1D
    all_th1d = []
    keys = tfile.GetListOfKeys()
    for ikey in keys:
        if  'TH1' in ikey.GetClassName()  : 
            hname = ikey.GetName()
            all_th1d.append(hname)
    print 'Getting these histograms'
    for item in all_th1d:
        print item
    return all_th1d

def getColors(name):
    colors = {'qqs':ROOT.kRed+1,'gg':ROOT.kRed-7,'other':ROOT.kMagenta,'wjets':ROOT.kGreen-3,'qcd':ROOT.kYellow}
    icolor = colors.get(name,0)
    return icolor
    
# This is specifically for comparing the stacked MC plots with data adding residule plots too
import math
def comparison_plot_v1(mc_,data_,legend,event_type='plots',draw_option = 'hist',logy=False):
    # Def axis title
    xaxis_name = data_.GetXaxis().GetTitle()
    canvas_title = data_.GetTitle()
    # create canvas
    c1 = ROOT.TCanvas(data_.GetName()+'_compare')
    if logy == "log" : 
        c1.SetLogy()
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
    h_res.SetDirectory(0)
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
    mc_.GetYaxis().SetTitle("Events")
    mc_.GetYaxis().SetTitleOffset(1.0)
    mc_.GetXaxis().SetLabelOffset(999)   
    mc_.SetTitle(canvas_title)
    data_.Draw('SAME PE1X0'); 

    # Draw data stat box
    ROOT.gStyle.SetOptStat("i");
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
    c1.SaveAs('%s.png'%event_type)
    return c1