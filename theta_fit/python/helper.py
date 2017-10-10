import ROOT
import sys
import numpy
from array import array
import template
"""
common helper functions
"""
XBINS = template.XBINS
YBINS = template.YBINS
ZBINS = template.ZBINS

ROOT.TH1.SetDefaultSumw2(True)
print '(info) ROOT.TH1.SetDefaultSumw2(True)'

def GetValidWeight(ttree,weights):
    """
    input: ttree and str ( "w1*w2*w3" )
    output: string
    """
    weights = weights.split('*')
    new_ws = [ w for w in weights if ttree.FindBranch(w)]
    return '*'.join(new_ws)
    

def GetListBanches(ttree):
    """
    input: a ttree
    output: List[br_names]
    """
    brs = ttree.GetListOfBranches()
    return  [item.GetName() for item in brs]

def Check_if_nominal_weight(w_name,tag):
    """
    input: a string
    output: bool, if it is a norminal weight
    """
    if ('w_' in tag and  w_name.startswith('w_')) or ('w_' not in tag and 'reweight' in w_name):
        if w_name.endswith('_down') or w_name.endswith('_up') or w_name.endswith('_low') or w_name.endswith('_hi'):
            return False
        else:
            return True
    else:
        return False 

def GetWeightNames(ttree):
    """
    input: a ttree
    output: a string contains all weights for quick TTree.Draw
    """
    blacklist = ['top_pT_reweight','GJR_reweight','CT10_reweight','cteq_reweight','w_all_corr']
    all_brs = GetListBanches(ttree)
    if [item for item in all_brs if 'reweight' in item]!=[] : tag = 'reweight' 
    else: tag = 'w_'
    w_brs = [item for item in all_brs if Check_if_nominal_weight(item,tag) and item not in blacklist]
    return '*'.join(w_brs)

def GetTTreeName(tfile):
    # Find the name of the ttree
    keys = tfile.GetListOfKeys()
    for ikey in keys:
        if ikey.GetClassName() == 'TTree' : treename = ikey.GetName()
    #print 'Getting ttree',treename
    return treename

def multiply(_list):
    """
    input: list
    output: multiply all elements in the list
    """
    return reduce(lambda x,y:x*y,_list)

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

# def GetQualityPlots(hist,verbose=False):
#     """
#     input: a TH1D 
#     output: a new TH1D contains the distribution of BinContents for spotting empty bins 
#     """
#     hname = hist.GetName()
#     N_neg_bin = 0
#     hist_quality = ROOT.TH1D('quality_%s'%hname,'quality_%s'%hname,61,-1,60)
#     hist_quality.SetXTitle("evts per bin")
#     hist_quality.SetYTitle("num of bins")
#     for k in range(hist.GetSize()):
#         if not hist.IsBinUnderflow(k) and not hist.IsBinOverflow(k) :
#             binCounts = hist.GetBinContent(k)
#             if binCounts<0:
#                 N_neg_bin += 1
#             hist_quality.Fill(hist.GetBinContent(k))
#     hist_quality.SetDirectory(0)
#     if verbose:
#         print '(verbose) %s has %i negative bins.'%(hname,N_neg_bin)
#     return hist_quality

def GetQualityPlots_data(hist,verbose=True):
    """
    input: a TH1D 
    output: a new TH1D contains the distribution of fractional error for each bin
    Asumming unweighted histogram
    sigma_yi/yi = 1/sqrt(nentries_i)
    """
    hname = hist.GetName()
    N_abnormal = 0
    hist_quality = ROOT.TH1D('err_%s'%hname,'err_%s'%hname,40,-0.1,2)
    hist_quality.SetXTitle("sigma_yi/yi")
    hist_quality.SetYTitle("num of bins")

    for k in range(hist.GetSize()):
        if not hist.IsBinUnderflow(k) and not hist.IsBinOverflow(k) :
            yi = hist.GetBinContent(k)
            if yi>0:
                frac_err = 1/numpy.sqrt(yi)
            else:
                    frac_err = 1.99
            if frac_err>1 and frac_err!=1.99:
                N_abnormal += 1
            hist_quality.Fill(frac_err)
    hist_quality.SetDirectory(0)
    if verbose:
        print '(verbose) %s has %i bins with sigma_i/yi>1.'%(hname,N_abnormal)
    return hist_quality

def GetQualityPlots_MC(hist,verbose=True):
    """
    input: a TH1D 
    output: a new TH1D contains the distribution of fractional error for each bin
    Assumming sum of weights are conserved
    sigma_yi^2 = Sum_over(w_i^2) i=1,...,Number_of_entries in this bin
    yi = Sum_over(w_i), which is number of events in the scaled hist
    so sigma_yi/yi = sqrt(sumw2_i)/yi
    """
    hname = hist.GetName()
    N_abnormal = 0
    hist_quality = ROOT.TH1D('err_%s'%hname,'err_%s'%hname,40,-0.1,2)
    hist_quality.SetXTitle("sigma_yi/yi")
    hist_quality.SetYTitle("num of bins")
    sumofw2 = hist.GetSumw2()

    if sumofw2.GetSize()==0:
        print '(error) No Sumw2 stored for %s. Will exit.'%hist.GetName()
        sys.exit(1)

    for k in range(hist.GetSize()):
        if not hist.IsBinUnderflow(k) and not hist.IsBinOverflow(k) :
            yi = hist.GetBinContent(k)
            sigma_i = sumofw2.At(k)
            if yi>0:
                frac_err = numpy.sqrt(sigma_i)/yi
            else:
                frac_err = 1.99
            if frac_err>1 and frac_err!=1.99:
                N_abnormal += 1
            hist_quality.Fill(frac_err)
    hist_quality.SetDirectory(0)
    if verbose:
        print '(verbose) %s has %i bins with sigma_i/yi>1.'%(hname,N_abnormal)
    return hist_quality

def getColors(name):
    colors = {'qq':ROOT.kRed+1,'gg':ROOT.kRed-7,'ttbar':ROOT.kRed-7,'signal':ROOT.kRed-7,'tt_bkg':ROOT.kBlue+3,'other_bkg':ROOT.kMagenta,'singleT':ROOT.kMagenta,'WJets':ROOT.kGreen-3,'qcd':ROOT.kYellow-3,'zjets':ROOT.kAzure-2 }
    icolor = colors.get(name,0)
    return icolor
    
# This is specifically for comparing the stacked MC plots with data adding residule plots too
import math
def comparison_plot_v1(mc_,data_,legend,event_type='plots',bin_type = 'fixed',logy=False,draw_option = 'hist',outputdir='',lep_type = '#mu + Jets'):
    # reset some tdrStyle
    ROOT.tdrStyle.SetTitleY(0.9);
    ROOT.tdrStyle.SetStatY(0.97);

    # Def axis title
    xaxis_name = data_.GetXaxis().GetTitle()
#    canvas_title = data_.GetTitle()
    canvas_title = 'CMS Private Work, 19.7 fb^{-1} at #sqrt{s} = 8 TeV'
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
    if bin_type == 'fixed':
        # print 'Making fixed bin res plot'
        h_res = ROOT.TH1D(event_type+'_residuals',';; Data/MC',h_data.GetNbinsX(),h_data.GetXaxis().GetXmin(),h_data.GetXaxis().GetXmax())
        h_ref = ROOT.TH1D(event_type+'_reference',';; Data/MC',h_data.GetNbinsX(),h_data.GetXaxis().GetXmin(),h_data.GetXaxis().GetXmax())

    else:
        data_hname = data_.GetName()
        # print 'Making variable bins res plot for %s'%data_hname
        if '_x' in data_hname:
            bins = XBINS
        elif '_y' in data_hname:
            bins = YBINS
        elif '_z' in data_hname:
            bins = ZBINS
        else:
            print 'cannot determin which projection!'
            sys.exit(1)
        h_res = ROOT.TH1D(event_type+'_residuals',';; Data/MC',len(bins)-1,bins)
        h_ref = ROOT.TH1D(event_type+'_reference',';; Data/MC',len(bins)-1,bins)

    h_res.SetDirectory(0)
    h_ref.SetDirectory(0)

    # h_res.GetXaxis().SetName(h_data.GetXaxis().GetName())
    # h_res.SetXTitle(h_data.GetXaxis().GetName())

    maxxdeviations = 0.0
    for ibin in range(1,h_data.GetNbinsX()+1):
        databin = max(0,h_data.GetBinContent(ibin))
        mcbin = max(0,h_stack.GetBinContent(ibin))
        # Calculate residual
        if mcbin != 0 and databin != 0:
            res = databin*1.0/mcbin
            # Calculate error of residual, delta(res) = residual*sqrt(1/data+1/mc)
            try:
                res_err = res*math.sqrt(1.0/databin+1.0/mcbin)
            except ValueError:
                res_err = 0
                print 'res %f, databin %f, mcbin %f'%(res,databin,mcbin)
	    ref_err = math.sqrt(1.0/mcbin)
            # Find maximum residual
            maxxdeviations = max(maxxdeviations,max(abs(res+res_err-1.0),abs(res-res_err-1.0)))
        else :
            res = 0 ; res_err = 0 ; ref_err = 0;
        # Set residual histograms
        h_res.SetBinContent(ibin,res)
        h_res.SetBinError(ibin,res_err)
	h_ref.SetBinContent(ibin,1.0)
	h_ref.SetBinError(ibin,ref_err)
    # print 'maxxdeviations',maxxdeviations
    # Setup residual histograms
    h_res.SetStats(0)
    h_ref.SetStats(0)
    h_res.GetXaxis().SetLabelSize((0.05*0.72)/0.28); h_res.GetXaxis().SetTitleOffset(0.8)
    h_res.GetYaxis().SetLabelSize((0.05*0.72)/0.28); h_res.GetYaxis().SetTitleOffset(0.4)
    h_res.GetXaxis().SetTitleSize((0.72/0.28)*h_res.GetXaxis().GetTitleSize())
    h_res.GetYaxis().SetTitleSize((0.72/0.28)*h_res.GetYaxis().GetTitleSize())
    maxx = 1.0+1.1*maxxdeviations
    minx = 1.0-1.1*maxxdeviations
    h_res.GetYaxis().SetRangeUser(minx,maxx)
    h_res.GetYaxis().SetNdivisions(503)
    h_res.GetYaxis().SetDecimals(1)
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
    x_histo_pad.SetGridx(); x_histo_pad.SetGridy();

    if logy == 'log': x_histo_pad.SetLogy()
    x_resid_pad=ROOT.TPad(channame+'_x_residuals_pad',channame+'_x_residuals_pad',0,0,1.,0.25)
    x_histo_pad.SetCanvas(c1); x_resid_pad.SetCanvas(c1)
    x_histo_pad.SetLeftMargin(0.16); x_histo_pad.SetRightMargin(0.05) 
    x_histo_pad.SetTopMargin(0.11);  x_histo_pad.SetBottomMargin(0.02)
    x_histo_pad.SetBorderMode(0)
    x_resid_pad.SetLeftMargin(0.16); x_resid_pad.SetRightMargin(0.05)
    x_resid_pad.SetTopMargin(0.0);   x_resid_pad.SetBottomMargin(0.3)
    x_resid_pad.SetBorderMode(0)
    x_resid_pad.SetGridx(); x_resid_pad.SetGridy();
    x_resid_pad.Draw(); x_histo_pad.Draw()
    x_histo_pad.cd(); 
    mc_.Draw(draw_option);
    mc_.GetYaxis().SetTitle("Events")
    mc_.GetYaxis().SetTitleOffset(1.0)
    mc_.GetXaxis().SetTitle("")
    mc_.GetXaxis().SetLabelOffset(999)   
    mc_.GetXaxis().SetLabelSize(0)
    mc_.SetTitle(canvas_title)
    # Draw mc stack error 
    final_hist = mc_.GetStack().Last().Clone()
    final_hist.SetName('%s_err'%event_type)
    final_hist.SetFillColor(ROOT.kBlue)
    final_hist.SetMarkerColor(final_hist.GetFillColor())
    final_hist.SetFillStyle(3002)
    #final_hist.Draw('SAMEs E2')
    # Draw data points 
    data_.Draw('SAME PE1X0'); 

    # Draw data stat box
    ROOT.gStyle.SetOptStat("i");
    ROOT.gStyle.SetStatFormat("5.5g");

    data_.SetStats(1)
    data_.Draw('SAMEs PE1X0'); 
    x_histo_pad.Update()
    statbox1 = data_.FindObject("stats")
    statbox1.SetLineColor(0)
    statbox1.Draw('sames')

    legend.SetShadowColor(0)
    legend.SetLineColor(0)
    legend.Draw()
    legend.SetFillColor(0)   

    x_histo_pad.cd();
    pt = ROOT.TLatex(.2,.84,lep_type);  

    pt.SetNDC(ROOT.kTRUE);
    pt.Draw();

    latex2 = ROOT.TLatex()
    # latex2.SetNDC()
    # latex2.SetTextSize(0.5*c.GetTopMargin())
    # latex2.SetTextFont(42)
    # latex2.SetTextAlign(31) # align right
    # latex2.DrawLatex(0.87, 0.95,"12.9 fb^{-1} (13 TeV)")
    latex2.SetTextFont(22)
    latex2.SetTextSize(0.05587301);
    latex2.SetLineWidth(2);
#    latex2.DrawLatex(-0.8250709,2726.119, "e+jets")

    x_resid_pad.cd(); 
    h_res.Draw('PE1X0 '); xline.Draw()
    h_res.GetXaxis().SetTitle(xaxis_name)
    h_res.SetMinimum(0.7);h_res.SetMaximum(1.3);
    h_ref.SetFillColor(ROOT.kBlue)
    h_ref.SetMarkerColor(h_ref.GetFillColor())
    h_ref.SetFillStyle(3002)
    h_ref.Draw('SAMEs E2')

    print 'xaxis_name',xaxis_name

    obj_title = c1.FindObject("title")
    obj_title.SetShadowColor(0)
    obj_title.SetLineColor(0)
    c1.Update()    
    c1.SaveAs('%s.pdf'%event_type)
    c1.SaveAs('%s_source.root'%(event_type))
    return c1,final_hist
