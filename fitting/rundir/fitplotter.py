# Make data/MC comparison plots of cos_theta, mtt, xf, charge ratio for template fit results
# Usage: python ../../fitplotter.py --dir fix_xi_delta_mtt_300
# Output: will makedir a plots in --dir, with cos_theta.png etc
# 12-10-2015

import ROOT
from optparse import OptionParser
import os
import glob
import math
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
# Job steering

# Input inputFiles to use. This is in "glob" format, so you can use wildcards.
# If you get a "cannot find file" type of error, be sure to use "\*" instead
# of "*" to make sure you don't confuse the shell. 

parser = OptionParser()

parser.add_option('--name', metavar='F', type='string', action='store',
              default = "blank",
                  dest='name',
                  help='')


parser.add_option('--dir', metavar='F', type='string', action='store',
              default = "",
                  dest='dir',
                  help='')

(options, args) = parser.parse_args()

argv = []


canvas_title = 'CMS Private Work, 19.7 fb^{-1} at #sqrt{s} = 8 TeV'


def main():
    global rundir

    rundir = options.dir

    # Get inputfile   
    fdata = ROOT.TFile(rundir+'final_stack.root')
    x_stack = fdata.Get('x_stack')
    y_stack = fdata.Get('y_stack')
    z_stack = fdata.Get('z_stack')
    event_numbers_stack = fdata.Get('event_numbers_stack')
    data_x = fdata.Get('data_x')
    data_y = fdata.Get('data_y')
    data_z = fdata.Get('data_z')
    event_numbers_data = fdata.Get('event_numbers_data') 
    leg = fdata.Get('legend')
    all_plots = [(x_stack,data_x,'cs'),(y_stack,data_y,'fx'),(z_stack,data_z,'mtt')]

    leg.SetX1NDC(0.6187291);leg.SetY1NDC(0.5104762);leg.SetX2NDC(0.9214047);leg.SetY2NDC(0.855873);  
    leg.SetFillColor(0)   

    # Make data/MC comparison plot
    for item in all_plots:
        mc_stack=item[0]
        data_hist = item[1]
        hname = item[2]
        draw_option = ''
        comparison_plot_v1(mc_stack,data_hist,leg,hname,draw_option)
    # plot charge ratio
    c1 = ROOT.TCanvas()
    event_numbers_stack.Draw('bar1')
    event_numbers_data.Draw('same PE1X0') 
    # leg.Draw('same')
    
    c1.SaveAs(rundir+'plots/charge_compare.png')       

def GetSampleColor(itype):
    if itype=='gg' : return 38
    if itype=='qq' or itype=='signal' : return 46
    if itype=='bck' : return 41
    if itype=='WJets': return ROOT.kGreen-3
    if itype=='tt_bkg': return ROOT.kMagenta
    return 0

# This is specifically for comparing the stacked MC plots with data adding residule plots too

def comparison_plot_v1(mc_,data_,legend,event_type='plots',draw_option = 'h',logy=False):
    global fout,var
    prefix = rundir+'plots'
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
    h_res.GetXaxis().SetTitle(mc_.GetXaxis().GetTitle())
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
    mc_.Draw(draw_option);
    mc_.GetYaxis().SetTitle("events")
    mc_.GetYaxis().SetTitleOffset(1.0)
    mc_.GetXaxis().SetLabelOffset(999)   
    data_.Draw('SAME PE1X0'); 
    mc_.SetTitle(canvas_title)
    # Draw data stat box
    ROOT.gStyle.SetOptStat("e");
    data_.SetStats(1)
    data_.Draw('SAMEs PE1X0'); 
    x_histo_pad.Update()
    statbox1 = data_.FindObject("stats")
    statbox1.SetLineColor(0)
    statbox1.Draw('sames')

    leg = legend.Clone()
    # leg.SetX1NDC(0.7550143);leg.SetY1NDC(0.4991304);leg.SetX2NDC(0.9054441);leg.SetY2NDC(0.8469565);            
    leg.SetShadowColor(0)
    leg.SetLineColor(0)
    leg.Draw()

    x_resid_pad.cd(); 
    h_res.Draw('PE1X0 '); xline.Draw()
    # h_res.GetXaxis().SetTitle(xaxis_name)

    obj_title = c1.FindObject("title")
    obj_title.SetShadowColor(0)
    obj_title.SetLineColor(0)    
    c1.Update()    
    # Saving
    c1.SaveAs(name)
    c1.Write()
    # fout.Close()
    return (c1)


# This is specifically for comparing the stacked MC plots with data
def comparison_plot(mc_,data_,legend,event_type='MC',draw_option = 'h',logy=False):
    prefix = rundir+'plots'
    # check if plotting dir is made. If not , make it now
    if not os.path.exists(prefix):
        os.mkdir(prefix)
        print 'Making '+prefix
    # Set the dir to put all plots
    plotdir = prefix+event_type+'/'
    if not os.path.exists(plotdir):
        os.mkdir(plotdir)
        os.system('cp ~/index.php '+plotdir)
        print 'Creating new dir '+plotdir

    fout = ROOT.TFile(plotdir+event_type+'_plots.root','recreate')

    # plotting
    if logy == "log" : 
        c1 = ROOT.TCanvas(data_.GetName()+'_log')
        c1.SetLogy()
        name = plotdir+event_type+'_compare_log.png'
    else :
        c1 = ROOT.TCanvas(data_.GetName())
        name = plotdir+event_type+'_compare.png'

    # Draw two histgrams
    mc_.Draw('draw_option')
    data_.Draw('SAME PE1X0')
    leg.Draw('same')

    # Saving
    c1.SaveAs(name)
    c1.Write()
    # file closure
    fout.Close()



main()

