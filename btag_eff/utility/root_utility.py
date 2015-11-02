import ROOT
import os

def normalized_compare_plots(histlist):
    """Make a compared plots with different color and renormalized"""
    if not len(histlist) > 0 : 
        print 'No histgrams in the input'
        #continue
    i = 1 
    # Make a new canvas
    #global  tmp_canvas = make_canvas()
    for hist in histlist:
        hist.SetLineColor(i)
        if i==1 : hist.Draw() 
        else : hist.Draw("same")
        i+=1
     
# Make plots and upload to webpage
def plotting(histlist,event_type='MC',upload = False,logy=False,legend = None,options_ = '',createmode='update'):
    prefix = './plots/'
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

    fout = ROOT.TFile(plotdir+event_type+'_plots.root',createmode)

    # plotting
    for ihist in histlist:
        c1 = ROOT.TCanvas(ihist.GetName())
        if logy == "log" : 
            c1.SetLogy()
            name = plotdir+event_type+'_'+ihist.GetName()+'_log.png'
        else :
            name = plotdir+event_type+'_'+ihist.GetName()+'.png'
        ihist.Draw(options_)
        if legend : legend.Draw()
        c1.SaveAs(name)
        c1.Write()

    fout.Write()
    # dump to webpage
    if upload == "dump":
        os.system('scp -r '+plotdir+'  ~/index.php pha:/home/lfeng/public_html/research/Dump/')
    # file closure
    fout.Close()

# This is specifically for comparing the stacked MC plots with data
def comparison_plot(mc_,data_,legend,event_type='MC',upload = False,logy=False,options_ = 'elp',createmode='update'):
    prefix = './plots/'
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

    fout = ROOT.TFile(plotdir+event_type+'_plots.root',createmode)

    # plotting
    if logy == "log" : 
        c1 = ROOT.TCanvas(data_.GetName()+'_log')
        c1.SetLogy()
        name = plotdir+event_type+'_'+mc_.GetName()+'_compare_log.png'
    else :
        c1 = ROOT.TCanvas(data_.GetName())
        name = plotdir+event_type+'_'+mc_.GetName()+'_compare.png'
    # Find the max of both histograms
    max_mc = mc_.GetMaximum()
    max_data = data_.GetMaximum()
    max_ = max(max_mc,max_data)*1.1
    data_.SetMaximum(max_)
    data_.SetTitle("")
    # Draw two histgrams
    data_.Draw('PE')
    mc_.Draw('same')
    data_.Draw('SAME PE')
    legend.Draw()

    # Saving
    c1.SaveAs(name)
    c1.Write()
        
    # dump to webpage
    if upload == "dump":
        os.system('scp -r '+plotdir+'  ~/index.php pha:/home/lfeng/public_html/research/Dump/')
    # file closure
    fout.Close()


# Save histograms to root files. Use for interactive run only
def saving(histlist,event_type='MC',index = 0000,createmode='recreate'):
    prefix = './output_rootfiles/'
    # check if plotting dir is made. If not , make it now
    if not os.path.exists(prefix):
        os.mkdir(prefix)
        print 'Making '+prefix
    # Set the dir to put all output rootfiles
    savedir = prefix+event_type+'/'
    if not os.path.exists(savedir) : 
        os.mkdir(savedir)
        print 'Creating new dir '+savedir
    # Saving root files
    outname = savedir+event_type+'_selection_output_'+str(index)+'.root'  
    fout = ROOT.TFile(outname,createmode)
    print 'saving output into file: ',outname
    for ihist in histlist:
        ihist.SetDirectory(fout)
        ihist.Write()
    #fout.Write()
    # file closure
    fout.Close()
    # return output names  
    return outname

# A more general function that saves a list of objects into a single root file
def savetoroot(objects,outputdir='histograms',event_type='test',fname='',createmode='recreate'):
    prefix = outputdir
    # check if plotting dir is made. If not , make it now
    if not os.path.exists(prefix):
        os.mkdir(prefix)
        print 'Making '+prefix
    # Set the dir to put all output rootfiles
    savedir = prefix+'/'+event_type+'/'
    if not os.path.exists(savedir) : 
        os.mkdir(savedir)
        print 'Creating new dir '+savedir
    # Saving root files
    fout_name = savedir+fname+'.root'
    fout = ROOT.TFile(fout_name,createmode)
    print 'saving output into file: '+fout_name
    fout.cd()
    for iobj in objects:
        iobj.Write()
    #fout.Write()
    # file closure
    fout.Close()    
    return fout_name 

# Save histograms to root files in current directory. Can be used in grid jobs
def gridsaving(histlist,event_type='MC',index = 0000,createmode='recreate'):
    # Saving root files
    outname = event_type+'_selection_output_'+str(index)+'.root'
    fout = ROOT.TFile(outname,createmode)
    print 'saving output into file: ',outname
    for ihist in histlist:
        ihist.SetDirectory(fout)
        ihist.Write()
    fout.Write()
    # file closure
    fout.Close()
    return outname

# Normalize histogram according to first bin
def norm( hist):
    norm = hist.GetMaximum()
    h1 = hist.Clone()
    h1.SetName(h1.GetName()+'_norm')
    h1.Scale(1.0/norm)
    return h1

# Normalize a stack
def normstack(stack):
    norm = stack.GetMaximum()
    st_ = ROOT.THStack(stack.GetName()+'_norm',stack.GetName())
    hists = stack.GetHists()
    for ihist in hists:
        ihist.Scale(1.0/norm)
        st_.Add(ihist)
    return st_

# Rebin a stack
def RebinStack(stack,nbin): # means merge nbin into one
    st_ = ROOT.THStack(stack.GetName()+'_rebinned',stack.GetName())
    hists = stack.GetHists()
    for ihist in hists:
        ihist.Rebin(nbin)
        st_.Add(ihist)
    return st_

   
def vector( _list):
    _vector = ROOT.vector("string")()
    for ifile in _list:
        _vector.push_back(ifile)
    return _vector 

def loadroot(file_):
    f_ = ROOT.TFile(file_)
    return(f_)


# This is specifically for comparing the stacked MC plots with data adding residule plots too
def comparison_plot_v1(mc_,data_,legend,event_type='MC',upload = False,logy=False,options_ = '',createmode='update'):
    prefix = './plots/'
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

    fout = ROOT.TFile(plotdir+event_type+'_plots.root',createmode)

    # plotting
    c1 = ROOT.TCanvas(data_.GetName())
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
    # mc_.SetMaximum(max_)
    # Calculating the residual of each bin
    h_data = data_
    h_stack = mc_.GetStack().Last() # This is the combined histogram in stack
    # Make residual histogram
    h_res = ROOT.TH1D(event_type+'_residuals',';; Data/MC',h_data.GetNbinsX(),h_data.GetXaxis().GetXmin(),h_data.GetXaxis().GetXmax())
    h_res.GetXaxis().SetName(h_data.GetXaxis().GetName())

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
    h_res.GetXaxis().SetTitle(data_.GetXaxis().GetTitle())
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
    mc_.Draw(); data_.Draw('SAME PE1X0'); mc_.GetXaxis().SetLabelOffset(999)
    legend.Draw()
    x_resid_pad.cd(); 
    h_res.Draw('PE1X0 '); xline.Draw()
    c1.Update()    

    # Saving
    c1.SaveAs(name)
    c1.Write()
        
    # dump to webpage
    if upload == "dump":
        os.system('scp -r '+plotdir+'  ~/index.php pha:/home/lfeng/public_html/research/Dump/')
    # file closure
    fout.Close()
    