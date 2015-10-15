# import ROOT in batch mode
import ROOT,os
ROOT.gROOT.SetBatch(True)  # (on lxplus, the X-connection is opened anyways) 
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
 
# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()
 
# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

# load other useful libs
import math
import glob
 
# superimportant shortcuts!!
#DeltaR = ROOT.Math.VectorUtil.DeltaR
#DeltaPhi = ROOT.Math.VectorUtil.DeltaPhi
#DeltaR2 = lambda a, b: DeltaR(a.p4(), b.p4())  # for reco::Candidates
#DeltaR2 = lambda a, b: a.p4().DeltaR(b.p4())
DeltaR2Lorentz = lambda a, b: a.DeltaR(b)
#DeltaPhi2 = lambda a, b: DeltaPhi(a.p4(), b.p4())  # for reco::Candidates

# work-around for a bug in root: currently "+" does not work on a LorenzVector type
def addLVs(a, b):
    """add two Lorenz-vectors. a and b should be of the same type"""
    LV = type(a) # return the type of object a
    return LV(a.x()+b.x(), a.y()+b.y(), a.z()+b.z(), a.t() + b.t())

def DeltaR2(a,b):
    """ Calculate deltaR of two ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > """
    a_p4,b_p4= ROOT.TLorentzVector(),ROOT.TLorentzVector()
    a_p4.SetPtEtaPhiE(a.pt(),a.eta(),a.phi(),a.energy())
    b_p4.SetPtEtaPhiE(b.pt(),b.eta(),b.phi(),b.energy())
    return a_p4.DeltaR(b_p4)
    

def write_out (inputs,outfile):
    with open(outfile, 'w') as f:
        for line in inputs:
            f.write(line + '\n')

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
        c1 = ROOT.TCanvas()
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
    c1 = ROOT.TCanvas()
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
    # Draw two histgrams
    data_.Draw(options_)
    mc_.Draw('same')
    data_.Draw(options_+' same')
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
def saving(histlist,event_type='MC',index = 0000):
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
    fout = ROOT.TFile(savedir+event_type+'_selection_output_'+str(index)+'.root','recreate')
    print 'saving output into file: '+savedir+event_type+'_selection_output_'+str(index)+'.root'
    for ihist in histlist:
        ihist.Write()
    fout.Write()
    # file closure
    fout.Close()

# Save histograms to root files in current directory. Can be used in grid jobs
def gridsaving(histlist,event_type='MC',index = 0000):
    # Saving root files
    fout = ROOT.TFile(event_type+'_selection_output_'+str(index)+'.root','recreate')
    print 'saving output into file: '+event_type+'_selection_output_'+str(index)+'.root'
    for ihist in histlist:
        ihist.Write()
    fout.Write()
    # file closure
    fout.Close()

# Print pdgIds of all daughter particles
def check_daus( gen ):
    daus_pdgid = []
    daus_status = []
    for idau in range(gen.numberOfDaughters()):
        daus_pdgid.append(gen.daughter(idau).pdgId())
        daus_status.append(gen.daughter(idau).status())
    print daus_pdgid
    print daus_status
       
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
    c1 = ROOT.TCanvas()
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
    h_res = data_.Clone()
    h_data = data_
    h_stack = mc_.GetStack().Last() # This is the combined histogram in stack
    for ibin in range(len(h_data.GetNbinsX())):
        databin = h_data.GetBinContent(ibin)
        mcbin = h_stack.GetBinContent(ibin)
        # Calculate residual
        if mcbin != 0 and databin != 0:
            res = databin*1.0/mcbin
            # Calculate error of residual, delta(res) = residual*sqrt(1/data+1/mc)
            res_err = res*math.sqrt(1.0/databin+1.0/mcbin)
        else :
            res = 0 ; res_err = 0
        # Set residual histograms
        h_res.SetBinContent(ibin,res)
        h_res.SetBinError(ibin,res_err)
    # Setup residual histograms
    h_res.SetStats(0)
    h_res.GetYaxis().SetRangeUser(1.1*h_res.GetMaximum(),1.1*h_res.GetMinimum())
    h_res.GetYaxis().SetNdivisions(503) 
    #Build the lines that go at 1 on the residuals plots
    xline = TLine(h_res.GetXaxis().GetXmin(),1.0,h_res.GetXaxis().GetXmax(),1.0); line_.SetLineWidth(2); line_.SetLineStyle(2)
    #plot stacks with data overlaid and residuals. Totally stole from Nick :)
    c1.cd()
    channame = event_type
    # Make and adjust pads
    x_histo_pad=TPad(channame+'_x_histo_pad',channame+'_x_histo_pad',0,0.25,1,1)
    x_resid_pad=TPad(channame+'_x_residuals_pad',channame+'_x_residuals_pad',0,0,1.,0.25)
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
    mc_.Draw('PE1X0 '); xline.Draw()
    c1.Update()    

    # Saving
    c1.SaveAs(name)
    c1.Write()
        
    # dump to webpage
    if upload == "dump":
        os.system('scp -r '+plotdir+'  ~/index.php pha:/home/lfeng/public_html/research/Dump/')
    # file closure
    fout.Close()
