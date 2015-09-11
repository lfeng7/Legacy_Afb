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
def plotting(histlist,event_type='MC',upload = False,testing = 'testing',logy=False):
    if testing == "testing":
        plotdir = 'testing_plots/'
    else : plotdir = 'plots/'
    fout = ROOT.TFile(plotdir+event_type+'_plots.root','recreate')
    for ihist in histlist:
        name = plotdir+event_type+'_'+ihist.GetName()+'.png'
        c1 = ROOT.TCanvas()
        if logy == "setlogy" : c1.SetLogy()
        ihist.Draw()
        c1.SaveAs(name)
        c1.Write()
    if upload == "dump":
        import os
        if testing == "testing" : 
            os.system('source ./dump_testing.sh')
        else : os.system('source ./dump_plots.sh')

# Save histograms to root files
def saving(histlist,event_type='MC',testing = 'testing'):
    if testing == "testing":
        plotdir = 'testing_plots/'
    else : plotdir = 'plots/'
    fout = ROOT.TFile(plotdir+event_type+'_selection_output.root','recreate')
    for ihist in histlist:
        ihist.Write()

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
   
def vector( _list):
    _vector = ROOT.vector("string")()
    for ifile in _list:
        _vector.push_back(ifile)
    return _vector 
