import ROOT
import math

def M3(pt,eta,phi,mass):
    if not len(pt)>= 3: return 0
    m3_vec = ROOT.TLorentzVector()
    for i in range(3):
        ijet = ROOT.TLorentzVector()
        ijet.SetPtEtaPhiM(pt[i],eta[i],phi[i],mass[i])
        m3_vec += ijet
    return m3_vec.M()

def GetSampleColor(sample_type):
    if sample_type == 'singletop': return ROOT.kMagenta
    elif sample_type == 'wjets'    : return ROOT.kGreen-3
    elif sample_type == 'ttbar'    : return ROOT.kRed+1
    elif sample_type == 'zjets'    : return ROOT.kAzure-2
    else :
        return 0

def GetSampleType(color):
    if color == ROOT.kMagenta : return 'singletop'
    elif color == ROOT.kGreen-3 : return 'wjets'  
    elif color == ROOT.kRed+1 : return 'ttbar' 
    elif color == ROOT.kAzure-2 : return 'zjets'
    elif color == 0 : return 'data'
    else : return 'unknown'

# This will return a list contain the number of entries in each bin
def GetBinEntry(hist) :
    nbins = hist.GetSize()
    entry = []
    for i in range(1,nbins) : 
        entry.append(int(hist.GetArray()[i]))
    return entry

################################################################
#               Functions for MC correctons                    #
################################################################

# PU weights
def LoadPUfiles() :
    data_pufile = ROOT.TFile('/uscms_data/d3/lfeng7/Payloads/run1/kinfit/data_pileup_distribution.root')
    mc_pufile = ROOT.TFile('/uscms_data/d3/lfeng7/Payloads/run1/kinfit/dumped_Powheg_TT.root') 
    data_pu_dist = data_pufile.Get('pileup').Clone('npv_data')
    MC_pu_dist   = mc_pufile.Get('pileup').Clone('npv_mc')
    data_pu_dist.Scale(1.0/data_pu_dist.Integral())
    MC_pu_dist.Scale(1.0/MC_pu_dist.Integral())
    pu_dists = [data_pu_dist,MC_pu_dist]
    for item in pu_dists: 
        item.SetDirectory(0)
    print '(info) NPV distribution loaded!'
    return pu_dists 

def GetPUWeights(npvRealTrue,pu_dists):
    # Set up pileup distribution files used in pileup reweighting
    data_pu_dist = pu_dists[0]
    MC_pu_dist = pu_dists[1]
    # Calculate PU corrections based on npvRealTrue
    w_PU = data_pu_dist.GetBinContent(data_pu_dist.FindFixBin(1.0*npvRealTrue))/MC_pu_dist.GetBinContent(MC_pu_dist.FindFixBin(1.0*npvRealTrue))
    return w_PU

# Formula for top pT weights
def GetTopPtWeights(a,b,pt1,pt2):
    pt_w = 1.0
    if pt1<400 and pt2<400 :
        pt_w = math.exp(a+b*(pt1+pt2)/2)
    return pt_w

# Trigger efficiency SF for HLT_Ele27_WP80_v*.
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/KoPFAElectronTagAndProbe
def GetTriggerSFs(pt,eta):
    sf = (1.0,0,0)
    eta = abs(eta)    
    if 30 <= pt <= 40:
        if 0<=eta<=0.8       : sf = (0.987,0.012,0.017)
        if 0.8<=eta<=1.478   : sf = (0.964,0.002,0.001)
        if 1.478<=eta<=2.500 : sf = (1.004,0.006,0.006)
    if 40 <= pt <= 50:
        if 0<=eta<=0.8       : sf = (0.997,0.001,0.001)
        if 0.8<=eta<=1.478   : sf = (0.980,0.001,0.001)
        if 1.478<=eta<=2.500 : sf = (1.033,0.007,0.007)
    if 50 <= pt <= 200:
        if 0<=eta<=0.8       : sf = (0.998,0.002,0.002)
        if 0.8<=eta<=1.478   : sf = (0.988,0.002,0.002)
        if 1.478<=eta<=2.500 : sf = (0.976,0.015,0.012) 
    return sf       

# Electron cut-based ID efficiency SF
# https://twiki.cern.ch/twiki/bin/view/Main/EGammaScaleFactors2012#2012_8_TeV_Jan22_Re_recoed_data
def LoadEleSFs():
    with open('/uscms_data/d3/lfeng7/Payloads/run1/SFs/sFGsfIdTight.txt','r') as fin:
        data = fin.readlines()
    all_SFs = []
    for line in data:
        line = line.split('    ')
        if not len(line) == 8 : continue
        irow = []
        for item in line :
            if item.strip() != '':
                irow.append(float(item.strip()))
        all_SFs.append(irow)
    print '(info) Electron cut based ID efficiency SFs loaded!'
    return all_SFs

def GetEleSFs(pt,eta,SFtable):
    sf = [1.,0,0]
    eta = abs(eta)
    for irow in SFtable:
        pt1 = irow[0]
        pt2 = irow[1]
        eta1 = irow[2]
        eta2 = irow[3]
        sf_ = irow[4]
        err_up = irow[5]
        err_down = irow[6]
        if pt1 <=pt<=pt2 and eta1<=eta<=eta2 :
            sf = [sf_,err_up,err_down]
    return sf

#btagging efficiency

def GetTypeBtagging(sample_name):
    btag_type = 'ttbar'
    for item in ['T_s','T_t','T_tW']:
        if item in sample_name : btag_type = 'singletop'

    for item in ['Tbar_s','Tbar_t','Tbar_tW']:
        if item in sample_name : btag_type = 'singletopbar'

    for item in ['DY1Jets','DY2Jets','DY3Jets','DY4Jets','DYJets']:
        if item in sample_name : btag_type = 'zjets' 

    for item in ['W1Jets','W2Jets','W3Jets','W4Jets','WJets']:
        if item in sample_name : btag_type = 'wjets'

    for item in ['TT']:
        if item in sample_name : btag_type = 'ttbar'   

    for item in ['SingleEl','SingleMu','data','Data']:
        if item in sample_name : btag_type = 'data' 

    return btag_type
     

def LoadBtagEfficiency(sampletype):
    #Set up btag efficiency files
    prepend = '/uscms_data/d3/lfeng7/Payloads/run1/btagging_efficiency/regular/'
    eff_files = []
    eff_files += [('ttbar',prepend+'TT_CT10_TuneZ2star_8TeV-powheg-tauola_CSVM_bTaggingEfficiencyMap.root')]
    eff_files += [('wjets',prepend+'WnJetsToLNu_TuneZ2Star_8TeV-madgraph_CSVM_bTaggingEfficiencyMap.root')]
    eff_files += [('zjets',prepend+'DYnJetsToLL_M-50_TuneZ2Star_8TeV-madgraph_CSVM_bTaggingEfficiencyMap.root')]
    eff_files += [('singletop',prepend+'T_star-channel_TuneZ2star_8TeV-powheg-tauola_CSVM_bTaggingEfficiencyMap.root')]
    eff_files += [('singletopbar',prepend+'Tbar_star-channel_TuneZ2star_8TeV-powheg-tauola_CSVM_bTaggingEfficiencyMap.root')]

    payload_found = 0
    for ifile in eff_files:
        if sampletype == ifile[0] :
            F_eff = ifile[1]
            payload_found = 1

    if payload_found == 0 :
        print ' (debug) No b tagging efficiency payload found for this sample!' 
        return 'None'

    file_tmp = ROOT.TFile(F_eff)
    efficiency_b = file_tmp.Get('efficiency_b').Clone()
    efficiency_c = file_tmp.Get('efficiency_c').Clone()
    efficiency_udsg = file_tmp.Get('efficiency_udsg').Clone()  
    print '(info) Btagging efficiency loaded for type',sampletype
    eff_hists = [efficiency_b,efficiency_c,efficiency_udsg]    
    for item in eff_hists:
        item.SetDirectory(0)
    return eff_hists


def get_btag_eff (pt,eta,jet_flavor,eff_hists):    
     #Debug only
#    print 'eff hist name',efficiency_udsg.GetName()
    # x,y of TH2F of efficiency are pt and eta
    efficiency_b = eff_hists[0]
    efficiency_c = eff_hists[1]
    efficiency_udsg = eff_hists[2]
    if jet_flavor == 5 :
        binx = efficiency_b.GetXaxis().FindBin(pt)
        biny = efficiency_b.GetYaxis().FindBin(eta)
        bins = efficiency_b.GetBin(binx,biny)
        return efficiency_b.GetBinContent(bins)
    elif jet_flavor == 4 :
        binx = efficiency_c.GetXaxis().FindBin(pt)
        biny = efficiency_c.GetYaxis().FindBin(eta)
        bins = efficiency_c.GetBin(binx,biny)
        return efficiency_c.GetBinContent(bins)  
    else :
        binx = efficiency_udsg.GetXaxis().FindBin(pt)
        biny = efficiency_udsg.GetYaxis().FindBin(eta)
        bins = efficiency_udsg.GetBin(binx,biny)
        return efficiency_udsg.GetBinContent(bins)   

def getptbin_for_btag(pt):  
    if(pt<30) : pt_bin = 0;
    elif(pt<40) : pt_bin = 1;
    elif(pt<50) : pt_bin = 2;
    elif(pt<60) : pt_bin = 3;
    elif(pt<70) : pt_bin = 4;
    elif(pt<80) : pt_bin = 5;
    elif(pt<100) : pt_bin = 6;
    elif(pt<120) : pt_bin = 7;
    elif(pt<160) : pt_bin = 8;
    elif(pt<210) : pt_bin = 9;
    elif(pt<260) : pt_bin = 10;
    elif(pt<320) : pt_bin = 11;
    elif(pt<400) : pt_bin = 12;
    elif(pt<500) : pt_bin = 13;
    elif(pt<600) : pt_bin = 14;
    else : pt_bin = 15;
    return pt_bin

def get_eta_bin_jet(eta) : 
    eta = math.fabs(eta);
    # float etaNBins[5] = {0., 0.9, 1.2, 2.1, 2.4};
    if(eta<0.9) : return 0;
    elif(eta<1.2) : return 1;
    elif(eta<2.1) : return 2;
    elif(eta<2.4) : return 3;
    else : return -1;

# //only for CSVM point 
def get_SF_btag(ptJet,etaJet,flavJet):      # checked

    # btagging efficiency constants for CSVM
    SFb_error = [
        0.0415694,
        0.023429,
        0.0261074,
        0.0239251,
        0.0232416,
        0.0197251,
        0.0217319,
        0.0198108,
        0.0193,
        0.0276144,
        0.0205839,
        0.026915,
        0.0312739,
        0.0415054,
        0.0740561,
        0.0598311]

    x = ptJet # the pt of the jet 
    eta = abs(etaJet); # abs(eta) 

    result = []         # the first is SF, second is SF error
    if(eta>2.4) :
        # print 'warning SF_btag_eta>2.4 ??  ' +str(eta)+'' 
        result.append(1)
        result.append(0)
        return result

    if(x<20) : x=20; 
    if(x>800) : x= 800;
    # SF for b or c. flavJet refers to the MC-true flavor of the jet    
    if(abs(flavJet)==5 or abs(flavJet) == 4) : 
        SF  = (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x)) 
        ptbin = getptbin_for_btag( ptJet ) #--> ptJet[indj] refers to the pt of this jet
        SFerr = SFb_error[ptbin]
        if(x>800 or x<20 ) : SFerr *= 2;
        if(abs(flavJet) == 4 ) : SFerr *= 2; 
    # SF for light jets
    else : 
        if(eta<=0.8) :
            SF = ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)))
            SF_up = ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)))
            SF_down = ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)))      
        elif(eta<=1.6) :
            SF = ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)))
            SF_up = ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)))
            SF_down = ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)))
        elif (eta>1.6 and eta<=2.4) :
            SF = ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)))
            SF_up = ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)))
            SF_down = ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)))
        # Calculate error for light jets
        SFerr = max( abs(SF_up-SF),abs(SF_down-SF) )    

    result.append(SF)
    result.append(SFerr)
    return result

# ///this is the functin to call to get the event-by-event b-tag weight
# ///as implemented, it only works for CSVM. 
# Input : given a list of selected jets info such as jets_info = [(pt,eta,flavor,csv)]
def get_weight_btag(jets_info, eff_hists) :
    # return default weight if no eff_hists is given
    if eff_hists == 'None': 
        return [1,0,0]
    # Initiate probability and its uncertainty
    # Probability
    mcTag = 1.;
    mcNoTag = 1.;
    dataTag = 1.;
    dataNoTag = 1.;
    # uncertainty
    err1 = 0; 
    err2 = 0; 
    err3 = 0; 
    err4 = 0; 
    # Loop over all jet candidates
    for ijet in jets_info:
        jet_Pt = ijet[0]
        jet_eta = abs(ijet[1])
        jet_flavor = ijet[2]
        jet_csv = ijet[3]
        # Basic cuts to make sure things go right.
        if(jet_eta>2.4): 
            continue
        if(jet_flavor<=0) :
            continue    #for jets with flavor 0, we ignore. 
        etabin = get_eta_bin_jet(jet_eta);
        # Get b-tagging efficiency using pt,eta and jet_flavor infor
        eff = get_btag_eff(jet_Pt,jet_eta,jet_flavor,eff_hists)
        # Get SF for this jet
        SF_result = get_SF_btag(jet_Pt,jet_eta,jet_flavor)
        jet_SF = SF_result[0]
        jet_SFerr = SF_result[1]
        # Get probability
        istag = jet_csv > 0.679 and math.fabs(jet_eta)<2.4 ;
        if istag :
            mcTag *= eff; 
            dataTag *= eff*jet_SF; 
            # debug
            if eff == 0 or eff*jet_SF == 0:
                print 'jet pt,eta,flavor,csv,eff,SF: '+str(jet_Pt)+','+str(jet_eta)+','+str(jet_flavor)+','+str(jet_csv)+','+str(eff)+','+str(jet_SF)
            # Get error
            if(jet_flavor==5 or jet_flavor ==4) : 
                err1 += jet_SFerr/jet_SF; # correlated for b/c
            else : 
                err3 += jet_SFerr/jet_SF; # correlated for light                                
        else :
            mcNoTag *= (1- eff); 
            dataNoTag *= (1- eff*jet_SF); 
            #debug
            if (1-eff) == 0 or (1-eff*jet_SF)==0 :
                print 'jet pt,eta,flavor,csv,eff,SF: '+str(jet_Pt)+','+str(jet_eta)+','+str(jet_flavor)+','+str(jet_csv)+','+str(eff)+','+str(jet_SF)
            # Get error
            if(jet_flavor==5 or jet_flavor ==4 ) : 
                err2 += (-eff*jet_SFerr)/(1-eff*jet_SF); # correlated for b/c
            else : 
                err4 += (-eff*jet_SFerr)/(1-eff*jet_SF);  # correlated for light

    # Check if any of the probability is zero. If so, then set weight to be 1, aka do nothing to this event
    if dataTag*dataNoTag == 0  or mcTag*mcNoTag == 0:
        print 'one of the probability is zero!'
        wtbtag = 1
        wtbtagErr = 0
    else :
        # Get event weight for this event
        wtbtag = (dataNoTag * dataTag ) / ( mcNoTag * mcTag ); 
        wtbtagErr = math.sqrt( pow(err1+err2,2) + pow( err3 + err4,2)) * wtbtag; # un-correlated for b/c and light
    # # debug
    # print 'wtbtag = '+str(wtbtag)+'+/-'+str(wtbtagErr)
    # Set return
    to_return = []
    to_return.append(wtbtag)
    to_return.append(wtbtagErr)
    return to_return

