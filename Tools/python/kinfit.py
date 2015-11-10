import ROOT
import math
from array import array

#Global constants
MW = 80.4 
MT = 173.3
# if options.mc == 'yes' :
#   MT = 172.5

MMUON     = 0.105658
MELECTRON = 0.000511
MLEP = MELECTRON
# if options.lep == 'mu' :
#   MLEP = MMUON
# elif options.lep == 'el' :
#   MLEP = MELECTRON

# Global payloads

#Set up file and read histograms used for CSV info in kinematic fit.
f1 = ROOT.TFile('/uscms_data/d3/lfeng7/Payloads/run1/kinfit/dumped_Powheg_TT.root')  # The CSV distribution input files for Powheg signal
# f1 = TFile('/uscms_data/d3/lfeng7/CMSSW_5_3_11/src/Analysis/codes/tagging/tagged_files_9_14/dump_TTJets_SemiLep.root') # The CSV input for Madgraph signals
bdisc = f1.Get('bDisc')
Wsubdisc = f1.Get('WsubDisc')
otherdisc = f1.Get('otherDisc')
bdiscInt = bdisc.Integral()
WsubdiscInt = Wsubdisc.Integral()
otherdiscInt = otherdisc.Integral()

def GetNvPz(lepton,metPt,metPhi):
    met1 = ROOT.TLorentzVector()
    #find the best first guess of the Pz based on the W mass and the lepton we found
    met1.SetPtEtaPhiM(metPt,0.0,metPhi,0.0)
    met2   = met1.Clone()
    pTv    = metPt
    phivec = [math.cos(metPhi),math.sin(metPhi)]
    Elep   = lepton.E()
    plep   = math.sqrt(lepton.Px()*lepton.Px()+lepton.Py()*lepton.Py()+lepton.Pz()*lepton.Pz())
    pZlep  = lepton.Pz()
    pPhi   = lepton.Px()*phivec[0]+lepton.Py()*phivec[1]
    arg0   = MW*MW+plep*plep-Elep*Elep+2.*pTv*pPhi
    arg    = Elep*Elep*(4.*pTv*pTv*(pZlep*pZlep-Elep*Elep)+arg0*arg0) #discriminant in the quadratic equation solution
    if not arg > 0 : #If discriminant is imaginary
        pzv1 = pZlep*arg0/(2.*(Elep*Elep-pZlep*pZlep))
        met1.SetPz(pzv1)
        met1.SetE(math.sqrt(met1.Px()*met1.Px()+met1.Py()*met1.Py()+met1.Pz()*met1.Pz()))
        met2 = met1.Clone()
    else : #have two choices for the neutrino Pz from the quadratic equation
        pzv1 = (pZlep*arg0+math.sqrt(arg))/(2.*(Elep*Elep-pZlep*pZlep))
        met1.SetPz(pzv1)
        met1.SetE(math.sqrt(met1.Px()*met1.Px()+met1.Py()*met1.Py()+met1.Pz()*met1.Pz()))
        pzv2 = (pZlep*arg0-math.sqrt(arg))/(2.*(Elep*Elep-pZlep*pZlep))
        met2.SetPz(pzv2)
        met2.SetE(math.sqrt(met2.Px()*met2.Px()+met2.Py()*met2.Py()+met2.Pz()*met2.Pz()))
    return(met1,met2)

#four vector constrained rescaling function for use in kinematic fit
def pscale(scalefactor, vector) :
    p2 = scalefactor*scalefactor*(vector.Px()*vector.Px()+vector.Py()*vector.Py()+vector.Pz()*vector.Pz())
    m2 = vector.M2()
    newE = math.sqrt(p2+m2)
    return ROOT.TLorentzVector(scalefactor*vector.Px(),scalefactor*vector.Py(),scalefactor*vector.Pz(),newE)

#Minimization function for kinematic fitting
minf = 1000000000.
def fcn(npar, deriv, f, par, flag) :
    global minf
    # Some constants
    #for debug purpose, see if MT is set correctly
    #print 'MT = ',MT

    QW = (MW*MW)/(MT*MT)
    ZW = (2.0*2.0)/(MT*MT)
    ZT = (1.4*1.4)/(MT*MT)
    CDFT  = (math.acos(0.)+math.atan(1./math.sqrt(ZT)))/math.sqrt(ZT)
    CDFW  = 0.5+2.*QW+(1.5*QW*QW-0.5*ZW*QW-1.5)*math.log(((1.-QW)*(1.-QW)+ZW*QW)/(QW*QW+QW*ZW))
    CDFW += ((QW*QW*QW-3.*ZW*QW*QW-3.*QW+2.)/math.sqrt(ZW*QW))*(math.atan((1.-QW)/math.sqrt(ZW*QW))+math.atan(QW/math.sqrt(ZW*QW)))
    SIGMAJ = 0.10
    SIGMAL = 0.03   
    # Set up fit function
    lnL = 0.0
    l = pscale(par[1], lepton)
    bl = pscale(par[2], blep)
    bh = pscale(par[3], bhad)
    whs1 = pscale(par[4],Wsub1)
    whs2 = pscale(par[5],Wsub2)
    wh = whs1+whs2
    newmetx = met.Px()+(1.0-par[1])*lepton.Px()+(1.0-par[2])*blep.Px()+(1.0-par[3])*bhad.Px()+(1.0-par[4])*Wsub1.Px()+(1.0-par[5])*Wsub2.Px()
    newmety = met.Py()+(1.0-par[1])*lepton.Py()+(1.0-par[2])*blep.Py()+(1.0-par[3])*bhad.Py()+(1.0-par[4])*Wsub1.Py()+(1.0-par[5])*Wsub2.Py()
    v = pscale(1.0,met)
    v.SetPx(newmetx)
    v.SetPy(newmety)
    v.SetPz(par[0])
    v.SetE(math.sqrt(v.Px()*v.Px()+v.Py()*v.Py()+v.Pz()*v.Pz()))
    #print ' par0='+str(par[0])+' par1='+str(par[1])+' par2='+str(par[2])+' par3='+str(par[3])+' par4='+str(par[4])+' par5='+str(par[5])+''
    wl = v + l
    tl = wl + bl
    th = wh + bh
    mwl2 = wl.M2()
    mtl2 = tl.M2()
    mwh2 = wh.M2()
    mth2 = th.M2()
    ql = mwl2/(MT*MT)
    xl = mtl2/(MT*MT)
    qh = mwh2/(MT*MT)
    xh = mth2/(MT*MT)
    pdftl = 1./((xl - 1.)*(xl - 1.) + ZT)
    pdfth = 1./((xh - 1.)*(xh - 1.) + ZT)
    pdfwl = (1. - ql)*(1. - ql)*(2. + ql)/((ql-QW)*(ql-QW)+ZW*QW)
    pdfwh = (1. - qh)*(1. - qh)*(2. + qh)/((qh-QW)*(qh-QW)+ZW*QW)
    pdf = pdftl*pdfth*pdfwl*pdfwh/(CDFT*CDFT*CDFW*CDFW)
    if pdf > 0.0 :
        lnL += math.log(pdf)    #need positive f
    else :
        print('WARNING -- pdf is negative!!!')
        pdf = 1.e-50
        lnL += math.log(pdf)
    pdiscblep = bdisc.GetBinContent(bdisc.FindFixBin(blepCSV))/bdiscInt
    pdiscbhad = bdisc.GetBinContent(bdisc.FindFixBin(bhadCSV))/bdiscInt
    if WCSV1 == -1.0 :
        pdiscw1 = 1.0
    else :
        pdiscw1 = Wsubdisc.GetBinContent(Wsubdisc.FindFixBin(WCSV1))/WsubdiscInt
    if WCSV2 == -1.0 :
        pdiscw2 = 1.0
    else :
        pdiscw2 = Wsubdisc.GetBinContent(Wsubdisc.FindFixBin(WCSV2))/WsubdiscInt
    if extraCSV == -1.0 :
        pdiscextra = 1.0
    else :
        pdiscextra = otherdisc.GetBinContent(otherdisc.FindFixBin(extraCSV))/otherdiscInt
    if (pdiscblep*pdiscbhad*pdiscw1*pdiscw2*pdiscextra) <= 0 :
        print('WARNING -- problem with discriminant values!!!')
        print 'Values are: '+str(blepCSV)+' '+str(bhadCSV)+' '+str(WCSV1)+' '+str(WCSV2)+''
        pdf = 1.e-50
        f[0] += math.log(pdf)
    else :
        f[0] = -2.0*lnL + (par[1]-1.)*(par[1]-1.)/(SIGMAL*SIGMAL) + (par[2]-1.)*(par[2]-1.)/(SIGMAJ*SIGMAJ)
        f[0] = f[0] + (par[3]-1.)*(par[3]-1.)/(SIGMAJ*SIGMAJ) + (par[4]-1.)*(par[4]-1.)/(SIGMAJ*SIGMAJ) + (par[5]-1.)*(par[5]-1.)/(SIGMAJ*SIGMAJ)
        f[0] = f[0] - 2.0*math.log(pdiscblep*pdiscbhad*pdiscw1*pdiscw2*pdiscextra)
    if f[0] < minf :
        minf = f[0]

############################################################################################################
#At this point, we have selected all the leptons, met, and jets for the event. The fourvector of the lepton#
#is called 'lepton', the 2 cursory (unfitted) fourvectors of the met are in 'met1/2', and the list of jet  #
#fourvectors is called 'jetCands'. Next we want to reconstruct the t and tbar from the event.              #
############################################################################################################
def DoReco(jetCands,jetCandCSVs,lep_p4,metPt,metPhi,lep_type,mcordata):
    global minf
    global lepton,blep,bhad,Wsub1,Wsub2,met
    global blepCSV,bhadCSV,WCSV1,WCSV2,extraCSV
    global MT,MLEP
    ######################################################
    ##              EVENT RECONSTRUCTION                ##
    ######################################################
    # Set some mass constant
    if lep_type in ['el','electron'] : MLEP = MELECTRON
    if lep_type in ['mu','muon'] : MLEP = MMUON
    if mcordata == 'mc': MT = 172.5
    if mcordata == 'data': MT = 173.3

    # Book leptons and neutrinos
    lepton = lep_p4
    nv_solutions = GetNvPz(lepton,metPt,metPhi)
    met1 = nv_solutions[0]
    met2 = nv_solutions[1]
    #make a list of all combinations of the jet candidates
    combos    = []
    comboCSVs = []
    icombos   = [] #this just checks the combinatorics. Print it if you want.
    extraJetCSVs = []
    if len(jetCands) == 4 :
        for i in range(len(jetCands)) :
            iblep = i
            iothers = []
            for j in range(len(jetCands)) :
                if j!=iblep :
                    iothers.append(j)
            for j in range(len(iothers)) :
                combos.append([])
                combos[len(combos)-1].append(jetCands[iblep])
                comboCSVs.append([])
                comboCSVs[len(comboCSVs)-1].append(jetCandCSVs[iblep])
                icombos.append([])
                icombos[len(icombos)-1].append(iblep)
                extraJetCSVs.append(-1.0)
                for k in range(len(iothers)) :
                    combos[len(combos)-1].append(jetCands[iothers[(k-j+len(iothers))%len(iothers)]])
                    comboCSVs[len(comboCSVs)-1].append(jetCandCSVs[iothers[(k-j+len(iothers))%len(iothers)]])
                    icombos[len(icombos)-1].append(iothers[(k-j+len(iothers))%len(iothers)])
    elif len(jetCands) == 5 :
        for i in range(len(jetCands)) :
            iblep = i
            for j in range(len(jetCands)) :
                if j == iblep :
                    continue
                iextraJet = j
                iothers   = []
                for k in range(len(jetCands)) :
                    if k != iextraJet and k != iblep :
                        iothers.append(k)
                for k in range(len(iothers)) :
                    combos.append([])
                    combos[len(combos)-1].append(jetCands[iblep])
                    comboCSVs.append([])
                    comboCSVs[len(comboCSVs)-1].append(jetCandCSVs[iblep])
                    icombos.append([])
                    icombos[len(icombos)-1].append(str(iblep))
                    extraJetCSVs.append(jetCandCSVs[iextraJet])
                    for n in range(len(iothers)) :
                        combos[len(combos)-1].append(jetCands[iothers[(n-k+len(iothers))%len(iothers)]])
                        comboCSVs[len(comboCSVs)-1].append(jetCandCSVs[iothers[(n-k+len(iothers))%len(iothers)]])
                        icombos[len(icombos)-1].append(str(iothers[(n-k+len(iothers))%len(iothers)]))
    #Combos now includes all possible orderings of the jets (order: blep, bhad, Wsub, Wsub)
    #For each combo, the extra jet goes in comboExtraJet
    #Remove any combos that have neither supposed b tagged as such, or any combo where the CSV value excludes one of the bs from being a b
    nbtags = len([ icsv for icsv in jetCandCSVs if icsv>0.679])
    i = 0
    while True :
        twountagged = comboCSVs[i][0] < 0.679 and comboCSVs[i][1] < 0.679
        oneuntagged = comboCSVs[i][0] < 0.679 or comboCSVs[i][1] < 0.679
        if (nbtags == 1 and twountagged) or (nbtags >= 2 and oneuntagged) or comboCSVs[i][0] < 0 or comboCSVs[i][1] < 0  :
            combos.pop(i)
            comboCSVs.pop(i)
            icombos.pop(i)
            extraJetCSVs.pop(i)
        else :
            i += 1
        if i >= len(combos) :
            break
    i = 0
    if not len(combos) > 0 :
        print 'no valid combos'
        return 'none'
    #Variables for the kinematic fits
    bestParValues1 = []
    bestParValues2 = []
    bestParValues = []
    Chis1 = []
    Chis2 = []
    parErrs = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    parNames = ['pZv','scaleLep','scaleblep','scalebhad','scaleWsub1','scaleWsub2']
    for i in range(len(combos)) :
        bestParValues1.append([0.0,0.0,0.0,0.0,0.0,0.0])
        bestParValues2.append([0.0,0.0,0.0,0.0,0.0,0.0])
    #Now perform kinematic fits for each combination of jetsand both neutrino solutions.
    nFits = 1
    if met1.Pz() != met2.Pz() :
        nFits = 2
    for iFit in range(nFits) : #one for each neutrino solution
        #set which neutrino to use for this iteration
        if iFit == 0 :
            met = met1.Clone()
        elif iFit == 1:
            met = met2.Clone()
        #Stuff common to every fit
        minuit = ROOT.TMinuit(6)
        minuit.SetFCN(fcn)
        ierflag = ROOT.Long(1)
        arglist = array( 'd', 1*[0.] )
        arglist[0] = -1.
        minuit.mnexcm('SET PRINT', arglist, 1,ierflag)
        minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
        arglist[0] = 100000.
        #do fits for every comination in combos
        for i in range(len(combos)) :
            #set the jet particle variables
            blep    = combos[i][0]
            bhad    = combos[i][1]
            Wsub1   = combos[i][2]
            Wsub2   = combos[i][3]
            blepCSV = comboCSVs[i][0]
            bhadCSV = comboCSVs[i][1]
            WCSV1   = comboCSVs[i][2]
            WCSV2   = comboCSVs[i][3]
            extraCSV = extraJetCSVs[i]
            #set the parameters in minuit
            minuit.mnparm(0,parNames[0],ROOT.Double(met.Pz()),1.0,0,0,ierflag)
            for j in range(1,6) :
                minuit.mnparm(j,parNames[j],1.0,1.0,0,0,ierflag)
            #Fix the scale 
            for j in range(0,6) :
#                minuit.FixParameter(j)
                pass
            #minimize
            minuit.mnexcm('MIGRAD', arglist, 1,ierflag)
            if ierflag != 0 :
                print 'PROBLEM IN FIT: ierflag = '+str(ierflag)+''
                continue
            #Set fit Chi of this particular combination
            if iFit == 0 :
                Chis1.append((minf,i))
            elif iFit == 1:
                Chis2.append((minf,i))
            minf = 1000000000.
            #Get the best parameters back from minuit
            for j in range(6) :
                tmp = ROOT.Double(1.0)
                minuit.GetParameter(j,tmp,ROOT.Double(parErrs[j]))
                if iFit == 0 :
                    bestParValues1[i][j] = tmp
                elif iFit == 1:
                    bestParValues2[i][j] = tmp
    #Find the best fit for this event and record it
    Chis1.sort()
    Chis2.sort()
    plot_final_chi = 1.0
    if len(Chis2) > 0 and Chis2[0][0] < Chis1[0][0] :
        j = Chis2[0][1]
        met = met2.Clone()
        blep = combos[j][0]
        bhad = combos[j][1]
        Wsub1 = combos[j][2]
        Wsub2 = combos[j][3]
        for i in range(6) :
            bestParValues.append(bestParValues2[j][i])
        finalChi = Chis2[0][0]
        plot_final_chi = Chis2[0][0]
    else :
        j = Chis1[0][1]
        met = met1.Clone()
        blep = combos[j][0]
        bhad = combos[j][1]
        Wsub1 = combos[j][2]
        Wsub2 = combos[j][3]
        for i in range(6) :
            bestParValues.append(bestParValues1[j][i])
        finalChi = Chis1[0][0]
        plot_final_chi = Chis1[0][0]
    #Rescale the particle fourvectors based on the optimal parameters
    lepton = pscale(bestParValues[1], lepton)
    blep = pscale(bestParValues[2], blep)
    bhad = pscale(bestParValues[3], bhad)
    Wtag = pscale(bestParValues[4],Wsub1)+pscale(bestParValues[5],Wsub2)
    newmetx = met.Px()+ (1.0-bestParValues[1])*lepton.Px()+(1.0-bestParValues[2])*blep.Px()+(1.0-bestParValues[3])*bhad.Px()
    newmetx += (1.0-bestParValues[4])*Wsub1.Px()+(1.0-bestParValues[5])*Wsub2.Px()
    newmety = met.Py()+ (1.0-bestParValues[1])*lepton.Py()+(1.0-bestParValues[2])*blep.Py()+(1.0-bestParValues[3])*bhad.Py()
    newmety += (1.0-bestParValues[4])*Wsub1.Py()+(1.0-bestParValues[5])*Wsub2.Py()
    met.SetPx(newmetx)
    met.SetPy(newmety)
    met.SetPz(bestParValues[0])
    met.SetE(math.sqrt(met.Px()*met.Px()+met.Py()*met.Py()+met.Pz()*met.Pz()))
    # Return top pair, Whad, Wlep, p4 and momentum scales
    tlep_p4 = (lepton+met+blep).Clone()
    thad_p4 = (Wtag+bhad).Clone()
    wlep_p4 = (lepton + met).Clone()
    whad_p4 = Wtag.Clone()

    toreturn = [ plot_final_chi,(tlep_p4,thad_p4,wlep_p4,whad_p4),(bestParValues) ]
    return toreturn

