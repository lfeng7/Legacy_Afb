#jetHelper contains helper functions for jet selection
#NICK EMINIZER JOHNS HOPKINS UNIVERSITY JANUARY 2015 nick.eminizer@gmail.com
#This code available on github at https://github.com/eminizer/TTBar_FB_Asym

import ROOT
from math import *
import numpy as np

#Global variables
#constants
TOP_MASS = 172.5
CSVL_WORKING_POINT = 0.244
#JEC values
JES_NOM  = 0.00
JES_UP	 = 0.03
JES_DOWN = -0.03
JER_NOM  = 0.10
JER_UP   = 0.20
JER_DOWN = 0.00
#cut values
MAX_LEPB_MASS = 50. #GeV
MIN_HAD_TOP_MASS = 100. #GeV

#selectJets makes the first element of the jets list the top candidate and the second element the leptonic b candidate
#hadronic top: hardest jet with appropriate jet mass
#leptonic b: hardest jet with low jet mass
def selectJets(jetlist) :
	#start with just the hardest jets
	jetlist.sort(key = lambda x: x.vec.Pt(),reverse=True)
	newjets = []
	if len(jetlist)>1 :
		newjets.append(jetlist[0]); newjets.append(jetlist[1])
	else :
		return jetlist
	#make some mass requirements
	top_cands = []; b_cands = []
	for i in range(0,2) :
		mass = jetlist[i].vec.M()
		if mass > MIN_HAD_TOP_MASS :
			top_cands.append(jetlist[i])
		elif mass < MAX_LEPB_MASS :
			b_cands.append(jetlist[i])
	if len(top_cands)>0 and len(b_cands)>0 :
		newjets[0] = top_cands[0]
		newjets[1] = b_cands[0]
	return newjets

#adjustJEC corrects the jets for scaling the JES or JER up or down
#this is a fix so that I don't need to rerun all the nTuples five times
def adjustJEC(jet,jecuncpos,jecuncneg,corr,ptsmear,etascale,phiscale,matchedJetEta,jec_type=0) :
	returnJet = ROOT.TLorentzVector(jet.X(),jet.Y(),jet.Z(),jet.T())
	newCorr = corr; newPtSmear = ptsmear; newEtaScale = etascale; newPhiScale = phiscale
	if jec_type == 1 :
		newCorr *= (1+sqrt(jecuncpos*jecuncpos+JES_UP*JES_UP))
	elif jec_type == 2 :
		newCorr *= (1-sqrt(jecuncneg*jecuncneg+JES_DOWN*JES_DOWN))
	elif jec_type == 3 :
		if etascale != 0. :
			newEtaScale = max(0.0,1.+JER_UP*(etascale-1.))
		if phiscale != 0. :
			newPhiScale = max(0.0,1.+JER_UP*(phiscale-1.))
		if ptsmear != 0. :
			abseta = abs(matchedJetEta)
			if abseta <= 0.5 :
				newPtSmear = max(0.0,1.+(0.105/0.079)*(ptsmear-1.))
			elif abseta <= 1.1 :
				newPtSmear = max(0.0,1.+(0.127/0.099)*(ptsmear-1.))
			elif abseta <= 1.7 :
				newPtSmear = max(0.0,1.+(0.150/0.121)*(ptsmear-1.))
			elif abseta <= 2.3 :
				newPtSmear = max(0.0,1.+(0.254/0.208)*(ptsmear-1.))
			elif abseta <= 2.8 :
				newPtSmear = max(0.0,1.+(0.316/0.254)*(ptsmear-1.))
			elif abseta <= 3.2 :
				newPtSmear = max(0.0,1.+(0.458/0.395)*(ptsmear-1.))
			elif abseta < 5.0 :
				newPtSmear = max(0.0,1.+(0.247/0.056)*(ptsmear-1.))
	elif jec_type == 4 :
		if etascale != 0. :
			newEtaScale = max(0.0,1.+JER_DOWN*(etascale-1.))
		if phiscale != 0. :
			newPhiScale = max(0.0,1.+JER_DOWN*(phiscale-1.))
		if ptsmear != 0. :
			abseta = abs(matchedJetEta)
			if abseta <= 0.5 :
				newPtSmear = max(0.0,1.+(0.053/0.079)*(ptsmear-1.))
			elif abseta <= 1.1 :
				newPtSmear = max(0.0,1.+(0.071/0.099)*(ptsmear-1.))
			elif abseta <= 1.7 :
				newPtSmear = max(0.0,1.+(0.092/0.121)*(ptsmear-1.))
			elif abseta <= 2.3 :
				newPtSmear = max(0.0,1.+(0.162/0.208)*(ptsmear-1.))
			elif abseta <= 2.8 :
				newPtSmear = max(0.0,1.+(0.192/0.254)*(ptsmear-1.))
			elif abseta <= 3.2 :
				newPtSmear = max(0.0,1.+(0.332/0.395)*(ptsmear-1.))
			elif abseta < 5.0 :
				newPtSmear = max(0.0,1.+(-0.135/0.056)*(ptsmear-1.))
	if phiscale != 0. :
		returnJet.SetPhi(returnJet.Phi()/phiscale)
	if etascale != 0. :
		returnJet.SetPtEtaPhiM(returnJet.Pt(),returnJet.Eta()/etascale,returnJet.Phi(),returnJet.M())
	if corr != 0. and ptsmear != 0. :
		returnJet*=(newCorr*newPtSmear)/(corr*ptsmear)
	returnJet.SetPtEtaPhiM(returnJet.Pt(),returnJet.Eta()*newEtaScale,returnJet.Phi(),returnJet.M())
	returnJet.SetPhi(returnJet.Phi()*newPhiScale)
	return returnJet

#matchUnprunedVec matches a pruned jet to an unpruned jet and returns the index of said unpruned jet in the list
def matchUnprunedVec(pruned_jet_vec,unpruned_fourvecs) :
	if len(unpruned_fourvecs)<1 :
		return -1
	thisVec = ROOT.TLorentzVector(unpruned_fourvecs[0].X(),unpruned_fourvecs[0].Y(),unpruned_fourvecs[0].Z(),unpruned_fourvecs[0].T())
	closestDR = thisVec.DeltaR(pruned_jet_vec)
	matchedJetIndex = 0
	for i in range(1,len(unpruned_fourvecs)) :
		thisVec = ROOT.TLorentzVector(unpruned_fourvecs[i].X(),unpruned_fourvecs[i].Y(),unpruned_fourvecs[i].Z(),unpruned_fourvecs[i].T())
		checkDR = thisVec.DeltaR(pruned_jet_vec)
		if checkDR < closestDR :
			closestDR = checkDR
			matchedJetIndex = i
	return matchedJetIndex

##################################################################################################
##########						b-tagging efficiency stuff 								##########
##################################################################################################

#btagging efficiency stuff
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation53XReReco
#https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods

#getSF takes a jet 4vector, flavour and a CSV working point 
#and returns a tuple of (SF, SF low, SF hi)
def getSF(jetvec,jetflav,CSV_wp) :
	if CSV_wp == CSVL_WORKING_POINT :
		sf_func_b = get_CSVL_SFb
		sf_func_c = get_CSVL_SFc
		sf_func_light = get_CSVL_SFlight
	else :
		print 'WARNING: CSV WORKING POINT NOT RECOGNIZED FOR BTAG EFFICIENCY SCALEFACTOR CALCULATION'
		return (1.0,1.0,1.0)
	#get SF, low, hi from functions for jets depending on flavor
	if abs(jetflav) == 5 : #b jets
		new_SF_tuple = sf_func_b(jetvec.Pt())
	elif abs(jetflav) == 4 : #c jets
		new_SF_tuple = sf_func_c(jetvec.Pt())
	elif abs(jetflav) == 0 or abs(jetflav) == 1 or abs(jetflav) == 2 or abs(jetflav) == 3 or abs(jetflav) == 21 : #light
		new_SF_tuple = sf_func_light(jetvec.Pt(),jetvec.Eta())
	elif abs(jetflav) == 100 : #skip jets without flavor info
		return (1.0,1.0,1.0)
	else :
		print 'WARNING: UNRECOGNIZED JET FLAVOUR IN BTAG EFFICIENCY SF CALCULATION: flavour = '+str(jetflav)+''
	return_SF = new_SF_tuple[0]
	return_SF_low = new_SF_tuple[1]
	return_SF_hi  = new_SF_tuple[2]
	return (return_SF,return_SF_low,return_SF_hi)
	#return (1.0,1.0,1.0) #DEBUG RETURN

#payload here: https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation53XReReco/SFb-pt_WITHttbar_payload_EPS13.txt
btag_eff_pt_bins = [20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800]
def get_CSVL_SFb(pT) : 
	SF = 1.00572*((1.+(0.013676*pT))/(1.+(0.0143279*pT)))
	CSVL_SFb_error = ( [0.033408,0.015446,0.0146992,0.0183964,0.0185363,0.0145547,0.0176743,0.0203609,0.0143342,
					0.0148771,0.0157936,0.0176496,0.0209156,0.0278529,0.0346877,0.0350101] )
	SF_err = 1.0
	if pT < btag_eff_pt_bins[0] :
		SF_err = 2.0*CSVL_SFb_error[0]
	elif pT > btag_eff_pt_bins[len(btag_eff_pt_bins)-1] :
		SF_err = 2.0*CSVL_SFb_error[len(CSVL_SFb_error)-1]
	else:
		for i in range(len(btag_eff_pt_bins)-1) :
			if pT > btag_eff_pt_bins[i] and pT < btag_eff_pt_bins[i+1] :
				SF_err = CSVL_SFb_error[i] 
	return SF,SF-SF_err,SF+SF_err

def get_CSVL_SFc(pT) : 
	SFb_tuple = get_CSVL_SFb(pT)
	return SFb_tuple[0],SFb_tuple[0]-(2.0*(SFb_tuple[0]-SFb_tuple[1])),SFb_tuple[0]+(2.0*(SFb_tuple[2]-SFb_tuple[0]))

#payload here: https://twiki.cern.ch/twiki/pub/CMS/BtagRecommendation53XReReco/SFlightFuncs_EPS2013.C
def get_CSVL_SFlight(pT,eta) : 
	eta = abs(eta)
	SF = 1.0; SF_low = 1.0; SF_hi = 1.0
	if eta < 0.5 :
		SF = ((1.01177+(0.0023066*pT))+(-4.56052e-06*(pT*pT)))+(2.57917e-09*(pT*(pT*pT)))
		SF_low = ((0.977761+(0.00170704*pT))+(-3.2197e-06*(pT*pT)))+(1.78139e-09*(pT*(pT*pT)))
		SF_hi  = ((1.04582+(0.00290226*pT))+(-5.89124e-06*(pT*pT)))+(3.37128e-09*(pT*(pT*pT)))
	elif eta > 0.5 and eta < 1.0 :
		SF = ((0.975966+(0.00196354*pT))+(-3.83768e-06*(pT*pT)))+(2.17466e-09*(pT*(pT*pT)))
		SF_low = ((0.945135+(0.00146006*pT))+(-2.70048e-06*(pT*pT)))+(1.4883e-09*(pT*(pT*pT)))
		SF_hi  = ((1.00683+(0.00246404*pT))+(-4.96729e-06*(pT*pT)))+(2.85697e-09*(pT*(pT*pT)))
	elif eta > 1.0 and eta < 1.5 :
		SF = ((0.93821+(0.00180935*pT))+(-3.86937e-06*(pT*pT)))+(2.43222e-09*(pT*(pT*pT)))
		SF_low = ((0.911657+(0.00142008*pT))+(-2.87569e-06*(pT*pT)))+(1.76619e-09*(pT*(pT*pT)))
		SF_hi  = ((0.964787+(0.00219574*pT))+(-4.85552e-06*(pT*pT)))+(3.09457e-09*(pT*(pT*pT)))
	elif eta > 1.5 and eta < 2.4 :
		SF = ((1.00022+(0.0010998*pT))+(-3.10672e-06*(pT*pT)))+(2.35006e-09*(pT*(pT*pT)))
		SF_low = ((0.970045+(0.000862284*pT))+(-2.31714e-06*(pT*pT)))+(1.68866e-09*(pT*(pT*pT)))
		SF_hi  = ((1.03039+(0.0013358*pT))+(-3.89284e-06*(pT*pT)))+(3.01155e-09*(pT*(pT*pT)))
#	else :
#		print 'WARNING: EVENT WITH JET OUTSIDE OF ETA RANGE! CANNOT GET BTAG EFFICIENCY! |eta| = '+str(eta)+''
	return SF,SF_low,SF_hi

# get maximal variation of PDF
def get_PDF_SF(pdf_w):
	"""
	input: a vector 
	output: a tuple (w_up,w_down)
	"""
	pdf = [ item for item in pdf_w]
	pdf = np.array(pdf)
	pdf_nom = pdf[0]
	pdf = pdf/pdf_nom
	return (pdf.max(),pdf.min())

	