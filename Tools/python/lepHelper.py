import pickle
import math

# for all lepton related helper functions
prefix = '/uscms_data/d3/lfeng7/Payloads/run1/SFs'
class muon_helper(object):
	"""docstring for muon_helper
       Muon helper object
       include all muon corrections here
	"""
	def __init__(self):
		super(muon_helper, self).__init__()
		self.reset_weights()
		self.load_payloads()


	def reset_weights(self):
		self.weight_tracking = 1.0
		self.weight_tracking_low = 1.0
		self.weight_tracking_hi = 1.0

		self.weight_lep_ID = 1.0
		self.weight_lep_ID_low = 1.0
		self.weight_lep_ID_hi = 1.0

		self.weight_lep_iso = 1.0
		self.weight_lep_iso_low = 1.0
		self.weight_lep_iso_hi = 1.0

		self.weight_trig_eff = 1.0
		self.weight_trig_eff_low = 1.0
		self.weight_trig_eff_hi = 1.0	

	def load_payloads(self):
		print '(info) Muon SF payloads loaded.'
		#tracking efficiency constants and errors
		self.tracking_eff_consts = [[-2.4,-2.1,0.9869,0.07],
							   [-2.1,-1.6,0.9948,0.02],
							   [-1.6,-1.2,0.9967,0.02],
							   [-1.2,-0.9,0.9974,0.02],
							   [-0.9,-0.6,0.9980,0.01],
							   [-0.6,-0.3,0.9980,0.01],
							   [-0.3,-0.2,0.9972,0.02],
							   [-0.2, 0.2,0.9963,0.01],
							   [ 0.2, 0.3,0.9978,0.02],
							   [ 0.3, 0.6,0.9977,0.01],
							   [ 0.6, 0.9,0.9976,0.01],
							   [ 0.9, 1.2,0.9968,0.02],
							   [ 1.2, 1.6,0.9959,0.03],
							   [ 1.6, 2.1,0.9970,0.02],
							   [ 2.1, 2.4,0.9836,0.08]]
		#unpickle the lepton ID, isolation, and trigger files
		#grid filenames
		muon_id_filename  = '%s/MuonEfficiencies_Run2012ReReco_53X.pkl'%prefix
		muon_iso_filename = '%s/MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl'%prefix
		muon_trigger_filename = '%s/SingleMuonTriggerEfficiencies_eta2p1_Run2012ABCD_v5trees.pkl'%prefix
		muon_id_file  = open(muon_id_filename)
		muon_iso_file = open(muon_iso_filename)
		muon_trigger_file = open(muon_trigger_filename)
		self.muon_id_dict  = pickle.load(muon_id_file)
		self.muon_iso_dict = pickle.load(muon_iso_file)
		self.muon_trigger_dict = pickle.load(muon_trigger_file)


	def getSF(self,prewiggle_leta,prewiggle_lpt,pileup_events):
		tracking_eff_consts = self.tracking_eff_consts
		#tracking efficiency SFs
		for i in range(len(tracking_eff_consts)) :
			if prewiggle_leta > tracking_eff_consts[i][0] and prewiggle_leta < tracking_eff_consts[i][1] :
				self.weight_tracking     = tracking_eff_consts[i][2]
				self.weight_tracking_hi  = self.weight_tracking+tracking_eff_consts[i][3]
				self.weight_tracking_low = self.weight_tracking-tracking_eff_consts[i][3]

		#Lepton ID, isolation, and trigger scale factors
		muon_id_dict = self.muon_id_dict
		muon_iso_dict = self.muon_iso_dict		
		muon_trigger_dict = self.muon_trigger_dict
		#vertex corrections
		for vtx_key in muon_id_dict['Tight']['vtxpt20-500'].keys() :
			#print 'vtx_key = '+vtx_key+'comp 1 = '+vtx_key.split('_')[0]+' comp 2 = '+vtx_key.split('_')[1]+' pileup = '+str(pileup_events[0])
			if (pileup_events > float(vtx_key.split('_')[0]) and pileup_events < float(vtx_key.split('_')[1])) or (vtx_key == '28.5_30.5' and pileup_events > 30) :
				muon_id_vtx_weight = muon_id_dict['Tight']['vtxpt20-500'][vtx_key]['data/mc']['efficiency_ratio']
				muon_id_vtx_weight_low = muon_id_dict['Tight']['vtxpt20-500'][vtx_key]['data/mc']['err_low']
				muon_id_vtx_weight_hi = muon_id_dict['Tight']['vtxpt20-500'][vtx_key]['data/mc']['err_hi']
				muon_iso_vtx_weight = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['vtxpt20-500'][vtx_key]['data/mc']['efficiency_ratio']
				muon_iso_vtx_weight_low = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['vtxpt20-500'][vtx_key]['data/mc']['err_low']
				muon_iso_vtx_weight_hi = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['vtxpt20-500'][vtx_key]['data/mc']['err_hi']
				muon_trigger_vtx_weight = muon_trigger_dict['IsoMu24']['TightID_IsodB']['VTX'][vtx_key]['data']['efficiency']
				muon_trigger_vtx_weight_low = muon_trigger_dict['IsoMu24']['TightID_IsodB']['VTX'][vtx_key]['data']['err_low']
				muon_trigger_vtx_weight_hi = muon_trigger_dict['IsoMu24']['TightID_IsodB']['VTX'][vtx_key]['data']['err_hi']
		#eta corrections
		for eta_key in muon_id_dict['Tight']['etapt20-500'].keys() :
			if prewiggle_leta > float(eta_key.split('_')[0]) and prewiggle_leta < float(eta_key.split('_')[1]) :
				muon_id_eta_weight = muon_id_dict['Tight']['etapt20-500'][eta_key]['data/mc']['efficiency_ratio']
				muon_id_eta_weight_low = muon_id_dict['Tight']['etapt20-500'][eta_key]['data/mc']['err_low']
				muon_id_eta_weight_hi = muon_id_dict['Tight']['etapt20-500'][eta_key]['data/mc']['err_hi']
				muon_iso_eta_weight = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['etapt20-500'][eta_key]['data/mc']['efficiency_ratio']
				muon_iso_eta_weight_low = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['etapt20-500'][eta_key]['data/mc']['err_low']
				muon_iso_eta_weight_hi = muon_iso_dict['combRelIsoPF04dBeta<012_Tight']['etapt20-500'][eta_key]['data/mc']['err_hi']
				muon_trigger_eta_weight = muon_trigger_dict['IsoMu24']['TightID_IsodB']['ETA'][eta_key]['data']['efficiency']
				muon_trigger_eta_weight_low = muon_trigger_dict['IsoMu24']['TightID_IsodB']['ETA'][eta_key]['data']['err_low']
				muon_trigger_eta_weight_hi = muon_trigger_dict['IsoMu24']['TightID_IsodB']['ETA'][eta_key]['data']['err_hi']
		#pt corrections
		nextKey = ''
		nextKey_trig = ''
		if abs(prewiggle_leta) < 0.9 :
			nextKey = 'ptabseta<0.9'
			nextKey_trig = 'PT_ABSETA_Barrel_0to0p9'
		elif abs(prewiggle_leta) < 1.2 :
			nextKey = 'ptabseta0.9-1.2'
			nextKey_trig = 'PT_ABSETA_Transition_0p9to1p2'
		elif abs(prewiggle_leta) < 2.1 :
			nextKey = 'ptabseta1.2-2.1'
			nextKey_trig = 'PT_ABSETA_Endcaps_1p2to2p1'
		elif abs(prewiggle_leta) < 2.4 :
			nextKey = 'ptabseta2.1-2.4'
			print 'WARNING: lepton found with abs(eta)>2.4! Cannot apply pt-based trigger SF!'
		else :
			print 'WARNING: lepton found with abs(eta) > 2.4! Cannot apply pt-based ID, iso, or trigger SFs!'
		for pt_key in muon_id_dict['Tight'][nextKey].keys() :
			if (prewiggle_lpt > float(pt_key.split('_')[0]) and prewiggle_lpt < float(pt_key.split('_')[1])) or (pt_key == '140_300' and prewiggle_lpt > 300) :
				muon_id_pt_weight = muon_id_dict['Tight'][nextKey][pt_key]['data/mc']['efficiency_ratio']
				muon_id_pt_weight_low = muon_id_dict['Tight'][nextKey][pt_key]['data/mc']['err_low']
				muon_id_pt_weight_hi = muon_id_dict['Tight'][nextKey][pt_key]['data/mc']['err_hi']
				muon_iso_pt_weight = muon_iso_dict['combRelIsoPF04dBeta<012_Tight'][nextKey][pt_key]['data/mc']['efficiency_ratio']
				muon_iso_pt_weight_low = muon_iso_dict['combRelIsoPF04dBeta<012_Tight'][nextKey][pt_key]['data/mc']['err_low']
				muon_iso_pt_weight_hi = muon_iso_dict['combRelIsoPF04dBeta<012_Tight'][nextKey][pt_key]['data/mc']['err_hi']
		for pt_key in muon_trigger_dict['IsoMu24']['TightID_IsodB'][nextKey_trig].keys() :
			if prewiggle_lpt > float(pt_key.split('_')[0]) and prewiggle_lpt < float(pt_key.split('_')[1]) :
				muon_trigger_pt_weight = muon_trigger_dict['IsoMu24']['TightID_IsodB'][nextKey_trig][pt_key]['data']['efficiency']
				muon_trigger_pt_weight_low = muon_trigger_dict['IsoMu24']['TightID_IsodB'][nextKey_trig][pt_key]['data']['err_low']
				muon_trigger_pt_weight_hi = muon_trigger_dict['IsoMu24']['TightID_IsodB'][nextKey_trig][pt_key]['data']['err_hi']
		self.weight_lep_ID       = muon_id_vtx_weight*muon_id_eta_weight*muon_id_pt_weight
		self.weight_lep_ID_low   = self.weight_lep_ID-self.weight_lep_ID*math.sqrt((muon_id_vtx_weight_low/muon_id_vtx_weight)**2+(muon_id_eta_weight_low/muon_id_eta_weight)**2+(muon_id_pt_weight_low/muon_id_pt_weight)**2)
		self.weight_lep_ID_hi    = self.weight_lep_ID+self.weight_lep_ID*math.sqrt((muon_id_vtx_weight_hi/muon_id_vtx_weight)**2+(muon_id_eta_weight_hi/muon_id_eta_weight)**2+(muon_id_pt_weight_hi/muon_id_pt_weight)**2)

		self.weight_lep_iso      = muon_iso_vtx_weight*muon_iso_eta_weight*muon_iso_pt_weight
		self.weight_lep_iso_low  = self.weight_lep_iso-self.weight_lep_iso*math.sqrt((muon_iso_vtx_weight_low/muon_iso_vtx_weight)**2+(muon_iso_eta_weight_low/muon_iso_eta_weight)**2+(muon_iso_pt_weight_low/muon_iso_pt_weight)**2)
		self.weight_lep_iso_hi   = self.weight_lep_iso+self.weight_lep_iso*math.sqrt((muon_iso_vtx_weight_hi/muon_iso_vtx_weight)**2+(muon_iso_eta_weight_hi/muon_iso_eta_weight)**2+(muon_iso_pt_weight_hi/muon_iso_pt_weight)**2)

		self.weight_trig_eff     = muon_trigger_vtx_weight*muon_trigger_eta_weight*muon_trigger_pt_weight
		self.weight_trig_eff_low = self.weight_trig_eff-self.weight_trig_eff*math.sqrt((muon_trigger_vtx_weight_low/muon_trigger_vtx_weight)**2+(muon_trigger_eta_weight_low/muon_trigger_eta_weight)**2+(muon_trigger_pt_weight_low/muon_trigger_pt_weight)**2)
		self.weight_trig_eff_hi  = self.weight_trig_eff+self.weight_trig_eff*math.sqrt((muon_trigger_vtx_weight_hi/muon_trigger_vtx_weight)**2+(muon_trigger_eta_weight_hi/muon_trigger_eta_weight)**2+(muon_trigger_pt_weight_hi/muon_trigger_pt_weight)**2)

