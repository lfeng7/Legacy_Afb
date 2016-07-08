# 
import ROOT
import template
import helper
import samples
import copy
from array import array

ROOT.gROOT.SetBatch(True)

class thetaTemp(object):
	"""
	Take a root file w/ ttree (for Data) or smoothed 3D templates (for MC) 
	and create a single root file with all templates 
	"""
	def __init__(self, outputName, inputFile, isTTree=False, txtfile=None, verbose = False, bin_type = 'fixed'):
		"""
		inputFile is a list of path of input template root files (for ttree) or a single file for MC
		outputName is the name for output thetaTemp.root file
		"""
		super(thetaTemp, self).__init__()
		self.process = {}
		self.outputName = outputName
		self.outfile = ROOT.TFile('templates/thetaTemplates_%s.root'%outputName,'recreate')
		self.thetaHistList = []
		self.template_file = inputFile
		self.outfile_aux = ROOT.TFile('templates/thetaTemplates_%s_aux.root'%self.outputName,'recreate')
		self.outfile_aux.mkdir('hists/')
		self.outfile_aux.mkdir('plots/')
		self.isTTree = isTTree
		self.QCD_SF = 0.06 # this number is from the ABCD method estimation
		self.verbose = verbose
		self.bin_type = bin_type
		self.AFB_sigma=1.0 # one sigma deviation of AFB from zero
		if txtfile is None:
			txtfile = 'MC_input_with_bkg.txt'
		self.samples_obj = samples.samples(txtfile)

	def main(self):
		self.define_process()
		self.define_sys()
		if self.isTTree:
			print '(info) Process template file w/ TTree.'
			self.assignFiles()
			self.GetWeights()
			# Loop over all processes, e.g., wjets, gg etc
			for key,value in self.process_files.iteritems():
				print '(info) Begin process %s samples.'%key
				self.TTtreeToTemplates(ttree_file=value, process_name=key)
		else:
			print '(info) Process template file w/ TH3D. '
			self.importTemp(self.template_file)
		self.makeControlPlots()
		print '(info) Done thetaTemp.main()'

	def assignFiles(self):
		"""
		Assign the files into a hashtable of processes
		Input: a list of root files as self.template_file
		output: a hashtable of files belong to which process
		"""
		process_files = {}
		for ifile in self.template_file:
			# if it's data , hard code it
			if 'Run' in ifile :
				if process_files.get('DATA',1):
					process_files['DATA'] = [ifile]
				else:
					process_files['DATA'].append(ifile)
			# if not, it's an MC and will look up it's process type and add into hashtable	
			else:					
				# load info about samples from txt file
				sample_info_obj = self.samples_obj.get_sample_info(ifile)
				# if cannot find info just skip this file for good
                                if sample_info_obj is None: continue
				sample_type = sample_info_obj.type
				# Add file path to dictionary
				if process_files.get(sample_type,0)==0:
					process_files[sample_type] = [ifile]
		#			if self.verbose: print 'Add type %s and file %s'%(sample_type,ifile)
				else:
					process_files[sample_type].append(ifile)
		#			if self.verbose: print 'Add sample %s'%ifile
		self.process_files = process_files
		if self.verbose:
			print '(DEBUG) %s'%process_files
		print '(info) Done assignFiles.'


	def define_process(self):
		"""
		define template processes in a dict. [0] is string correspond to legacy template th3f name, [1] is 
		more detailed describtion of what's in there.
		"""
		self.process['wjets']=['WJets','wjets']
		self.process['other']  =['bck','single_top/non-semilep_ttbar']
		self.process['gg']   =['gg','gg/qg_ttbar']
		self.process['qqs']   =['qqs','symetric qq_ttbar']
		self.process['qcd'] =['ntmj','qcd']
		print '(info) Done define_process'


	def define_sys(self):
		"""
		define types of systematic variation
		"""
		self.systematics_vary = ['btag_eff_reweight','lepID_reweight','trigger_reweight']
		self.systematics_fixed = ['top_pT_reweight','tracking_reweight','lepIso_reweight','pileup_reweight']
		self.systematics_all = self.systematics_fixed+self.systematics_vary
		self.sys_versions = ['plus','minus'] # naming conventions for theta
		self.sys_names =['up','down'] # naming coventions for legacy ttree
		print '(info) Done define_sys'


	def getTempKeys(self,old_key_string,new_key_string,title):
		# old_keys are original name of th3f in templates.root, new_keys are names of TH1D in thetaTemp.root
		old_keys,new_keys,new_titles = [],[],[]
		old_keys.append('f_%s_plus'%old_key_string) 
		old_keys.append('f_%s_minus'%old_key_string) 
		new_keys.append('f_plus__%s'%new_key_string)
		new_keys.append('f_minus__%s'%new_key_string) 
		new_titles.append('%s lep charge plus'%title)
		new_titles.append('%s lep charge minus'%title)
		return old_keys,new_keys,new_titles
	
	def Add_1D_temp(self,template,tempName,tempTitle):
		"""
		input: an object of template class
		output: add 1D template into thetaTemp.root 
		"""
		self.outfile.cd()
		hist_1D = template.convertTo1D()
		hist_1D.SetName(tempName)
		hist_1D.SetTitle(tempTitle)
		hist_1D.SetDirectory(0)
		hist_1D.Write()
		self.thetaHistList.append(hist_1D)

	def output_original_templates(self,template):
		"""
		input: an object of template class
		output: add original 3D and 3 1D projections into aux file
		"""
		original_templates=template.getOriginalTemps() # 3D,x,y,z
		self.write_templates_to_auxfile(hist_list=original_templates)


	def importTemp(self,tmp_file):
		"""
		process smoothed templates.root file
		input: original root file with TH3D templates
		output: new root file for theta feedin file with an unrolled TH1D

		The templates.root file has keys in the form below, with seperate th3f of plus and minus lep charge
		KEY: TH3F	f_WJets_minus;1	f_WJets_minus
		KEY: TH3F	f_WJets_plus;1	f_WJets_plus
		KEY: TH3F	f_bck_minus;1	f_bck_minus
		KEY: TH3F	f_bck_plus;1	f_bck_plus
		"""
		tfile = ROOT.TFile(tmp_file)
		self.outfile.cd()
		# loop over all processes,e.g., wjets, gg , and convert 3D templates into 1D and write into new template file
		for key,value in self.process.iteritems():
			old_key_string = value[0]
			new_key_string = key
			hist_title = value[1]
			# get keys of histograms in original templates and in new thetaTemplates
			old_keys,new_keys,new_titles = self.getTempKeys(old_key_string=old_key_string,new_key_string=new_key_string,title=hist_title)
			# Loop over all TH3D for this process, such as f_wjets_plus,f_wjets_minus,f_wjets_minus_btag_up etc
			for i,value in enumerate(old_keys):
				i_oldkey = old_keys[i]
				i_newkey =  new_keys[i]
				i_newtitle = new_titles[i]
				tmp_hist = tfile.Get(i_oldkey)
				if not tmp_hist:
					print '(debug) histogram %s not found in %s!'%(i_oldkey,tfile)
				else:
					# rescale Data-driven QCD templates using extrapolation factors get from ABCD
					if 'qcd' in i_newkey:
						print '(info) QCD templates rescaling to %.3f'%self.QCD_SF
						tmp_hist.Scale(self.QCD_SF)
					# use template class to convert the TH3D to TH1D
					# 	def __init__(self,name,formatted_name,bin_type='fixed') :
					template_obj = template.template(i_newkey,i_newkey,self.bin_type)
					template_obj.histo_3D = tmp_hist
					self.Add_1D_temp(template=template_obj,tempName=i_newkey,tempTitle=i_newtitle)
					# also include original template TH3D and projections in aux file for later cross check
					# original_temps
					# self.write_templates_to_auxfile()
		print '(info) Done importTemp'
		tfile.Close()

	def TTtreeToTemplates(self,ttree_file,process_name = None):
		"""
		Input: a root with ttree
		Output: add 1D templates into template file, and control plots to aux file
		for data, just one template
		for MC, will include nominal, and sys_up and down templates for each sys
		"""
		print '(info) Process TTtreeToTemplates with file %s'%(ttree_file)
		# set up template name and weight
		if process_name=='DATA':
			self.makeTemplates(file_list=ttree_file,process_name=process_name)
			print '(info) Done TTtreeToTemplates for %s.'%process_name
		else:
			# Loop over files in the list
			weight_lists = {} # e.g. wjets__btag_eff_reweight__up:[w1,w2,..w_up]
			norm_weights = []
			# Get normalization weights for each file in the list
			for ifile in ttree_file:					
				# load info about samples from txt file
				sample_info_obj = self.samples_obj.get_sample_info(ifile)
				norm_weight = sample_info_obj.weight
				norm_weights.append(norm_weight)
			# prepare for sys weights of each variation 
			# first make the nominal templates
			weight = []
			weight += self.fixed_weight+self.varied_weight_nominal
			template_name = process_name # e.g. wjets
			weight_lists[template_name] = weight
			# then make up and down templates for each sys for all MC templates BUT QCD
			if process_name != 'qcd':
				for i,isys in enumerate(self.systematics_vary):
					for j,jversion in enumerate(self.sys_versions): # ['plus','minus']
						weight = []
						# append fixed weights
						weight += self.fixed_weight
						# add a variant of isys weight
						weight.append(self.varied_weight[jversion][i])
						# add the rest of variable sys using the nominal values
						nominal_weights = copy.copy(self.varied_weight_nominal)
						nominal_weights.pop(i)
						weight += nominal_weights
						# set up names for this variant template: e.g. wjets__btag_eff_reweight__up
						template_name = '%s__%s__%s'%(process_name,isys,jversion)
						# finally, add a new template
						weight_lists[template_name]=weight
			# finally make a template for each version of weight list
			for key,value in weight_lists.iteritems():
				self.makeTemplates(file_list=ttree_file,process_name=key,weight_list=value,norm_weights=norm_weights)
			print '(info) Done making all templates for MC.'
	

	def GetWeights(self) :
		"""
		Make lists contains strings of the names of the weights in the ttree
		self.fixed_weight = []
		self.varied_weight_nominal = []
		self.varied_weight = {}
		"""
		# initiate
		self.fixed_weight = []
		self.varied_weight_nominal = []
		self.varied_weight = {'plus':[],'minus':[]}
		# set up fixed sys reweight 
		for iweight in self.systematics_fixed:
			self.fixed_weight.append(iweight)
		# set up varying weight
		for iweight in self.systematics_vary:
			# first set nominal ones
			self.varied_weight_nominal.append(iweight)
			# set up high version
			weight_key_up = '%s_hi'%iweight # like btag_eff_reweight_high
			self.varied_weight['plus'].append(weight_key_up)
			# set up low version
			weight_key_down = '%s_low'%iweight
			self.varied_weight['minus'].append(weight_key_down)
		print '(info) Done GetWeights'

	def makeTemplates(self,file_list,process_name,weight_list=[1.0],norm_weights = [1.0]):
		"""
		input: list of ttree, list of weight ( a list of 1D arrays ), tmp_name
		output: a 1D template and control plots of 3d projections
		"""
		if self.verbose:
			print '(DEBUG) makeTemplates process %s'%process_name

		tmpNames=[]
		tmpNames.append('f_plus__%s'%process_name)
		tmpNames.append('f_minus__%s'%process_name)
		# accomodate fqq,fgg specially
		# Determing if it is gg/qq process, if so , need to addtwice
		add_twice=False
		if 'gg' in process_name or 'qq' in process_name:
			add_twice = True
		# Add two additional templates for qq, qq__AFB__plus and qq__AFB__minus, but do this only during making nominal qq templates
		is_qq = False
		if process_name == 'qq':
			tmpNames.append('f_plus__%s__AFB__plus'%process_name)
			tmpNames.append('f_plus__%s__AFB__minus'%process_name)
			tmpNames.append('f_minus__%s__AFB__plus'%process_name)			
			tmpNames.append('f_minus__%s__AFB__minus'%process_name)	
			is_qq = True		
		# create template objects
		tmp_objects=[]
		for item in tmpNames:
			tmp_objects.append(template.template(name=item,formatted_name=item,bin_type=self.bin_type))
		if self.verbose: print '(DEBUG) all weights ',weight_list
		# Loop over ttree in list of ttrees
		for i,ifile in enumerate(file_list):
			# loading ttree and load branches
			tfile = ROOT.TFile(ifile)
			self.outfile.cd()
			tree_name = helper.GetTTreeName(tfile)
			ttree = tfile.Get(tree_name)
			weights = weight_list
			# check if norm weight is loaded correctly
			if self.verbose:
				print '(DEBUG) %s norm weight = %.3f'%(ifile,norm_weights[i])
			# Loop over entries in ttree and fill templates
			n_entries = ttree.GetEntries()
			for iev in range(n_entries):
				ttree.GetEntry(iev)
				# load observables
				cs = ttree.cos_theta_cs
				xf = ttree.Feynman_x
				mtt = ttree.ttbar_mass
				lep_charge = ttree.Q_l
				# get the weight right by looping over a list of arrays(or float)
				if process_name=='DATA': # for data, with no weights
					total_weight = 1
				elif process_name=='qcd':
					# QCD ttree is a sum of data and MC events in sideband, with MC events having negative weights for substraction purpose
					norm_weight = ttree.normalization_weight*self.QCD_SF
					total_weight = norm_weight
				else:
					norm_weight = norm_weights[i]
					total_weight = [getattr(ttree,item) for item in weights]
					total_weight.append(norm_weight)
					if self.verbose and iev<1: print '(DEBUG) total_weight=',total_weight
					total_weight = helper.multiply(total_weight)
					if self.verbose and iev<1: print '(DEBUG) total_weight=%.3f'%total_weight
				# for added twice case:
				if add_twice:
					motherPIDs = ttree.motherPIDs
					if 21 in motherPIDs or -21 in motherPIDs : add_twice=False
					else: total_weight *= 0.5
					w_a = ttree.w_a
				# fill templates
				# special note for qq sample: tmp_obj[2,3,4,5] are f_plus_up,down, f_minus_up,down
				if lep_charge>0:
					tmp_objects[0].Fill(cs,abs(xf),mtt,total_weight)
					if add_twice:
						tmp_objects[1].Fill(-cs,abs(xf),mtt,total_weight)
					if is_qq:
						tmp_objects[2].Fill(cs,abs(xf),mtt,total_weight*(1+self.AFB_sigma*w_a)) # Q>0, fwd
						tmp_objects[3].Fill(cs,abs(xf),mtt,total_weight*(1-self.AFB_sigma*w_a)) # Q>0, bwd
						tmp_objects[4].Fill(-cs,abs(xf),mtt,total_weight*(1-self.AFB_sigma*w_a)) # Q<0, fwd
						tmp_objects[5].Fill(-cs,abs(xf),mtt,total_weight*(1+self.AFB_sigma*w_a)) # Q<0, bwd
				elif lep_charge<0:
					tmp_objects[1].Fill(cs,abs(xf),mtt,total_weight)
					if add_twice:
						tmp_objects[0].Fill(-cs,abs(xf),mtt,total_weight)
					if is_qq:
						tmp_objects[4].Fill(cs,abs(xf),mtt,total_weight*(1+self.AFB_sigma*w_a)) # Q<0, fwd
						tmp_objects[5].Fill(cs,abs(xf),mtt,total_weight*(1-self.AFB_sigma*w_a)) # Q<0, bwd
						tmp_objects[2].Fill(-cs,abs(xf),mtt,total_weight*(1-self.AFB_sigma*w_a)) # Q>0, fwd
						tmp_objects[3].Fill(-cs,abs(xf),mtt,total_weight*(1+self.AFB_sigma*w_a)) # Q>0, bwd
				else:
					print '(debug) lep_charge==0. Something is wrong!'
					sys.exit(1)
				if self.verbose:
					# print 'cs %.2f,xf %.2f,mtt %.2f,lep_charge %i,total_weight %.2f'%(cs,xf,mtt,lep_charge,total_weight)
					pass
		# Write proper unrolled 1D templates into thetaTemp file
		for i,itemp in enumerate(tmp_objects):
			self.Add_1D_temp(template=tmp_objects[i],tempName=tmpNames[i],tempTitle=tmpNames[i])
			# further, write original projections into aux file for later conparisons
			self.output_original_templates(template=tmp_objects[i])



	def makeControlPlots(self):
		"""
		Take 1D templates and create all control plots, such as 3D and x,y,z projections
		Write all control stuff in an aux.root file
		"""	
		# for every 1D templates, get 3D original and 3 1D projection histograms
		for ihist in self.thetaHistList:
			hname = ihist.GetName()+'_proj'
			template_obj =  template.template(hname,hname+' projected back from 1D hist', self.bin_type)
			# 	getTemplateProjections	return [self.histo_3D,self.histo_x,self.histo_y,self.histo_z]
			hist_proj = template_obj.getTemplateProjections(ihist)
			self.write_templates_to_auxfile(hist_list=hist_proj)
			del(template_obj)
		print '(info) Done makeControlPlots.'


	def write_templates_to_auxfile(self,hist_list):
		"""
		input: a list of [TH3D,x_proj,y_proj,z_proj]
		output: write histograms and nicer plots into aux.root file
		"""
		for i,item in enumerate(hist_list):
			if i!=0:
				item.SetMinimum(0)
			self.outfile_aux.cd('hists/')
			item.Write()
			# make nicer plots and write into another dir
			c = ROOT.TCanvas()
			c.SetName(item.GetName())
			item.Draw('hist E1')
			self.outfile_aux.cd('plots')
			c.Write()		

	def __del__(self):
		self.outfile.Close()
		self.outfile_aux.Close()
		print '(info) Closeup thetaTemp object.'



