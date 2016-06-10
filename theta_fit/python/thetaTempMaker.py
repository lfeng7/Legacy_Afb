# 
import ROOT
import template
import helper

ROOT.gROOT.SetBatch(True)

class thetaTemp(object):
	"""
	Take a root file w/ ttree (for Data) or smoothed 3D templates (for MC) 
	and create a single root file with all templates 
	"""
	def __init__(self, outputName, inputFile, isTTree=False, weight=1):
		"""
		inputFile is the path of input template root file
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
		self.weight = weight
		self.QCD_SF = 0.2

	def main(self):
		self.define_process()
		if self.isTTree:
			print '#### Process template file w/ TTree. ####'
			self.TTtreeToTemplates(ttree_file=self.template_file,weight=self.weight)
		else:
			print '#### Process template file w/ TH3D. ####'
			self.importTemp(self.template_file)
		self.makeControlPlots()
		print 'Done thetaTemp.main()'

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
					print 'histogram %s not found in %s!'%(i_oldkey,tfile)
				else:
					# rescale Data-driven QCD templates using extrapolation factors get from ABCD
					if 'qcd' in i_newkey:
						print 'QCD templates rescaling to %.3f'%self.QCD_SF
						tmp_hist.Scale(self.QCD_SF)
					# use template class to convert the TH3D to TH1D
					# 	def __init__(self,name,formatted_name,bin_type='fixed') :
					template_obj = template.template(i_newkey,i_newkey)
					template_obj.histo_3D = tmp_hist
					self.Add_1D_temp(template=template_obj,tempName=i_newkey,tempTitle=i_newtitle)
					# also include original template TH3D and projections in aux file for later cross check
					# original_temps
					# self.write_templates_to_auxfile()
		print 'Done importTemp'
		tfile.Close()

	def TTtreeToTemplates(self,ttree_file,weight=1,isData=True):
		"""
		Input: a root with ttree
		Output: add 1D templates into template file, and control plots to aux file
		"""
		print 'Process TTtreeToTemplates with file %s, weight %s, isData=%s'%(ttree_file,weight,isData)

		tfile = ROOT.TFile(ttree_file)
		self.outfile.cd()
		tree_name = helper.GetTTreeName(tfile)
		ttree = tfile.Get(tree_name)
		# set up template name
		if isData:
			tmpProcess='DATA'
		tmpNames=[]
		tmpNames.append('f_plus__%s'%tmpProcess)
		tmpNames.append('f_minus__%s'%tmpProcess)
		# create template objects
		tmp_objects=[]
		for item in tmpNames:
			tmp_objects.append(template.template(name=item,formatted_name=item))
		# Loop over entries in 'data' ttree and fill templates
		n_entries = ttree.GetEntries()
		for iev in range(n_entries):
			ttree.GetEntry(iev)
			# load observables
			cs = ttree.cos_theta_cs
			xf = ttree.Feynman_x
			mtt = ttree.ttbar_mass
			lep_charge = ttree.Q_l
			# fill templates
			if lep_charge>0:
				tmp_objects[0].Fill(cs,abs(xf),mtt,weight)
			elif lep_charge<0:
				tmp_objects[1].Fill(cs,abs(xf),mtt,weight)
			else:
				print 'lep_charge==0. Something is wrong!'
				sys.exit(1)
		# Write proper unrolled 1D templates into thetaTemp file
		for i,itemp in enumerate(tmp_objects):
			self.Add_1D_temp(template=tmp_objects[i],tempName=tmpNames[i],tempTitle=tmpNames[i])
			# further, write original projections into aux file for later conparisons
			self.output_original_templates(template=tmp_objects[i])
		print 'Done TTtreeToTemplates.'


	def makeControlPlots(self):
		"""
		Take 1D templates and create all control plots, such as 3D and x,y,z projections
		Write all control stuff in an aux.root file
		"""	
		# for every 1D templates, get 3D original and 3 1D projection histograms
		for ihist in self.thetaHistList:
			hname = ihist.GetName()+'_proj'
			template_obj =  template.template(hname,hname+' projected back from 1D hist')
			# 	getTemplateProjections	return [self.histo_3D,self.histo_x,self.histo_y,self.histo_z]
			hist_proj = template_obj.getTemplateProjections(ihist)
			self.write_templates_to_auxfile(hist_list=hist_proj)
			del(template_obj)
		print 'Done makeControlPlots.'


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
		print 'Closeup thetaTemp object.'



