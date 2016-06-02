# 
import ROOT
import template

class thetaTemp(object):
	"""
	Take a root file w/ ttree (for Data) or smoothed 3D templates (for MC) 
	and create a single root file with all templates 
	"""
	def __init__(self, outputName, inputFile):
		"""
		inputFile is the path of input template root file
		outputName is the name for output thetaTemp.root file
		"""
		super(thetaTemp, self).__init__()
		self.process = {}
		self.outfile = ROOT.TFile('thetaTemplate_%s.root'%outputName,'recreate')
		self.thetaHistList = []
		self.template_file = inputFile

	def main(self):
		self.define_process()
		self.importTemp(self.template_file)

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
	

	def importTemp(self,tmp_file):
		"""
		process smoothed templates.root file
		The templates.root file has keys in the form below, with seperate th3f of plus and minus lep charge
		KEY: TH3F	f_WJets_minus;1	f_WJets_minus
		KEY: TH3F	f_WJets_plus;1	f_WJets_plus
		KEY: TH3F	f_bck_minus;1	f_bck_minus
		KEY: TH3F	f_bck_plus;1	f_bck_plus
		"""
		tfile = ROOT.TFile(tmp_file)
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
					# use template class to convert the TH3D to TH1D
					template_obj = template.template(i_newkey,i_newkey)
					template_obj.histo_3D = tmp_hist
					hist_1D = template_obj.convertTo1D()
					hist_1D.SetName(i_newkey)
					hist_1D.SetTitle(i_newtitle)
					hist_1D.SetDirectory(0)
					self.thetaHistList.append(hist_1D)
		tfile.Close()

	def __del__(self):
		self.outfile.cd()
		for item in self.thetaHistList:
			item.Write()
		self.outfile.Close()
		print 'Closeup thetaTemp object.'






		



