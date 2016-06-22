# Count number of events in ABCD regions , and estimate number in signal region
import ROOT
import sys
import glob
import numpy as np
import re
import os

ROOT.gROOT.SetBatch(True)

argv = sys.argv[1:]
if len(argv)< 1 or ('help' in argv):
	print """
	Usage: python  simple_counting.py dir outputname "cuts"
	"""
	sys.exit(1)

# inputdir = argv.pop(0)
# outname = argv.pop(0)

# if len(argv)>0 :
#     cut = argv.pop(0)
# else:
#     cut = ''

class counter():
	"""docstring for counter"""
	def __init__(self, inputdir,outname,verbose=False):
		self.inputdir = inputdir

		self.verbose = verbose
		self.weight = 'weight_norm'
		self.outname=outname
		self.regions = {}
		self.printQCDcounts = 0
		self.progress = 0


	def main(self):
		print 'process findSamples '
		self.findSamples()
		print 'process getCounts_table'
		self.getCounts_table()
		print 'process printCounts'
		self.printCounts()
	    
		
	# Combine samples of the same type together
	def findSamples(self):
		# samples to plot
		self.samples = ['DY.*Jets.*','W.*Jets.*','TT_.*','T(bar)?_.*','QCD.*','SingleEl.*']
		self.sample_name = ['DYJets','WJets','TTbar','SingleTop','QCD MC','Data']
		self.all_files = glob.glob('%s/*.root'%self.inputdir)
		# Make a hashtable of form {'DYJets':['DY1.root','DY2.root']}
		self.sample_table = {}
		for i,item in enumerate(self.sample_name):
			ipattern = self.samples[i]
			iname = self.sample_name[i]
			ifiles = [item for item in self.all_files if re.match('.*/%s'%ipattern,item)]
			self.sample_table[iname] = ifiles


	# get normalized counts of several samples
	def getCounts(self,tchain,cuts,sample_type):
		self.show_progress()
		cuts = '%s&&%s'%(self.basecuts,cuts)
		if 'QCD' in sample_type: # do not scare QCD MC to data Lumi
			weight = 1.0
		else:
			weight = self.weight
		h1 = ROOT.TH1D('h1','h1',10,-1,5)
		tchain.Draw('n_btags>>h1','(%s)*(%s)'%(cuts,weight))
		if self.verbose:
			print 'cuts and weights: %s | %s'%(cuts,weight)
		return int(h1.Integral())


	# Loop over types of samples to get numbers
	def getCounts_table(self):
		self.count_table = {}
		# struc of table: {A,B,C,...}
		# A: [['DYJets',counts,fraction=counts/data in region A],['WJets',...]...]
		for isample,ifiles in self.sample_table.iteritems():
			chain = ROOT.TChain("selected")
			print '\nProcessing %s files.'%isample
			for item in ifiles:
				if self.verbose: print item
				chain.Add(item)
			# Loop over regions and fill the table, i.e. loop over ABCD... where iregion = A,B,..
			for iregion,icut in self.regions.iteritems():
				icounts = self.getCounts(tchain = chain, cuts=icut,sample_type=isample)
				ierr = np.sqrt(icounts)
				if self.count_table.get(iregion,0)==0 :
					self.count_table[iregion] = []
				self.count_table[iregion].append([isample,icounts,ierr])
		# Add QCD entry
		for iregion,value in self.count_table.iteritems():
			# each value is a column in the table, like A:[entry1,entry2...]
			# calculate total Data and MC counts, exculding QCD MC
			c_data = [item[1] for item in value if 'Data' in item[0]][0]
			c_MC = [item[1] for item in value if 'Data' not in item[0] and 'QCD' not in item[0]]
			total_c_MC = np.sum(c_MC)
			err_c_totalMC = np.sqrt(total_c_MC)
			# add an entry for total MC
			value.append(['all MC',total_c_MC,err_c_totalMC])
			# add an entry of QCD model for this column
			c_qcd_model = c_data - total_c_MC
			err_c_qcd_model = np.sqrt(total_c_MC+c_data)
			value.append(['QCD model',c_qcd_model,err_c_qcd_model])
			# add fraction for MC and data entries
			# each entry looks like ['DYJets',counts,error]
			for ientry in value:
				icount = ientry[1]
				ifrac = icount*1.0/c_data
				if icount == 0:
					err_frac = 0
				else:
					err_frac = ifrac*np.sqrt(1.0/abs(c_data)+1.0/abs(icount))
				ientry.extend([ifrac,err_frac])
		# final structure of count table
		# A:[['DYJets',counts,fraction],..]
		print 'Finish getting counts.'

	def lookupCount(self,region,sample):
		col = self.count_table.get(region,0)
		if col==0:
			print 'Cannot find %s %s in count table!'%(region,sample)
			sys.exit(1)
		sample_name = sample
		entry = [item for item in col if sample_name in item[0]][0]
		count = entry[1]
		count_err = entry[2]
		return count,count_err

	def QCDprojection(self,regA,regB,regC,regD):
		# Get data-QCD counts and error in region A,B and C
		count_A,err_count_A = self.lookupCount(region=regA,sample='QCD model')
		count_B,err_count_B = self.lookupCount(region=regB,sample='QCD model')
		count_C,err_count_C = self.lookupCount(region=regC,sample='QCD model')
		# calculate QCD projection in signal region and error
		count_qcd_sig = count_B*1.0/count_A*count_C
		err_count_qcd_sig = count_qcd_sig*np.sqrt(pow(err_count_A*1.0/count_A,2)+pow(err_count_B*1.0/count_B,2)+pow(err_count_C*1.0/count_C,2))
		# QCD projection fraction , Nqcd/N_mc_total in signal region
		count_all_mc_sig,err_count_all_mc_sig = self.lookupCount(regD,sample='all MC')
		R_qcd_sig = count_qcd_sig*1.0/count_all_mc_sig
		err_R_qcd_sig = R_qcd_sig*np.sqrt(pow(err_count_qcd_sig/count_qcd_sig,2)+pow(err_count_all_mc_sig/count_all_mc_sig,2))
		# return 
		return count_qcd_sig,err_count_qcd_sig,R_qcd_sig,err_R_qcd_sig

	def printCounts(self):
		self.fout = open('%s.txt'%self.outname,'w')
		self.fout.write('basecuts: %s\n\n'%self.basecuts)
		# print out counts table
		table_print = []
		for iregion,icut in self.regions.iteritems():
			to_print = '%s : %s\n'%(iregion,icut)
			# print header
			to_print += '%15s %10s %10s %12s %12s\n'%('sample','counts','err_c','fraction','err_frac')
			icol = self.count_table[iregion]
			for ientry in icol:
				# sample_name, counts, error, fraction
				to_print += '%15s %10i %10i %10.3f %% %10.3f %%\n'%(ientry[0],ientry[1],ientry[2],ientry[3]*100,ientry[4]*100)
			to_print += '\n'
			table_print.append([iregion,to_print])
		table_print.sort()
		for item in table_print:
			print item[1]
			self.fout.write(item[1])

	def printQCDprojection(self,regA='A',regB='B',regC='C',regD='D'):
		# print out QCD projection in signal region
		proj_print =''
		if self.printQCDcounts==0:
			proj_print += '\n'
			proj_print += 'QCD projection in Signal region\n'
			proj_print += '%15s %15s %15s %17s %17s\n'%('Type','QCD counts,','counts err,','R_qcd,','err R_qcd')
		proj_type = 'QCD in %s=%s/%s*%s'%(regD,regB,regA,regC)
		count_qcd_sig,err_count_qcd_sig,R_qcd_sig,err_R_qcd_sig = self.QCDprojection(regA=regA,regB=regB,regC=regC,regD=regD)
		proj_print += '%15s %15i %15i %17.3f %% %17.3f %%\n'%(proj_type,count_qcd_sig,err_count_qcd_sig,R_qcd_sig*100,err_R_qcd_sig*100)

		self.fout.write(proj_print)
		print proj_print
		self.printQCDcounts +=1

	def makeControlPlots(self,cut,cutNum):
		# RelIso 
		cmd = './masterplot.py --var lep_iso --Min 0 --Max 1.2 --xaxis RelIso --yaxis Events'
		cmd+= ' --title "%s && %s" --dir %s --stack yes --plotname %s_iso_%s'%(self.basecuts,cutNum,self.inputdir,self.outname,cutNum)
		cmd+= ' --cut "%s && %s"\n'%(self.basecuts,cut)
		# MET
		cmd+= './masterplot.py --var met_pt_vec --Min 0 --Max 150 --xaxis "MET(GeV)" --yaxis Events'
		cmd+= ' --title "%s && %s" --dir %s --stack yes --plotname %s_met_%s'%(self.basecuts,cutNum,self.inputdir,self.outname,cutNum)
		cmd+= ' --cut "%s && %s"\n'%(self.basecuts,cut)
		# lep_pT
		cmd+= './masterplot.py --var lep_pt --Min 25 --Max 150 --xaxis "ele pT(GeV)" --yaxis Events'
		cmd+= ' --title "%s && %s" --dir %s --stack yes --plotname %s_ele_pT_%s'%(self.basecuts,cutNum,self.inputdir,self.outname,cutNum)
		cmd+= ' --cut "%s && %s"\n'%(self.basecuts,cut)
		os.system(cmd)
		print cmd


	def show_progress(self, freq=2):
		"""Print a dot to stderr every 2 calls of getCounts (frequency can be changed)."""
		self.progress += 1
		if self.progress % freq == 1:
			sys.stderr.write('.')


	# close memory resident objects
	def save(self):
		self.fout.close()

# execute main function
if __name__ == '__main__':
	counter(inputdir=inputdir,verbose=False,outname=outname).main()
