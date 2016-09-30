#imports
from ROOT import *
from array import array
import numpy

#global variables
#histogram limits
# x is c*, y is xf, z is mtt
XBINS = numpy.arange(-1,1.1,0.1)
#YBINS = array('d',[0.,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.7]) # fine binning
#YBINS = array('d',[0.,0.05,0.15,0.3,0.7]) # coarse binning
YBINS = array('d',[0.,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.275,0.3,0.7]) # more fine binning
ZBINS = array('d',[350.,400,450,500,550,600,650,700,750,800,850,900,950,1000])#,1750])

binx = [20,-1,1]
biny = [30,0,0.6]
binz = [40,350,1750]

#TDR Style
#gROOT.Macro('rootlogon.C')

##############################		   Template Class  		##############################

class template :
	#docstring
	"""template class"""
	
	#__init__function
	def __init__(self,name,formatted_name,bin_type='fixed') :
		print '				Adding template with name '+name
		self.name = name
		self.formatted_name = formatted_name
		self.bin_type = bin_type
		self.createHists()

	def createHists(self):
		if self.bin_type !='fixed':
			self.histo_3D = TH3D(self.name     ,self.formatted_name+'; c*; |x_{F}|; M (GeV)',len(XBINS)-1,XBINS,len(YBINS)-1,YBINS,len(ZBINS)-1,ZBINS)
			self.histo_x  = TH1D(self.name+'_x',self.formatted_name+' X Projection; c*',len(XBINS)-1,XBINS)
			self.histo_y  = TH1D(self.name+'_y',self.formatted_name+' Y Projection; |x_{F}|',len(YBINS)-1,YBINS)
			self.histo_z  = TH1D(self.name+'_z',self.formatted_name+' Z Projection; M (GeV)',len(ZBINS)-1,ZBINS)
		else:
			self.histo_3D = TH3D(self.name     ,self.formatted_name+'; c*; |x_{F}|; M (GeV)',binx[0],binx[1],binx[2],biny[0],biny[1],biny[2],binz[0],binz[1],binz[2])
			self.histo_x  = TH1D(self.name+'_x',self.formatted_name+' X Projection; c*',binx[0],binx[1],binx[2])
			self.histo_y  = TH1D(self.name+'_y',self.formatted_name+' Y Projection; |x_{F}|',biny[0],biny[1],biny[2])
			self.histo_z  = TH1D(self.name+'_z',self.formatted_name+' Z Projection; M (GeV)',binz[0],binz[1],binz[2])
		self.histo_3D.SetDirectory(0); self.histo_x.SetDirectory(0); self.histo_y.SetDirectory(0); self.histo_z.SetDirectory(0)

	def getTemplateProjections(self,histo_1D):
		"""
		Convert a unrolled 1D template into the original 3D hist and three projections
		"""
		self.make_from_1D_histo(histo_1D)
		return [self.histo_3D,self.histo_x,self.histo_y,self.histo_z]

	def getOriginalTemps(self):
		if self.histo_x.Integral==0:
			print 'No original templates filled! Will exit.'
			sys.exit(1)
		return [self.histo_3D,self.histo_x,self.histo_y,self.histo_z]
		
	# def Fill(self,c,x,m,w) :
	# 	self.histo_3D.Fill(c,x,m,w)
	# 	self.histo_x.Fill(c,w)
	# 	self.histo_y.Fill(x,w)
	# 	self.histo_z.Fill(m,w)

	def Fill(self,c,x,m,w) :
		# determine if this event is in the range of all three variables
		if self.bin_type != 'fixed':
			inxbounds = c>=XBINS[0] and c<XBINS[len(XBINS)-1]
			inybounds = x>=YBINS[0] and x<YBINS[len(YBINS)-1]
			inzbounds = m>=ZBINS[0] and m<ZBINS[len(ZBINS)-1]
		else:
			inxbounds = c>=binx[1] and c<binx[2]
			inybounds = x>=biny[1] and x<biny[2]
			inzbounds = m>=binz[1] and m<binz[2]
		# Fill 3D and 3 1D projections
		if inxbounds and inybounds and inzbounds :
			self.histo_3D.Fill(c,x,m,w)
			self.histo_x.Fill(c,w)
			self.histo_y.Fill(x,w)
			self.histo_z.Fill(m,w)

	#convertTo1D takes a 3D distribution and makes it 1D for use with theta
	def convertTo1D(self) :
		nBins = self.histo_3D.GetNbinsX()*self.histo_3D.GetNbinsY()*self.histo_3D.GetNbinsZ()
		newHisto = TH1F(self.histo_3D.GetName()+'_1D',self.histo_3D.GetTitle(),nBins,0.,nBins-1.)
		newHisto.SetDirectory(0)
		realbincounter = 1
		nglobalbins = self.histo_3D.GetSize()
		for k in range(nglobalbins) :
			if not self.histo_3D.IsBinOverflow(k) and not self.histo_3D.IsBinUnderflow(k) :
				if not self.histo_3D.GetBinContent(k) < 0. :
					newHisto.SetBinContent(realbincounter,self.histo_3D.GetBinContent(k))
					newHisto.SetBinError(realbincounter,self.histo_3D.GetBinError(k))
				realbincounter+=1
		return newHisto

	#make_from_1D_histo takes a 1D distribution and makes a template out of it!
	def make_from_1D_histo(self,histo_1D) :
		nglobalbins = self.histo_3D.GetSize()
		global1Dbincounter = 1
		for k in range(nglobalbins) :
			if not self.histo_3D.IsBinOverflow(k) and not self.histo_3D.IsBinUnderflow(k) :
				content = histo_1D.GetBinContent(global1Dbincounter)
				error   = histo_1D.GetBinError(global1Dbincounter)
				self.histo_3D.SetBinContent(k,content)
				self.histo_3D.SetBinError(k,error)
				binx = array('i',[0]); biny = array('i',[0]); binz = array('i',[0])
				self.histo_3D.GetBinXYZ(k,binx,biny,binz)
				self.histo_x.SetBinContent(binx[0],self.histo_x.GetBinContent(binx[0])+content)
				self.histo_y.SetBinContent(biny[0],self.histo_y.GetBinContent(biny[0])+content)
				self.histo_z.SetBinContent(binz[0],self.histo_z.GetBinContent(binz[0])+content)
				self.histo_x.SetBinError(binx[0],self.histo_x.GetBinError(binx[0])+error*error)
				self.histo_y.SetBinError(biny[0],self.histo_y.GetBinError(biny[0])+error*error)
				self.histo_z.SetBinError(binz[0],self.histo_z.GetBinError(binz[0])+error*error)
				global1Dbincounter+=1
		hs = [self.histo_x,self.histo_y,self.histo_z]
		for h in hs :
			for k in range(h.GetSize()) :
				if not h.IsBinUnderflow(k) and not h.IsBinOverflow(k) :
					h.SetBinError(k,sqrt(h.GetBinError(k)))

