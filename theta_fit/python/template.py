#imports
from ROOT import *
from array import array

#global variables
#histogram limits
XBINS = array('d',[-1.0,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.0])
YBINS = array('d',[0.,0.05,0.15,0.3,0.7])
ZBINS = array('d',[700.,900.,1100.,1300.,1500.,2500.])

#TDR Style
#gROOT.Macro('rootlogon.C')

##############################		   Template Class  		##############################

class template :
	#docstring
	"""template class"""
	
	#__init__function
	def __init__(self,name,formatted_name) :
		print '				Adding template with name '+name
		self.name = name
		self.formatted_name = formatted_name
		self.histo_3D = TH3D(name,	   formatted_name+'; c*; |x_{F}|; M (GeV)',len(XBINS)-1,XBINS,len(YBINS)-1,YBINS,len(ZBINS)-1,ZBINS)
		self.histo_x  = TH1D(name+'_x',formatted_name+' X Projection; c*',len(XBINS)-1,XBINS)
		self.histo_y  = TH1D(name+'_y',formatted_name+' Y Projection; |x_{F}|',len(YBINS)-1,YBINS)
		self.histo_z  = TH1D(name+'_z',formatted_name+' Z Projection; M (GeV)',len(ZBINS)-1,ZBINS)
		self.histo_3D.SetDirectory(0); self.histo_x.SetDirectory(0); self.histo_y.SetDirectory(0); self.histo_z.SetDirectory(0)

	def Fill(self,c,x,m,w) :
		self.histo_3D.Fill(c,x,m,w)
		self.histo_x.Fill(c,w)
		self.histo_y.Fill(x,w)
		self.histo_z.Fill(m,w)

	#convertTo1D takes a 3D distribution and makes it 1D for use with theta
	def convertTo1D(self) :
		nBins = self.histo_3D.GetNbinsX()*self.histo_3D.GetNbinsY()*self.histo_3D.GetNbinsZ()
		newHisto = TH1F(self.histo_3D.GetName(),self.histo_3D.GetTitle(),nBins,0.,nBins-1.)
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