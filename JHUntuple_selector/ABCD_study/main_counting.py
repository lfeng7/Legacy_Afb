# main code to do the ABCD counting
from simple_counting import *

argv = sys.argv[1:]
if len(argv)< 1 or ('help' in argv):
	print """
	Usage: python  main_counting.py dir outputname makeplots (yes or no)
	"""
	sys.exit(1)

inputdir = argv.pop(0)
outname = argv.pop(0)

makeplots = False
if len(argv)>0:
	makeplots = argv.pop(0)
	if makeplots=='yes':
		makeplots = True
	else:
		makeplots = False

verbose = True

def setABCD1():
	# Define ABCD here 
	self.basecuts = 'njets>=4&&lep_isTight&&n_btags==2'
	self.outname = 'met_iso_ABCD'

	x1 = 'met_pt_vec<20'
	x2 = 'met_pt_vec>=20'
	y1 = 'lep_iso<0.1'
	y2 = 'lep_iso>=0.1 && lep_iso<0.2'
	y3 = 'lep_iso>=0.2 && lep_iso<1.2'
	self.regions['A'] = '%s && %s'%(x1,y3)
	self.regions['B'] = '%s && %s'%(x1,y1)
	self.regions['C'] = '%s && %s'%(x2,y3)
	self.regions['D'] = '%s && %s'%(x2,y1)
	self.regions['E'] = '%s && %s'%(x1,y2)
	self.regions['F'] = '%s && %s'%(x2,y2)

	self.main()
	self.printQCDprojection(regA='A',regB='B',regC='C',regD='D')
	self.printQCDprojection(regA='E',regB='B',regC='F',regD='D')
	# make control plots
	xcuts = [x1,x2]
	if makeplots:
		for i,icut in enumerate(xcuts):
			self.makeControlPlots(cut=icut,cutNum=i)
	self.save()


def setABCD2():
	# Define ABCD here 
	self.basecuts = 'njets>=4 && n_btags==2'
	self.outname = 'ID_iso_ABCD'

	x1 = '!lep_isLoose'
	x2 = 'lep_isLoose && !lep_isTight'
	x3 = 'lep_isTight'
	y1 = 'lep_iso<0.1'
	y2 = 'lep_iso>=0.2 && lep_iso<1.2'
	self.regions['A'] = '%s && %s'%(x1,y2)
	self.regions['B'] = '%s && %s'%(x1,y1)
	self.regions['C'] = '%s && %s'%(x3,y2)
	self.regions['D'] = '%s && %s'%(x3,y1)
	self.regions['E'] = '%s && %s'%(x2,y2)
	self.regions['F'] = '%s && %s'%(x2,y1)

	self.main()
	self.printQCDprojection(regA='A',regB='B',regC='C',regD='D')
	self.printQCDprojection(regA='E',regB='F',regC='C',regD='D')
	# make control plots
	xcuts = [x1,x2,x3]
	if makeplots:
		for i,icut in enumerate(xcuts):
			self.makeControlPlots(cut=icut,cutNum=i)
	self.save()

def setABCD3():
	# Define ABCD here 
	self.basecuts = 'lep_isTight && n_btags==2'
	self.outname = 'njets_iso_ABCD'
	x1 = 'njets==2'
	x2 = 'njets==3'
	x3 = 'njets>=4'
	y1 = 'lep_iso<0.1'
	y2 = 'lep_iso>=0.2 && lep_iso<1.2'
	self.regions['A'] = '%s && %s'%(x1,y2)
	self.regions['B'] = '%s && %s'%(x1,y1)
	self.regions['C'] = '%s && %s'%(x3,y2)
	self.regions['D'] = '%s && %s'%(x3,y1)
	self.regions['E'] = '%s && %s'%(x2,y2)
	self.regions['F'] = '%s && %s'%(x2,y1)

	self.main()
	self.printQCDprojection(regA='A',regB='B',regC='C',regD='D')
	self.printQCDprojection(regA='E',regB='F',regC='C',regD='D')
	# make control plots
	xcuts = [x1,x2,x3]
	if makeplots:
		for i,icut in enumerate(xcuts):
			self.makeControlPlots(cut=icut,cutNum=i)
	self.save()

# main program

self = counter(inputdir=inputdir,verbose=verbose,outname=outname)
setABCD1()	
del(self)

self = counter(inputdir=inputdir,verbose=verbose,outname=outname)
setABCD2()
del(self)

self = counter(inputdir=inputdir,verbose=verbose,outname=outname)
setABCD3()
del(self)


