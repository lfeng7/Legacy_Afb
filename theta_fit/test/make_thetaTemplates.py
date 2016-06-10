"""
main code that produce templates for theta
"""
import Legacy_Afb.theta_fit.thetaTempMaker as thetaTempMaker 
import sys

argv=sys.argv[1:]
if len(argv)==0:
	print """
	Usage: 
	python python make_thetaTemplates.py temp_angles/angles_data_SingleEl_Run2012ABCD.root data ttree
	or
	python make_thetaTemplates.py combined_tempMerge/templates.root test
	"""
	sys.exit(1)

inputfile = argv.pop(0)
output_name = argv.pop(0)

argv = ' '.join(argv)
if 'ttree' in argv:
	isTTree = True
else:
	isTTree = False

print 'create a thetaTemp object.'
self = thetaTempMaker.thetaTemp(outputName=output_name, inputFile=inputfile, isTTree=isTTree)
self.main()
print 'All finished.'


