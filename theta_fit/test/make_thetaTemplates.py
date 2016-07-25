"""
main code that produce templates for theta
"""
import Legacy_Afb.theta_fit.thetaTempMaker as thetaTempMaker 
import sys
import glob

argv=sys.argv[1:]
if len(argv)==0:
	print """
	Usage: 
	python python make_thetaTemplates.py temp_angles/angles_data_SingleEl_Run2012ABCD.root outputname  (ttree) (variable) (verbose)
	or
	python make_thetaTemplates.py combined_tempMerge/templates.root 
	"""
	sys.exit(1)

inputfile = argv.pop(0)
output_name = argv.pop(0)
#output_name = inputfile.split('/')[-1]
#output_name = 'thetaTemp_%s'%output_name

# some other options
isTTree = False
bin_type = 'variable'
isVerbose = False

argv = ' '.join(argv)
if 'ttree' in argv:
	isTTree = True
if 'fixed' in argv:
	bin_type = 'fixed'
if 'verbose' in argv:
	isVerbose = True

if isTTree:
	"Processing ttree files"
	inputfile = glob.glob(inputfile)
        output_name = 'thetaTemplates'

print 'create a thetaTemp object.'
print 'options: thetaTemp(outputName=%s, inputFile=%s, isTTree=%s, bin_type=%s, verbose=%s)'%(output_name,inputfile,isTTree,bin_type,isVerbose)
self = thetaTempMaker.thetaTemp(outputName=output_name, inputFile=inputfile, isTTree=isTTree, bin_type=bin_type, verbose = isVerbose)
self.main()
print 'All finished.'


