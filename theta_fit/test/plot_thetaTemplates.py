"""
main code that produce templates for theta
"""
import Legacy_Afb.theta_fit.tempPlotter as tempPlotter 
import sys

argv=sys.argv[1:]
if len(argv)==0:
	print """
	Usage: 
	python python plot_thetaTemplates.py temp_angles/angles_data_SingleEl_Run2012ABCD.root  (verbose)
	"""
	sys.exit(1)

inputfile = argv.pop(0)
argv = ' '.join(argv)

is_verbose = False
bin_type = 'variable'
if 'verbose' in argv:
	is_verbose = True
if 'fixed' in argv:
	bin_type = 'fixed'
	

self = tempPlotter.plotter(template_file=inputfile, verbose=is_verbose, bin_type=bin_type)
self.main()
print '(MAIN) All finished.'


