"""
main code that produce templates for theta
"""
from python import thetaTempMaker
import os

argv=os.argv[1:]
if len(argv)==0:
	print """
	Usage: 
	"""
inputfile = argv.pop(0)
output_name = argv.pop(0)

print 'create a thetaTemp object.'
self = thetaTemp.template(outputName=output_name, inputFile=inputfile)
self.main()
print 'All finished.'


