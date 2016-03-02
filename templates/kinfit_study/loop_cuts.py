#! /usr/bin/python 
import sys
import os

argv = sys.argv[1:]

if len(argv)<2:
	print """
	usage of the code
	./loop_cuts.py  chi2_new 10 30 5 output.txt
	chi2_new < from 10 to 30 every 5 and output to output.txt
	if not output is given, use 
	"""
	sys.exit(1)

cut_name = argv.pop(0)
cut_low = float(argv.pop(0))
cut_high = float(argv.pop(0))
cut_step = float(argv.pop(0))
if len(argv) == 1:
	output_file = argv.pop(0)
else :
	output_file = ''

cut_value = cut_low
while cut_value < cut_high :
    cmd = 'python efficiency_calculator.py --dir template_files/pscale_removed_v2 --cut \'%s < %.2f\' '%(cut_name,cut_value) 
    if output_file != '':
    	cmd += ' --output %s'%output_file
    os.system(cmd)
    print cmd
    cut_value += cut_step
