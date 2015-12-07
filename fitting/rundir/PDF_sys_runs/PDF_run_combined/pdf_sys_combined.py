################### PDF uncertainty calculator #################
## Edited by Raymond
## Edited on 4/16/14
## A better version!

#!/usr/bin/python
import re
import os
import math
import ROOT
import numpy
from ROOT import *
from array import array
from optparse import OptionParser
# MC truth of fitted vars
afb_mc = '0.036'
rqqbar_mc = '0.075'
rbck_mc = '0.33'
delta_mc = '0'
nevents = '2M'
# Declare global vars
num_pdf = 53 
Rqqbar= num_pdf*[0.0]
Rbck= num_pdf*[0.0]
delta= num_pdf*[0.0]
Afb= num_pdf*[0.0]
pdf = 0
Rqqbar_err=Rbck_err=delta_err=Afb_err=0
# rqq4_0=rqq5_0=rbck4_0=rbck5_0=del4_0=del5_0=xi4_0=xi5_0=afb_0=None
output = open('PDF_uncertainty.txt','w')

# COMMAND LINE OPTIONS
parser = OptionParser()
parser.add_option('--plot', metavar='F', type='string', action='store',
	default='no',
	dest='make_plots',
                	help='') ## Sets which files to run on
(options, args) = parser.parse_args()


def main():
	# global num_pdf,Rqqbar_4jet,Rbck_4jet,Rqqbar_5jet,Rbck_5jet,delta_4jet,delta_5jet,Afb
	# Parse fit results and write into files
	parsing()
	# Calculate PDF systematics
	output.write('Number of Events: '+nevents+'\n')
	# Prepare for inputs into error function
	list_results = []
	list_results.append(Afb)
	list_results.append(Rqqbar)
	list_results.append(Rbck)
	list_results.append(delta)
	list_para_names = []
	list_para_names.append('Afb   ')
	list_para_names.append('Rqqbar')
	list_para_names.append('Rbck  ')
	list_para_names.append('delta ')
	list_stat_err = []
	list_stat_err.append(Afb_err)
	list_stat_err.append(Rqqbar_err)
	list_stat_err.append(Rbck_err)
	list_stat_err.append(delta_err)
	list_mc_truth = []
	list_mc_truth.append(afb_mc)
	list_mc_truth.append(rqqbar_mc)
	list_mc_truth.append(rbck_mc)
	list_mc_truth.append(delta_mc)
	# Run error functions
	for i in range(len(list_results)):
		get_errors(list_results[i],list_para_names[i],list_stat_err[i],list_mc_truth[i])
	# Make plots of vars
	if options.make_plots == 'yes':
		print 'Making plots...'
		for i  in range(len(list_results)):
			plot(list_results[i],list_para_names[i])
	else:
		print 'No plots are made!'
	output.close()

########################################### Parse the output txt and write into lists of vars ############################
def parsing():
	global num_pdf,Rqqbar,Rbck,delta,Afb
	global Rqqbar_err,Rbck_err,delta_err,Afb_err
	# Open and read the list of result.txt
	f = open("result_list.txt")
	result_list = list(f)
	result_list = [word.strip() for word in result_list]	#word.strip() get rid of \n in this line
	f.close()
	# Loop over result.txt
	for fname in result_list:
		# Read in results
		f = open(fname)
		result = list(f)
		result = [word.strip() for word in result]
		# Find the pdf_index
		pdf_index1 = re.search(r'^\D*(\d+)\D*',fname)
		if pdf_index1 :
			pdf = int(pdf_index1.group(1))
			# Find the numbers
			numbers = result[0].split(',')
			Rqqbar[pdf] = numbers[1]
			Rbck[pdf] = numbers[3]
			delta[pdf] = numbers[7]
			Afb[pdf] = numbers[9]
			# Find errors for only nominal PDF
			if pdf == 0:
				Rqqbar_err = numbers[2]
				Rbck_err = numbers[4]
				delta_err = numbers[8]
				Afb_err = numbers[10]
				# Validation
				print str(Rqqbar[pdf])+','+str(Rbck[pdf])+','+str(delta[pdf])+','+str(Afb[pdf])
		else :
			print 'wrong file name! will stop.'
			break
		f.close()

def get_errors(var,var_name,stat_err,mc_truth):	# var_name here should be the string of nominal var results,e.g. rqq4_0
	global output
	mid = float(var[0])
	sum_up = 0.
	sum_down = 0.
	for i in range (1,num_pdf):
		if i%2 == 1:
			# print str(i)
			up = float(var[i])
			down = float(var[i+1])
			temp_up = max_up(up,down,mid)
			temp_down = max_down(up,down,mid)
			sum_up += pow(temp_up,2)
			sum_down += pow(temp_down,2)
			error_up   = math.sqrt(sum_up)
			error_down = math.sqrt(sum_down)
		#Print out result and write into files
	error_up = format(error_up,'.4f')
	error_down = format(error_down,'.4f')
	printout = var_name+' = '+str(var[0])+' +/- '+str(stat_err)+' (statistics) +/- '+ str(error_up)+'/' +str(error_down)+' (PDF systematics)'+' | MC truth '+str(mc_truth)
	#convert the array into float type
	var_float = []
	for j in range(len(var)):
		var_float.append(float(var[j]))
	printout += ' stdev '+format((numpy.std(var_float)),'.4f')

	output.write(printout+'\n')


def max_up(up,down,mid):
	return max(up-mid,down-mid,0)

def max_down(up,down,mid):
	return max(mid-up,mid-down,0)

def plot(var,var_name):
	fig_name = var_name+'.eps'
	dim = len(var)
	# Make x and y array for TGraph
	y = array('f',dim*[0.0])
	index = array('f',dim*[0])
	for i in range(dim):
		index[i] = i
		y[i] = float(var[i])
	# Make a reference formula
	nominal = TF1('func1',var[0],0,dim)
	# Set canvas and graph
	c1 = TCanvas('c1','var vs pdf_member',200,10,700,500)
	c1.SetGrid()
	gr1 = TGraph(dim,index,y)
	gr1.SetMarkerStyle( 21 )
	gr1.SetTitle(var_name)
	gr1.GetXaxis().SetTitle('pdf index')
	gr1.GetYaxis().SetTitle(var_name)
	gr1.Draw('ACP')
	# Draw reference line
	nominal.Draw('same')
	c1.Update()
	# Save to eps
	c1.SaveAs(fig_name)
	# #Set up output file
	# outputname = "fit_result.root"
	# f = TFile( outputname, "Recreate" )  
	# f.cd()
	# c1.Write()
	# f.Write()
	# f.Close()  

main()



