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
afb_mc_4j = '0.0900'
rqqbar_mc_4j = '0.087'
rbck_mc_4j = '0.375'
delta_mc_4j = '0'
afb_mc_5j = '-0.07'
rqqbar_mc_5j = '0.059'
rbck_mc_5j = '0.257'
delta_mc_5j = '0'
nevents = '2M'
# Declare global vars
pdf = 0
num_pdf = 53

Rqqbar_4j= num_pdf*[0.0]
Rbck_4j= num_pdf*[0.0]
delta_4j= num_pdf*[0.0]
Afb_4j= num_pdf*[0.0]
Rqqbar_5j= num_pdf*[0.0]
Rbck_5j= num_pdf*[0.0]
delta_5j= num_pdf*[0.0]
Afb_5j= num_pdf*[0.0]
Rqqbar_err_4j=Rbck_err_4j=delta_err_4j=Afb_err_4j=0
Rqqbar_err_5j=Rbck_err_5j=delta_err_5j=Afb_err_5j=0
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
	output.write('Number of Events: '+nevents+'\n\n')
	# Prepare for inputs into error function
	list_results = []
	list_para_names = []
	list_stat_err = []
	list_mc_truth = []

	list_results.append(Afb_4j)
	list_results.append(Rqqbar_4j)
	list_results.append(Rbck_4j)
	list_results.append(delta_4j)
	list_para_names.append('Afb_4j   ')
	list_para_names.append('Rqqbar_4j')
	list_para_names.append('Rbck_4j  ')
	list_para_names.append('delta_4j ')
	list_stat_err.append(Afb_err_4j)
	list_stat_err.append(Rqqbar_err_4j)
	list_stat_err.append(Rbck_err_4j)
	list_stat_err.append(delta_err_4j)
	list_mc_truth.append(afb_mc_4j)
	list_mc_truth.append(rqqbar_mc_4j)
	list_mc_truth.append(rbck_mc_4j)
	list_mc_truth.append(delta_mc_4j)

	list_results.append(Afb_5j)
	list_results.append(Rqqbar_5j)
	list_results.append(Rbck_5j)
	list_results.append(delta_5j)
	list_para_names.append('Afb_5j   ')
	list_para_names.append('Rqqbar_5j')
	list_para_names.append('Rbck_5j  ')
	list_para_names.append('delta_5j ')
	list_stat_err.append(Afb_err_5j)
	list_stat_err.append(Rqqbar_err_5j)
	list_stat_err.append(Rbck_err_5j)
	list_stat_err.append(delta_err_5j)
	list_mc_truth.append(afb_mc_5j)
	list_mc_truth.append(rqqbar_mc_5j)
	list_mc_truth.append(rbck_mc_5j)
	list_mc_truth.append(delta_mc_5j)

	# Run error functions
	for i in range(len(list_results)):
		get_errors(list_results[i],list_para_names[i],list_stat_err[i],list_mc_truth[i])
		if i == 3:
			output.write('\n')
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
	global num_pdf,Rqqbar_4j,Rbck_4j,delta_4j,Afb_4j,Rqqbar_5j,Rbck_5j,delta_5j,Afb_5j
	global Rqqbar_err_4j,Rbck_err_4j,delta_err_4j,Afb_err_4j,Rqqbar_err_5j,Rbck_err_5j,delta_err_5j,Afb_err_5j
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
			Rqqbar_4j[pdf] = numbers[1]
			Rbck_4j[pdf] = numbers[3]
			Afb_4j[pdf] = numbers[5]
			delta_4j[pdf] = numbers[15]

			Rqqbar_5j[pdf] = numbers[7]
			Rbck_5j[pdf] = numbers[9]
			Afb_5j[pdf] = numbers[11]
			delta_5j[pdf] = numbers[19]			

			# Find errors for only nominal PDF
			if pdf == 0:
				Rqqbar_err_4j = numbers[2]
				Rbck_err_4j = numbers[4]
				Afb_err_4j = numbers[6]
				delta_err_4j = numbers[16]
				Rqqbar_err_5j = numbers[8]
				Rbck_err_5j = numbers[10]
				Afb_err_5j = numbers[12]
				delta_err_5j = numbers[20]
				# Validation
				print str(Rqqbar_4j[pdf])+','+str(Rbck_4j[pdf])+','+str(delta_4j[pdf])+','+str(Afb_4j[pdf])
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



