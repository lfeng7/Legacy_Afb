#! /usr/bin/python
# This is the master program that call different plotting codes and make various kind of plots

import sys
import os

dir_plotting_code = '/uscms_data/d3/lfeng7/B2G_FW/CMSSW_7_2_0/src/Legacy_Afb/plotter/'

argv = sys.argv[1:]

if len(argv)< 1 or ('help' in argv):
	print """
	Making various kind of plots. 
	Have to modify the hard coded dir_plotting_code to the path where all plotting codes locate
	Example usage:
	(a) Make simple 1D plot
	./masterplot.py file TT_CT10_TuneZ2star_8TeV_reco_angles_0.root, var mtt, bin 50, min 350, max 1500, plotname Mtt_a , xaxis "M_{tt}(GeV)", cut "final_chi2<25" , title "example (a)", saveroot yes
	
	(b) Make simple 1D plot with two different cuts and scaled
	./masterplot.py file TT_CT10_TuneZ2star_8TeV_reco_angles_0.root, var mtt, bin 50, min 350, max 1500, plotname Mtt_b , xaxis "M_{tt}(GeV)", cut "final_chi2<25" , title "example (b)", saveroot yes , cut2 "final_chi2>25" ,scale yes
	
	(c) Make a plot comparing the same var in several different files
	./masterplot.py file TT_CT10_TuneZ2star_8TeV_reco_angles_0.root  W4JetsToLNu_TuneZ2Star_8TeV_reco_angles_0.root T_t_reco_angles_0.root , var mtt, bin 50, min 350, max 1500, xaxis "M_{tt}(GeV)" , title "example (c)" , plotname Mtt_c ,scale yes , label "TT W4Jets T_t"	
	
	(d) A 2D plot of a single file with var1 vs var2
	./masterplot.py --file TT_CT10_TuneZ2star_8TeV_reco_angles_0.root --varx mtt --vary mtt_mc --binx 30 --biny 30 --Minx 350 --Maxx 1000 --Miny 350 --Maxy 1000 --xaxis "M_{tt}(reco)(GeV)" --yaxis "M_{tt}(mc)(GeV)" --title "example plot (d)" --plotname Mtt_d

	(e) Compare two vars from same root file in a 1D plot
	masterplot.py file TT_CT10_TuneZ2star_8TeV_reco_angles_0.root , var1 mtt, var2 mtt_mc , bin 50, min 350 , max 1500, xaxis "M_{tt}(GeV)", title "example (e)", scale yes ,plotname Mtt_e, label1 "mtt", label2 "mtt_mc"

	(f) Compare a var from different type of events, e.g., signal, bkg ..  using a txt file with cross section and Nevts_gen, and a dir with input files to plot
	./masterplot.py dir . , var1 mtt , bin 50, min 350 , max 1500, xaxis "M_{tt}(GeV)", title "example (f)" ,plotname Mtt_f, scale yes

	(g) Make a stack plot of a var from different type of events, e.g., signal, bkg ..  using a txt file with cross section and Nevts_gen, and a dir with input files to plot
	./masterplot.py dir signal , var1 mtt , bin 50, min 350 , max 1500, xaxis "M_{tt}(GeV)", title "example (g)" ,plotname Mtt_g, stack yes
	"""
	sys.exit(1)

argv = ' '.join(argv)

if '--' in argv:
	argv = argv.split('--')
elif ',' in argv:
	argv = argv.split(',')
else:
	print 'arguments should be seperated by either -- or ,'
	sys.exit(1)
# python ../multi_sample_plotter_v2.py --var cos_theta_cs --Min -1 --Max 1 --name cs_sideband --bin 30 --dir template_files/sideband/ --title "cos_theta_cs sideband" --xaxis "cos#theta*" --yaxis "events" --weight correction_weight

# initiate the values of arguments
file1 = ""
file2 = ""
var1 = ""
var2 = ""
Minx = ""
Maxx = ""
Miny = ""
Maxy = ""
plotname = ""
binx = ""
biny = ""
cut1 = ""
cut2 = ""
log = ""
scale = ""
title = ""
xaxis = ""
yaxis = ""
label1 = ""
label2 = ""
saveroot = ""
stack = ""
dir_path = ""
plot = ""
weight = ""
yields = ""
mode = ""
verbose = ''

# find the input arguments
for i,item in enumerate(argv):
	item = item.split()
	if len(item)<2:
		continue
	word = item.pop(0).lower( ) 
	if word in ['file1','file','files']: file1 = '\"%s\"'%' '.join(item)
	if word in ['file2']: file2 = item[0]
	if word in ['var1','var','varx']: var1 = '\"%s\"'%' '.join(item)
	if word in ['var2','vary']: var2 = '\"%s\"'%' '.join(item)
	if word in ['min','minx']: Minx = item[0]
	if word in ['max','maxx']: Maxx = item[0]
	if word in ['miny']: Miny = item[0]
	if word in ['maxy']: Maxy = item[0]	
	if word in ['name','fname','plotname']: plotname = item[0]
	if word in ['bin','binx']: binx = item[0]
	if word in ['bin','biny']: biny = item[0]
	if word in ['dir']: dir_path = item[0]
	if word in ['log']: log = item[0]
	if word in ['scale']: scale = item[0]
	if word in ['xaxis']: xaxis = ' '.join(item)
	if word in ['yaxis']: yaxis = ' '.join(item)
	if word in ['label','label1']: label1 = ' '.join(item)
	if word in ['label2']: label2 = ' '.join(item)
	if word in ['verbose']: verbose = item[0]
	if word in ['cut','cut1']: cut1 = ' '.join(item)
	if word in ['cut2']: cut2 = ' '.join(item)	
	if word in ['title']: title = ' '.join(item)
	if word in ['saveroot','root']: saveroot = item[0]
	if word in ['stack']: stack = item[0]
	if word in ['plot']: plot = item[0]
	if word in ['weight']: weight = item[0]
	if word in ['yields']: yields = item[0]
	if word in ['mode']: mode = item[0]
	if word in ['mcinfo']: MCinfo = item[0]

# Determine what kind of plot to make
plotting_code = ''
if file1 != "":
	if file2 == "": # single file to plot
		if var2 == "": # single var to plot
			if len(file1.split())>1 or '*' in file1:	# more than one file to plot 
				if stack == 'yes':
					plotting_code = 'simple_stacker.py'
				else:		
					plotting_code = 'multi_sample_plotter.py'
			else:
				plotting_code = 'plotter.py'
		else:
			plotting_code = 'double_var_plotter.py'
		if Miny != "": # 2D plots
			plotting_code = '2D_plotter.py'
	else: # double files
		plotting_code = 'double_sample_plotter.py'
if dir_path != "":
	if stack == 'yes':
		plotting_code = 'stack_plotter.py'
	else:
		plotting_code = 'multi_sample_plotter_v2.py'
# Overwrite plotter type by plotting mode parameter
if mode == 'plotter':
	plotting_code = 'plotter.py'

# Make actual plotting cmd
cmd = 'python %s/%s '%(dir_plotting_code,plotting_code)

if plotname!="": cmd+=' --name %s'%plotname
if log!="": cmd+=' --log %s'%log
if scale!="": cmd+=' --scale %s'%scale

if binx!="": 
	if plotting_code == '2D_plotter.py':
		cmd+=' --binx %s'%binx
	else:
		cmd+=' --bin %s'%binx

if title!="": cmd+=' --title \"%s\"'%title
if xaxis!="": cmd+=' --xaxis \"%s\"'%xaxis
if yaxis!="": cmd+=' --yaxis \"%s\"'%yaxis
if plot!="": cmd+=' --plot %s'%plot
if cut1!="": cmd+=' --cut \"%s\"'%cut1
if cut2!="": cmd+=' --cut2 \"%s\"'%cut2
if saveroot!="": cmd+=' --save %s'%saveroot
if label1!="": cmd+=' --label \"%s\"'%label1
if label2!="": cmd+=' --label2 \"%s\"'%label2
if weight!="": cmd+=' --weight \"%s\"'%weight
if yields!="": cmd+=' --yields %s'%yields
if verbose!="": cmd+=' --verbose '


# check if some key arguments are provided
if Minx == "":
	print 'x Min is not provided. Please type in now.'
	Minx = raw_input()
if Maxx == "":
	print 'x Max is not provided. Please type in now.'
	Maxx = raw_input()
if file1 == "" and dir_path == "":
	print 'input file is not provided. Please type in and run again!'
	sys.exit(1)
if var1 == "":
	print 'plotting var is not provided. Please type in now.'
	var1 = raw_input()
if plotname == "":
	print 'plotname var is not provided. Please type in now.'
	plotname = raw_input()

# Make plots
if plotting_code == 'plotter.py':
	cmd += ' --Min %s --Max %s --file %s --var %s'%(Minx,Maxx,file1,var1)

if plotting_code == 'double_sample_plotter.py':
	cmd += ' --Min %s --Max %s --file1 %s --file2 %s --var %s'%(Minx,Maxx,file1,file2,var1)

if plotting_code == 'multi_sample_plotter.py':
	cmd += ' --Min %s --Max %s --files %s --var %s'%(Minx,Maxx,file1,var1)

if plotting_code == 'simple_stacker.py':
	cmd += ' --Min %s --Max %s --files %s --var %s'%(Minx,Maxx,file1,var1)

if plotting_code == 'multi_sample_plotter_v2.py':
	cmd += ' --Min %s --Max %s --dir %s --var %s'%(Minx,Maxx,dir_path,var1)

if plotting_code == 'double_sample_plotter.py':
	cmd += ' --Min %s --Max %s --file1 %s --file2 %s --var %s'%(Minx,Maxx,file1,file2,var1)

if plotting_code == 'double_var_plotter.py':
	cmd += ' --Min %s --Max %s --file %s --var1 %s --var2 %s'%(Minx,Maxx,file1,var1,var2)

if plotting_code == 'stack_plotter.py':
	cmd += ' --Min %s --Max %s --dir %s  --var %s --MCinfo %s'%(Minx,Maxx,dir_path,var1,MCinfo)

if plotting_code == '2D_plotter.py':
	cmd += ' --Minx %s --Maxx %s --Miny %s --Maxy %s --biny %s --file %s --varx %s --vary %s '%(Minx,Maxx,Miny,Maxy,biny,file1,var1,var2)

# decide if we want to print out the cmd only

if ('printout' in ''.join(argv)) or ('debug' in ''.join(argv)):
	print '\n%s\n'%cmd
else:
	print '\ncall: python %s/%s\n '%(dir_plotting_code,plotting_code)
	os.system(cmd)








