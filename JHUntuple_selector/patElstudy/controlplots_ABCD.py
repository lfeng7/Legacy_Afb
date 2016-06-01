# make control plots for ABCD study

import glob
import sys
import re
# ./masterplot.py --file SingleEl_Run2012A_all_v4.root --varx lep_iso --vary met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 1.2 --Miny 0 --Maxy 150 --xaxis "RelIso" --yaxis "MET" --title "Data" --plotname iso_met_data

argv = sys.argv[1:]
if len(argv)< 1 or ('help' in argv):
	print """
	Usage:
	"""
	sys.exit(1)
inputdir = argv.pop(0)

# cuts
loose_cut = 'lep_isLoose'
signal_cut = 'njets>=4&&n_btags==2&&lep_isTight'
norm = 'weight_norm'
# samples to plot
samples = ['DY.*Jets.*','W.*Jets.*','TT_.*','T(bar)?_.*','QCD.*','SingleEl.*']
sample_name = ['DYJets','WJets','TTbar','SingleTop','QCD','Data']
all_files = glob.glob('%s/*.root'%inputdir)

sample_table = {}
# make plotting cmds
fout = open('plot_controlABCD.sh','w')
for i,item in enumerate(sample_name):
	ipattern = samples[i]
	iname = sample_name[i]
	ifiles = [item for item in all_files if re.match('.*/%s'%ipattern,item)]
	ifiles = ' '.join(ifiles)
	# actual plots to make
	# iso vs MET , signal cuts
	ititle = '%s %s'%(iname,signal_cut)
	icut = ' "(%s)*(%s)" '%(signal_cut,norm)
	icmd = './masterplot.py --file %s --title "%s" --plotname iso_met_signal_cut%s --cut %s '%(ifiles,ititle,iname,icut)
	icmd += ' --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" \n\n' 
	# iso vs MET , very loose cuts
	ititle = '%s %s'%(iname,loose_cut)
	icut = ' "(%s)*(%s)" '%(loose_cut,norm)
	icmd += './masterplot.py --file %s --title "%s" --plotname iso_met_loose_cut%s --cut %s '%(ifiles,ititle,iname,icut)
	icmd += ' --vary lep_iso --varx met_pt_vec --binx 20 --biny 20 --Minx 0 --Maxx 150 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "MET" \n\n' 
	# iso vs impact parameters
	ititle = '%s %s'%(iname,loose_cut)
	icut = ' "(%s)*(%s)" '%(loose_cut,norm)
	icmd += './masterplot.py --file %s --title "%s" --plotname iso_transIP_loose_cut%s --cut %s '%(ifiles,ititle,iname,icut)
	icmd += ' --vary lep_iso --varx TransverseIP --binx 20 --biny 20 --Minx 0 --Maxx 0.1 --Miny 0 --Maxy 1.2 --yaxis "RelIso" --xaxis "TransverseIP" \n\n' 


	fout.write(icmd)




