import Legacy_Afb.theta_fit.samples as samples 
import ROOT
import glob
import sys
from array import array

import argparse

parser = argparse.ArgumentParser(description='merge MC templates to form a pseudo-data file')
parser.add_argument('-input',type=str, help='input files')
parser.add_argument('-output',type=str, default='test', help='input is smoothed histograms')
parser.add_argument('-MC_info',type=str, default='MC_input_with_bkg.txt', help='path of MC_info.txt')
parser.add_argument('-nevts', type=int, default=1000, help='nevts to run on each sample. default=1000')
parser.add_argument('-verbose', help='if verbose', action='store_true')

args = parser.parse_args()


input_files = glob.glob(args.input)
evt_end = int(args.nevts)
txtfile = args.MC_info


newtree_name = 'angles_data'
oldtree_name = 'angles'

QCD_SF = 0.06

all_weights = ['pileup_reweight','top_pT_reweight','btag_eff_reweight','tracking_reweight']
all_weights+= ['lepID_reweight','lepIso_reweight','trigger_reweight']
# set up output root file
fout = ROOT.TFile('%s_mc_data.root'%args.output,'recreate')
newtree = ROOT.TTree(newtree_name,newtree_name)
# Add new branches to the output tree
br_defs = []
# vars
ttbar_mass = array('f',[0.])
cos_theta_cs = array('f',[0.])
Feynman_x = array('f',[0.])
Q_l = array('i',[0])
total_w = array('f',[0.])
n_bTags = array('i',[-1])
sample_type = ROOT.vector('string')()
# names 
br_defs += [('ttbar_mass',ttbar_mass,'ttbar_mass/F')]
br_defs += [('cos_theta_cs',cos_theta_cs,'cos_theta_cs/F')]
br_defs += [('Feynman_x',Feynman_x,'Feynman_x/F')]
br_defs += [('Q_l',Q_l,'Q_l/I')]
br_defs += [('total_w',total_w,'total_w/F')]
br_defs += [('n_bTags',n_bTags,'n_bTags/I')]

# Add branches to the tree
for ibr in br_defs:
    newtree.Branch(ibr[0],ibr[1],ibr[2])
newtree.Branch('sample_type',sample_type)

# def mulitplication of all elements in a list
multiply = lambda l:reduce(lambda x, y: x*y, l)
# loop over input root file
samples_obj = samples.samples(txtfile)
for ifile in input_files:
	sample_info_obj = samples_obj.get_sample_info(ifile)
	if not sample_info_obj: continue
	print '(info) Processing...  %s'%sample_info_obj.print_out()
	root_file = ROOT.TFile(ifile)
	tmptree = root_file.Get(oldtree_name)
	norm_weight = sample_info_obj.weight
	sample_type.clear()
	sample_type.push_back(sample_info_obj.type)
	# QCD need another SF to match with data
	if 'QCD' in sample_info_obj.key:
		norm_weight *= QCD_SF
        # check if gen_w is in the tree
        if tmptree.FindBranch('gen_weight'):
	        has_genW = True
		print '(info) Has_genW'
        else:
        	has_genW = False
	# loop over entries and fill
	n_evt = 0
	for iev in xrange(0,tmptree.GetEntries()):
		# Progress report
		if iev%5000 == 0 : print 'processing event',iev 
		# Break at the given event number
		if iev == evt_end : 
			print 'Finish processing. Quit at event',evt_end
			break 
		tmptree.GetEntry(iev)
		ttbar_mass[0] = tmptree.ttbar_mass
		cos_theta_cs[0] = tmptree.cos_theta_cs
		Feynman_x[0] = tmptree.Feynman_x
		Q_l[0] = tmptree.Q_l
		n_bTags[0] = tmptree.n_bTags
		# Calculate total weight
		total_weight = [getattr(tmptree,item) for item in all_weights]
		total_weight.append(norm_weight)
		total_w[0] = multiply(total_weight)
                # load genW
                if has_genW:
	                gen_weight = tmptree.gen_weight
                else:
                	gen_weight = 1.0
		total_w[0] = total_w[0]*gen_weight
		# Fill current entry into ttree
		newtree.Fill()
	root_file.Close()
fout.Write()
fout.Close()
# closing
