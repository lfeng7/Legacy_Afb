# input: regular root file, with both MC and data
# output: new root file, with cloning of the original , plus an additional branch called w_norm, that contains the correct normalization information for MC, multiply by -1 , and 1 for data

import ROOT as rt
import samples
import sys
import argparse
from array import array
import glob

parser = argparse.ArgumentParser(description='merge MC templates to form a pseudo-data file')
parser.add_argument('-input',type=str, help='input files')
parser.add_argument('-MC_info',type=str, default='MC_input_with_bkg.txt', help='path of MC_info.txt')
parser.add_argument('-nevts', type=int, default=1000, help='nevts to run on each sample. default=1000')
parser.add_argument('-verbose', help='if verbose', action='store_true')

args = parser.parse_args()

input_files = glob.glob(args.input)
txtfile = args.MC_info
samples_obj = samples.samples(txtfile)

for ifile in input_files:
    sample_info_obj = samples_obj.get_sample_info(ifile)
    try:
        norm_w = sample_info_obj.weight
    except AttributeError:
        print '(info) %s is Data file!'%ifile
        isData = True
    else:
        isData = False
    
    old_file = rt.TFile(ifile)
    old_tree = old_file.Get('selected')
    old_fname = ifile.split('/')[-1].split('.root')[0]
    new_file = rt.TFile('%s_sideband.root'%old_fname,'recreate')
    newtree = old_tree.CloneTree(0)
    # add new branches
    w_norm = array('f',[1.0])
    newtree.Branch('w_norm',w_norm,'w_norm/F')
    # for data only add some w_corr branches
    new_brs_to_add = []
    new_brs_to_add += [('weight_top_pT','f'),('w_PU','f'),('w_PDF','f'),('w_btag','f')\
                       ,('w_lepID','f'),('w_trigger','f'),('w_tracking','f'),('w_lepIso','f') ]
    # add non-vector new branches
    added_brs = {}
    br_defs = []
    for item in new_brs_to_add:
        name,br_type = item
        if br_type == 'f':
            tmp_array = array('f',[1.])
            br_defs+=[(name,tmp_array,'%s/F'%name)]
            added_brs[name] = tmp_array            
        elif br_type == 'i':
            tmp_array = array('i',[-10])
            br_defs+=[(name,tmp_array,'%s/I'%name)]
            added_brs[name] = tmp_array

    if isData:
        # Add float branch into ttree
        print '(info) Adding these branches:',[item[2] for item in br_defs]
        for ibr in br_defs:
            newtree.Branch(ibr[0],ibr[1],ibr[2])

    # loop over evts and fill
    for iev in range(old_tree.GetEntries()):
        if iev == args.nevts: 
            print 'Stop at %i evt'%iev
            break
        if isData:
            w_norm[0] = 1.0
#            for ikey in added_brs:
#                added_brs[ikey][0]=1.0
        else:
            w_norm[0] = -1.0*sample_info_obj.weight
        old_tree.GetEntry(iev)
        newtree.Fill()
    new_file.Write()
    new_file.Close()

    
    
    

