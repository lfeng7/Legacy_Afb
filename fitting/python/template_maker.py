import ROOT
import copy
from array import array
import sys
import numpy

import template
import helper



class tempMaker(object):
    """
    take in a ttree in legacy format
    make TH3D templates that is compatible with legacy fitter
    """
    def __init__(self, inputFile='aggregated_distributions.root',outputname='template',verbose=False,nevts_run=-1):
        super(tempMaker, self).__init__()
        self.inputFile = inputFile
        self.outputname = outputname
        self.verbose=verbose
        self.nevts_run=nevts_run
        self.versions = ['plus','minus'] # Ql=+/- version of templates
        self.binType = 'variable'
        self.all_templates = {} # a hash table with key=template_name, value=template.object

    def main(self):
        """
        Main function for the class
        """
        self.initOutput()
        self.initTemplates()
        self.makeTemplates(self.inputFile)
        self.writeTemplates()
        print '(info) Finish main!'     

    def initOutput(self):
        """
        create new root files for templates and aux hists
        """
        self.outfile = ROOT.TFile('%s.root'%self.outputname,'recreate')
        self.outfile_aux = ROOT.TFile('%s_aux.root'%self.outputname,'recreate')


    def initTemplates(self):
        """
        Intialize all templates, create a template object for each template
        Then organize these templates in the form of {'dist_type1':[[temp_obj1,w1],[temp_obj2,w2],...]],'dist_type2':[],...}
        """
        self.processes = []
        self.processes.extend([['WJets',1],['bck',1],['gg',1],['ntmj',1],['qqa','wa'],['qqs',1]])
        self.processes.extend([['qqa_delta','wadelta'],['qqa_xi','waxi'],['qqs_delta','wsdelta'],['qqs_xi','wsxi']])
        prefix = 'f'
        all_templates = {}
        for iversion in self.versions:
            for iprocess in self.processes:
                #   def __init__(self,name,formatted_name,bin_type='fixed') :
                iprocess_name = iprocess[0]
                iname = '%s_%s_%s'%(prefix,iprocess_name,iversion)
                itemp = template.template(iname,iname,self.binType)
                all_templates[iname]=[itemp,iprocess[1]]
        # organize templates in the form of process_plus/minus for quick look up and filling
        for key,value in all_templates.iteritems():
            # all keys are in the form of 'f_qqs_delta_minus'
            temp_keys = key.split('_')
            new_key = '%s_%s'%(temp_keys[1],temp_keys[-1])
            if self.all_templates.get(new_key,0)==0:
                self.all_templates[new_key]=[value]
            else:
                self.all_templates[new_key].append(value)
        print '(info) Finished initTemplates!'

    def writeTemplates(self):
        """
        write all templates in self.all_templates into output file
        """
        for key,value in self.all_templates.iteritems():
            for item in value:
                itemp = item[0]
                self.outfile.cd()
                itemp.histo_3D.Write()
                self.outfile_aux.cd()
                itemp.histo_x.Write()
                itemp.histo_y.Write()
                itemp.histo_z.Write()
                print '(info) Writing %s.'%itemp.name

    def makeTemplates(self,tree_file):
        """
        input: a single tree file with all processes , list of weight
        ( a list of 1D arrays, for the purpose of sys templates )
        output: a 1D template and control plots of 3d projections
        """

        # loading ttree and load branches
        tfile = ROOT.TFile(tree_file)
        tree_name = helper.GetTTreeName(tfile)
        ttree = tfile.Get(tree_name)

        # Loop over entries in ttree and fill templates
        n_entries = ttree.GetEntries()
        for iev in range(n_entries):
            # progress report
            if iev%5000==1:
                print '(info) Processing %i evt.'%iev
            if iev==self.nevts_run:
                print '(info) Reaching %i evts. Will quit here.'%iev
                break

            ttree.GetEntry(iev)
            # load observables
            cs = ttree.costheta
            xf = ttree.xF
            mtt = ttree.tt_M
            lep_charge = ttree.Ql
            total_weight = ttree.weight*ttree.normalization_weight
            dist_type = ttree.dist_type
            # first assign process type, e.g., qqs_plus,qqa_plus
            if lep_charge>0:
                # is_for_distribution = (1,qq) (2,gg) (3,bkg) (4,WJets) (-1,sb)
                if dist_type==-1:
                    proc_type = ['ntmj_plus']
                elif dist_type==1:
                    proc_type = ['qqs_plus','qqa_plus']
                elif dist_type==2:
                    proc_type = ['gg_plus']
                elif dist_type==3:
                    proc_type = ['bck_plus']
                elif dist_type==4:
                    proc_type = ['WJets_plus']
            elif lep_charge<0:
                # is_for_distribution = (1,qq) (2,gg) (3,bkg) (4,WJets) (-1,sb)
                if dist_type==-1:
                    proc_type = ['ntmj_minus']
                elif dist_type==1:
                    proc_type = ['qqs_minus','qqa_minus']
                elif dist_type==2:
                    proc_type = ['gg_minus']
                elif dist_type==3:
                    proc_type = ['bck_minus']
                elif dist_type==4:
                    proc_type = ['WJets_minus']
            else:
                print '(error) lep_charge = 0. This event is corrupted! Will skip this.'
                continue
            # then fill the corresponding templates
            # proc_type = ['qqs_minus','qqa_minus']
            for itype in proc_type:
                itemplates = self.all_templates.get(itype,0)
                # itemplates = [['qqs_xi_minus',wxi],['qqs_minus','ws'],['qqs_delta_minus','wsdelta']]
                if itemplates==0:
                    print '(error) Cannot find current templates %s.'%itype
                    sys.exit(1)
                else:
                    for itemp in itemplates:
                        # itemp=[template_object,'wsxi']
                        itemp_obj = itemp[0]
                        iweight = itemp[1]
                        if iweight!=1:
                            additional_weight = getattr(ttree,iweight)
                        else:
                            additional_weight=1
                        new_w = total_weight*additional_weight
                        itemp_obj.Fill(cs,abs(xf),mtt,new_w)
                        if self.verbose and iev%5000==1:
                            print '(debug) Filling into %s with weight %s.'%(itemp_obj.name,iweight)
        print '(info) Finish loop over all entries!'
        tfile.Close()

    def __del__(self):
        self.outfile.Close()
        self.outfile_aux.Close()
        print '(info) Closeup thetaTemp object.'

if __name__=='__main__':
    argv = sys.argv[1:]
    # help message
    if len(argv)==0:
        print """
        Usage: 

        """
        sys.exit(1)
    infile = argv.pop(0)
    outname = argv.pop(0)
    nevts_run = int(argv.pop(0))
    if 'verbose' in  argv:
        verbose=True
    else:
        verbose=False
    # def __init__(self, inputFile='aggregated_distributions.root',outputname='template',verbose=False,nevts_run=-1):
    self = tempMaker(inputFile=infile,outputname=outname,verbose=verbose,nevts_run=nevts_run)
    self.main()