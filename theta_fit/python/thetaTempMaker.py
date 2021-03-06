# 
import ROOT
import template
import helper
import samples
import copy
from array import array
import sys
import numpy as np

ROOT.gROOT.SetBatch(True)
rate_sigma=0.8

class thetaTemp(object):
    """
    Take a root file w/ ttree (for Data) or smoothed 3D templates (for MC) 
    and create a single root file with all templates 
    """
    def __init__(self, outputName, inputFile, lep_type, isTTree=False, txtfile=None, use_MC_DATA=False, \
        top_reweight=True,JEC_sys='nominal', use_sys = True, verbose = False, \
        bin_type = 'fixed', nevts = -1, lep_combo = False):
        """
        inputFile is a list of path of input template root files (for ttree) or a single file for MC
        outputName is the name for output thetaTemp.root file
        """
        super(thetaTemp, self).__init__()
        self.nevts = nevts
        self.process = {}
        self.outputName = outputName
        self.outfile = ROOT.TFile('templates/%s_template.root'%outputName,'recreate')
        self.thetaHistList = []
        self.template_file = inputFile
        self.outfile_aux = ROOT.TFile('templates/%s_template_aux.root'%self.outputName,'recreate')
        self.outfile_aux.mkdir('hists/')
        self.outfile_aux.mkdir('plots/')
        self.isTTree = isTTree
        self.QCD_SF = 0.06 # this number is from the ABCD method estimation
        self.btag_cut = 2
        self.verbose = verbose
        self.bin_type = bin_type
        self.use_MC_DATA =use_MC_DATA
        self.top_reweight = top_reweight
        self.lep_type = lep_type
        self.use_sys = use_sys
        self.JEC_sys = JEC_sys # nominal,JES__minus,JES__plus, JER__minus,JER__plus
        self.AFB_sigma=1.0 # one sigma deviation of AFB from zero
        self.lep_combo = lep_combo
        if txtfile is None:
            txtfile = 'MC_input_with_bkg.txt'
        self.samples_obj = samples.samples(txtfile,self.verbose)
        # Added features
        # rate_change of each process in terms of 1sigma deviation of SF from 1.
        self.rate_change = (('WJets',rate_sigma),('other_bkg',rate_sigma),('qq',rate_sigma))
        ###### very special ######
        ROOT.TH1.SetDefaultSumw2(True)

    def main(self):
        self.define_process()
        self.define_sys()
        if self.isTTree:
            print '\n(info) Process template file w/ TTree.'
            self.assignFiles()
            self.GetWeights()
            # Loop over all processes, e.g., wjets, gg etc
            for key,value in self.process_files.iteritems():
                print '###############################################################'
                print '\n(info) Begin process %s samples.'%key
                self.TTtreeToTemplates(ttree_file=value, process_name=key)
            if self.JEC_sys == 'nominal':
                self.get_rate_templates()
        else:
            print '\n(info) Process template file w/ TH3D. '
            self.importTemp(self.template_file)
        self.makeControlPlots()
        print '(info) Done thetaTemp.main()'

    def assignFiles(self):
        """
        Assign the files into a hashtable of processes
        Input: a list of root files as self.template_file
        output: a hashtable of files belong to which process
        """
        if self.use_MC_DATA:
            self.Data_name = 'mc_data'
        else:
            self.Data_name = 'Run'

        process_files = {}
        for ifile in self.template_file:
            # if it's data , hard code it
            if self.Data_name in ifile :
                if process_files.get('DATA',1):
                    process_files['DATA'] = [ifile]
                else:
                    process_files['DATA'].append(ifile)
            # if not, it's an MC and will look up it's process type and add into hashtable  
            else:                   
                # load info about samples from txt file
                sample_info_obj = self.samples_obj.get_sample_info(ifile)
                # if cannot find info just skip this file for good
                if sample_info_obj is None: continue
                sample_type = sample_info_obj.type
                # Add file path to dictionary
                if process_files.get(sample_type,0)==0:
                    process_files[sample_type] = [ifile]
        #           if self.verbose: print 'Add type %s and file %s'%(sample_type,ifile)
                else:
                    process_files[sample_type].append(ifile)
        #           if self.verbose: print 'Add sample %s'%ifile
        self.process_files = process_files
        if self.verbose:
            print '(DEBUG) %s'%process_files
        print '(info) Done assignFiles.'


    def define_process(self):
        """
        define template processes in a dict. [0] is string correspond to legacy template th3f name, [1] is 
        more detailed describtion of what's in there.
        """
        self.process['wjets']=['WJets','wjets']
        self.process['other']  =['bck','single_top/non-semilep_ttbar']
        self.process['gg']   =['gg','gg/qg_ttbar']
        self.process['qqs']   =['qqs','symetric qq_ttbar']
        self.process['qcd'] =['ntmj','qcd']
        print '(info) Done define_process'


    def define_sys(self):
        """
        define types of systematic variation
        """
        if self.lep_type == 'mu':
            if self.use_sys:
                self.systematics_vary = ['btag_eff_reweight','lepID_reweight','trigger_reweight','tracking_reweight','lepIso_reweight','Pdf_weights']
                self.systematics_fixed = ['pileup_reweight']
            else:
                self.systematics_vary = []
                self.systematics_fixed = ['btag_eff_reweight','lepID_reweight','trigger_reweight','tracking_reweight','lepIso_reweight','pileup_reweight']
        elif self.lep_type == 'el':
            if self.use_sys:
                    self.systematics_vary = ['btag_eff_reweight','lepID_reweight','trigger_reweight','Pdf_weights']
                    self.systematics_fixed = ['tracking_reweight','lepIso_reweight','pileup_reweight']
            else:
                    self.systematics_vary = []
                    self.systematics_fixed = ['btag_eff_reweight','lepID_reweight','trigger_reweight','tracking_reweight','lepIso_reweight','pileup_reweight']
        if self.top_reweight:
            self.systematics_fixed.append('top_pT_reweight')

        self.systematics_all = self.systematics_fixed+self.systematics_vary
        self.sys_versions = ['plus','minus'] # naming conventions for theta
        self.sys_names =['up','down'] # naming coventions for legacy ttree
        print '(info) Done define_sys'


    def getTempKeys(self,old_key_string,new_key_string,title):
        # old_keys are original name of th3f in templates.root, new_keys are names of TH1D in thetaTemp.root
        old_keys,new_keys,new_titles = [],[],[]
        old_keys.append('f_%s_plus'%old_key_string) 
        old_keys.append('f_%s_minus'%old_key_string) 
        new_keys.append('f_plus__%s'%new_key_string)
        new_keys.append('f_minus__%s'%new_key_string) 
        new_titles.append('%s lep charge plus'%title)
        new_titles.append('%s lep charge minus'%title)
        return old_keys,new_keys,new_titles
    
    def Add_1D_temp(self,template,tempName,tempTitle):
        """
        input: an object of template class
        output: add 1D template into thetaTemp.root 
        """
        self.outfile.cd()
        hist_1D = template.convertTo1D()
        hist_1D.SetName(tempName)
        hist_1D.SetTitle(tempTitle)
        hist_1D.SetDirectory(0)
        hist_1D.Write()
        self.thetaHistList.append(hist_1D)

    def output_original_templates(self,template):
        """
        input: an object of template class
        output: add original 3D and 3 1D projections into aux file
        """
        original_templates=template.getOriginalTemps() # 3D,x,y,z
        self.write_templates_to_auxfile(hist_list=original_templates)


    def importTemp(self,tmp_file):
        """
        process smoothed templates.root file
        input: original root file with TH3D templates
        output: new root file for theta feedin file with an unrolled TH1D

        The templates.root file has keys in the form below, with seperate th3f of plus and minus lep charge
        KEY: TH3F   f_WJets_minus;1 f_WJets_minus
        KEY: TH3F   f_WJets_plus;1  f_WJets_plus
        KEY: TH3F   f_bck_minus;1   f_bck_minus
        KEY: TH3F   f_bck_plus;1    f_bck_plus
        """
        tfile = ROOT.TFile(tmp_file)
        self.outfile.cd()
        # loop over all processes,e.g., wjets, gg , and convert 3D templates into 1D and write into new template file
        for key,value in self.process.iteritems():
            old_key_string = value[0]
            new_key_string = key
            hist_title = value[1]
            # get keys of histograms in original templates and in new thetaTemplates
            old_keys,new_keys,new_titles = self.getTempKeys(old_key_string=old_key_string,new_key_string=new_key_string,title=hist_title)
            # Loop over all TH3D for this process, such as f_wjets_plus,f_wjets_minus,f_wjets_minus_btag_up etc
            for i,value in enumerate(old_keys):
                i_oldkey = old_keys[i]
                i_newkey =  new_keys[i]
                i_newtitle = new_titles[i]
                tmp_hist = tfile.Get(i_oldkey)
                if not tmp_hist:
                    print '(debug) histogram %s not found in %s!'%(i_oldkey,tfile)
                else:
                    # rescale Data-driven QCD templates using extrapolation factors get from ABCD
                    if 'qcd' in i_newkey:
                        print '(info) QCD templates rescaling to %.3f'%self.QCD_SF
                        tmp_hist.Scale(self.QCD_SF)
                    # use template class to convert the TH3D to TH1D
                    #   def __init__(self,name,formatted_name,bin_type='fixed') :
                    template_obj = template.template(i_newkey,i_newkey,self.bin_type)
                    template_obj.histo_3D = tmp_hist
                    self.Add_1D_temp(template=template_obj,tempName=i_newkey,tempTitle=i_newtitle)
                    # also include original template TH3D and projections in aux file for later cross check
                    # original_temps
                    # self.write_templates_to_auxfile()
        print '(info) Done importTemp'
        tfile.Close()

    def TTtreeToTemplates(self,ttree_file,process_name = None):
        """
        Input: a root with ttree
        Output: add 1D templates into template file, and control plots to aux file
        for data, just one template
        for MC, will include nominal, and sys_up and down templates for each sys
        JEC/JER sys is handled differently. They come from different versions of root files
        Controled by self.JEC_sys, append __JES__plus/minus to process_name, and ignore other correction SFs
        """
        print '(info) Process TTtreeToTemplates with file %s'%(ttree_file)
        # set up template name and weight
        tmp_type = 'NA'
        if process_name=='DATA':
            self.makeTemplates(file_list=ttree_file,process_name=process_name)
            print '(info) Done TTtreeToTemplates for %s.'%process_name
        else:
            # Loop over files in the list
            weight_lists = {} # e.g. wjets__btag_eff_reweight__up:[w1,w2,..w_up]
            norm_weights = []
            # Get normalization weights for each file in the list
            for ifile in ttree_file:                    
                # load info about samples from txt file
                sample_info_obj = self.samples_obj.get_sample_info(ifile)
                norm_weight = sample_info_obj.weight
                norm_weights.append(norm_weight)
                if tmp_type == 'NA':
                    tmp_type = sample_info_obj.type
            print '(Debug) Process type %s'%tmp_type
            # prepare for sys weights of each variation 
            # first make the nominal templates
            weight = []
            weight += self.fixed_weight+self.varied_weight_nominal
            # template_name is what comes into the theta template naming convention
            if self.JEC_sys == 'nominal':
                template_name = process_name # e.g. wjets
            else:
                template_name = '%s__%s'%(process_name,self.JEC_sys)
            weight_lists[template_name] = weight
            # then make up and down templates for each sys for all MC templates BUT QCD. 
            # also skip this if JEC_sys is not None
            if process_name != 'qcd' and self.JEC_sys == 'nominal':
                for i,isys in enumerate(self.systematics_vary):
                    for j,jversion in enumerate(self.sys_versions): # ['plus','minus']
                        weight = []
                        # append fixed weights
                        weight += self.fixed_weight
                        # add a variant of isys weight
                        weight.append(self.varied_weight[jversion][i])
                        # add the rest of variable sys using the nominal values
                        nominal_weights = copy.copy(self.varied_weight_nominal)
                        nominal_weights.pop(i)
                        weight += nominal_weights
                        # set up names for this variant template: e.g. wjets__btag_eff_reweight__up
                        template_name = '%s__%s__%s'%(process_name,isys,jversion)
                        # finally, add a new template
                        weight_lists[template_name]=weight
            # finally make a template for each version of weight list
            for key,value in weight_lists.iteritems():
                self.makeTemplates(file_list=ttree_file,process_name=key,weight_list=value,norm_weights=norm_weights,tmp_type = tmp_type)
            print '(info) Done making all templates for MC.'
    

    def GetWeights(self) :
        """
        Make lists contains strings of the names of the weights in the ttree
        self.fixed_weight = []
        self.varied_weight_nominal = []
        self.varied_weight = {}
        """
        # initiate
        self.fixed_weight = []
        self.varied_weight_nominal = []
        self.varied_weight = {'plus':[],'minus':[]}
        # set up fixed sys reweight 
        for iweight in self.systematics_fixed:
            self.fixed_weight.append(iweight)
        # set up varying weight
        for iweight in self.systematics_vary:
            # first set nominal ones
            self.varied_weight_nominal.append(iweight)
            # set up high version
            weight_key_up = '%s_hi'%iweight # like btag_eff_reweight_high
            self.varied_weight['plus'].append(weight_key_up)
            # set up low version
            weight_key_down = '%s_low'%iweight
            self.varied_weight['minus'].append(weight_key_down)
        print '(info) Done GetWeights'

    def makeTemplates(self,file_list,process_name,weight_list=[1.0],norm_weights = [1.0],tmp_type='NA'):
        """
        input: list of ttree, list of weight ( a list of 1D arrays ), tmp_name
        output: a 1D template and control plots of 3d projections
        """
        if self.verbose:
            print '(DEBUG) makeTemplates process %s'%process_name

        tmpNames=[]
        tmpNames.append('f_plus__%s'%process_name)
        tmpNames.append('f_minus__%s'%process_name)
        # accomodate fqq,fgg specially
        # Determing if it is gg/qq process, if so , need to addtwice
        add_twice=False
        if 'gg' in tmp_type or 'qq' in tmp_type:
            add_twice = True
            print '(Debug) Add twice for %s, type %s'%(process_name,tmp_type)
        else:
            print '(Debug) Not add twice for %s, type %s'%(process_name,tmp_type)
        # Add two additional templates for qq, qq__AFB__plus and qq__AFB__minus, but do this only during making nominal qq templates
        is_qq = False
        if process_name == 'qq':
            tmpNames.append('f_plus__%s__AFB__plus'%process_name)
            tmpNames.append('f_plus__%s__AFB__minus'%process_name)
            tmpNames.append('f_minus__%s__AFB__plus'%process_name)          
            tmpNames.append('f_minus__%s__AFB__minus'%process_name) 
            is_qq = True        
        # create template objects
        tmp_objects=[]
        for item in tmpNames:
            tmp_objects.append(template.template(name=item,formatted_name=item,bin_type=self.bin_type))
        if self.verbose: print '(DEBUG) all weights ',weight_list
  
        # Loop over ttree in list of ttrees
        for i,ifile in enumerate(file_list):
            print '\n(DEBUG) Processing %s'%ifile
            # loading ttree and load branches
            tfile = ROOT.TFile(ifile)
            self.outfile.cd()
            tree_name = helper.GetTTreeName(tfile)
            ttree = tfile.Get(tree_name)
            weights = weight_list
            # check if norm weight is loaded correctly
            if self.verbose:
                print '(DEBUG) %s norm weight = %.3f'%(ifile,norm_weights[i])
            # check if gen_w is in the tree
            if ttree.FindBranch('gen_weight'):
                has_genW = True
            else:
                has_genW = False
            # Loop over entries in ttree and fill templates
            final_total_w = []
            n_entries = ttree.GetEntries()
            n_fills_plus,n_fills_minus = 0,0
            for iev in range(n_entries):
                if iev == self.nevts:
                    print '(info) Reach evts %i. Move on next sample!'%iev
                    break
                ttree.GetEntry(iev)
                # load observables
                cs = ttree.cos_theta_cs
                xf = ttree.Feynman_x
                mtt = ttree.ttbar_mass
                lep_charge = ttree.Q_l
                n_bTags = ttree.n_bTags
                # require nbtags>=2. Need to discuss if should include eventes with more than 2 bs
                if n_bTags < self.btag_cut: continue
                # load genW
                if has_genW:
                    gen_weight = ttree.gen_weight
                else:
                    gen_weight = 1.0
                # get the weight right by looping over a list of arrays(or float)
                if process_name=='DATA': # for data, with no weights
                    if self.use_MC_DATA:
                        total_weight = ttree.total_w
                    else:
                        total_weight = 1
                elif process_name=='qcd':
                    # QCD ttree is a sum of data and MC events in sideband, with MC events having negative weights for substraction purpose
                    norm_weight = ttree.normalization_weight*self.QCD_SF
                    total_weight = norm_weight*gen_weight
                else:
                    norm_weight = norm_weights[i]*gen_weight
                    total_weight = [getattr(ttree,item) for item in weights]
                    total_weight.append(norm_weight)
                    if self.verbose and iev<1: print '(DEBUG) total_weight=',total_weight
                    total_weight = helper.multiply(total_weight)
                    if self.verbose and iev<1: print '(DEBUG) total_weight=%.3f'%total_weight
                # keep final weight for every event for sanity checks later
                final_total_w.append(total_weight)
                # for added twice case, which means templates is qq_* or gg_* temp:
                tmp_add_twice = False
                if add_twice:
                    motherPIDs = ttree.motherPIDs
                    # Do symmetrization only for gg or qq process. Not symmetrize qg process
                    if (motherPIDs[0]==21 and motherPIDs[1]==21) or (motherPIDs[0]+motherPIDs[1]==0):
                                            total_weight *= 0.5
                                            tmp_add_twice = True
                    else: 
                        tmp_add_twice=False 
                #       print '(DEBUG) Set add_twice False!'
                    w_a = ttree.w_a
                # fill templates
                # special note for qq sample: tmp_obj[2,3,4,5] are f_plus_up,down, f_minus_up,down
                if lep_charge>0:
                    tmp_objects[0].Fill(cs,abs(xf),mtt,total_weight)
                    n_fills_plus += 1
                    if tmp_add_twice:
                        tmp_objects[1].Fill(-cs,abs(xf),mtt,total_weight)
                        n_fills_minus +=1
                    if is_qq:
                        tmp_objects[2].Fill(cs,abs(xf),mtt,total_weight*(1+self.AFB_sigma*w_a)) # Q>0, fwd
                        tmp_objects[3].Fill(cs,abs(xf),mtt,total_weight*(1-self.AFB_sigma*w_a)) # Q>0, bwd
                        tmp_objects[4].Fill(-cs,abs(xf),mtt,total_weight*(1-self.AFB_sigma*w_a)) # Q<0, fwd
                        tmp_objects[5].Fill(-cs,abs(xf),mtt,total_weight*(1+self.AFB_sigma*w_a)) # Q<0, bwd
                elif lep_charge<0:
                    tmp_objects[1].Fill(cs,abs(xf),mtt,total_weight)
                    n_fills_minus +=1
                    if tmp_add_twice:
                        tmp_objects[0].Fill(-cs,abs(xf),mtt,total_weight)
                        n_fills_plus += 1
                    if is_qq:
                        tmp_objects[4].Fill(cs,abs(xf),mtt,total_weight*(1+self.AFB_sigma*w_a)) # Q<0, fwd
                        tmp_objects[5].Fill(cs,abs(xf),mtt,total_weight*(1-self.AFB_sigma*w_a)) # Q<0, bwd
                        tmp_objects[2].Fill(-cs,abs(xf),mtt,total_weight*(1-self.AFB_sigma*w_a)) # Q>0, fwd
                        tmp_objects[3].Fill(-cs,abs(xf),mtt,total_weight*(1+self.AFB_sigma*w_a)) # Q>0, bwd
                else:
                    print '(debug) lep_charge==0. Something is wrong!'
                    sys.exit(1)
                if self.verbose:
                    # print 'cs %.2f,xf %.2f,mtt %.2f,lep_charge %i,total_weight %.2f'%(cs,xf,mtt,lep_charge,total_weight)
                    pass
            # print out mean and stdev of final total weight
            final_total_w = np.array(final_total_w)
            print '(Info) %s total_weight mean = %.3f, stdev = %.3f'%(tfile.GetName(),final_total_w.mean(),final_total_w.std()) 
            print '(Debug) Fill plus/minus template %i/%i times'%(n_fills_plus,n_fills_minus)

            if self.verbose : print '(Debug) Add twice is %s\n'%str(add_twice)
        # Write proper unrolled 1D templates into thetaTemp file
        for i,itemp in enumerate(tmp_objects):
            # where theta template name is given
            if self.lep_combo:
                template_name = '%s_%s'%(self.lep_type,tmpNames[i])
            else:
                template_name = tmpNames[i]
            self.Add_1D_temp(template=tmp_objects[i],tempName=template_name,tempTitle=template_name)
            # further, write original projections into aux file for later conparisons
            self.output_original_templates(template=tmp_objects[i])

    def get_rate_templates(self):
        """
        input: self.thetaHistList contains all 1D unrolled templates for theta
        output: build, on top of nominal templates, a rate up and down templates for each bkg as well fqq
        """
        # case1: SF_bkgs, in this case other bkg remains the same, only Fqq and Fgg changes
        # case2: change SFqq, in this case all bkgs remains the same and only Fqq and Fgg changes

        # case1
        # rate_change = (('Wjets',0.2),('other_bkg',0.2),('qq',0.2))

        self.outfile.cd()
        # First get Ntt,Nbkg,Nqq
        self.nominal_hist_map={}
        self.nominal_hist_map['qq'] = ('%s_f_plus__qq'%self.lep_type,'%s_f_minus__qq'%self.lep_type)
        self.nominal_hist_map['gg'] = ('%s_f_plus__gg'%self.lep_type,'%s_f_minus__gg'%self.lep_type)
        self.nominal_hist_map['WJets'] = ('%s_f_plus__WJets'%self.lep_type,'%s_f_minus__WJets'%self.lep_type)
        self.nominal_hist_map['other_bkg'] = ('%s_f_plus__other_bkg'%self.lep_type,'%s_f_minus__other_bkg'%self.lep_type)
        hist_map = self.nominal_hist_map

        Nqq = 0
        Ngg = 0
        # calculate Ntt etc and also make a new hashmap for 1D templates for quick look up using hist.GetName() as key
        self.thetaHistList_map={}
        for hist in self.thetaHistList:
            iname = hist.GetName()
            self.thetaHistList_map[iname] = hist
            if iname in hist_map['qq']: # iname is like f_plus__qq or el_f_plus__qq, depend on if is lep combined
                Nqq += hist.Integral()
            elif iname in hist_map['gg']:
                Ngg += hist.Integral()
        Ntt = Nqq+Ngg

        for item in self.rate_change:
            process_name = item[0]
            SF_sigma_up = item[1]
            SF_sigma_down = -item[1]
            new_templates=[]
            # calculate SF for each process in R_process templates making
            # case 1
            if process_name == 'qq':
                SF_qq_up = 1+SF_sigma_up
                SF_gg_up = 1-SF_sigma_up*Nqq/Ngg
                SF_qq_down = 1+SF_sigma_down
                SF_gg_down = 1-SF_sigma_down*Nqq/Ngg
                new_templates+=[[hist_map['qq'],SF_qq_up,SF_qq_down],[hist_map['gg'],SF_gg_up,SF_gg_down]]
            #case 2
            else:
                # get Nbkg_i
                Nbkg_i=0
                for iname in hist_map[process_name]:
                    ihist = self.thetaHistList_map.get(iname,0)
                    if ihist==0:
                        print '(Error) No %s 1D temp found!'%iname
                        sys.exit(1)
                    else:
                        Nbkg_i += ihist.Integral()
                # calculate SFs
                # SF_bkg = 1+sigma_bkg
                # SF_tt = 1-sigma_bkg*N_bkg/N_tt
                SF_bkg_up = 1+SF_sigma_up
                SF_bkg_down = 1+SF_sigma_down
                SF_qq_up = 1-SF_sigma_up*Nbkg_i/Ntt
                SF_qq_down = 1-SF_sigma_down*Nbkg_i/Ntt
                SF_gg_up = SF_qq_up
                SF_gg_down = SF_qq_down
                new_templates += [[hist_map['qq'],SF_qq_up,SF_qq_down],[hist_map['gg'],SF_gg_up,SF_gg_down],[hist_map[process_name],SF_bkg_up,SF_bkg_down]]

            # Apply SFs on nominal templates and add the new ones into template.root files
            # Loop over all templates that are affected by R_process
            for item in new_templates:
                itemplate_names = item[0] # e.g. ('f_plus__qq','f_minus__qq')
                iSF_plus = item[1]
                iSF_minus = item[2]
                # loop over all subversions of templates and create scaled ones
                for iname in itemplate_names:
                    ihist_nominal = self.thetaHistList_map.get(iname,0)
                    if ihist_nominal==0:
                        print '(Error) No %s 1D temp found!'%iname
                        sys.exit(1)
                    if process_name != 'qq':
                        print 'process_name = ',process_name
                        new_plus_name = '%s__R_%s_%s__plus'%(iname,process_name,self.lep_type)
                        new_minus_name = '%s__R_%s_%s__minus'%(iname,process_name,self.lep_type)
                    else:                        
                        new_plus_name = '%s__R_%s__plus'%(iname,process_name)
                        new_minus_name = '%s__R_%s__minus'%(iname,process_name)
                    ihist_scaled_plus = ihist_nominal.Clone(new_plus_name)
                    ihist_scaled_minus = ihist_nominal.Clone(new_minus_name)
                    ihist_scaled_plus.Scale(iSF_plus)
                    ihist_scaled_minus.Scale(iSF_minus)
                    # Add into all thetaTemp list
                    self.thetaHistList.append(ihist_scaled_plus)
                    self.thetaHistList.append(ihist_scaled_minus)
                    # also write new templates into temp root file
                    ihist_scaled_plus.Write()
                    ihist_scaled_minus.Write()
                    # print out info message
                    print '             Adding template with name %s'%new_plus_name
                    print '             Adding template with name %s'%new_minus_name

        # finish making all new templates
        print '(info) Finish get_rate_templates.'


    def makeControlPlots(self):
        """
        Take 1D templates and create all control plots, such as 3D and x,y,z projections
        Write all control stuff in an aux.root file
        """ 
        # for every 1D templates, get 3D original and 3 1D projection histograms
        for ihist in self.thetaHistList:
            hname = ihist.GetName()+'_proj'
            template_obj =  template.template(hname,hname+' projected back from 1D hist', self.bin_type)
            #   getTemplateProjections  return [self.histo_3D,self.histo_x,self.histo_y,self.histo_z]
            hist_proj = template_obj.getTemplateProjections(ihist)
            self.write_templates_to_auxfile(hist_list=hist_proj)
            del(template_obj)
        print '(info) Done makeControlPlots.'


    def write_templates_to_auxfile(self,hist_list):
        """
        input: a list of [TH3D,x_proj,y_proj,z_proj]
        output: write histograms and nicer plots into aux.root file
        """
        for i,item in enumerate(hist_list):
            if i!=0:
                item.SetMinimum(0)
            self.outfile_aux.cd('hists/')
            item.Write()
            # make nicer plots and write into another dir
            c = ROOT.TCanvas()
            c.SetName(item.GetName())
            item.Draw('hist E1')
            self.outfile_aux.cd('plots')
            c.Write()       

    def __del__(self):
        self.outfile.Close()
        self.outfile_aux.Close()
        print '(info) Closeup thetaTemp object.'



