"""
plotter class for 1D unrolled histogram in templates
"""
import thetaTempMaker
import ROOT
import template
import helper
import os
import numpy as np
from collections import OrderedDict

class plotter(object):
    """docstring for plotter"""
    def __init__(self, template_file , verbose = False,bin_type='fixed'):
        super(plotter, self).__init__()
        self.input_file = template_file
        self.verbose = verbose
        self.template_file = ROOT.TFile(template_file)
        self.all_hist = []
        self.projections = OrderedDict() # {'wjets_plus':[hist_cs,hist_xf,hist_mass],,,}
        self.process = OrderedDict()
        self.observables = ['plus','minus']
        self.stack_lists = ['x','y','z']
        self.stack_xaxis = ['cos#theta*','|x_{F}|','M_{tt}(GeV)']
        self.stacks = {} # keys: 'plus_x','minus_y' etc
        self.DATA_proj = {} # contains data projection th1 with same key as self.stacks
        self.process_counts = OrderedDict() # for total number of events given the 1D templates , for R_process calculation
        self.bin_type = bin_type

    def main(self):
        """
        Main function
        """
        self.define_process()
        self.defineIO()
        self.makeControlPlots() # make 3 projections for each 1D templates
        self.stackMaker() # add projections into proper stacks
        self.make_data_mc_comparison()
        self.write_counts_to_file()
        self.makeQualityPlots()

    def define_process(self):
        """
        define template processes in a dict. 
        Order of entries adding will decide the order of stack fill process.
        """
        self.process['zjets'] ='DY+Jets'
        self.process['qcd']   ='QCD'
        self.process['WJets'] ='W+Jets'
        self.process['singleT'] ='Single Top'
        self.process['tt_bkg'] ='None-semilep TT'
        self.process['gg']    ='gg/qg TT'
        self.process['qq']   ='qq TT'
        self.process['DATA']   ='DATA'
        # initiate a process count table
        for ikey in self.process:
            self.process_counts[ikey]=0
        print '(info) Done define_process.'


    def stackMaker(self):
        """
        input: self.projections
        output: stackplots
        """
        for iobs in self.observables: # three projections for plus and minus charge templates 
            # set up stacks
            for i,value in enumerate(self.stack_lists):
                ikey = '%s_%s'%(iobs,value)
                istack = ROOT.THStack('stack_%s_%s'%(iobs,value),'%s %s projection comparison'%(iobs,value))
                self.stacks[ikey] = istack
                # keys of self.stack = plus_x

        self.legend = ROOT.TLegend(0.5,0.5,0.66,0.86)
        self.legend.SetName('legend')
        # Loop over projections of processes and add proper hist into stacks
        for key,value in self.projections.iteritems():
            # key is wjets_plus
            iprocess = key.split('__')[0] #wjets
            iprocess_title = self.process[iprocess]
            iobs = key.split('__')[1] #plus
            icolor = helper.getColors(iprocess)
            # Loop over x,y,z projected hists
            for i in range(len(self.stack_lists)):
                var = self.stack_lists[i] # like x
                ihist = value[i]
                stack_key = '%s_%s'%(iobs,var)
                if iprocess != 'DATA':
                    # Add to proper THStack if it is not DATA
                    ihist.SetFillColor(icolor)
                    self.stacks[stack_key].Add(ihist)
                else:
                    # Keep data projections in another hashtable
                    self.DATA_proj[stack_key] = ihist
                # add to legend
                if stack_key=='plus_x':
                    print '(info) Adding %s into stack'%iprocess_title
                    self.legend.AddEntry(ihist,iprocess_title,"F")
                # add total integral of hists into a list for later calculation of R_process
                if 'x' in stack_key:
                    self.process_counts[iprocess] += ihist.Integral()
                    # print '(info) Add counts of process %s from hist %s into count table'%(iprocess,stack_key)
        self.write_stack_to_auxfile(hist_list=self.stacks.values(),legend=self.legend)
        print '(info) Done stackMaker.'


    def defineIO(self):
        input_name = self.input_file.split('/')[-1]
        input_name = input_name.split('.root')[0]
        self.input_name = input_name
        self.output_dir = input_name+'_plots'
        # make sure output dir exists
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
            print '(info) Making new dir %s'%self.output_dir
        output_rootfile = '%s/control_%s.root'%(self.output_dir,input_name)
        self.outfile_aux = ROOT.TFile(output_rootfile,'recreate')
        self.outfile_aux.mkdir('hists/')
        self.outfile_aux.mkdir('plots/')
        print '(info) Making aux output file %s.'%output_rootfile
        # Create a txt file for counts
        self.txt_file = open('%s/counts_%s.txt'%(self.output_dir,input_name),'w')


    def makeControlPlots(self):
        """
        Take 1D templates and create all control plots, such as 3D and x,y,z projections
        Write all control stuff in an aux.root file
        """ 
        all_templates_names = helper.GetListTH1D(self.template_file)
        self.all_templates = []
        print '(info) Keeping these hists.'
        for name in all_templates_names:
            # only keep nominal templates for simplicity
            if len(name.split('__'))>2: continue
            self.all_templates.append(self.template_file.Get(name))
            print name
        # for every 1D templates, get 3D original and 3 1D projection histograms
        tmp_projections = {}
        for ihist in self.all_templates:
            # get 3 projections from 1D templates and write in aux file
            hname = ihist.GetName()+'_proj'
            template_obj =  template.template(name=hname,formatted_name=hname+' projected back from 1D hist',bin_type=self.bin_type)
            #   getTemplateProjections  return [self.histo_3D,self.histo_x,self.histo_y,self.histo_z]
            hist_proj = template_obj.getTemplateProjections(ihist)
            self.write_templates_to_auxfile(hist_list=hist_proj)
            del(template_obj)
            # write 1D templates into aux file
            self.write_templates_to_auxfile([ihist])
            # assign projections to corresponding MC processes
            for key,value in self.process.iteritems():
                for iobs in self.observables:
                    newkey = '%s__%s'%(key,iobs) # wjets_plus etc
                    if key in hname and iobs in hname:
                        tmp_projections[newkey] = hist_proj[1:]

        # re-arrange projections with the same order or self.process
        for iprocess in self.process:
            for key,value in tmp_projections.iteritems():
                if iprocess in key:
                    self.projections[key]=value

        print '(info) Done makeControlPlots.'

    def GetQualityPlots(self,hist):
        """
        input: a TH1D 
        output: a new TH1D contains the distribution of BinContents for spotting empty bins 
        """
        hname = hist.GetName()
        N_neg_bin = 0
        hist_quality = ROOT.TH1D('quality_%s'%hname,'quality_%s'%hname,61,-1,60)
        for k in range(hist.GetSize()):
            if not hist.IsBinUnderflow(k) and not hist.IsBinOverflow(k) :
                binCounts = hist.GetBinContent(k)
                if binCounts<0:
                    N_neg_bin += 1
                hist_quality.Fill(hist.GetBinContent(k))
        hist_quality.SetDirectory(0)
        if self.verbose:
            print '(verbose) %s has %i negative bins.'%(hname,N_neg_bin)
        return hist_quality


    def makeQualityPlots(self):
        """
        take the 1D templates and make a 1D hist with number of events in each bin
        input: self.all_templates
        output: add hists h_quality_plus and h_quality_minus to aux file
        """
        # Get all MC templates
        template_quality_hists = []
        data_quality_hists = []

        # loop over plus or minus templates
        for iobs in self.observables:
            total_temp = None
            # Loop over all templates
            for ihist in self.all_templates:
                hname = ihist.GetName()
                # pick exactly the process and observable from the templates
                if iobs not in hname: continue
                # make quality hists for Data
                if 'DATA' in hname:
                    hist_quality = self.GetQualityPlots(ihist) 
                    hist_quality.SetLineColor(ROOT.kRed)                    
                    data_quality_hists.append(hist_quality)
                else:
                    # if it's not data, start to build up total templates
                    # if no total template available already, initiate it by clone it.
                    if total_temp is None:
                        total_temp = ihist.Clone('%s_total_template'%iobs)
                        if self.verbose:
                            print '(info) Create a new total_template hist %s'%total_temp.GetName()
                    # if existed, add directly on it
                    else:
                        total_temp.Add(ihist)
                        if self.verbose:
                            print '(info) Add %s in %s'%(hname,total_temp.GetName())
            # Get template_quality_hists
            quality_total_temp = self.GetQualityPlots(total_temp)
            quality_total_temp.SetLineColor(ROOT.kBlue)
            template_quality_hists.append(quality_total_temp)

        # legends
        tmp_legend = ROOT.TLegend(0.7,0.7,0.9,0.9)
        tmp_legend.AddEntry(template_quality_hists[0],'templates','LP')
        tmp_legend.AddEntry(data_quality_hists[0],'DATA','LP')
        # write to aux file
        for i,item in enumerate(template_quality_hists):
            self.outfile_aux.cd('hists/')
            item.Write()
            if len(data_quality_hists)>0:
                data_quality_hists[i].Write()                
            # make nicer plots and write into another dir
            self.outfile_aux.cd('plots')
            c = ROOT.TCanvas()
            c.SetLogy()
            c.SetName(item.GetName())
            item.Draw('hist')
            # set xaxis title
            item.GetXaxis().SetTitle('sum of weights')
            c.Update()
            if len(data_quality_hists)>0:
                data_quality_hists[i].Draw('same hist')
            tmp_legend.Draw('Same')
            c.Write()
        print '(info) Done makeQualityPlots.'


    def write_stack_to_auxfile(self,hist_list,legend):
        """
        input: a list of [TH3D,x_proj,y_proj,z_proj]
        output: write histograms and nicer plots into aux.root file
        """
        # write stack and plots
        for i,item in enumerate(hist_list):
            self.outfile_aux.cd('hists/')
            item.Write()
            # make nicer plots and write into another dir
            c = ROOT.TCanvas()
            c.SetName(item.GetName())
            item.Draw('hist')
            # set xaxis title
            ititle = item.GetHists()[0].GetXaxis().GetTitle()
            item.GetXaxis().SetTitle(ititle)
            c.Update()
            legend.Draw('Same')
            self.outfile_aux.cd('plots')
            c.Write()
        # write legend
        self.outfile_aux.cd('hists/')
        legend.Write()  

    def make_data_mc_comparison(self):
        """
        if self.DATA_proj is not empty, then make Data_MC comparison plots
        """
        if len(self.DATA_proj)==0:
            print '(info) Data is not available. Will not make Data/MC comparison plots!'
            pass
        # Loop over data projections and make plots
        for key,value in self.DATA_proj.iteritems():
            istack = self.stacks[key]
            data_hist = value
            # def comparison_plot_v1(mc_,data_,legend,event_type='plots',draw_option = 'h',logy=False):
            c_compare = helper.comparison_plot_v1(mc_=istack,data_=data_hist,legend=self.legend,event_type='%s/%s'%(self.output_dir,key),bin_type = self.bin_type)
            self.outfile_aux.cd('plots')
            c_compare.Write()
        print '(info) Done making data/mc comparison plots.'


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

    def write_counts_to_file(self):
        """
        Write counts of processes into a txt file
        """
        towrite = '%15s,%15s,%15s\n'%('process','Nevts','R_process')
        # Calculate total number of events in post fit templates
        temp_counts = [self.process_counts[key] for key in self.process_counts if key!='DATA']
        total_temp_counts = np.sum(temp_counts)
        self.process_counts['Total']=total_temp_counts
        # Loop over process to write into txt file
        for key,value in self.process_counts.iteritems():
            towrite += '%15s,%15i,%15.1f%%\n'%(key,value,value*1.0/total_temp_counts*100)
        # write into file
        self.txt_file.write(towrite)
        print '(info) Done write_counts_to_file .'


    def __del__(self):
        self.template_file.Close()
        self.outfile_aux.Close()
        self.txt_file.close()
        print '(info) Closeup output file.'



