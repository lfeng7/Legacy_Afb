"""
plotter class for 1D unrolled histogram in templates
"""
import thetaTempMaker
import ROOT
import template
import helper
import os

class plotter(object):
    """docstring for plotter"""
    def __init__(self, template_file):
        super(plotter, self).__init__()
        self.input_file = template_file
        self.template_file = ROOT.TFile(template_file)
        self.all_hist = []
        self.projections = {} # {'wjets_plus':[hist_cs,hist_xf,hist_mass],,,}
        self.process = {}
        self.observables = ['plus','minus']
        self.stack_lists = ['x','y','z']
        self.stack_xaxis = ['cos#theta*','|x_{F}|','M_{tt}(GeV)']
        self.stacks = {} # keys: 'plus_x','minus_y' etc
        self.DATA_proj = {} # contains data projection th1 with same key as self.stacks

    def main(self):
        """
        Main function
        """
        self.define_process()
        self.defineIO()
        self.makeControlPlots() # make 3 projections for each 1D templates
        self.stackMaker() # add projections into proper stacks
        self.make_data_mc_comparison()



    def define_process(self):
        """
        define template processes in a dict. [0] is string correspond to legacy template th3f name, [1] is 
        more detailed describtion of what's in there.
        """
        self.process['wjets'] ='WJets'
        self.process['other'] ='s_t/tt_other'
        self.process['gg']    ='gg/qg_ttbar'
        self.process['qqs']   ='qqs_ttbar'
        self.process['qcd']   ='QCD'
        self.process['DATA']   ='DATA'


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
            iprocess = key.split('_')[0] #wjets
            iprocess_title = self.process[iprocess]
            iobs = key.split('_')[1] #plus
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
                    self.legend.AddEntry(ihist,iprocess_title,"F")
        self.write_stack_to_auxfile(hist_list=self.stacks.values(),legend=self.legend)
        print 'Done stackMaker.'

    def defineIO(self):
        input_name = self.input_file.split('/')[-1]
        self.output_dir = input_name.split('.root')[0]+'_plots'
        # make sure output dir exists
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
            print 'Making new dir %s'%self.output_dir
        output_rootfile = '%s/control_%s'%(self.output_dir,input_name)
        self.outfile_aux = ROOT.TFile(output_rootfile,'recreate')
        self.outfile_aux.mkdir('hists/')
        self.outfile_aux.mkdir('plots/')
        print 'Making aux output file %s.'%output_rootfile

    def makeControlPlots(self):
        """
        Take 1D templates and create all control plots, such as 3D and x,y,z projections
        Write all control stuff in an aux.root file
        """ 
        all_templates_names = helper.GetListTH1D(self.template_file)
        self.all_templates = []
        for name in all_templates_names:
            self.all_templates.append(self.template_file.Get(name))
        # for every 1D templates, get 3D original and 3 1D projection histograms
        for ihist in self.all_templates:
            # get 3 projections from 1D templates and write in aux file
            hname = ihist.GetName()+'_proj'
            template_obj =  template.template(hname,hname+' projected back from 1D hist')
            #   getTemplateProjections  return [self.histo_3D,self.histo_x,self.histo_y,self.histo_z]
            hist_proj = template_obj.getTemplateProjections(ihist)
            self.write_templates_to_auxfile(hist_list=hist_proj)
            del(template_obj)

            # write 1D templates into aux file
            self.write_templates_to_auxfile([ihist])

            # assign projections to corresponding MC processes
            for key,value in self.process.iteritems():
                for iobs in self.observables:
                    newkey = '%s_%s'%(key,iobs) # wjets_plus etc
                    if key in hname and iobs in hname:
                        self.projections[newkey] = hist_proj[1:]

        print 'Done makeControlPlots.'

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
            print 'Data is not available. Will not make Data/MC comparison plots!'
            pass
        # Loop over data projections and make plots
        for key,value in self.DATA_proj.iteritems():
            istack = self.stacks[key]
            data_hist = value
            # def comparison_plot_v1(mc_,data_,legend,event_type='plots',draw_option = 'h',logy=False):
            c_compare = helper.comparison_plot_v1(mc_=istack,data_=data_hist,legend=self.legend,event_type='%s/%s'%(self.output_dir,key))
            self.outfile_aux.cd('plots')
            c_compare.Write()
        print 'Done making data/mc comparison plots.'


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
        self.template_file.Close()
        self.outfile_aux.Close()
        print 'Closeup output file.'



