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

DEBUG_RUNS = -1
SYS_debug = 'all'

#DEBUG_RUNS = 1
#SYS_debug = 'JER'

class plotter(object):
    """docstring for plotter"""
    def __init__(self, template_file , lep_type, verbose = False,bin_type='fixed'):
        super(plotter, self).__init__()
        self.input_file = template_file
        self.verbose = verbose
        self.is_postfit = False
        if 'postfit' in template_file:
            self.is_postfit = True
        self.template_file = ROOT.TFile(template_file)
        self.all_hist = []
        self.projections = OrderedDict() # {'wjets_plus':[hist_cs,hist_xf,hist_mass],,,}
        self.process = OrderedDict()
        self.observables = ['comb']
        self.stack_lists = ['x','y','z']
        self.stack_xaxis = ['cos#theta*','|x_{F}|','M_{tt}(GeV)']
        self.stacks = {} # keys: 'plus_x','minus_y' etc
        self.temp_shapes = {} # same key as self.stacks, to organize projections of all templates
        self.DATA_proj = {} # contains data projection th1 with same key as self.stacks
        self.process_counts = OrderedDict() # for total number of events given the 1D templates , for R_process calculation
        self.bin_type = bin_type
        self.shapes_to_compare = ['qq','gg','WJets','other_bkg']
        self.lep_type = lep_type

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
        if not self.is_postfit:
            self.makeQualityPlots()
            self.plotShapes()

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
        self.process['other_bkg'] ='TT/SingleT/DY'
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
                # also set up temp_shapes ={'plus_x':[h1,h2,...],...} for shape comparison control plots as well
                self.temp_shapes[ikey] =[]

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
                # set xaxis title
                ihist.GetXaxis().SetTitle(self.stack_xaxis[i])
                stack_key = '%s_%s'%(iobs,var)
                if iprocess != 'DATA':
                    # Add to proper THStack if it is not DATA
                    ihist.SetFillColor(icolor)
                    self.stacks[stack_key].Add(ihist)
                    self.temp_shapes[stack_key].append(ihist)
                else:
                    # Keep data projections in another hashtable
                    self.DATA_proj[stack_key] = ihist
                # add to legend
                if stack_key in ['plus_x','el_f_plus_x']:
                    print '(info) Adding %s into stack'%iprocess_title
                    if iprocess != 'DATA':
                        self.legend.AddEntry(ihist,iprocess_title,"F")
                    else:
                        self.legend.AddEntry(ihist,iprocess_title,"lep")
                # add total integral of hists into a list for later calculation of R_process
                if 'comb_x' in stack_key:
                    self.process_counts[iprocess] += ihist.Integral()
                #    print '(info) Add counts of process %s from hist %s into count table'%(iprocess,stack_key)
        self.write_stack_to_auxfile(hist_list=self.stacks.values(),legend=self.legend)
        print '(info) Done stackMaker.'


    def defineIO(self):
        input_name = self.input_file.split('/')[-1]
        input_name = input_name.split('.root')[0]
        self.input_name = input_name
        self.output_dir = 'plots_'+input_name
        self.output_src_dir = self.output_dir+'/src'
        self.output_extra_dir = self.output_dir+'/extra' # store other formats of plots, say, pdf
        print '(info) self.output_dir:',self.output_dir
        print '(info) self.output_src_dir:',self.output_src_dir
        print '(info) self.output_extra_dir:',self.output_extra_dir

        # make sure output dir exists
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
            print '(info) Making new dir %s'%self.output_dir
        if not os.path.exists(self.output_src_dir):
            os.mkdir(self.output_src_dir)
            print '(info) Making new dir',self.output_src_dir
        if not os.path.exists(self.output_extra_dir):
            os.mkdir(self.output_extra_dir)
            print '(info) Making new dir',self.output_extra_dir
        # create output files for hists and plots
        output_rootfile = '%s/control_%s.root'%(self.output_src_dir,input_name)
        self.outfile_aux = ROOT.TFile(output_rootfile,'recreate')
        self.outfile_aux.mkdir('hists/')
        self.outfile_aux.mkdir('plots/')
        print '(info) Making aux output file %s.'%output_rootfile
        plots_file = '%s/plots.root'%self.output_src_dir
        self.plots_file = ROOT.TFile(plots_file,'recreate')
        # Create a txt file for counts
        self.txt_file = open('%s/counts_%s.txt'%(self.output_dir,input_name),'w')

    def combineCharge(self):
        """
        input: self.all_templates ( all 1D nominal or postfit templates )
        output: add charge combined templates into self.all_templates
        """
        comb_temps = {}
        for itemp in self.all_templates:
            # find process names 
            # template name: f_plus__DATA
            temp_name = itemp.GetName()
            proc_name = '__'.join(temp_name.split('__')[1:])
            combined_key = 'f_comb__%s'%proc_name
            if comb_temps.get(combined_key,0)==0:
                print '(info) Making combined temp %s'%combined_key
                tmp_hist = itemp.Clone(combined_key)
                tmp_hist.SetDirectory(0)
                comb_temps[combined_key]=tmp_hist
                # add nominal temp combined to the corresponding list
                if len(temp_name.split('__'))<=2:
                    self.nominal_templates.append(tmp_hist)
            else:
                comb_temps[combined_key].Add(itemp)
        # Append comb templates into self.templates
        tmp_hists = [val for key,val in comb_temps.iteritems()]
        self.all_templates.extend(tmp_hists)
        print '(info) Done combineCharge!'

    def makeControlPlots(self):
        """
        Take 1D templates and create all control plots, such as 3D and x,y,z projections
        Write all control stuff in an aux.root file
        """ 
        all_templates_names = helper.GetListTH1D(self.template_file)
        # get all observable names
        obs_names = list(set([item.split('__')[0] for item in all_templates_names]))
        self.observables += obs_names
        # add templates into nominal or systematic list
        self.all_templates = []
        self.nominal_templates = []
        # print '(info) Keeping these hists.'
        for name in all_templates_names:
            # only keep nominal templates for simplicity
            if len(name.split('__'))<=2:
                 self.nominal_templates.append(self.template_file.Get(name))
            if SYS_debug in name or len(name.split('__'))<=2 or SYS_debug == 'all' : # just for debugging purpose 
                self.all_templates.append(self.template_file.Get(name))
            # print name

        # make charge combined templates here
        self.combineCharge()

        # for every 1D templates, get 3D original and 3 1D projection histograms, and put into a new dict
        tmp_proj_all = {} # has all template projections
        tmp_projections = {} # only has nominal template projections
        for ihist in self.all_templates:
            # get 3 projections from 1D templates and write in aux file
            hname = ihist.GetName()+'_proj'
            template_obj =  template.template(name=hname,formatted_name=hname+' projected back from 1D hist',bin_type=self.bin_type)
            #   getTemplateProjections  return [self.histo_3D,self.histo_x,self.histo_y,self.histo_z]
            hist_proj = template_obj.getTemplateProjections(ihist)
            self.write_templates_to_auxfile(hist_list=hist_proj[1:])
            tmp_proj_all[hname] = hist_proj[1:]
            del(template_obj)
            # write 1D templates into aux file
            self.write_templates_to_auxfile([ihist])
            # assign projections to corresponding MC processes, for nominal templates only
            if ihist not in self.nominal_templates: continue
            for key,value in self.process.iteritems():
                for iobs in self.observables:
                    newkey = '%s__%s'%(key,iobs) # wjets_plus etc
                    if key in hname and iobs in hname:
                        tmp_projections[newkey] = hist_proj[1:]

        # re-arrange projections with the same order or self.process, only add nominal templates for plot stacks later
        for iprocess in self.process:
            for key,value in tmp_projections.iteritems():
                if iprocess in key:
                    self.projections[key]=value

        # loop over all templates, assign [up,down,nominal] 3d hist triples and plot templates comparison plots
        # key: f_plus__gg__JES__minus_proj, val:[proj_x,proj_y,proj_z]
        # key form: channel__process__SYS__plus/minus_proj
        tmp_index = 0
        for key in tmp_proj_all:
            try:
                channel,process,SYS,pm = key.split('__')
            except ValueError:
                continue
            if 'plus' not in pm or 'comb' not in channel: continue
            if tmp_index==DEBUG_RUNS:
                print 'End debug run for sys plotter!'
                break
            tmp_index += 1
            proj_plus = tmp_proj_all['%s__%s__%s__plus_proj'%(channel,process,SYS)]
            proj_nom = tmp_proj_all['%s__%s_proj'%(channel,process)]
            proj_minus = tmp_proj_all['%s__%s__%s__minus_proj'%(channel,process,SYS)]
            for i in range(len(tmp_proj_all[key])):
                # loop over proj x,y,z
                hlist = [ proj_plus[i],proj_nom[i],proj_minus[i] ]
                self.plot_sys_temp(hlist,'__'.join([channel, process, SYS]),axis=i)


    def plot_sys_temp(self,hlist,title,axis):
        """
        input: hlist=[h_sys_up,h_nom,h_sys_down]
        output: plot systematics templates comparison plots
        """
        self.shape_hists = []
        tmp_legend = ROOT.TLegend(0.67,0.25,0.83,0.40)
        axis_name = ['cstar','xf','mtt']
        title = title+'__'+axis_name[axis]
        #  ['cos#theta*','|x_{F}|','M_{tt}(GeV)']
        c = ROOT.TCanvas()
        c.SetName(title)
        colors = [ROOT.kRed,ROOT.kBlack,ROOT.kBlue]
        leg_entry = ['up','nominal','down']
        # first get ymax
        ymax = 1.1*max(item.GetMaximum() for item in hlist)
        for i in range(len(hlist)):
            h1 = hlist[i]
            h1.SetMaximum(ymax)
            h1.SetMinimum(0)
            # draw histogram
            h1.SetTitle(' '.join(title.split('__')))
            h1.SetFillColor(0)
            h1.SetLineColor(colors[i])
            tmp_legend.AddEntry(h1,leg_entry[i],"L")
            if i==0:
                h1.Draw("hist e0")
                h1.GetXaxis().SetTitle(self.stack_xaxis[axis])
            else:
                h1.Draw("hist e0 same")  
        tmp_legend.SetFillColor(0)          
        tmp_legend.Draw("same")
        c.SetTitle(' '.join(title.split('__')))
        # save
        # c.SaveAs('%s/%s_sys.png'%(self.output_dir,title))
        c.SaveAs('%s/%s_sys.pdf'%(self.output_extra_dir,title))
        self.plots_file.cd()
        c.Write()       
        del(c)

    def plotShapes(self):
        """
        input: self.temp_shapes
        output: plot hist and save in output dir
        """
        self.shape_hists = []
        tmp_legend = ROOT.TLegend(0.7,0.7,0.9,0.9)

        has_leg = False
        for key,hist_list in self.temp_shapes.iteritems():
            c = ROOT.TCanvas()
            isFirstHist=True
            # first get ymax
            ymax= 0
            for i,ihist_0 in enumerate(hist_list):
                ihist = ihist_0.Clone()
                ihist.Scale(1.0/ihist.Integral())
                ymax = max(ymax,ihist.GetMaximum())
            ymax *= 1.1
            # draw histogram
            for i,ihist_0 in enumerate(hist_list):
                if 'other_bkg' in ihist_0.GetName():
#                    continue
                    pass
                # Need to find correct color first....
                # f_plus__qcd_proj_x 
                proc_name = ihist_0.GetName().split('__')[-1].split('_proj')[0]
                icolor = helper.getColors(proc_name)
                ihist = ihist_0.Clone()
#                ihist.SetDirectory(0)
                ihist.SetFillColor(0)
                ihist.SetLineColor(icolor)
                # normalize to 1
                ihist.Scale(1.0/ihist.Integral())
                ihist.SetMaximum(ymax)
                ihist.SetMinimum(0)
                self.shape_hists.append(ihist)
                if not has_leg:
                    tmp_legend.AddEntry(ihist,proc_name,"L")
                if isFirstHist:
                    ihist.Draw("hist e0")
                    isFirstHist=False
                else:
                    ihist.Draw("hist same e0")
            print '%s done'%key
            has_leg = True
            self.legend.Draw("same")
            # c.SaveAs('%s/%s_shapes.png'%(self.output_dir,key))
            c.SaveAs('%s/%s_shapes.pdf'%(self.output_dir,key))

            del(c)
        print '(info) Done plotShapes!'


    def makeQualityPlots(self):
        """
        take the 1D templates and make a 1D hist with number of events in each bin
        input: self.nominal_templates
        output: add hists h_quality_plus and h_quality_minus to aux file
        """
        # Get all MC templates
        template_quality_hists = []
        data_quality_hists = []

        # loop over plus or minus templates
        for iobs in self.observables:
            total_temp = None
            # Loop over all templates
            for ihist in self.nominal_templates:
                hname = ihist.GetName()
                # pick exactly the process and observable from the templates
                if iobs not in hname: continue
                # make quality hists for Data
                if 'DATA' in hname:
                    hist_quality = helper.GetQualityPlots_data(ihist) 
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
            quality_total_temp = helper.GetQualityPlots_MC(total_temp)
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
            # c.SetLogy()
            c.SetName(item.GetName())
            item.Draw('hist')
            # set xaxis title
            item.GetXaxis().SetTitle('frac error')
            c.Update()
            if len(data_quality_hists)>0:
                data_quality_hists[i].Draw('same hist')
            tmp_legend.Draw('Same')
            c.Write()
            c.SaveAs('%s/%s.pdf'%(self.output_dir,c.GetName()))
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
            print item.GetName()
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
            # print '[data_proj] keys:',key
            istack = self.stacks[key]
            data_hist = value
            # def comparison_plot_v1(mc_,data_,legend,event_type='plots',draw_option = 'h',logy=False):
            lep_type_ = self.lep_type
            # find lep charge type, plus or minus
            if 'plus' in key:
                charge_type = 'Q>0'
            elif 'minus' in key:
                charge_type = 'Q<0'
            else:
                charge_type = 'Charge combined'
            # assign lep type, to be e+jets, mu+jets or l+jets
            if self.lep_type=='combo' :
                if 'el' in key:
                    lep_type_ = 'e+Jets'
                elif 'mu' in key:
                    lep_type_ = '#mu + Jets'
                else:
                    lep_type_ = 'l+Jets'
            c_compare,h_err = helper.comparison_plot_v1( mc_=istack,data_=data_hist,outputdir=self.output_src_dir,\
                legend=self.legend,event_type='%s/%s'%(self.output_dir,key),bin_type = self.bin_type,\
                lep_type='%s,%s'%(lep_type_,charge_type) )

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
            # item.Write()
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
            if value==0:continue
            towrite += '%15s,%15i,%15.1f%%\n'%(key,value,value*1.0/total_temp_counts*100)
        # # Add AFB from data
        # data_cs_proj = [self.projections['DATA__minus'][0],self.projections['DATA__plus'][0]]
        # bin_edge1 = data_cs_proj[0].FindFixBin(-1)
        # bin_edge2 = data_cs_proj[0].FindFixBin(-0.001)
        # bin_edge3 = data_cs_proj[0].FindFixBin(0.001)
        # bin_edge4 = data_cs_proj[0].FindFixBin(1)
        # # check if bin_edge2 and 3 correspond to 1 bin difference
        # if bin_edge3-bin_edge2!=1:
        #     print '(debug) cs=-0.001 and 0.001 are not two adjecant bins!'
        # N_fwd,N_bwd = 0,0
        # for ihist in data_cs_proj:
        #     N_bwd += ihist.Integral(bin_edge1,bin_edge2)
        #     N_fwd += ihist.Integral(bin_edge3,bin_edge4)
        # AFB = (N_fwd-N_bwd)/(N_fwd+N_bwd)
        # towrite += '                Data AFB=%.3f\n'%AFB
        # write into file
        self.txt_file.write(towrite)
        print '(info) Done write_counts_to_file .'


    def __del__(self):
        self.template_file.Close()
        self.outfile_aux.Close()
        self.txt_file.close()
        self.plots_file.Close()
        print '(info) Closeup output file.'



