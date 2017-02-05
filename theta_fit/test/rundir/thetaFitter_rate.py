execfile("extras.py")
import sys
import math
import os
import copy
import ROOT
import numpy as np
from collections import OrderedDict
execfile("helper.py")
execfile("common.py")
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )

"""
theta fitting code
"""
epsilon = 1E-8

class thetaFitter(object):
    """docstring for thetaFitter"""
    def __init__(self, AFB_sigma=1.0):
        super(thetaFitter, self).__init__()
        
        self.AFB_CENTRAL_VALUE = 0
        self.AFB_sigma = AFB_sigma

        self.shape_sys_gauss = ['btag_eff_reweight','trigger_reweight','lepID_reweight']
	self.shape_sys_gauss_white = []
        #self.shape_sys_gauss_white = 'all'	
        self.flat_param = ['AFB','R_qq','R_WJets','R_other_bkg','qcd_rate','lumi']
        self.obs = 'f_minus'
        self.pois = ['AFB','R_qq','R_other_bkg','R_WJets','qcd_rate']

        # define sigma value here for all parameters
        self.sigma_values = {}
        self.sigma_values['AFB']=AFB_sigma
        self.sigma_values['qcd_rate']=0.2
        self.sigma_values['R_qq']=0.8
        self.sigma_values['R_other_bkg']=0.8
        self.sigma_values['R_WJets']=0.8
        self.sigma_values['lumi']=0.045

        # Toy experiments settings
        self.ntoys = 1000
        self.Toys_per_thread = 50
        # define AFB toy params
        AFB_list = np.arange(-1,1.01,0.5).tolist()
        Rqq_list = np.arange(-0.8,1.21,0.4).tolist()

        #AFB_list=[-1.0,-0.5,-0.2,0,0.2,0.5,1.0]
#        AFB_list = [-0.4,0]
        #AFB_list = [-50,-20,-10,0,10,20,50]#-5,-2,-0.5,0,0.5,2,5,10]#,0,0.3,0.7]
        #Rqq_list = [0,0.5,1.0,1.5]

        self.toy_params = {'AFB':AFB_list,'R_qq':Rqq_list}

        # # counter setup
        # self.process = OrderedDict()
        # self.process_counts = OrderedDict() # for total number of events given the 1D templates , for R_process calculation

    def defineIO(self,template_file_path):
        """
        Handle all IO stuff
        input: template_file_path
        """
        # output
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
            print '(info) Creating new dir %s.'%self.outdir
        self.txtfile = open('%s/results.txt'%self.outdir,'w')
        self.fout = ROOT.TFile('%s/plots.root'%self.outdir,'recreate')

    # def define_process(self):
    #     """
    #     define template processes in a dict. 
    #     Order of entries adding will decide the order of stack fill process.
    #     """
    #     self.process['zjets'] ='DY+Jets'
    #     self.process['qcd']   ='QCD'
    #     self.process['WJets'] ='W+Jets'
    #     self.process['singleT'] ='Single Top'
    #     self.process['tt_bkg'] ='None-semilep TT'
    #     self.process['other_bkg'] ='TT/SingleT/DY'
    #     self.process['gg']    ='gg/qg TT'
    #     self.process['qq']   ='qq TT'
    #     self.process['DATA']   ='DATA'
    #     # initiate a process count table
    #     for ikey in self.process:
    #         self.process_counts[ikey]=0
    #     print '(info) Done define_process.'

    def simple_counter(self,theta_histo):
        """
        input: theta-template histo objects?
        output: R_process table , a dictionary {'Rqq':float,'Rgg':float,etc}
        """
        ## histos object is like this 
        ## {'f_minus': {'qq': <theta_auto.Model.Histogram object at 0xa78f90>, 'gg': <theta_auto.Model.Histogram object at 0x9b4050>, 
        ## 'WJets': <theta_auto.Model.Histogram object at 0x28b20d0>, 'other_bkg': <theta_auto.Model.Histogram object at 0x28b7cd0>, 
        ##'qcd': <theta_auto.Model.Histogram object at 0x28b7e50>}, 'f_plus': { ... } } 
        R_process = OrderedDict()
        counter = {'total':0} 
        for obs in theta_histo: # 'f_minus','f_plus'
            for proc in theta_histo[obs]: # 'qq','gg','WJets','other_bkg','qcd'
                if proc not in R_process:
                    R_process[proc],counter[proc] = 0.,0.
                count = theta_histo[obs][proc].get_value_sum()
                counter[proc] += count
                counter['total'] += count
        total_count = counter['total']
        for proc in R_process:
            # print '(debug) %s = %i'%(proc,counter[proc])
            if 'qq' not in proc:
                R_proc = counter[proc]/total_count
            else:
                R_proc = counter['qq']/(counter['gg']+counter['qq'])
            R_process[proc] = R_proc
            #print 'R_%s = %.3f'%(proc,counter[proc]/total_count)
        R_process['tt'] = 1-sum([val for key,val in R_process.iteritems() if key not in ['qq','gg']])
        return R_process 


    # main function starts here
    def main(self,template_file_path,outdir,useToys=False):

        # IO
        self.template_file = template_file_path
        self.outdir = outdir
        self.defineIO(template_file_path)

        # Define model
        print '(info) Begin defining model.'
        self.model = self.get_model(self.template_file)
        self.model_pars =  self.model.get_parameters('')
        self.model_bins = self.model.get_range_nbins(self.obs)[-1]
        # change range and width for some distributions
        self.resetModel()
        # add model info into results
        self.report_model(self.model,self.txtfile)
        # Report the model in an html file
        par_summary = model_summary(self.model)

        # maximum likelihood fit on data without signal (background-only).
        print '(info) Begin MLE fit.'
        parVals,theta_options = self.mleFit(self.model)
        self.toy_options = theta_options.copy()
        # keep a copy of parVals and theta_options in global var
        self.parVals,self.theta_options = parVals,theta_options
        # print fit result
        self.mle_result_print(parVals)

        # Toy experiments for determination of error of POI
        if useToys:
            self.doToys()

        # Save post fit hists to root file
        self.savePostFit(self.postfit_histos)

        # Write as html+
        report.write_html('%s/htmlout'%self.outdir)
        # close files
        self.txtfile.close()
        self.fout.Close()

    def doToys(self): # checked
        """
        inputs: self.AFB_list
        outputs: fitted toy experiments and plots
        """

        print '(info) Begin toy experiments on %s'%self.toy_params
        nthreads = str(self.ntoys/self.Toys_per_thread)
        # debug
        #nthreads = '20'

        self.toy_options.set('main','n_threads',nthreads)
        #toy_dist =  model.distribution.copy()
        for key,value in self.toy_params.iteritems():
            toy_dist = self.getToyDist() 

            toy_sigma = self.sigma_values[key]
            toy_name = key

            fit_hists=[] # for containing fit results of different input param values
            toy_input,toy_fit_mean,toy_fit_sigma = [],[],[]
            for toy_val in value:
                """ def fitToys(self,toy_param_val,toy_param_sigma,toy_param_name): """
                fit_hist, fit_results, Rqq_input = self.fitToys(toy_dist = toy_dist,toy_param_val=toy_val,toy_param_sigma=toy_sigma,toy_param_name=toy_name)
                # fit_results=[fit_mean,fit_sigma]
                if 'qq' in toy_name:
                    toy_input_val = Rqq_input
                    toy_input.append(Rqq_input)
                else:
                    toy_input_val = toy_val*toy_sigma
                    toy_input.append(toy_input_val)
            
                toy_fit_mean.append(fit_results[0])
                toy_fit_sigma.append(fit_results[1])
                fit_hists.append(fit_hist)
                print '(info) toy exp fit for %s = %.3f is %.3f +/- %.3f'%(toy_name,toy_input_val,fit_results[0],fit_results[1])

            # plot Neyman bands
            # from helper.py
            # def makeTGraphErrors(x,y,y_err,x_err=None,x_title='x',y_title='y',title='TGraph'):
            # makeTGraphErrors is from helper.py
            neyman_plot = makeTGraphErrors(x=toy_input,y=toy_fit_mean,y_err=toy_fit_sigma,x_title='param_input',y_title='param_fit',title='%s Based on %i toy experiments'%(toy_name,self.ntoys))
            neyman_plot.Fit("pol1")
            canv_neyman = ROOT.TCanvas()
            canv_neyman.SetName('%s_neyman'%toy_name)
            self.fout.cd()
            neyman_plot.Draw()
            canv_neyman.SaveAs('%s/%s_neyman.png'%(self.outdir,toy_name))
            canv_neyman.Write()
            neyman_plot.Write()
            # plot pull plots
            self.plotPull(self.fout)
            print '(info) Done all toy experiments!'

    def histogram_filter(self,hname):
        """
        Filter out any histgrams in black list
        """
	if self.shape_sys_gauss_white == 'all': return True
        blacklist = [item for item in self.shape_sys_gauss if item not in self.shape_sys_gauss_white]
        toReturn = True
        for item in blacklist:
            if item in hname:
                toReturn = False
        return toReturn

    def report_model(self,_model,txt_output):
        dist = _model.distribution.distributions
        str_out = ''
        for key in dist:
            str_out += '%s:%s\n'%(key,str(dist[key]))
        print str_out
        txt_output.write(str_out)
        print '(info) Done report_model.'
        

    def get_model(self,template):
        """
        theta statistical model defined here
        input: template.root
        output: model 
        """
        print 'Building theta model from %s'%template
        model = build_model_from_rootfile(template, self.histogram_filter)
        
        # add a minimum MC stat. uncertainty corresponding to +-1 MC event in each bins (esp. empty bins)
        # Need to study how the smoothing affect likelihood
        model.fill_histogram_zerobins()
        
        # Specifying all uncertainties. Internally, this adds a factor exp(lambda * p)
        # where p is the parameter specified as first argument and lambda is the constant
        # in the second argument:
        # check if has qcd temp, if not , not adding qcd_rate sys
        if 'qcd' in model.get_processes('f_minus'):
            self.has_qcd = True
            model.add_lognormal_uncertainty('qcd_rate', math.log(1.2), 'qcd') # gauss, sigma=50%
        else:
            self.has_qcd = False

        # Set parameter ranges
        for p in model.distribution.get_parameters() :
            if p=='AFB' :
                pass
    #            model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=[-1.0,1.0])
            elif p=='R_qq' or p=='wjets_rate' or p=='gg_rate' :
                pass
    #            model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=[-5.0,inf])
            elif p=='qcd_rate':
    #       pass
                model.distribution.set_distribution_parameters(p, range = [-1.0, 1.0]) 
            else :
                d = model.distribution.get_distribution(p)
                if d['typ'] == 'gauss' :
    #       pass
                    model.distribution.set_distribution_parameters(p, range = [-1.0, 1.0])        
        
        # the qcd is derived from data, so do not apply a lumi uncertainty on that:
        for p in model.processes:
            if p == 'qcd': continue
            model.add_lognormal_uncertainty('lumi', math.log(1.045), p)

        # get MC R_process 
        histos = evaluate_prediction(model,model.distribution.get_means(),include_signal = False)
        self.R_process_MC = self.simple_counter(histos)
        self.Rqq_MC = self.R_process_MC['qq']

        print '(info) Done getmodel.'
        return model

    def arrayToStr(self,m):
        """
        Convert a matrix into a nice str representation
        """
        str_rep = ''
        for row in m:
            for item in row:
                str_rep += '%15.3f,'%item
            str_rep += '\n'
        return str_rep

    # Not used anymore  
    def setRange(self):
        """
        partially reset range of some parameters
        """
        for p in ('R_qq','wjets_rate','gg_rate'):
            self.model.distribution.set_distribution_parameters(p,range=[-inf,inf])
        self.model.distribution.set_distribution_parameters('qcd_rate',range=[-inf,inf]) 
        self.model.distribution.set_distribution_parameters('AFB',range=[-AFB_range,AFB_range])

        # for sys
        print '(info) reset range for sys.'
        for p in self.model.distribution.get_parameters():
            d = self.model.distribution.get_distribution(p)
            if d['typ'] == 'gauss' and d['mean'] == 0.0 and d['width'] == 1.0:
                self.model.distribution.set_distribution_parameters(p, range = [-inf, inf])
        print '(info) Done setRange.'

    def resetModel(self):
        """
        partially reset some prior distribution for later use
        """
        # set flat prior for some parameters
        for p in self.flat_param:
            if (not self.has_qcd) and p=='qcd_rate': continue     
            self.model.distribution.set_distribution_parameters(p,width=inf,range=[-100.0,100.0])

        # for Gaussian prior params, set the range to inf
        for p in self.model.distribution.get_parameters():
            d = self.model.distribution.get_distribution(p)
            if d['width'] == 1.0: # i.e, non-flat prior parameters
                self.model.distribution.set_distribution_parameters(p, range = [-inf, inf])
        print '(info) Done resetModel for finer control of model parameter priors.'

    def mle_result_print(self,result,sp='',n=None):
        """
        input: a dict from mle output, with all mle fit param results
        output: write all kinds of mle fit results into output txt file
        """
        str_result = ''
        n = len(result[sp]['__nll'])
        # fit parameters of interests
        for p in self.pois:
            if '__' in p or result[sp].get(p,0)==0: continue
            # n is number of toys to print
            if n is None: n = len(result[sp][p])
            str_result += "%20s =" % p
            sigma_p = self.sigma_values.get(p,1.0)
            # for parameter results
            for i in range(min([n, 10])):
                str_result +=  " %5.3f +- %5.3f \n" % (result[sp][p][i][0]*sigma_p, result[sp][p][i][1]*sigma_p)
        str_result += '\n'
        # other nuisance params
        for p in result[sp]:
            if '__' in p or p in self.pois: continue
            # n is number of toys to print
            if n is None: n = len(result[sp][p])
            str_result += "%20s =" % p
            sigma_p = self.sigma_values.get(p,1.0)
            # for parameter results
            for i in range(min([n, 10])):
                str_result +=  " %5.3f +- %5.3f \n" % (result[sp][p][i][0]*sigma_p, result[sp][p][i][1]*sigma_p)
        # stdev of pars for each experiment
        sigmas = []
        for i in range(min([n, 10])):
            isigma = []
            for par in self.model_pars:
                isigma.append(result[sp][par][i][1])
            sigmas.append(isigma)
        # get Sigma_ij = sigma_i*sigma_j
        sig_matrices = []
        for isigma in sigmas:
            sig_arr = np.array(isigma)
            isig_matrix = np.outer(sig_arr,sig_arr)
            sig_matrices.append(isig_matrix)
        # chi2,nll and cov
        str_result += '-------------------------\n'
        for p in ('__nll','__chi2','__cov'):
            if not p in result[sp]: continue      
            if n is None: n = len(result[sp][p])
            str_result += "%20s =" % p
            for i in range(min([n, 10])):
                if 'cov' not in p:
                    if 'chi2' in p:
                        res = result[sp][p][i]/(self.model_bins*2)
                    else:
                        res = result[sp][p][i]
                    str_result += " %5.3f\n"%res
                else:
                    # This is the cov matrix
                    cov_matrix = result[sp][p][i]
                    str_result += '\nCovariant matrix\n'
                    for item in self.model_pars:
                        str_result += '%15s,'%item
                    str_result += '\n%s\n'%arrayToStr(cov_matrix)
                    # Add correlation matrix, rho_matrix
                    rho_matrix = cov_matrix/sig_matrices[i]
                    str_result += 'Correlation matrix\n'
                    for item in self.model_pars:
                        str_result += '%15s,'%item
                    str_result += '\n%s\n'%arrayToStr(rho_matrix)    

        # add R_process before and after fit with error
        str_result += '---------------------------------------------------------------------------\n\n'
        # first add MC_input values
        prefit_R = '(info) R_process prefit\n'
        prefit_R += self.R_proc_to_str(self.R_process_MC)

        # Next add post fit R_process from mle fit results
        # >>> self.parVals['']['R_qq'][0] = (0.805501721611904, 0.12242412236195577)
        postfit_R = '\n(info) R_process postfit from MLE output\n'
        fit_AFB_mean = self.parVals['']['AFB'][0][0]*self.sigma_values['AFB']
        fit_AFB_err = self.parVals['']['AFB'][0][1]*self.sigma_values['AFB']
        postfit_R += 'AFB = %.3f +/- %.3f\n\n'%(fit_AFB_mean,fit_AFB_err)
        R_tt,R_tt_err_sq = 1.,0.
        fit_R = {}
        temp_str = {}
        for key,val in self.R_process_MC.iteritems():
            param_name = [par for par in self.parVals[''] if key in par]
            if param_name : 
                param_name = param_name[0]
            else:
                continue
            if 'reweight' in param_name: continue
            # print '(debug) param_name %s.'%param_name
            # R_process = R_process_MC*(1+N_sigma*sigma)
            # e.g., R_qq_MC = 0.067, fit 'R_qq' = 0.8 +/- 0.12, sigma_Rqq = 0.8
            # fianl fit R_qq_fit = 0.067(1+0.8*0.8) +/- 0.067*(0.12*0.8)
            fit_R_mean = self.R_postfit[key]
            fit_R_err  = val*(self.parVals[''][param_name][0][1]*self.sigma_values[param_name])
            temp_str[key] = '%s = %.3f +/- %.3f\n'%(param_name,fit_R_mean,fit_R_err)
            fit_R[key] = [fit_R_mean,fit_R_err]
            if 'qq' not in key:
                R_tt -= fit_R_mean
                R_tt_err_sq += fit_R_err*fit_R_err

        R_tt_err = numpy.sqrt(R_tt_err_sq)
        R_qq = fit_R['qq'][0]; R_qq_err = fit_R['qq'][1]
        R_gg = R_tt*(1-R_qq)
        R_gg_err = numpy.sqrt(  numpy.power(R_gg*R_tt_err/R_tt,2)+numpy.power(R_tt*R_qq_err,2) )
        temp_str['gg'] = 'R_gg = %.3f +/- %.3f\n'%( R_gg,R_gg_err )
        temp_str['tt'] = 'R_tt = %.3f +/- %.3f\n'%( R_tt,R_tt_err )
        # Insert fit results with same order as counter
        for key in self.R_process_MC:
            postfit_R += temp_str[key]

        str_result += prefit_R+postfit_R
        # save result to txt
        print str_result
        self.txtfile.write(str_result)
        print '(info) Done mle_result_print'

    def mleFit(self,model):
        """
        input: an theta model object
        output: mle fit param estimations
        """
        # Set some options
        options = Options()
        options.set('global','debug','True')
        options.set('minimizer','strategy','robust')
        options.set('minimizer','minuit_tolerance_factor','10')

        # Do mle fit here
        """
        MLE fit and output postfit templates to a root file
        """
        parVals = mle(model, 'data', 1,with_covariance=True, signal_process_groups = {'': [] },chi2=True,options = options)
        mle_LL = parVals['']['__nll'][0] 
        mle_chi2 = parVals['']['__chi2'][0]

        # Get 1sigma and 2 sigma interval of AFB
        afb_interval = pl_interval(model, 'data', 1,signal_process_groups = {'': [] }, parameter='AFB')
        afb_interval = afb_interval['']
        str_write = 'profile likelihood AFB significance interval:\n'
        for key,value in afb_interval.iteritems():
            str_write += 'confidence level = %.3f\n'%key 
            for item in value: 
                if key==0.0:
                    str_write += '[%.3f]\n'%item
                else:
                    str_write += '[%.3f,%.3f]\n'%(item[0],item[1])

        self.txtfile.write(str_write)
        print str_write

        # NLL scan for AFB
        mle_nllscan = nll_scan(model, 'data', 1, npoints=100, range=[-AFB_range, AFB_range], signal_process_groups = {'': [] }, parameter='AFB',adaptive_startvalues=False)
        mle_nllscan = mle_nllscan[''][0]
        print mle_nllscan 
        # plot nll_scan result
        plotutil.plot(mle_nllscan,'AFB(sigma=%.1f)'%self.AFB_sigma,'NLL','%s/AFB_nll.png'%self.outdir)

        # Get postfit histo object for later plotting
        parameter_values = {}
        for p in self.model.get_parameters([]):
            parameter_values[p] = parVals[''][p][0][0]
        # parameter_values['beta_signal'] = 1.0 # do not scale signal
        # Save postfit templates
        histos = evaluate_prediction(self.model,parameter_values,include_signal = False)
        self.postfit_histos = histos
        # Write R_process postfit into txt file
        R_proc_post = self.simple_counter(histos)
        self.R_postfit = R_proc_post
        R_proc_post = '\n(info) R_process postfit from counting\n'+self.R_proc_to_str(R_proc_post)
        self.txtfile.write(R_proc_post)
        print R_proc_post
        print '(info) Done mleFit.' 
        return parVals,options

    def reject_outliers(self,data, m = 10.):
        data=np.array(data)
        d = np.abs(data - np.median(data))
        mdev = np.median(d)
        s = d/mdev if mdev else 0.
        #print 'mdev=%.2f,median=%.2f'%(mdev,np.median(data))
        return data[s<m]

    def plot_hist(self,data,name='test',xtitle='x',ytitle='Events',title='Histogram'):
        """
        input: a list of data
        output: a canvas
        """
        print '(info) plot_hist for %s'%name
        # first exclude outliers
        data = self.reject_outliers(data)
        data = np.array(data)
        data_mean = data.mean()
        data_std = data.std()
        #print data
        xmin = data_mean-6*data_std
        xmax = data_mean+6*data_std 
        nbins = 50
        hist = ROOT.TH1D(name,'%s;%s;%s'%(title,xtitle,ytitle),nbins,xmin,xmax)
        print "name=%s,mean=%.1f,stdev=%.2f,min=%.1f,max=%.1f"%(name,data_mean,data_std,xmin,xmax)
        for item in data:
            hist.Fill(item)
        # fit with gaussian
        hist_fit = hist.Fit('gaus','sq')
        fit = hist.GetFunction('gaus')
        fit_mean = fit.GetParameter(1)
        fit_sigma = fit.GetParameter(2)
        canv = ROOT.TCanvas()
        canv.SetName('c_%s'%name)
        hist.Draw('hist e')
       
        # output 
        self.fout.cd()
        hist.Write()
        canv.Write()
        return hist,canv,[fit_mean,fit_sigma]

    def getToyDist(self):
        """ 
        get prior dist for toy generations
        fix all nuisance parameters by default
        """
        toy_dist =  self.model.distribution.copy()
        for p in toy_dist.get_parameters() :
            toy_dist.set_distribution(p,typ='gauss',mean=0.0,width=0.0,range=[0.0,0.0])
        return toy_dist

    def dict_to_str(self,idict):
        str_out = ''
        for key,val in idict.iteritems():
            str_out+='%s: %.3f, '%(key,val)
        return str_out

    def R_proc_to_str(self,idict):
        str_out = ''
        for key,val in idict.iteritems():
            str_out+='R_%s = %.3f\n'%(key,val)
        return str_out

############### still working on this!!!! ################
    def fitToys(self,toy_dist,toy_param_val,toy_param_sigma,toy_param_name): # checked
        """
        Fit to toy experiment with input AFB fixed
        """
        #toy_dist =  self.model.distribution.copy()
        # calculate the AFB input in terms of toy_param_sigma of templates. e.g, if template is toy_param_sigma=0.1, AFB_input=-0.6 => new_input = -7

        new_mean = toy_param_val/toy_param_sigma
        toy_dist.set_distribution_parameters(toy_param_name,mean=new_mean,width=0.0,range=[new_mean,new_mean])

        # calculate R_proc before fit
        histos = evaluate_prediction(self.model,toy_dist.get_means(),include_signal = False)
        R_proc_toy = self.simple_counter(histos)
        Rqq_input = R_proc_toy['qq']

        print '\n(info) Toy experiments with %s = %s'%(toy_param_name,toy_param_val)
        print '(info) Input R_process is { %s }'%self.dict_to_str(R_proc_toy)

        # do mle fit for generated toys
        toy_fit = mle(self.model, 'toys:0.0', self.ntoys,with_covariance=False, signal_process_groups = {'': [] },chi2=True,options = self.toy_options,nuisance_prior_toys=toy_dist)
        fit_AFB,fit_chi2 = [],[]
        # self.ntoys = len(fit_AFB)
        all_AFB = toy_fit[''][toy_param_name]
        all_chi2 = toy_fit['']['__chi2']

        ###### need to implement getting post fit R_proc for toy experiments on R_qq only! ######
        # toy_fit['']['lumi'] = [(-0.9563586273699287, 0.8486084490392977), (-0.8390765758997597, 0.89665414376131)] 
        if 'qq' in toy_param_name:
            fit_Rqq = []
            for i in range(len(all_AFB)):
                continue    
                if i%50==1: print '(progress) simple_counter at %i toy'%i
                parameter_values = {}
                for param in self.model_pars:
                    parameter_values[param] = toy_fit[''][param][i][0]
                histos = evaluate_prediction(self.model,parameter_values,include_signal = False)
                R_proc_fit = self.simple_counter(histos)
                fit_Rqq.append(R_proc_fit['qq'])
            # calculate mean and stdev for all toys given current input val
            fit_Rqq = numpy.array(fit_Rqq)
            #mean_and_std = [fit_Rqq.mean(),fit_Rqq.std()]

        # all_AFB is in the form of [(-0.9563586273699287, 0.8486084490392977), (-0.8390765758997597, 0.89665414376131)]
        # which is self.ntoys number of tuples with first as central, second as error
        converted_fit_val = []
        for i in range(len(all_AFB)):
            # get post fit Rqq
            fit_AFB_value = all_AFB[i][0]

            # remember to convert fit AFB central value back to actual AFB, which is independent of choice of sigma
            actual_AFB = all_AFB[i][0]*toy_param_sigma
            converted_fit_val.append(actual_AFB)
            if 'qq' in toy_param_name:
                actual_AFB = self.Rqq_MC*(1+actual_AFB)
            # fit_AFB contains all fit results that has converted to physically meanningful val, independent of sigma
            fit_AFB.append(actual_AFB)
            # convert chi2 with ndof to chi2/ndof
            fit_chi2.append(all_chi2[i]*1.0/self.model_bins)

        converted_fit_val = numpy.array(converted_fit_val)
        print '(debug)  converted_fit_val = %.3f +/- %.3f'%(converted_fit_val.mean(),converted_fit_val.std())

        # make histograms
        if toy_param_val<0:
            postfix = 'minus%ipct'%(abs(toy_param_val)*100)
        else:
            postfix = 'plus%ipct'%(toy_param_val*100)
        hist_AFB,canv_AFB,fit_results_AFB = self.plot_hist(data=fit_AFB,name='%s_%s'%(toy_param_name,postfix),xtitle='%s(fit)'%toy_param_name,title='%s for %i toys, input = %.2f'%(toy_param_name,self.ntoys,toy_param_val))
        hist_chi2,canv_chi2,fit_results_chi2 = self.plot_hist(data=fit_chi2,name='chi2_%s'%postfix,xtitle='chi2/%i'%self.model_bins,title='chi2 for %i toys, AFB_input = %.2f'%(self.ntoys,toy_param_val))

        # compare Theta-fitter output for fit param
        print '(debug) Theta output %s = %.3f +/- %.3f'%(toy_param_name,fit_results_AFB[0],fit_results_AFB[1])

        # finish
        print '(info) Done fitToys for Param %s = %.2f'%(toy_param_name,toy_param_val)
        #    canv_AFB.SaveAs('%s/AFB_toys_%s.png'%(self.outdir,postfix))
        #    canv_chi2.SaveAs('%s/chi2_toys_%s.png'%(self.outdir,postfix))
        if 'qq' in toy_param_name:
            mean_and_std = fit_results_AFB
            return [hist_AFB,hist_chi2],mean_and_std,Rqq_input
        else:
            return [hist_AFB,hist_chi2],fit_results_AFB,Rqq_input


    def plotPull(self,tfile):
        """
        plot all TH1s in the file
        """
        keys = tfile.GetListOfKeys()
        hist_sets = set()
        for ikey in keys:
            if 'TH1' in  ikey.GetClassName() : histname = ikey.GetName()
            if histname in hist_sets or not [ item for item in self.toy_params.keys() if item in histname]:
                continue
            else:
                hist_sets.add(histname)
    #        print '(debug) plotting %s'%histname
            ihist = tfile.Get(histname)
            canv = ROOT.TCanvas()
            ihist.Draw()
            canv.SaveAs('%s/pull_%s.png'%(self.outdir,histname))
        print '(info) Done plotPull.'

    def savePostFit(self,histos):
        # Save original data histograms
        for o in self.model.get_observables():
            histos[o]['DATA'] = self.model.get_data_histogram(o)
        # Output to root file
        write_histograms_to_rootfile(histos,'%s/postfit_histos.root'%self.outdir)#,template_file.split('/')[-1]))
        print '(info) Done savePostFit'

if __name__ == '__main__':
    """
    Parse inputs to get all parameters needed for thetaFitter class
    params: template_file,outdir,AFB_sigma,useToys
    """
    # Input
    argv = sys.argv[2:]
    if len(argv)==0:
        print """
        Usage:
        python ../../theta-auto.py thetaFitter_rate.py templates/final_test.root test 1.0 toy
        """
        sys.exit(1)
    template_file = argv.pop(0)
    # output dir name
    if len(argv)==0:
        outdir = 'runs/test'
    else:
        outdir = 'runs/%s'%argv.pop(0)    
    # AFB_signma def
    if len(argv)==0:
        AFB_sigma = 1.0
    else:
        AFB_sigma=float(argv.pop(0))
    AFB_range = 2.0/AFB_sigma
    argv = ' '.join(argv)
    if 'toy' in argv:
        useToys = True
    else:
        useToys = False

    # create a new class instance
    self = thetaFitter(AFB_sigma=AFB_sigma)
    self.main(template_file_path = template_file,useToys=useToys,outdir='%s_postfit'%outdir)
