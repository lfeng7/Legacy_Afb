execfile("extras.py")
import sys
import math
import os
import copy
import ROOT
execfile("helper.py")
execfile("common.py")
ROOT.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )

"""
theta fitting code
"""

epsilon = 1E-8
AFB_CENTRAL_VALUE = 0

#shape_sys_gauss = []
shape_sys_gauss = ['btag_eff_reweight','lepID_reweight']
#shape_sys_gauss = ['btag_eff_reweight','trigger_reweight','lepID_reweight']
flat_param = ['AFB','R_qq','R_WJets']
obs = 'f_minus'
pois = ['AFB','R_qq','R_other_bkg','R_WJets','qcd_rate']

# define sigma value here for all parameters
sigma_values = {}
sigma_values['AFB']=1.0
sigma_values['qcd_rate']=0.2
sigma_values['R_qq']=0.8
sigma_values['R_other_bkg']=0.8
sigma_values['R_WJets']=0.8
sigma_values['lumi']=0.045


ntoys = 2000
Toys_per_thread = 100

# define AFB toy params
toy_param = 'AFB'
range_toy_param = [-1,1.01]
AFB_toy_step = 0.5

# define Rqq toy params
#toy_param = 'R_qq'
#range_toy_param = [-4.0,4.01]
#AFB_toy_step = 2.0

AFB_list = []
afb_tmp = range_toy_param[0]
while afb_tmp<range_toy_param[1]:
    AFB_list.append(afb_tmp)
    afb_tmp += AFB_toy_step
#AFB_list=[-1.0,-0.5,-0.2,0,0.2,0.5,1.0]
#AFB_list = [-0.4,0,0.4]
#AFB_list = [-50,-20,-10,0,10,20,50]#-5,-2,-0.5,0,0.5,2,5,10]#,0,0.3,0.7]


def histogram_filter(hname):
    """
    Filter out any histgrams in black list
    """
    blacklist = []
    blacklist += shape_sys_gauss
    toReturn = True
    for item in blacklist:
        if item in hname:
            toReturn = False
    return toReturn

def report_model(_model,txt_output):
    dist = _model.distribution.distributions
    str_out = ''
    for key in dist:
        str_out += '%s:%s\n'%(key,str(dist[key]))
    print str_out
    txt_output.write(str_out)
    print '(info) Done report_model.'
    

def get_model(template):
    """
    theta statistical model defined here
    input: template.root
    output: model 
    """
    print 'Building theta model from %s'%template
    model = build_model_from_rootfile(template, histogram_filter)
    
    # add a minimum MC stat. uncertainty corresponding to +-1 MC event in each bins (esp. empty bins)
    # Need to study how the smoothing affect likelihood
    model.fill_histogram_zerobins()
    
    # Specifying all uncertainties. Internally, this adds a factor exp(lambda * p)
    # where p is the parameter specified as first argument and lambda is the constant
    # in the second argument:
    model.add_lognormal_uncertainty('qcd_rate', math.log(1.2), 'qcd') # gauss, sigma=50%

    # Set parameter ranges
    for p in model.distribution.get_parameters() :
        if p=='AFB' :
            pass
#            model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=[-1.0,1.0])
        elif p=='R_qq' or p=='wjets_rate' or p=='gg_rate' :
            pass
#            model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=[-5.0,inf])
        elif p=='qcd_rate':
#	    pass
            model.distribution.set_distribution_parameters(p, range = [-1.0, 1.0]) 
        else :
            d = model.distribution.get_distribution(p)
            if d['typ'] == 'gauss' :
#		pass
                model.distribution.set_distribution_parameters(p, range = [-1.0, 1.0])        
    
    # the qcd is derived from data, so do not apply a lumi uncertainty on that:
    for p in model.processes:
        if p == 'qcd': continue
        model.add_lognormal_uncertainty('lumi', math.log(1.045), p)
    print '(info) Done getmodel.'
    return model

def arrayToStr(m):
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
def setRange():
    """
    partially reset range of some parameters
    """
    for p in ('R_qq','wjets_rate','gg_rate'):
        model.distribution.set_distribution_parameters(p,range=[-inf,inf])
    model.distribution.set_distribution_parameters('qcd_rate',range=[-inf,inf]) 
    model.distribution.set_distribution_parameters('AFB',range=[-AFB_range,AFB_range])

    # for sys
    print '(info) reset range for sys.'
    for p in model.distribution.get_parameters():
        d = model.distribution.get_distribution(p)
        if d['typ'] == 'gauss' and d['mean'] == 0.0 and d['width'] == 1.0:
            model.distribution.set_distribution_parameters(p, range = [-inf, inf])
    print '(info) Done setRange.'

def resetModel():
    """
    partially reset some prior distribution for later use
    """
    # set flat prior for some parameters
    for p in flat_param:
        model.distribution.set_distribution_parameters(p,width=inf,range=[-100.0,100.0])

    # for Gaussian prior params, set the range to inf
    for p in model.distribution.get_parameters():
        d = model.distribution.get_distribution(p)
        if d['width'] == 1.0: # i.e, non-flat prior parameters
            model.distribution.set_distribution_parameters(p, range = [-inf, inf])
    print '(info) Done resetModel for finer control of model parameter priors.'

def mle_result_print(result,sp='',n=None):
    str_result = ''
    n = len(result[sp]['__nll'])
    # fit parameters of interests
    for p in pois:
        if '__' in p or result[sp].get(p,0)==0: continue
        # n is number of toys to print
        if n is None: n = len(result[sp][p])
        str_result += "%20s =" % p
        sigma_p = sigma_values.get(p,1.0)
        # for parameter results
        for i in range(min([n, 10])):
            str_result +=  " %5.3f +- %5.3f \n" % (result[sp][p][i][0]*sigma_p, result[sp][p][i][1]*sigma_p)
    str_result += '\n'
    # other nuisance params
    for p in result[sp]:
        if '__' in p or p in pois: continue
        # n is number of toys to print
        if n is None: n = len(result[sp][p])
        str_result += "%20s =" % p
        sigma_p = sigma_values.get(p,1.0)
        # for parameter results
        for i in range(min([n, 10])):
            str_result +=  " %5.3f +- %5.3f \n" % (result[sp][p][i][0]*sigma_p, result[sp][p][i][1]*sigma_p)
    # stdev of pars for each experiment
    sigmas = []
    for i in range(min([n, 10])):
        isigma = []
        for par in model_pars:
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
		    res = result[sp][p][i]/(model_bins*2)
		else:
		    res = result[sp][p][i]
                str_result += " %5.3f\n"%res
            else:
                # This is the cov matrix
                cov_matrix = result[sp][p][i]
                str_result += '\nCovariant matrix\n'
                for item in model_pars:
                    str_result += '%15s,'%item
                str_result += '\n%s\n'%arrayToStr(cov_matrix)
                # Add correlation matrix, rho_matrix
                rho_matrix = cov_matrix/sig_matrices[i]
                str_result += 'Correlation matrix\n'
                for item in model_pars:
                    str_result += '%15s,'%item
                str_result += '\n%s\n'%arrayToStr(rho_matrix)                
    # save result to txt
    print str_result
    txtfile.write(str_result)
    print '(info) Done mle_result_print'

def mleFit(theta_model):
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

    txtfile.write(str_write)
    print str_write

    # NLL scan for AFB
    mle_nllscan = nll_scan(model, 'data', 1, npoints=100, range=[-AFB_range, AFB_range], signal_process_groups = {'': [] }, parameter='AFB',adaptive_startvalues=False)
    mle_nllscan = mle_nllscan[''][0]
    # print mle_nllscan 
    # plot nll_scan result
    plotutil.plot(mle_nllscan,'AFB(sigma=%.1f)'%AFB_sigma,'NLL','%s/AFB_nll.png'%outdir)

    print '(info) Done mleFit.' 
    return parVals,options

def reject_outliers(data, m = 10.):
    data=np.array(data)
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    #print 'mdev=%.2f,median=%.2f'%(mdev,np.median(data))
    return data[s<m]

def plot_hist(data,name='test',xtitle='x',ytitle='Events',title='Histogram'):
    """
    input: a list of data
    output: a canvas
    """
    print '(info) plot_hist for %s'%name
    # first exclude outliers
    data = reject_outliers(data)
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
    fout.cd()
    hist.Write()
    canv.Write()
    return hist,canv,[fit_mean,fit_sigma]

def getToyDist():
    """ 
    get prior dist for toy generations
    fix all nuisance parameters by default
    """
    toy_dist =  model.distribution.copy()
    for p in toy_dist.get_parameters() :
        toy_dist.set_distribution(p,typ='gauss',mean=0.0,width=0.0,range=[0.0,0.0])
    return toy_dist
    

def fitToys(AFB_toy):
    """
    Fit to toy experiment with input AFB fixed
    """
    #toy_dist =  model.distribution.copy()
    # calculate the AFB input in terms of AFB_sigma of templates. e.g, if template is AFB_sigma=0.1, AFB_input=-0.6 => new_input = -7
    new_mean = AFB_toy/AFB_sigma
    toy_dist.set_distribution_parameters(toy_param,mean=new_mean,width=0.0,range=[new_mean,new_mean])
    toy_fit = mle(model, 'toys:0.0', ntoys,with_covariance=False, signal_process_groups = {'': [] },chi2=True,options = toy_options,nuisance_prior_toys=toy_dist)
    fit_AFB,fit_chi2 = [],[]
    # ntoys = len(fit_AFB)
    all_AFB = toy_fit[''][toy_param]
    all_chi2 = toy_fit['']['__chi2']
    # all_AFB is in the form of [(-0.9563586273699287, 0.8486084490392977), (-0.8390765758997597, 0.89665414376131)]
    # which is ntoys number of tuples with first as central, second as error
    for i in range(len(all_AFB)):
        # remember to convert fit AFB central value back to actual AFB
        actual_AFB = all_AFB[i][0]*AFB_sigma
        fit_AFB.append(actual_AFB)
        # convert chi2 with ndof to chi2/ndof
        fit_chi2.append(all_chi2[i]*1.0/model_bins)
    # make histograms
    if AFB_toy<0:
        postfix = 'minus%ipct'%(abs(AFB_toy)*100)
    else:
        postfix = 'plus%ipct'%(AFB_toy*100)
    hist_AFB,canv_AFB,fit_results_AFB = plot_hist(data=fit_AFB,name='%s_%s'%(toy_param,postfix),xtitle='%s(fit)'%toy_param,title='%s for %i toys, input = %.2f'%(toy_param,ntoys,new_mean))
    hist_chi2,canv_chi2,fit_results_chi2 = plot_hist(data=fit_chi2,name='chi2_%s'%postfix,xtitle='chi2/%i'%model_bins,title='chi2 for %i toys, AFB_input = %.2f'%(ntoys,new_mean))

    # finish
    print '(info) Done fitToys for Param %s = %.2f'%(toy_param,AFB_toy)
#    canv_AFB.SaveAs('%s/AFB_toys_%s.png'%(outdir,postfix))
#    canv_chi2.SaveAs('%s/chi2_toys_%s.png'%(outdir,postfix))
    return [hist_AFB,hist_chi2],fit_results_AFB


def plotPull(tfile):
    """
    plot all TH1s in the file
    """
    keys = tfile.GetListOfKeys()
    hist_sets = set()
    for ikey in keys:
        if 'TH1' in  ikey.GetClassName() : histname = ikey.GetName()
        if histname in hist_sets or toy_param not in histname:
            continue
        else:
            hist_sets.add(histname)
#        print '(debug) plotting %s'%histname
        ihist = tfile.Get(histname)
        canv = ROOT.TCanvas()
        ihist.Draw()
        canv.SaveAs('%s/pull_%s.png'%(outdir,histname))
    print '(info) Done plotPull.'

def savePostFit(parVals):
    parameter_values = {}
    for p in model.get_parameters([]):
        parameter_values[p] = parVals[''][p][0][0]
    # parameter_values['beta_signal'] = 1.0 # do not scale signal
    # Save postfit templates
    histos = evaluate_prediction(model,parameter_values,include_signal = False)
    # Save original data histograms
    for o in model.get_observables():
        histos[o]['DATA'] = model.get_data_histogram(o)
    # Output to root file
    write_histograms_to_rootfile(histos,'%s/postfit_histos_%s'%(outdir,template_file.split('/')[-1]))
    print '(info) Done savePostFit'

# main function starts here

# Input
argv = sys.argv[2:]
if len(argv)==0:
    print """
    Usage:
    python ../../theta-auto.py thetaFitter.py templates/afb_one.root linear_interpolate_AFB_one 1.0 toy
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
    AFB_sigma=0.5
else:
    AFB_sigma=float(argv.pop(0))
AFB_range = 2.0/AFB_sigma
argv = ' '.join(argv)
if 'toy' in argv:
    useToys = True
else:
    useToys = False

# output
if not os.path.exists(outdir):
    os.mkdir(outdir)
    print '(info) Creating new dir %s.'%outdir
txtfile = open('%s/results.txt'%outdir,'w')
fout = ROOT.TFile('%s/plots.root'%outdir,'recreate')

# Define model
print '(info) Begin defining model.'
model = get_model(template_file)
model_pars =  model.get_parameters('')
model_bins = model.get_range_nbins(obs)[-1]
# debug
resetModel()
# add model info into results
report_model(model,txtfile)
# Report the model in an html file
par_summary = model_summary(model)

# maximum likelihood fit on data without signal (background-only).
print '(info) Begin MLE fit.'
parVals,theta_options = mleFit(model)
toy_options = theta_options.copy()
# print fit result
mle_result_print(parVals)

# Toy experiments for determination of error of POI
if useToys:
    print '(info) Begin toy experiments on AFB.'
    nthreads = str(ntoys/Toys_per_thread)
    toy_options.set('main','n_threads',nthreads)
    #toy_dist =  model.distribution.copy()
    toy_dist = getToyDist() 
    fit_hists=[]
    toy_AFB_input,toy_AFB_fit_mean,toy_AFB_fit_sigma = [],[],[]
    for AFB in AFB_list:
        fit_hist, fit_results = fitToys(AFB)
        # fit_results=[fit_mean,fit_sigma]
        toy_AFB_input.append(AFB)
        toy_AFB_fit_mean.append(fit_results[0])
        toy_AFB_fit_sigma.append(fit_results[1])
        fit_hists.append(fit_hist)
    # plot Neyman bands
    # from helper.py
    # def makeTGraphErrors(x,y,y_err,x_err=None,x_title='x',y_title='y',title='TGraph'):
    neyman_plot = makeTGraphErrors(x=toy_AFB_input,y=toy_AFB_fit_mean,y_err=toy_AFB_fit_sigma,x_title='param_input',y_title='param_fit',title='%s Based on %i toy experiments'%(toy_param,ntoys))
    neyman_plot.Fit("pol1")
    canv_neyman = ROOT.TCanvas()
    canv_neyman.SetName('%s_neyman'%toy_param)
    fout.cd()
    neyman_plot.Draw()
    canv_neyman.SaveAs('%s/%s_neyman.png'%(outdir,toy_param))
    canv_neyman.Write()
    neyman_plot.Write()
    # plot pull plots
    plotPull(fout)
    print '(info) Done all toy experiments!'

# Save post fit hists to root file
savePostFit(parVals)

# Write as html
report.write_html('%s/htmlout'%outdir)
# close files
txtfile.close()
fout.Close()
