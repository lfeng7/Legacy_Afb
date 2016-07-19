
execfile("extras.py")
import sys
import math
execfile("common.py")

"""
theta fitting code
"""

epsilon = 1E-4
AFB_sigma = 0.05
AFB_CENTRAL_VALUE = 0

shape_sys_gauss = ['btag_eff_reweight','trigger_reweight','lepID_reweight']

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
    model.add_lognormal_uncertainty('wjets_rate', math.log(1.2), 'WJets') # flat
    # model.add_lognormal_uncertainty('singleT_rate', 0.05, 'singleT')
    # model.add_lognormal_uncertainty('zjets_rate', 0.05, 'zjets')
    model.add_lognormal_uncertainty('other_bkg_rate', math.log(1.05), 'other_bkg') # gauss, sigma=5%
    model.add_lognormal_uncertainty('qq_rate', math.log(1.2), 'qq') # flat
    model.add_lognormal_uncertainty('qcd_rate', math.log(1.2), 'qcd') # gauss, sigma=50%
    model.add_lognormal_uncertainty('gg_rate', math.log(1.2), 'gg') # flat

    # Set parameter ranges
    for p in model.distribution.get_parameters() :
        if p=='AFB' :
            model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=[-1.0,1.0])
        elif p=='qq_rate' or p=='wjets_rate' or p=='gg_rate' :
            model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=[-5.0,inf])
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
#        if p == 'qcd': continue
        model.add_lognormal_uncertainty('lumi', 0.045, p)
    return model

def arrayToStr(m):
    """
    Convert a matrix into a nice str representation
    """
    str_rep = ''
    for row in m:
        for item in row:
            str_rep += '%10.5f,'%item
        str_rep += '\n'
    return str_rep

def setRange():
    """
    partially reset range of some parameters
    """
    for p in ('qq_rate','wjets_rate','gg_rate'):
        model.distribution.set_distribution_parameters(p,range=[-inf,inf])
    model.distribution.set_distribution_parameters('qcd_rate',range=[-inf,inf]) 
    # for sys
    print '(info) reset range for sys.'
    for p in model.distribution.get_parameters():
        d = model.distribution.get_distribution(p)
        if d['typ'] == 'gauss' and d['mean'] == 0.0 and d['width'] == 1.0:
            print p
            model.distribution.set_distribution_parameters(p, range = [-inf, inf])
    print '(info) Done setRange.'

def mle_result_print(result,sp='',n=None):
    str_result = ''
    for p in parVals[sp]:
        # n is number of toys to print
        if n is None: n = len(result[sp][p])
        str_result += "%20s =" % p
        # for parameter results
        for i in range(min([n, 10])):
            if '__' not in p:
                str_result +=  " %5.2f +- %5.2f \n" % (result[sp][p][i][0], result[sp][p][i][1])
            elif 'cov' not in p:
                str_result += " %5.2f\n"%result[sp][p][i][0]
            else:
                # This is the cov matrix
                cov_matrix = result[sp][p][i][0]
                str_result += '\n%s\n'%arrayToStr(cov_matrix)
    print str_result
    with open('results.txt','w') as txtfile:
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
    #mle_cov = parVals['']['__cov'][0]
    # print 'NLL=%.3f'%mle_LL
    # print 'chi2=%.1f'%mle_chi2
    # print parVals

    # Get 1sigma and 2 sigma interval of AFB
    afb_interval = pl_interval(model, 'data', 1,signal_process_groups = {'': [] }, parameter='AFB')
    afb_interval = afb_interval['']
    print afb_interval

    # NLL scan for AFB
    mle_nllscan = nll_scan(model, 'data', 1, npoints=100, range=[-1.0, 1.0], signal_process_groups = {'': [] }, parameter='AFB',adaptive_startvalues=False)
    mle_nllscan = mle_nllscan[''][0]
    # print mle_nllscan 
    # plot nll_scan result
    plotutil.plot(mle_nllscan,'AFB','NLL','AFB_nll.png')

    return parVals

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
    write_histograms_to_rootfile(histos,'postfit_histos_'+template_file.split('/')[-1])

# main function starts here

# Input
argv = sys.argv[2:]
template_file = argv.pop(0)
# Define model
model = get_model(template_file)
setRange()
# Report the model in an html file
par_summary = model_summary(model)
print 'Done get_summary.'
#maximum likelihood fit on data without signal (background-only).
parVals = mleFit(model)
print 'Done mle_fitter.'
# print fit result
mle_result_print(parVals)
# Save post fit hists to root file
savePostFit()

# Write as html
report.write_html('htmlout')



    


