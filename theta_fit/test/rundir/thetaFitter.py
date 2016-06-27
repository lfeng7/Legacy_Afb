
execfile("extras.py")
import sys
execfile("common.py")

"""
theta fitting code
"""

epsilon = 1E-4
AFB_sigma = 0.05
AFB_CENTRAL_VALUE = 0

shape_sys_gauss = ['btag_eff_reweight','trigger_reweight','lepID_reweight']

def get_model(template):
    """
    theta statistical model defined here
    input: template.root
    output: model 
    """
    print 'Building theta model from %s'%template
    model = build_model_from_rootfile(template, include_mc_uncertainties = True)
    
    # add a minimum MC stat. uncertainty corresponding to +-1 MC event in each bins (esp. empty bins)
    # Need to study how the smoothing affect likelihood
    model.fill_histogram_zerobins(epsilon)
    
    # Specifying all uncertainties. Internally, this adds a factor exp(lambda * p)
    # where p is the parameter specified as first argument and lambda is the constant
    # in the second argument:
    model.add_lognormal_uncertainty('wjets_rate', 0.05, 'WJets')
    model.add_lognormal_uncertainty('singleT_rate', 0.05, 'singleT')
    model.add_lognormal_uncertainty('zjets_rate', 0.05, 'zjets')
    model.add_lognormal_uncertainty('qq_rate', 0.05, 'qq')
    model.add_lognormal_uncertainty('gg_rate', 0.05, 'gg')

    # Set shape based morphine parameters
    for p in model.distribution.get_parameters() :
        if p=='AFB' :
            low_afb = (-0.7-AFB_CENTRAL_VALUE)/AFB_sigma
            hi_afb  = (0.7-AFB_CENTRAL_VALUE)/AFB_sigma
            model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=[low_afb,hi_afb])
        elif p=='qq_rate' :
            model.distribution.set_distribution_parameters(p,typ='flat_distribution',range=[-5.0,5.0])
        else :
            d = model.distribution.get_distribution(p)
            if d['typ'] == 'gauss' :
                model.distribution.set_distribution_parameters(p, range = [-5.0, 5.0])        
    
    # the qcd is derived from data, so do not apply a lumi uncertainty on that:
    for p in model.processes:
        if p == 'qcd': continue
        model.add_lognormal_uncertainty('lumi', 0.045, p)
    return model

def mleFit(theta_model):
    """
    Definition
    """
    # Set some options
    options = Options()
    options.set('global','debug','True')
    options.set('minimizer','strategy','robust')
    options.set('minimizer','minuit_tolerance_factor','10')
    # MLE Fit
    parVals = mle(theta_model, 'data', 1, signal_process_groups = {'': [] })
    parameter_values = {}
    for p in theta_model.get_parameters([]):
        parameter_values[p] = parVals[p][0][0]
    # parameter_values['beta_signal'] = 1.0 # do not scale signal
    #Save fit histograms
    histos = evaluate_prediction(theta_model,parameter_values,include_signal = False)
    write_histograms_to_rootfile(histos,'postfit_histos_'+template_file)
    print 'parameter_values: %s',parameter_values

# main function starts here
# Input
argv = sys.argv[2:]
template_file = argv.pop(0)
# Define model
model = get_model(template_file)
# Report the model in an html file
model_summary(model)
print 'Done get_summary.'

#maximum likelihood fit on data without signal (background-only).
# mleFit(model)
# print 'Done mle_fitter.'

# Do mle fit here
"""
MLE fit and output postfit templates to a root file
"""
parVals = mle_print(model, 'data', 1, signal_process_groups = {'': [] })
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

# Write as html
report.write_html('htmlout')



    


