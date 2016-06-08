
execfile("extras.py")
import sys
execfile("common.py")

"""
theta fitting code
"""

epsilon = 1E-10 

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
    model.add_lognormal_uncertainty('vjets_rate', 0.05, 'wjets')

    
    # the qcd is derived from data, so do not apply a lumi uncertainty on that:
    for p in model.processes:
        if p == 'qcd': continue
        model.add_lognormal_uncertainty('lumi', 0.045, p)
    return model

def mleFit(theta_model):
    """
    Definition
    """
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
# model_summary(model)
# report.write_html('htmlout')
# print 'Done get_summary.'

#maximum likelihood fit on data without signal (background-only).
# mleFit(model)
# print 'Done mle_fitter.'

# Do mle fit here
"""
Here's the testing code of mle fit to understand how theta works
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

# print 'parameter_values: ',parameter_values



    


