# Take in the selection output ttree, do reco reconstruction via kinfit

from Legacy_Afb.Tools.ttbar_utility import *
from Legacy_Afb.Tools.angles_tools import *
import glob
from optparse import OptionParser
import ROOT
from ROOT import TMinuit,Long,Double	
from array import array
import math
import numpy as np


# Some preset constants
data_lumi = 19748 
csvm = 0.679 
top_mass = 173.0


# from angles_tools import *


class gen_fitter(object):
    """docstring for gen_fitter"""
    def __init__(self, ttree,nevts,weighted):
        super(gen_fitter, self).__init__()
        self.ttree = ttree
        self.nevts = nevts
        self.weighted = weighted
    
    def fit(self):
        """
        Fit for alpha and AFB using gen fqq(c*,Mtt)
        """  
        #now we've got to do the fits for the "deweighting" alpha and epsilon values
        print 'fitting for alpha/epsilon values... '
        #Define the Minuit instance we'll use
        alpha_minuit = TMinuit(2); alpha_minuit.SetFCN(self.alpha_fcn)
        #miscellaneous minuit stuff
        ierflag = Long(1); arglist = array('d',2*[-1])
        arglist[0]=100000.
        arglist[1]=100000.
        start_val = 0.0
        error_up = 1.0

        #add parameter
        """
        void TMinuit::mnparm    (   Int_t   k1, TString     cnamj, Double_t    uk, Double_t    wk, 
        Double_t    a, Double_t    b, Int_t &     ierflg )   

        K (external) parameter number
        CNAMK parameter name
        UK starting value
        WK starting step size or uncertainty
        A, B lower and upper physical parameter limits and sets up (updates) the parameter lists. Output:
        IERFLG=0 if no problems
        >0 if MNPARM unable to implement definition
        """
        alpha_minuit.mnparm(0,'alpha',start_val,0.1,0,0,ierflag)
        alpha_minuit.mnparm(1,'AFB',start_val,0.1,0,0,ierflag)

        #minimize
        arglist[0] = 3.
        alpha_minuit.mnexcm('SET PRINT', arglist, 1,ierflag); alpha_minuit.mnexcm('SET NOWARNINGS',arglist,1,ierflag)
        # MIGRAD
        alpha_minuit.SetErrorDef(error_up)
        arglist[0] = 100000.
        alpha_minuit.mnexcm('MIGRAD',arglist,0,ierflag)
        # Minos
        alpha_minuit.SetErrorDef(error_up)
        arglist[0] = 0.
        alpha_minuit.mnexcm('MINOS',arglist,0,ierflag)

        #get back the fitted alpha value
        fitted_alpha=Double(0.0); fitted_alpha_err=Double(0.0)
        fitted_AFB=Double(0.0); fitted_AFB_err=Double(0.0)
        alpha_minuit.GetParameter(0,fitted_alpha,fitted_alpha_err)
        alpha_minuit.GetParameter(1,fitted_AFB,fitted_AFB_err)
        # print out
        result = 'Fit AFB = %.3f +/- %.3f\n'%(fitted_AFB,fitted_AFB_err)
        result += 'Fit alpha = %.3f +/- %.3f\n'%(fitted_alpha,fitted_alpha_err)
        print result
        return result


    def alpha_fcn(self,npar, deriv, f, par, flag):

        alpha = par[0]
        AFB = par[1]
        NLL=0. # the Negative log likelihood for minimization
        nevts = self.ttree.GetEntries()

        N_bad_evts = 0
        #loop over entries
        for i in range(nevts) :
            # stop at nevts
            if i==self.nevts: continue
            self.ttree.GetEntry(i)
            cstar = self.ttree.cos_theta_mc
            Mtt  = self.ttree.ttbar_mass_mc
            if self.weighted:
                gen_w = self.ttree.gen_weight
                if gen_w>0: 
                    gen_w = 1.
                else:
                    gen_w = -1.
            else:
                gen_w = 1.
            beta_tmp = 1-4*(top_mass/Mtt)**2
            if beta_tmp < 0 : continue
            beta = np.sqrt(beta_tmp)
            # calc likelihood , checked with formular (45)
            f_L = (2+beta**2*cstar**2-beta**2+alpha*(1-beta**2*cstar**2))/\
                    2*(2-2*beta**2/3+alpha*(1-beta**2/3)) + AFB*cstar
            if f_L>0 :
                NLL+= -2.0*np.log(f_L)*gen_w
            else :
                N_bad_evts += 1
                continue
                NLL+= 10000000000.
        f[0] = NLL

        print "f = %.5f, alpha = %.5f, AFB = %.5f, N_bad_evts=%i, flag = %i, nevts = %i"%\
                (f[0], alpha, AFB, N_bad_evts, flag, nevts);


if __name__=='__main__':


    # Job steering

    # Input inputFiles to use. This is in "glob" format, so you can use wildcards.
    # If you get a "cannot find file" type of error, be sure to use "\*" instead
    # of "*" to make sure you don't confuse the shell. 

    parser = OptionParser()

    parser.add_option('--inputfiles', metavar='F', type='string', action='store',
                      default = "",
                      dest='inputfiles',
                      help='Input files')

    parser.add_option('--nevts', metavar='F', type='int', action='store',
                      default = 1000,
                      dest='nevts',
                      help='number of events to run')


    parser.add_option('--weighted', metavar='F', action='store_true',
                      dest='weighted',
                      help='if weighted(amcnlo) or not (powheg)')


    parser.add_option('--verbose', metavar='F', type='string', action='store',
                      default = 'no',
                      dest='verbose',
                      help='If you want more information than you usually need.')


    (options, args) = parser.parse_args()

    argv = []

    tfile = ROOT.TFile(options.inputfiles)
    ttree = tfile.Get('angles')
    nevts = options.nevts
    weighted = options.weighted

    new_fitter = gen_fitter(ttree,nevts,weighted)
    result = new_fitter.fit()
 
    fname = options.inputFiles.split('/')[-1].split('.root')[0]
    outfile = open('alpha_fit_results/%s_result.txt'%fname,'w')
    outfile.write('Fit result for %s\n%s'%(fname,result))
    outfile.close()
