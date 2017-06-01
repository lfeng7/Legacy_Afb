import ROOT as rt #TFile, TGraphErrors, TGraphAsymmErrors, TH2D, TLine, TLegend, TCanvas
import numpy as np
import matplotlib.pyplot as plt
import os

rt.gROOT.Macro( os.path.expanduser( '~/rootlogon.C' ) )
rt.gROOT.SetBatch(True)

DEGREE = 1 # degree of polynomials
common_text = 'CMS Simulation, 19.7 fb^{-1} at #sqrt{s} = 8 TeV' 

def main():
    # def make_neyman_plots(x,y2_down,y1_down,y,y1_up,y2_up,par_name,fname,data_fit_central_value=0.1):
    x = np.linspace(-1,1,20)
    y = x
    sigma = 0.05
    y1_down = y-sigma
    y1_up = y+sigma
    y2_down = y-2*sigma
    y2_up = y+2*sigma
    make_neyman_plots(x,y2_down,y1_down,y,y1_up,y2_up,'A_{FB}','AFB',0.06)


#settings
def make_neyman_plots(x,y2_down,y1_down,y,y1_up,y2_up,par_name,fname,title='',data_fit_central_value=0.1):
    """
    input: lists. x_values: input val, y_values: fit vals
    output: root and pdf files
    """
    output_filename = '%s_Neyman.root'%fname

    #output file
    outfile = rt.TFile(output_filename,'recreate')

    #TGraphErrors setup
    n = len(x)
    x_errs   = np.zeros(n,'float')
    x = np.array(x)
    y = np.array(y)
    y1_up = np.array(y1_up)
    y2_up = np.array(y2_up)    
    y1_down = np.array(y1_down)
    y2_down = np.array(y2_down)
    #
    y1_err_up = abs(y1_up-y)
    y1_err_down = abs(y1_down-y)
    y2_err_up = abs(y2_up-y)
    y2_err_down = abs(y2_down-y)

    onesigma_graph = rt.TGraphAsymmErrors(n,x,y,x_errs,x_errs,y1_err_down,y1_err_up)
    onesigma_graph.SetName('one_sigma')
    twosigma_graph = rt.TGraphAsymmErrors(n,x,y,x_errs,x_errs,y2_err_down,y2_err_up)
    twosigma_graph.SetName('two_sigma')

    # debug
    to_print = [x,y2_down,y1_down,y,y1_up,y2_up,y1_err_up,y1_err_down,y2_err_up,y2_err_down]
    for item in to_print:
        continue
        print item
    #closure test plot setup

    #Set graph attributes
    onesigma_graph.SetFillColor(rt.kGreen)
    onesigma_graph.SetLineStyle(2)
    onesigma_graph.SetLineWidth(3)
    onesigma_graph.SetLineColor(rt.kBlack)
    onesigma_graph.SetMarkerStyle(20)
    onesigma_graph.SetMarkerColor(rt.kBlack)

    twosigma_graph.SetFillColor(rt.kYellow)
    twosigma_graph.SetTitle('Neyman Construction for '+par_name)
    twosigma_graph.GetXaxis().SetTitle('Input %s'%par_name)
    twosigma_graph.GetYaxis().SetTitle('Fit %s'%par_name)

    #Interpolate given the central value
    
    # first fit five arrays, corresponding to central, +/- sigma, +/- 2 sigma
    fit_central = np.polyfit(x, y, DEGREE)
    fit1_up = np.polyfit(x,y1_up,DEGREE)
    fit1_down = np.polyfit(x,y1_down,DEGREE)
    fit2_up = np.polyfit(x,y2_up,DEGREE)
    fit2_down = np.polyfit(x,y2_down,DEGREE)

    # make sample plot to show the effect of fit
    y_errs = [y1_up,y1_down,y2_up,y2_down]
    y_errs_labels = ['1 up','1 down','2 up','2 down']
    y_errs_fit = [fit1_up,fit1_down,fit2_up,fit2_down]

    x_range = max(x)-min(x)
    x_low = min(x)-x_range/20.
    x_high = max(x)+x_range/20.
    xp = np.linspace(min(x),max(x), 100)

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(x,y,'o',label='central')
    for i in range(len(y_errs)):
        ax1.plot(x,y_errs[i],'+',label=y_errs_labels[i])
        # make fit line plot
        ax1.plot(xp, np.poly1d(y_errs_fit[i])(xp),'--')
    plt.xlim(x_low,x_high)
    plt.xlabel('Input %s'%par_name)
    plt.ylabel('Fit %s'%par_name)
    plt.legend(loc='best')
    fig1.savefig('%s_Neyman_check.pdf'%fname)
    plt.gcf().clear()

    # print out debug info
    print '(debug) fit_central: %.2f,%.2f'%tuple(fit_central.tolist())

    # then inversely solve for the five values
    # x = (y-b)/k
    data_fit_mean = (data_fit_central_value-fit_central[1])/fit_central[0]
    data_fit_one_sigma_down = (data_fit_central_value - fit1_up[1])/fit1_up[0]
    data_fit_one_sigma_up = (data_fit_central_value - fit1_down[1])/fit1_down[0]
    data_fit_two_sigma_down = (data_fit_central_value - fit2_up[1])/fit2_up[0]
    data_fit_two_sigma_up = (data_fit_central_value - fit2_down[1])/fit2_down[0]

    print_res = '%s:\n'%par_name 
    print_res = 'Central value = %.5f\nminus one sigma = %.5f, plus one sigma = %.5f'%(data_fit_mean,data_fit_one_sigma_down,data_fit_one_sigma_up)
    print_res+='\nminus two sigma = %.5f, plus two sigma = %.5f'%(data_fit_two_sigma_down,data_fit_two_sigma_up)
    print_res+='\nFINAL RESULT: Parameter %s = %.5f + %.5f - %.5f (95%% CL = %.5f)\n'%(par_name,data_fit_mean,data_fit_one_sigma_up-data_fit_mean,abs(data_fit_one_sigma_down-data_fit_mean),data_fit_two_sigma_up)
    print print_res

    #Build the lines to indicate based on the data fit value
    lines = []
    x_low = twosigma_graph.GetXaxis().GetBinLowEdge(twosigma_graph.GetXaxis().GetFirst())
    x_high = twosigma_graph.GetXaxis().GetBinUpEdge(twosigma_graph.GetXaxis().GetLast())
    y_low = twosigma_graph.GetYaxis().GetBinLowEdge(twosigma_graph.GetYaxis().GetFirst())
    lines.append(rt.TLine(x_low,data_fit_central_value,data_fit_one_sigma_up,data_fit_central_value))
    lines.append(rt.TLine(data_fit_one_sigma_down,data_fit_central_value,data_fit_one_sigma_down,y_low))
    lines.append(rt.TLine(data_fit_one_sigma_up,data_fit_central_value,data_fit_one_sigma_up,y_low))
    lines.append(rt.TLine(data_fit_mean,data_fit_central_value,data_fit_mean,y_low))
    for line in lines :
        line.SetLineWidth(3); line.SetLineColor(rt.kRed); line.SetLineStyle(2)

    #Build a legend
    leg = rt.TLegend(0.23,0.62,0.43,0.77)
    leg.AddEntry(onesigma_graph,'Mean values','PL')
    leg.AddEntry(onesigma_graph,'#pm 1 #sigma','F')
    leg.AddEntry(twosigma_graph,'#pm 2 #sigma','F')
    leg.AddEntry(lines[0],'data fit','L')

    #The line to go on the closure test plots
    otherline = rt.TLine(x[0],x[0],x[n-1],x[n-1])
    otherline.SetLineColor(rt.kBlack)
    otherline.SetLineWidth(4)

    #Plot the neyman construction histogram
    neyman_canv = rt.TCanvas('neyman_canv','neyman_canv')
    neyman_canv.cd()
    twosigma_graph.Draw('A E3')
    onesigma_graph.Draw('SAME E3')
    onesigma_graph.Draw('SAME PLX')
    for line in lines :
        line.Draw()
    leg.Draw()

    pt = rt.TLatex(.23,.80,common_text);
    pt.SetNDC(rt.kTRUE);
    pt.SetTextSize(0.03)
    pt.Draw();

    neyman_canv.Write()
    twosigma_graph.Write()
    onesigma_graph.Write()
    neyman_canv.SaveAs('%s_Neyman.pdf'%fname)

    outfile.Close()
    return print_res


if __name__ == "__main__":
    main()
