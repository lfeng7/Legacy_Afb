#define angles_data_cxx
#include "angles_data.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TLegend.h>
#include <THStack.h>
#include <float.h> //needed for DBL_EPSILON

/////////////////////////////////////////////////////////////////////////
//     This file used during fitting to loop over all data events      //
//     Loop function has been modified to return a log likelihood      //
/////////////////////////////////////////////////////////////////////////

//Name of file holding complete chrage-separated distributions
const char* histo_filename = "aggregated_distributions.root";
const char* template_filename = "templates.root";
//Histograms to be gotten from inside the file
TH3D* fqqs_plus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_xi_plus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_delta_plus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_plus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_xi_plus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_delta_plus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_minus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_xi_minus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_delta_minus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_minus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_xi_minus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_delta_minus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_plus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_xi_plus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_delta_plus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_plus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_xi_plus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_delta_plus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_minus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_xi_minus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqs_delta_minus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_minus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_xi_minus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fqqa_delta_minus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fgg_plus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fbck_plus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* WJets_plus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fgg_minus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fbck_minus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* WJets_minus_4jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fgg_plus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fbck_plus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* WJets_plus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fgg_minus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* fbck_minus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH3D* WJets_minus_5jet = (TH3D*)malloc(sizeof(TH3D*));
TH1D* gtt_4jet = (TH1D*)malloc(sizeof(TH1D*));
TH1D* gbk_4jet = (TH1D*)malloc(sizeof(TH1D*));
TH1D* gtt_5jet = (TH1D*)malloc(sizeof(TH1D*));
TH1D* gbk_5jet = (TH1D*)malloc(sizeof(TH1D*));
//Histogram limits
int nbinsx, nbinsy, nbinsz, nbins_lnL;
int iterations = 0;
double x_low, x_high, y_low, y_high, z_low, z_high, lnL_low, lnL_high;
double xbinwidth, ybinwidth, zbinwidth;
//F constants
double F_xi_comb = 0.0;
double F_delta_comb = 0.0;
double F_xi_4jet = 0.0;
double F_xi_5jet = 0.0;
double F_delta_4jet = 0.0;
double F_delta_5jet = 0.0;
// Number of tagged bs required
int num_b_tag_cuts = 2;
// lnL cut
double lnL_cut = 100;

//Main loop function: returns a log likelihood for the function over all events
//COMBINED CASE
double angles_data::Loop(double Rqqbar, double Rbck, double RWJets, double xi, double delta, double Afb, int plot)
{
	//Set up loop
	if (fChain == 0) return 0.0;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	//ln(L) which will be returned
	double logL = 0.0;

	int count_fit = 0;
	int count_not_fit = 0;
	int count_almost_fit = 0;
	TH1D *event_likelihoods = (TH1D*)malloc(sizeof(TH1D*));
	if (iterations == 0 || plot==1)
		event_likelihoods = new TH1D("event_likelihoods","Event -2*Log(L); -2*ln(L)",100,0,50);
	double min_ev_lnL = 1000000.;
	double max_ev_lnL = 0.;
	//MAIN LOOP
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;

		//If Feynman x is negative, flip it.
		if (Feynman_x<0.0) 
			Feynman_x = -1.0*Feynman_x;
		//cut on whether the event is in range of the training set
		if (cos_theta_cs>=x_low && cos_theta_cs<=x_high && Feynman_x>=y_low && Feynman_x<=y_high && ttbar_mass>=z_low && ttbar_mass<=z_high && n_bTags >= num_b_tag_cuts && ln_L < lnL_cut) {
			//use parameters of event to look up values of functions from histograms (depending on lepton charge and jet number)
			double ev_fqqs = 0.0;
			double ev_fqqs_xi = 0.0;
			double ev_fqqs_delta = 0.0;
			double ev_fqqa = 0.0;
			double ev_fqqa_xi = 0.0;
			double ev_fqqa_delta = 0.0;
			double ev_fgg = 0.0;
			double ev_fbck = 0.0;
			double ev_WJets = 0.0;
			double ev_gtt = 0.0;
			double ev_gbk = 0.0;
			if (Q_l > 0) {
				if (n_valid_jets == 4) {
					ev_fqqs = fqqs_plus_4jet->GetBinContent(fqqs_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_xi = fqqs_xi_plus_4jet->GetBinContent(fqqs_xi_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_delta = fqqs_delta_plus_4jet->GetBinContent(fqqs_delta_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa = fqqa_plus_4jet->GetBinContent(fqqa_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_xi = fqqa_xi_plus_4jet->GetBinContent(fqqa_xi_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_delta = fqqa_delta_plus_4jet->GetBinContent(fqqa_delta_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fgg  = fgg_plus_4jet->GetBinContent(fgg_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fbck = fbck_plus_4jet->GetBinContent(fbck_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_WJets = WJets_plus_4jet->GetBinContent(WJets_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_gtt = gtt_4jet->GetBinContent(gtt_4jet->FindFixBin(ln_L));
					ev_gbk = gbk_4jet->GetBinContent(gbk_4jet->FindFixBin(ln_L));
				}
				else if (n_valid_jets == 5) {
					ev_fqqs = fqqs_plus_5jet->GetBinContent(fqqs_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_xi = fqqs_xi_plus_5jet->GetBinContent(fqqs_xi_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_delta = fqqs_delta_plus_5jet->GetBinContent(fqqs_delta_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa = fqqa_plus_5jet->GetBinContent(fqqa_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_xi = fqqa_xi_plus_5jet->GetBinContent(fqqa_xi_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_delta = fqqa_delta_plus_5jet->GetBinContent(fqqa_delta_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fgg  = fgg_plus_5jet->GetBinContent(fgg_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fbck = fbck_plus_5jet->GetBinContent(fbck_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_WJets = WJets_plus_5jet->GetBinContent(WJets_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_gtt = gtt_5jet->GetBinContent(gtt_5jet->FindFixBin(ln_L));
					ev_gbk = gbk_5jet->GetBinContent(gbk_5jet->FindFixBin(ln_L));
				}
				else {
					printf("JET NUMBER INFORMATION NOT AVAILABLE: EVENT SKIPPED!\n");
					continue;
				}
			}
			else if (Q_l < 0) {
				if (n_valid_jets == 4) {
					ev_fqqs = fqqs_minus_4jet->GetBinContent(fqqs_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_xi = fqqs_xi_minus_4jet->GetBinContent(fqqs_xi_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_delta = fqqs_delta_minus_4jet->GetBinContent(fqqs_delta_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa = fqqa_minus_4jet->GetBinContent(fqqa_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_xi = fqqa_xi_minus_4jet->GetBinContent(fqqa_xi_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_delta = fqqa_delta_minus_4jet->GetBinContent(fqqa_delta_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fgg  = fgg_minus_4jet->GetBinContent(fgg_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fbck = fbck_minus_4jet->GetBinContent(fbck_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_WJets = WJets_minus_4jet->GetBinContent(WJets_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_gtt = gtt_4jet->GetBinContent(gtt_4jet->FindFixBin(ln_L));
					ev_gbk = gbk_4jet->GetBinContent(gbk_4jet->FindFixBin(ln_L));
				}
				else if (n_valid_jets == 5) {
					ev_fqqs = fqqs_minus_5jet->GetBinContent(fqqs_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_xi = fqqs_xi_minus_5jet->GetBinContent(fqqs_xi_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_delta = fqqs_delta_minus_5jet->GetBinContent(fqqs_delta_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa = fqqa_minus_5jet->GetBinContent(fqqa_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_xi = fqqa_xi_minus_5jet->GetBinContent(fqqa_xi_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_delta = fqqa_delta_minus_5jet->GetBinContent(fqqa_delta_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fgg  = fgg_minus_5jet->GetBinContent(fgg_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fbck = fbck_minus_5jet->GetBinContent(fbck_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_WJets = WJets_minus_5jet->GetBinContent(WJets_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_gtt = gtt_5jet->GetBinContent(gtt_5jet->FindFixBin(ln_L));
					ev_gbk = gbk_5jet->GetBinContent(gbk_5jet->FindFixBin(ln_L));
				}
				else {
					printf("JET NUMBER INFORMATION NOT AVAILABLE: EVENT SKIPPED!\n");
					continue;
				}
			}
			else {
				printf("LEPTON CHARGE INFORMATION NOT AVAILABLE: EVENT SKIPPED!\n");
				continue;
			}
			//calculate L for this one event
			// double ev_L_bck = Rbck*ev_fbck*ev_gbk;
			// double ev_L_WJets = RWJets*ev_WJets*ev_gbk;
			// double ev_L_gg = (1.0-Rbck-RWJets)*ev_gtt*(1.0-Rqqbar)*ev_fgg;
			// double ev_L_qqs = (1.0-Rbck-RWJets)*ev_gtt*Rqqbar*((1.0)/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(ev_fqqs+xi*ev_fqqs_xi+delta*ev_fqqs_delta);
			// double ev_L_qqa = Afb*(1.0-Rbck-RWJets)*ev_gtt*Rqqbar*((1.0)/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(ev_fqqa+xi*ev_fqqa_xi+delta*ev_fqqa_delta);
			double ev_L_bck = Rbck*ev_fbck;
			double ev_L_WJets = RWJets*ev_WJets;
			double ev_L_gg = (1.0-Rbck-RWJets)*(1.0-Rqqbar)*ev_fgg;
			double ev_L_qqs = (1.0-Rbck-RWJets)*Rqqbar*((1.0)/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(ev_fqqs+xi*ev_fqqs_xi+delta*ev_fqqs_delta);
			double ev_L_qqa = Afb*(1.0-Rbck-RWJets)*Rqqbar*((1.0)/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(ev_fqqa+xi*ev_fqqa_xi+delta*ev_fqqa_delta);
			double ev_L = ev_L_bck+ev_L_WJets+ev_L_gg+ev_L_qqs+ev_L_qqa;
			//make sure that there were a nonzero number of training set events
			if (ev_fgg!=0 || ev_fqqs!=0 || ev_fbck!=0 || ev_WJets != 0) {
				if (ev_L<=0) {
					printf("Parameters are in a bad spot: event likelihood = %.12f\n",ev_L);
					printf("	costheta = %.4f,   x_F = %.4f,   M_tt = %.4f,   nJets = %d,   Q_l = %d\n", cos_theta_cs, Feynman_x, ttbar_mass, n_valid_jets, Q_l);
					printf("	ev_fqqs = %f,   ev_fqqs_delta = %f,   ev_fqqa = %f,   ev_fqqa_delta = %f,   ev_fgg = %f,   ev_fbck = %f,    ev_WJets = %f,	ev_gtt = %f,   ev_gbk = %f\n",
							ev_fqqs, ev_fqqs_delta, ev_fqqa, ev_fqqa_delta, ev_fgg, ev_fbck, ev_WJets, ev_gtt, ev_gbk);
					printf("	10^12*ev_L_bck = %f,   10^12*ev_L_bck = %f,    10^12*ev_L_gg = %f,   10^12*ev_L_qqs = %f,   10^12*ev_L_qqa = %f\n",
							ev_L_bck*1000000000000, ev_L_WJets*1000000000000, ev_L_gg*1000000000000, ev_L_qqs*1000000000000, ev_L_qqa*1000000000000);
					printf("	Rqqbar = %.4f,   Rbck = %.4f,   RWJets = %.4f,    delta = %.4f,   Afb = %.4f\n", Rqqbar, Rbck, RWJets, delta, Afb);
					ev_L = DBL_EPSILON;
				}
				++count_fit;
				//take -2*ln(L)
				double ev_lnL = -2.0*TMath::Log(ev_L);
				//add to running total of -2ln(L)
				logL += ev_lnL;
				if (iterations==0 || plot == 1)
					event_likelihoods->Fill(ev_lnL);
				if (ev_lnL<min_ev_lnL)
					min_ev_lnL=ev_lnL;
				if (ev_lnL>max_ev_lnL)
					max_ev_lnL=ev_lnL;
			}
			else
				++count_almost_fit;
		}
		else {
			//printf("cos_theta_cs = %f, Feynman_x = %f, ttbar_mass = %f \n", cos_theta_cs, Feynman_x, ttbar_mass);
			++count_not_fit;
		}
		nb = fChain->GetEntry(jentry);   nbytes += nb;
	}//end loop
	
	iterations+=1;
	if (iterations%5 == 0) {
		printf("Rqqbar = %f\n",Rbck);
		printf("Rbck = %f, RWJets  = %f\n",Rbck,RWJets);
		printf("xi     = %f, delta = %f\n",xi,delta);
		printf("		Afb = %f\n",Afb);
		printf("min event log(L) = %.4f\n",min_ev_lnL);
		printf("max event log(L) = %.4f\n",max_ev_lnL);
		printf("Fit: %d (%.2f%%), Almost Fit: %d (%.2f%%), Not Fit: %d (%.2f%%)\n", 
			count_fit, 100.0*count_fit/(count_fit+count_almost_fit+count_not_fit),
			count_almost_fit, 100.0*count_almost_fit/(count_fit+count_almost_fit+count_not_fit), 
			count_not_fit, 100.0*count_not_fit/(count_fit+count_almost_fit+count_not_fit));
	}
	//put out the first iteration event likelihood plot
	if (iterations == 1) {
		TCanvas *c = new TCanvas("c","Event -2*Log(L) at Initialized Values",900,900);
		c->cd();
		event_likelihoods->Draw();
		c->Print("ev_lnL_initial.pdf","pdf");
	}
	//put out the last iteration event likelihood plot
	if (iterations!=1 && plot == 1) {
		TCanvas *c = new TCanvas("c","Event -2*Log(L) at Final Values",900,900);
		c->cd();
		event_likelihoods->Draw();
		c->Print("ev_lnL_final.pdf","pdf");
	}
	//return the log likelihood
	return logL;
}

//Main loop function: returns a log likelihood for the function over all events
//SEPARATED CASE
double angles_data::Loop(double Rqqbar_4jet, double Rbck_4jet, double RW4Jets, double Rqqbar_5jet, double Rbck_5jet, double RW5Jets, 
						 double xi_4jet, double xi_5jet, double delta_4jet, double delta_5jet, double Afb_4jet, double Afb_5jet)
{
	//Set up loop
	if (fChain == 0) return 0.0;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;

	iterations+=1;
	if (iterations%5 == 0) {
		printf("----------------------------------------------------------\n");
		printf("Rqqbar_4jet = %f, Rqqbar_5jet = %f\n",Rqqbar_4jet,Rqqbar_5jet);
		printf("Rbck_4jet   = %f, Rbck_5jet   = %f\n",Rbck_4jet,Rbck_5jet);
		printf("RW4Jets     = %f, RW5Jets     = %f\n",RW4Jets,RW5Jets);
		printf("xi_4jet     = %f, xi_5jet     = %f\n",xi_4jet,xi_5jet);
		printf("delta_4jet  = %f, delta_5jet  = %f\n",delta_4jet,delta_5jet);
		printf("Afb_4jet    = %f, Afb_5jet    = %f\n",Afb_4jet,Afb_5jet); 
	}

	//create variable to hold ln(L) which will be returned
	double logL = 0.0;

	int count_fit = 0;
	int count_not_fit = 0;
	int count_almost_fit = 0;
	//MAIN LOOP
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;

		//If Feynman x is negative, flip it.
		if (Feynman_x<0.0) 
			Feynman_x = -1.0*Feynman_x;
		//Redefine the Rqqbar and Rbck
		double Rqqbar = 0.0;
		double Rbck = 0.0;
		double RWJets = 0.0;
		double xi = 0.0;
		double delta = 0.0;
		double Afb = 0.0;
		double F_xi = 0.0;
		double F_delta = 0.0;
		if (n_valid_jets==4) {
			Rqqbar = Rqqbar_4jet;
			Rbck = Rbck_4jet;
			RWJets = RW4Jets;
			xi = xi_4jet;
			delta = delta_4jet;
			F_xi = F_xi_4jet;
			F_delta = F_delta_4jet;
			Afb = Afb_4jet;
		}
		else if (n_valid_jets==5) {
			Rqqbar = Rqqbar_5jet;
			Rbck = Rbck_5jet;
			RWJets = RW5Jets;
			xi = xi_5jet;
			delta = delta_5jet;
			F_xi = F_xi_5jet;
			F_delta = F_delta_5jet;
			Afb = Afb_5jet;
		}
		//cut on whether the event is in range of the training set
		if (cos_theta_cs>=x_low && cos_theta_cs<=x_high && Feynman_x>=y_low && Feynman_x<=y_high && ttbar_mass>=z_low && ttbar_mass<=z_high && n_bTags >= num_b_tag_cuts && ln_L < lnL_cut) {
			//use parameters of event to look up values of four functions from histograms (depending on lepton charge)
			double ev_fqqs = 0.0;
			double ev_fqqs_xi = 0.0;
			double ev_fqqs_delta = 0.0;
			double ev_fqqa = 0.0;
			double ev_fqqa_xi = 0.0;
			double ev_fqqa_delta = 0.0;
			double ev_fgg = 0.0;
			double ev_fbck = 0.0;
			double ev_WJets = 0.0;
			double ev_gtt = 0.0;
			double ev_gbk = 0.0;
			if (Q_l > 0) {
				if (n_valid_jets == 4) {
					ev_fqqs = fqqs_plus_4jet->GetBinContent(fqqs_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_xi = fqqs_xi_plus_4jet->GetBinContent(fqqs_xi_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_delta = fqqs_delta_plus_4jet->GetBinContent(fqqs_delta_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa = fqqa_plus_4jet->GetBinContent(fqqa_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_xi = fqqa_xi_plus_4jet->GetBinContent(fqqa_xi_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_delta = fqqa_delta_plus_4jet->GetBinContent(fqqa_delta_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fgg  = fgg_plus_4jet->GetBinContent(fgg_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fbck = fbck_plus_4jet->GetBinContent(fbck_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_WJets = WJets_plus_4jet->GetBinContent(WJets_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_gtt = gtt_4jet->GetBinContent(gtt_4jet->FindFixBin(ln_L));
					ev_gbk = gbk_4jet->GetBinContent(gbk_4jet->FindFixBin(ln_L));
				}
				else if (n_valid_jets == 5) {
					ev_fqqs = fqqs_plus_5jet->GetBinContent(fqqs_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_xi = fqqs_xi_plus_5jet->GetBinContent(fqqs_xi_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_delta = fqqs_delta_plus_5jet->GetBinContent(fqqs_delta_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa = fqqa_plus_5jet->GetBinContent(fqqa_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_xi = fqqa_xi_plus_5jet->GetBinContent(fqqa_xi_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_delta = fqqa_delta_plus_5jet->GetBinContent(fqqa_delta_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fgg  = fgg_plus_5jet->GetBinContent(fgg_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fbck = fbck_plus_5jet->GetBinContent(fbck_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_WJets = WJets_plus_5jet->GetBinContent(WJets_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_gtt = gtt_5jet->GetBinContent(gtt_5jet->FindFixBin(ln_L));
					ev_gbk = gbk_5jet->GetBinContent(gbk_5jet->FindFixBin(ln_L));
				}
				else {
					printf("JET NUMBER INFORMATION NOT AVAILABLE: EVENT SKIPPED!\n");
					continue;
				}
			}
			else if (Q_l < 0) {
				if (n_valid_jets == 4) {
					ev_fqqs = fqqs_minus_4jet->GetBinContent(fqqs_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_xi = fqqs_xi_minus_4jet->GetBinContent(fqqs_xi_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_delta = fqqs_delta_minus_4jet->GetBinContent(fqqs_delta_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa = fqqa_minus_4jet->GetBinContent(fqqa_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_xi = fqqa_xi_minus_4jet->GetBinContent(fqqa_xi_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_delta = fqqa_delta_minus_4jet->GetBinContent(fqqa_delta_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fgg  = fgg_minus_4jet->GetBinContent(fgg_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fbck = fbck_minus_4jet->GetBinContent(fbck_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_WJets = WJets_minus_4jet->GetBinContent(WJets_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_gtt = gtt_4jet->GetBinContent(gtt_4jet->FindFixBin(ln_L));
					ev_gbk = gbk_4jet->GetBinContent(gbk_4jet->FindFixBin(ln_L));
				}
				else if (n_valid_jets == 5) {
					ev_fqqs = fqqs_minus_5jet->GetBinContent(fqqs_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_xi = fqqs_xi_minus_5jet->GetBinContent(fqqs_xi_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqs_delta = fqqs_delta_minus_5jet->GetBinContent(fqqs_delta_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa = fqqa_minus_5jet->GetBinContent(fqqa_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_xi = fqqa_xi_minus_5jet->GetBinContent(fqqa_xi_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fqqa_delta = fqqa_delta_minus_5jet->GetBinContent(fqqa_delta_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fgg  = fgg_minus_5jet->GetBinContent(fgg_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_fbck = fbck_minus_5jet->GetBinContent(fbck_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_WJets = WJets_minus_5jet->GetBinContent(WJets_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
					ev_gtt = gtt_5jet->GetBinContent(gtt_5jet->FindFixBin(ln_L));
					ev_gbk = gbk_5jet->GetBinContent(gbk_5jet->FindFixBin(ln_L));
				}
				else {
					printf("JET NUMBER INFORMATION NOT AVAILABLE: EVENT SKIPPED!\n");
					continue;
				}
			}
			else {
				printf("LEPTON CHARGE INFORMATION NOT AVAILABLE: EVENT SKIPPED!\n");
				continue;
			}
			//calculate L for this one event
			// double ev_L_bck = Rbck*ev_fbck*ev_gbk;
			// double ev_L_WJets = RWJets*ev_WJets*ev_gbk;
			// double ev_L_gg = (1.0-Rbck-RWJets)*ev_gtt*(1.0-Rqqbar)*ev_fgg;
			// double ev_L_qqs = (1.0-Rbck-RWJets)*ev_gtt*Rqqbar*((1.0)/(1.0+xi*F_xi+delta*F_delta))*(ev_fqqs+xi*ev_fqqs_xi+delta*ev_fqqs_delta);
			// double ev_L_qqa = Afb*(1.0-Rbck-RWJets)*ev_gtt*Rqqbar*((1.0)/(1.0+xi*F_xi+delta*F_delta))*(ev_fqqa+xi*ev_fqqa_xi+delta*ev_fqqa_delta);
			double ev_L_bck = Rbck*ev_fbck;
			double ev_L_WJets = RWJets*ev_WJets;
			double ev_L_gg = (1.0-Rbck-RWJets)*(1.0-Rqqbar)*ev_fgg;
			double ev_L_qqs = (1.0-Rbck-RWJets)*Rqqbar*((1.0)/(1.0+xi*F_xi+delta*F_delta))*(ev_fqqs+xi*ev_fqqs_xi+delta*ev_fqqs_delta);
			double ev_L_qqa = Afb*(1.0-Rbck-RWJets)*Rqqbar*((1.0)/(1.0+xi*F_xi+delta*F_delta))*(ev_fqqa+xi*ev_fqqa_xi+delta*ev_fqqa_delta);
			double ev_L = ev_L_bck + ev_L_WJets + ev_L_gg + ev_L_qqs + ev_L_qqa;
			//printf("ev_L = %f = %f * %f + %f * %f + %f * (%f + %f * %f), F_xi = %f, F_delta = %f\n",ev_L,(1.0-Rqqbar-Rbck),ev_fgg,Rbck,ev_fbck,Rqqbar,ev_fqqs,Afb,ev_fqqa,F_xi,F_delta);
			//printf("point is cos_theta_cs = %f, Feynman_x = %f, ttbar_mass = %f\n", cos_theta_cs, Feynman_x, ttbar_mass);
			//make sure that there just weren't any training set events at all
			if (ev_fgg!=0 || ev_fqqs!=0 || ev_fbck!=0 || ev_WJets!=0) {
				if (ev_L<=0) {
					printf("Parameters are in a bad spot: event likelihood = %.12f\n",ev_L);
					printf("	costheta = %.4f,   x_F = %.4f,   M_tt = %.4f,   nJets = %d,   Q_l = %d\n", cos_theta_cs, Feynman_x, ttbar_mass, n_valid_jets, Q_l);
					printf("	ev_fqqs = %f,   ev_fqqs_delta = %f,   ev_fqqa = %f,   ev_fqqa_delta = %f,   ev_fgg = %f,   ev_fbck = %f,    ev_WJets = %f,	ev_gtt = %f,   ev_gbk = %f\n",
							ev_fqqs, ev_fqqs_delta, ev_fqqa, ev_fqqa_delta, ev_fgg, ev_fbck, ev_WJets, ev_gtt, ev_gbk);
					printf("	10^12*ev_L_bck = %f,   10^12*ev_L_bck = %f,    10^12*ev_L_gg = %f,   10^12*ev_L_qqs = %f,   10^12*ev_L_qqa = %f\n",
							ev_L_bck*1000000000000, ev_L_WJets*1000000000000, ev_L_gg*1000000000000, ev_L_qqs*1000000000000, ev_L_qqa*1000000000000);
					printf("	Rqqbar = %.4f,   Rbck = %.4f,   RWJets = %.4f,    delta = %.4f,   Afb = %.4f\n", Rqqbar, Rbck, RWJets, delta, Afb);
					ev_L = DBL_EPSILON;
				}
				++count_fit;
				//take -2*ln(L)
				double ev_lnL = -2.0*TMath::Log(ev_L);
				//add to running total of -2ln(L)
				logL += ev_lnL;
			}
			else
				++count_almost_fit;
		}
		else {
			//printf("cos_theta_cs = %f, Feynman_x = %f, ttbar_mass = %f \n", cos_theta_cs, Feynman_x, ttbar_mass);
			++count_not_fit;
		}
		nb = fChain->GetEntry(jentry);   nbytes += nb;
	}//end loop
	
	//printf("Fit: %d (%f%%), Almost Fit: %d (%f%%), Not Fit: %d (%f%%)\n", 
	//		count_fit, 100.0*count_fit/(count_fit+count_almost_fit+count_not_fit),
	//		count_almost_fit, 100.0*count_almost_fit/(count_fit+count_almost_fit+count_not_fit), 
	//		count_not_fit, 100.0*count_not_fit/(count_fit+count_almost_fit+count_not_fit));
	//return the log likelihood
	//fprintf(test, "%.8f 		%.8f\n",delta_4jet,lnL);
	return logL;
}

//other loop function for making plots of the final fit compared to the data
//COMBINED CASE
double angles_data::Loop(double Rqqbar, double sigma_Rqqbar, double Rbck, double sigma_Rbck, double RWJets, double sigma_RWJets,
						 double xi, double sigma_xi, double delta, double sigma_delta, 
						 double Afb, double sigma_Afb, char* rname) {
	//Set up loop
	if (fChain == 0) return 0.0;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	
	//make histograms to hold the data and the fit and their projections
	TH3D* data_hist = new TH3D("data_hist","Data Distribution; cos(#theta *); Feynman x (x_{F}); M_{t #bar{t}} (GeV)",
							   nbinsx,x_low,x_high,nbinsy,y_low,y_high,nbinsz,z_low,z_high);
	TH1D* data_x = new TH1D("data_x","Data x Projection",nbinsx,x_low,x_high);
	TH1D* data_y = new TH1D("data_y","Data y Projection",nbinsy,y_low,y_high);
	TH1D* data_z = new TH1D("data_z","Data z Projection",nbinsz,z_low,z_high);
	TH1D* gg_x = new TH1D("gg_x","gg x Projection",nbinsx,x_low,x_high);
	TH1D* gg_y = new TH1D("gg_y","gg y Projection",nbinsy,y_low,y_high);
	TH1D* gg_z = new TH1D("gg_z","gg z Projection",nbinsz,z_low,z_high);
	TH1D* qq_x = new TH1D("qq_x","qqbar x Projection",nbinsx,x_low,x_high);
	TH1D* qq_y = new TH1D("qq_y","qqbar y Projection",nbinsy,y_low,y_high);
	TH1D* qq_z = new TH1D("qq_z","qqbar z Projection",nbinsz,z_low,z_high);
	TH1D* bg_x = new TH1D("bg_x","background x Projection",nbinsx,x_low,x_high);
	TH1D* bg_y = new TH1D("bg_y","background y Projection",nbinsy,y_low,y_high);
	TH1D* bg_z = new TH1D("bg_z","background z Projection",nbinsz,z_low,z_high);
	TH1D* wj_x = new TH1D("wj_x","WJets x Projection",nbinsx,x_low,x_high);
	TH1D* wj_y = new TH1D("wj_y","WJets y Projection",nbinsy,y_low,y_high);
	TH1D* wj_z = new TH1D("wj_z","WJets z Projection",nbinsz,z_low,z_high);
	TH1D* event_numbers_data = new TH1D("event_numbers_data","lepton charge and jet multiplicity in data",6,0.,6.);
	TH1D* event_numbers_bck = new TH1D("event_numbers_bck","lepton charge and jet multiplicity in background",6,0.,6.);
	TH1D* event_numbers_WJets = new TH1D("event_numbers_WJets","lepton charge and jet multiplicity in WJets",6,0.,6.);
	TH1D* event_numbers_gg = new TH1D("event_numbers_gg","lepton charge and jet multiplicity in gg",6,0.,6.);
	TH1D* event_numbers_qq = new TH1D("event_numbers_qq","lepton charge and jet multiplicity in qq",6,0.,6.);
	event_numbers_bck->SetBarWidth(0.70);
	event_numbers_WJets->SetBarWidth(0.70);
	event_numbers_gg->SetBarWidth(0.70);
	event_numbers_qq->SetBarWidth(0.70);
	event_numbers_bck->SetBarOffset(0.15);
	event_numbers_WJets->SetBarOffset(0.15);
	event_numbers_gg->SetBarOffset(0.15);
	event_numbers_qq->SetBarOffset(0.15);
	TH3D* sideband = new TH3D("sideband","Sideband Distribution; cos(#theta *); Feynman x (x_{F}); M_{t #bar{t}} (GeV)",
							   nbinsx,x_low,x_high,nbinsy,y_low,y_high,nbinsz,z_low,z_high);

	int count_added=0;
	//MAIN LOOP
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;

		//If Feynman x is negative, flip it.
		if (Feynman_x<0.0) 
			Feynman_x = -1.0*Feynman_x;
		
		//cut on whether the event is in range of the training set
		if (cos_theta_cs>=x_low && cos_theta_cs<=x_high && Feynman_x>=y_low && Feynman_x<=y_high && ttbar_mass>=z_low && ttbar_mass<=z_high && n_bTags >= num_b_tag_cuts && ln_L < lnL_cut) {
			//Fill the data histogram with the appropriate values
			double ev_fqqs_plus_4jet = fqqs_plus_4jet->GetBinContent(fqqs_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fqqs_minus_4jet = fqqs_minus_4jet->GetBinContent(fqqs_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fqqs_plus_5jet = fqqs_plus_5jet->GetBinContent(fqqs_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fqqs_minus_5jet = fqqs_minus_5jet->GetBinContent(fqqs_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fgg_plus_4jet  = fgg_plus_4jet->GetBinContent(fgg_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fgg_minus_4jet  = fgg_minus_4jet->GetBinContent(fgg_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fgg_plus_5jet  = fgg_plus_5jet->GetBinContent(fgg_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fgg_minus_5jet  = fgg_minus_5jet->GetBinContent(fgg_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fbck_plus_4jet = fbck_plus_4jet->GetBinContent(fbck_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fbck_minus_4jet = fbck_minus_4jet->GetBinContent(fbck_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fbck_plus_5jet = fbck_plus_5jet->GetBinContent(fbck_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fbck_minus_5jet = fbck_minus_5jet->GetBinContent(fbck_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_WJets_plus_4jet = WJets_plus_4jet->GetBinContent(WJets_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_WJets_minus_4jet = WJets_minus_4jet->GetBinContent(WJets_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_WJets_plus_5jet = WJets_plus_5jet->GetBinContent(WJets_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_WJets_minus_5jet = WJets_minus_5jet->GetBinContent(WJets_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			//double ev_gtt_4jet = gtt_4jet->GetBinContent(gtt_4jet->FindFixBin(ln_L));
			//double ev_gbk_4jet = gbk_4jet->GetBinContent(gbk_4jet->FindFixBin(ln_L));
			//double ev_gtt_5jet = gtt_5jet->GetBinContent(gtt_5jet->FindFixBin(ln_L));
			//double ev_gbk_5jet = gbk_5jet->GetBinContent(gbk_5jet->FindFixBin(ln_L));
			if ((Q_l>0 && n_valid_jets==4 && (ev_fgg_plus_4jet!=0 || ev_fqqs_plus_4jet!=0 || ev_fbck_plus_4jet!=0 || ev_WJets_plus_4jet)) ||
				(Q_l<0 && n_valid_jets==4 && (ev_fgg_minus_4jet!=0 || ev_fqqs_minus_4jet!=0 || ev_fbck_minus_4jet!=0 || ev_WJets_minus_4jet)) ||
				(Q_l>0 && n_valid_jets==5 && (ev_fgg_plus_5jet!=0 || ev_fqqs_plus_5jet!=0 || ev_fbck_plus_5jet!=0 || ev_WJets_plus_5jet)) ||
				(Q_l<0 && n_valid_jets==5 && (ev_fgg_minus_5jet!=0 || ev_fqqs_minus_5jet!=0 || ev_fbck_minus_5jet!=0 || ev_WJets_minus_5jet))) {
				data_hist->Fill(cos_theta_cs,Feynman_x,ttbar_mass);
				data_x->Fill(cos_theta_cs);
				data_y->Fill(Feynman_x);
				data_z->Fill(ttbar_mass);
				++count_added;
				if (n_valid_jets==4 && Q_l>0)
					event_numbers_data->Fill(1.5);
				else if (n_valid_jets==4 && Q_l<0)
					event_numbers_data->Fill(2.5);
				else if (n_valid_jets==5 && Q_l>0)
					event_numbers_data->Fill(3.5);
				else if (n_valid_jets==5 && Q_l<0)
					event_numbers_data->Fill(4.5);
			}
			else 
				sideband->Fill(cos_theta_cs,Feynman_x,ttbar_mass);
		}//end cut
		nb = fChain->GetEntry(jentry);   nbytes += nb;
	}//end loop

	printf("BUILDING FIT COMPARISON PLOT\n");
	printf("# of events in data = %d\n", count_added);

	//build the total histograms by summing over the lepton charge and jet number and rescaling
	TH3D *fqqs = (TH3D*)fqqs_plus_4jet->Clone("fqqs");
	TH3D *fqqa = (TH3D*)fqqs_plus_4jet->Clone("fqqa");
	TH3D *fgg = (TH3D*)fqqs_plus_4jet->Clone("fgg");
	TH3D *fbck = (TH3D*)fqqs_plus_4jet->Clone("fbck");
	TH3D *fWJets = (TH3D*)fqqs_plus_4jet->Clone("fWJets");
	for (double x = x_low+(0.5*xbinwidth); x<x_high; x=x+xbinwidth) {
		for (double y = y_low+(0.5*ybinwidth); y<y_high; y=y+ybinwidth) {
			for (double z = z_low+(0.5*zbinwidth); z<z_high; z=z+zbinwidth) {
				double bck_events = fbck_plus_4jet->GetBinContent(fbck_plus_4jet->FindFixBin(x,y,z)) + fbck_minus_4jet->GetBinContent(fbck_minus_4jet->FindFixBin(x,y,z))
								  + fbck_plus_5jet->GetBinContent(fbck_plus_5jet->FindFixBin(x,y,z)) + fbck_minus_5jet->GetBinContent(fbck_minus_4jet->FindFixBin(x,y,z));
				double WJets_events = WJets_plus_4jet->GetBinContent(WJets_plus_4jet->FindFixBin(x,y,z)) + WJets_minus_4jet->GetBinContent(WJets_minus_4jet->FindFixBin(x,y,z))
								  + WJets_plus_5jet->GetBinContent(WJets_plus_5jet->FindFixBin(x,y,z)) + WJets_minus_5jet->GetBinContent(WJets_minus_4jet->FindFixBin(x,y,z));
				double gg_events  = fgg_plus_4jet->GetBinContent(fgg_plus_4jet->FindFixBin(x,y,z)) + fgg_minus_4jet->GetBinContent(fgg_minus_4jet->FindFixBin(x,y,z))
								  + fgg_plus_5jet->GetBinContent(fgg_plus_5jet->FindFixBin(x,y,z)) + fgg_minus_5jet->GetBinContent(fgg_minus_5jet->FindFixBin(x,y,z));
				double qqs_events = fqqs_plus_4jet->GetBinContent(fqqs_plus_4jet->FindFixBin(x,y,z)) + fqqs_minus_4jet->GetBinContent(fqqs_minus_4jet->FindFixBin(x,y,z))
								  + fqqs_plus_5jet->GetBinContent(fqqs_plus_5jet->FindFixBin(x,y,z)) + fqqs_minus_5jet->GetBinContent(fqqs_minus_5jet->FindFixBin(x,y,z));
				double qqs_xi_events = fqqs_xi_plus_4jet->GetBinContent(fqqs_xi_plus_4jet->FindFixBin(x,y,z)) + fqqs_xi_minus_4jet->GetBinContent(fqqs_xi_minus_4jet->FindFixBin(x,y,z))
									 + fqqs_xi_plus_5jet->GetBinContent(fqqs_xi_plus_5jet->FindFixBin(x,y,z)) + fqqs_xi_minus_5jet->GetBinContent(fqqs_xi_minus_5jet->FindFixBin(x,y,z));
				double qqs_delta_events = fqqs_delta_plus_4jet->GetBinContent(fqqs_delta_plus_4jet->FindFixBin(x,y,z)) + 
										  fqqs_delta_minus_4jet->GetBinContent(fqqs_delta_minus_4jet->FindFixBin(x,y,z)) +
										  fqqs_delta_plus_5jet->GetBinContent(fqqs_delta_plus_5jet->FindFixBin(x,y,z)) + 
										  fqqs_delta_minus_5jet->GetBinContent(fqqs_delta_minus_5jet->FindFixBin(x,y,z));
				double qqa_events = fqqa_plus_4jet->GetBinContent(fqqa_plus_4jet->FindFixBin(x,y,z)) + fqqa_minus_4jet->GetBinContent(fqqa_minus_4jet->FindFixBin(x,y,z))
								  + fqqa_plus_5jet->GetBinContent(fqqa_plus_5jet->FindFixBin(x,y,z)) + fqqa_minus_5jet->GetBinContent(fqqa_minus_5jet->FindFixBin(x,y,z));
				double qqa_xi_events = fqqa_xi_plus_4jet->GetBinContent(fqqa_xi_plus_4jet->FindFixBin(x,y,z)) + fqqa_xi_minus_4jet->GetBinContent(fqqa_xi_minus_4jet->FindFixBin(x,y,z))
									 + fqqa_xi_plus_5jet->GetBinContent(fqqa_xi_plus_5jet->FindFixBin(x,y,z)) + fqqa_xi_minus_5jet->GetBinContent(fqqa_xi_minus_5jet->FindFixBin(x,y,z));
				double qqa_delta_events = fqqa_delta_plus_4jet->GetBinContent(fqqa_delta_plus_4jet->FindFixBin(x,y,z)) + 
										  fqqa_delta_minus_4jet->GetBinContent(fqqa_delta_minus_4jet->FindFixBin(x,y,z)) +
										  fqqa_delta_plus_5jet->GetBinContent(fqqa_delta_plus_5jet->FindFixBin(x,y,z)) + 
										  fqqa_delta_minus_5jet->GetBinContent(fqqa_delta_minus_5jet->FindFixBin(x,y,z));
				double bck = Rbck*bck_events;
				double WJets = RWJets*WJets_events;
				double gg  = (1.0-Rbck-RWJets)*(1.0-Rqqbar)*gg_events;
				double qqs = (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(qqs_events+xi*qqs_xi_events+delta*qqs_delta_events);
				double qqa = Afb*(1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(qqa_events+xi*qqa_xi_events+delta*qqa_delta_events);
				fqqs->SetBinContent(fqqs->FindFixBin(x,y,z),qqs);
				fqqa->SetBinContent(fqqa->FindFixBin(x,y,z),qqa);
				fgg->SetBinContent(fgg->FindFixBin(x,y,z),gg);
				fbck->SetBinContent(fbck->FindFixBin(x,y,z),bck);
				fWJets->SetBinContent(fWJets->FindFixBin(x,y,z),WJets);
				gg_x->Fill(x,count_added*gg);
				gg_y->Fill(y,count_added*gg);
				gg_z->Fill(z,count_added*gg);
				qq_x->Fill(x,count_added*(qqs+qqa));
				qq_y->Fill(y,count_added*(qqs+qqa));
				qq_z->Fill(z,count_added*(qqs+qqa));
				bg_x->Fill(x,count_added*bck);
				bg_y->Fill(y,count_added*bck);
				bg_z->Fill(z,count_added*bck);
				wj_x->Fill(x,count_added*WJets);
				wj_y->Fill(y,count_added*WJets);
				wj_z->Fill(z,count_added*WJets);
			}
		}
	}//end loops over bins
	
	//and a few more loops over the histograms to set the errors on the projections for the data
	for (double x = x_low+(0.5*xbinwidth); x<x_high; x=x+xbinwidth) {
		double x_bin_var_data=0.0;
		for (double y = y_low+(0.5*ybinwidth); y<y_high; y=y+ybinwidth) {
			double y_bin_var_data=0.0;
			for (double z = z_low+(0.5*zbinwidth); z<z_high; z=z+zbinwidth) {
				y_bin_var_data=y_bin_var_data+data_hist->GetBinError(data_hist->FindFixBin(x,y,z))*data_hist->GetBinError(data_hist->FindFixBin(x,y,z));
			}
			x_bin_var_data=x_bin_var_data+y_bin_var_data;
		}
		data_x->SetBinError(data_x->FindFixBin(x),TMath::Sqrt(x_bin_var_data));
	}
	for (double z = z_low+(0.5*zbinwidth); z<z_high; z=z+zbinwidth) {
		double z_bin_var_data=0.0;
		for (double x = x_low+(0.5*xbinwidth); x<x_high; x=x+xbinwidth) {
			double x_bin_var_data=0.0;
			for (double y = y_low+(0.5*ybinwidth); y<y_high; y=y+ybinwidth) {
				x_bin_var_data=x_bin_var_data+data_hist->GetBinError(data_hist->FindFixBin(x,y,z))*data_hist->GetBinError(data_hist->FindFixBin(x,y,z));
			}
			z_bin_var_data=z_bin_var_data+x_bin_var_data;
		}
		data_z->SetBinError(data_z->FindFixBin(z),TMath::Sqrt(z_bin_var_data));
	}
	for (double y = y_low+(0.5*ybinwidth); y<y_high; y=y+ybinwidth) {
		double y_bin_var_data=0.0;
		for (double z = z_low+(0.5*zbinwidth); z<z_high; z=z+zbinwidth) {
			double z_bin_var_data=0.0;
			for (double x = x_low+(0.5*xbinwidth); x<x_high; x=x+xbinwidth) {
				z_bin_var_data=z_bin_var_data+data_hist->GetBinError(data_hist->FindFixBin(x,y,z))*data_hist->GetBinError(data_hist->FindFixBin(x,y,z));
			}
			y_bin_var_data=y_bin_var_data+z_bin_var_data;
		}
		data_y->SetBinError(data_y->FindFixBin(y),TMath::Sqrt(y_bin_var_data));
	}

	//Build the event type distributions for the background, gg, and qq
	double n_fit_bck_plus_4jet = (Rbck*fbck_plus_4jet->Integral());
	double n_fit_WJets_plus_4jet = (RWJets*WJets_plus_4jet->Integral());  
	double n_fit_gg_plus_4jet  = (1.0-Rbck-RWJets)*((1.0-Rqqbar)*fgg_plus_4jet->Integral());
	double n_fit_qq_plus_4jet  = (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(fqqs_plus_4jet->Integral()+Afb*fqqa_plus_4jet->Integral());
		   n_fit_qq_plus_4jet += (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*xi*(fqqs_xi_plus_4jet->Integral()+Afb*fqqa_xi_plus_4jet->Integral());
		   n_fit_qq_plus_4jet += (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*delta*(fqqs_delta_plus_4jet->Integral()+Afb*fqqa_delta_plus_4jet->Integral());
	double n_fit_bck_minus_4jet = (Rbck*fbck_minus_4jet->Integral()); 
	double n_fit_WJets_minus_4jet = (RWJets*WJets_minus_4jet->Integral()); 
	double n_fit_gg_minus_4jet  = (1.0-Rbck-RWJets)*((1.0-Rqqbar)*fgg_minus_4jet->Integral());
	double n_fit_qq_minus_4jet  = (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(fqqs_minus_4jet->Integral()+Afb*fqqa_minus_4jet->Integral());
		   n_fit_qq_minus_4jet += (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*xi*(fqqs_xi_minus_4jet->Integral()+Afb*fqqa_xi_minus_4jet->Integral());
		   n_fit_qq_minus_4jet += (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*delta*(fqqs_delta_minus_4jet->Integral()+Afb*fqqa_delta_minus_4jet->Integral());
	double n_fit_bck_plus_5jet = (Rbck*fbck_plus_5jet->Integral()); 
	double n_fit_WJets_plus_5jet = (RWJets*WJets_plus_5jet->Integral()); 
	double n_fit_gg_plus_5jet  = (1.0-Rbck-RWJets)*((1.0-Rqqbar)*fgg_plus_5jet->Integral());
	double n_fit_qq_plus_5jet  = (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(fqqs_plus_5jet->Integral()+Afb*fqqa_plus_5jet->Integral());
		   n_fit_qq_plus_5jet += (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*xi*(fqqs_xi_plus_5jet->Integral()+Afb*fqqa_xi_plus_5jet->Integral());
		   n_fit_qq_plus_5jet += (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*delta*(fqqs_delta_plus_5jet->Integral()+Afb*fqqa_delta_plus_5jet->Integral());
	double n_fit_bck_minus_5jet = (Rbck*fbck_minus_5jet->Integral()); 
	double n_fit_WJets_minus_5jet = (RWJets*WJets_minus_5jet->Integral()); 
	double n_fit_gg_minus_5jet  = (1.0-Rbck-RWJets)*((1.0-Rqqbar)*fgg_minus_5jet->Integral());
	double n_fit_qq_minus_5jet  = (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*(fqqs_minus_5jet->Integral()+Afb*fqqa_minus_5jet->Integral());
		   n_fit_qq_minus_5jet += (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*xi*(fqqs_xi_minus_5jet->Integral()+Afb*fqqa_xi_minus_5jet->Integral());
		   n_fit_qq_minus_5jet += (1.0-Rbck-RWJets)*(Rqqbar/(1.0+xi*F_xi_comb+delta*F_delta_comb))*delta*(fqqs_delta_minus_5jet->Integral()+Afb*fqqa_delta_minus_5jet->Integral());
	event_numbers_bck->Fill(1.5,count_added*n_fit_bck_plus_4jet); event_numbers_bck->Fill(2.5,count_added*n_fit_bck_minus_4jet);
	event_numbers_bck->Fill(3.5,count_added*n_fit_bck_plus_5jet); event_numbers_bck->Fill(4.5,count_added*n_fit_bck_minus_5jet);
	event_numbers_WJets->Fill(1.5,count_added*n_fit_WJets_plus_4jet); event_numbers_WJets->Fill(2.5,count_added*n_fit_WJets_minus_4jet);
	event_numbers_WJets->Fill(3.5,count_added*n_fit_WJets_plus_5jet); event_numbers_WJets->Fill(4.5,count_added*n_fit_WJets_minus_5jet);
	event_numbers_gg->Fill(1.5,count_added*n_fit_gg_plus_4jet); event_numbers_gg->Fill(2.5,count_added*n_fit_gg_minus_4jet);
	event_numbers_gg->Fill(3.5,count_added*n_fit_gg_plus_5jet); event_numbers_gg->Fill(4.5,count_added*n_fit_gg_minus_5jet);
	event_numbers_qq->Fill(1.5,count_added*n_fit_qq_plus_4jet); event_numbers_qq->Fill(2.5,count_added*n_fit_qq_minus_4jet);
	event_numbers_qq->Fill(3.5,count_added*n_fit_qq_plus_5jet); event_numbers_qq->Fill(4.5,count_added*n_fit_qq_minus_5jet);
	
	//plot the new histograms on top of each other
	//build histogram stacks
	THStack *x_stack = new THStack("x_stack","cos(#theta *) Comparison; cos(#theta *)");
	THStack *y_stack = new THStack("y_stack","Feynman x Comparison; (x_{F})");
	THStack *z_stack = new THStack("z_stack","t #bar{t} Mass Comparison; M_{t #bar{t}}");
	THStack *event_numbers_stack = new THStack("event_numbers_stack","Lepton Charge and Jet Multiplicity Comparison");
	//add to stacks
	x_stack->Add(bg_x);
	y_stack->Add(bg_y);
	z_stack->Add(bg_z);
	x_stack->Add(wj_x);
	y_stack->Add(wj_y);
	z_stack->Add(wj_z);
	x_stack->Add(gg_x);
	y_stack->Add(gg_y);
	z_stack->Add(gg_z);
	x_stack->Add(qq_x);
	y_stack->Add(qq_y);
	z_stack->Add(qq_z);
	event_numbers_stack->Add(event_numbers_bck);
	event_numbers_stack->Add(event_numbers_WJets);
	event_numbers_stack->Add(event_numbers_gg);
	event_numbers_stack->Add(event_numbers_qq);
	//set drawing options for projections
	//glueglue
	gg_x->SetFillColor(38);
	gg_x->SetMarkerStyle(21);
	gg_x->SetMarkerColor(38);
	gg_y->SetFillColor(38);
	gg_y->SetMarkerStyle(21);
	gg_y->SetMarkerColor(38);
	gg_z->SetFillColor(38);
	gg_z->SetMarkerStyle(21);
	gg_z->SetMarkerColor(38);
	event_numbers_gg->SetFillColor(38);
	event_numbers_gg->SetMarkerStyle(21);
	event_numbers_gg->SetMarkerColor(38);
	//qqbar
	qq_x->SetFillColor(46);
	qq_x->SetMarkerStyle(21);
	qq_x->SetMarkerColor(46);
	qq_y->SetFillColor(46);
	qq_y->SetMarkerStyle(21);
	qq_y->SetMarkerColor(46);
	qq_z->SetFillColor(46);
	qq_z->SetMarkerStyle(21);
	qq_z->SetMarkerColor(46);
	event_numbers_qq->SetFillColor(46);
	event_numbers_qq->SetMarkerStyle(21);
	event_numbers_qq->SetMarkerColor(46);
	//background
	bg_x->SetFillColor(41);
	bg_x->SetMarkerStyle(21);
	bg_x->SetMarkerColor(41);
	bg_y->SetFillColor(41);
	bg_y->SetMarkerStyle(21);
	bg_y->SetMarkerColor(41);
	bg_z->SetFillColor(41);
	bg_z->SetMarkerStyle(21);
	bg_z->SetMarkerColor(41);
	event_numbers_bck->SetFillColor(41);
	event_numbers_bck->SetMarkerStyle(21);
	event_numbers_bck->SetMarkerColor(41);
	//WJets
	wj_x->SetFillColor(kGreen-3);
	wj_x->SetMarkerStyle(21);
	wj_x->SetMarkerColor(kGreen-3);
	wj_y->SetFillColor(kGreen-3);
	wj_y->SetMarkerStyle(21);
	wj_y->SetMarkerColor(kGreen-3);
	wj_z->SetFillColor(kGreen-3);
	wj_z->SetMarkerStyle(21);
	wj_z->SetMarkerColor(kGreen-3);
	event_numbers_WJets->SetFillColor(kGreen-3);
	event_numbers_WJets->SetMarkerStyle(21);
	event_numbers_WJets->SetMarkerColor(kGreen-3);
	//data
	data_x->SetMarkerColor(kBlack);
	data_y->SetMarkerColor(kBlack);
	data_z->SetMarkerColor(kBlack);
	data_x->SetLineWidth(2);
	data_y->SetLineWidth(2);
	data_z->SetLineWidth(2);
	event_numbers_data->SetMarkerColor(kBlack);
	event_numbers_data->SetLineWidth(2);
	//add a legend to the plot
	TLegend *leg = new TLegend(0.56,0.73,0.88,0.88);
	leg->SetName("legend");
	leg->AddEntry(data_x,"Data","LPEX0");
	leg->AddEntry(qq_x,"q#bar{q} #rightarrow t#bar{t}(j)","F");
	leg->AddEntry(gg_x,"gg(qg) #rightarrow t#bar{t}(j)","F");
	leg->AddEntry(wj_x,"WJets Background","F");
	leg->AddEntry(bg_x,"Other Backgrounds","F");
	TCanvas* c = new TCanvas("c","Fit Comparison",1200,1200);
	c->Divide(2,2);
	c->cd(1);
	event_numbers_stack->SetMaximum(1.1*event_numbers_data->GetMaximum());
	event_numbers_stack->Draw("bar1");
	event_numbers_stack->GetXaxis()->SetBinLabel(event_numbers_stack->GetXaxis()->FindFixBin(1.5),"4jets, l+");
	event_numbers_stack->GetXaxis()->SetBinLabel(event_numbers_stack->GetXaxis()->FindFixBin(2.5),"4jets, l-");
	event_numbers_stack->GetXaxis()->SetBinLabel(event_numbers_stack->GetXaxis()->FindFixBin(3.5),"5jets, l+");
	event_numbers_stack->GetXaxis()->SetBinLabel(event_numbers_stack->GetXaxis()->FindFixBin(4.5),"5jets, l-");
	event_numbers_data->Draw("SAME PE1X0");
	leg->Draw();
	c->cd(2);
	y_stack->SetMaximum(1.1*data_y->GetMaximum());
	y_stack->Draw();
	data_y->Draw("SAME PE1X0");
	c->cd(3);
	x_stack->SetMaximum(1.1*data_x->GetMaximum());
	x_stack->Draw();
	data_x->Draw("SAME PE1X0");	
	c->cd(4);
	z_stack->SetMaximum(1.1*data_z->GetMaximum());
	z_stack->Draw();
	data_z->Draw("SAME PE1X0");
	//save the plots
	c->Print("fit_comparison.png","png");	

	TCanvas* c1 = new TCanvas("c1","Sideband",900,900);
	c1->Divide(2,2,0.01,0.01,0);
	c1->cd(1);
	sideband->Draw("BOX");
	c1->cd(2);
	(sideband->ProjectionY())->Draw();
	c1->cd(3);
	(sideband->ProjectionX())->Draw();
	c1->cd(4);
	(sideband->ProjectionZ())->Draw();
	c1->Print("sideband.pdf","pdf");
	
	// Save final stacks and data hist and legend into a root file for post processing
	printf("Writing final stack and data hists into a root file!\n");
	TFile* stack_file = new TFile("final_stack.root","Recreate");
	stack_file->cd();
	x_stack->Write();data_x->Write();
	y_stack->Write();data_y->Write();
	z_stack->Write();data_z->Write();
	event_numbers_stack->Write();event_numbers_data->Write();
	leg->Write();
	stack_file->Write();
	stack_file->Close();

	//finally, save the results of the fit to a text file.
	FILE* o = fopen("fit_results.txt","w");
	fprintf(o,"*************		FIT       RESULTS		*************\n");
	fprintf(o,"Rqqbar = %-.4f +/- %-.4f\n",Rqqbar,sigma_Rqqbar);
	fprintf(o,"Rbck = %-.4f +/- %-.4f\n",Rbck,sigma_Rbck);
	fprintf(o,"RWJets = %-.4f +/- %-.4f\n",RWJets,sigma_RWJets);
	fprintf(o,"xi = %-.4f +/- %-.4f\n",xi,sigma_xi);
	fprintf(o,"delta = %-.4f +/- %-.4f\n",delta,sigma_delta);
	fprintf(o,"Afb    = %-.4f +/- %-.4f\n",Afb,sigma_Afb);
	fprintf(o,"\n\n");
	fprintf(o,"Runname is %s\n",rname);
	fprintf(o,"\n");
	fprintf(o,"Binning and limits:\n");
	fprintf(o,"    %-.2f < costheta* < %-.2f in %d bins\n",x_low,x_high,nbinsx);
	fprintf(o,"    %-.2f < Feynman x < %-.2f in %d bins\n",y_low,y_high,nbinsy);
	fprintf(o,"    %-.2f < ttbarMass < %-.2f in %d bins\n",z_low,z_high,nbinsz);
	fclose(o);
	//and save an abridged version to the summary file
	o = fopen("summary_combined.txt","a");
	fprintf(o,"%s	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%d	%.4f\n",
			rname,Rqqbar,sigma_Rqqbar,Rbck,sigma_Rbck,RWJets,sigma_RWJets,xi,sigma_xi,delta,sigma_delta,Afb,sigma_Afb,count_added,
			Loop(Rqqbar, Rbck, RWJets, xi, delta, Afb, 1));
	fclose(o);
	o = fopen("summary_combined.csv","a");
	fprintf(o,"%s,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%d,%.4f\n",
			rname,Rqqbar,sigma_Rqqbar,Rbck,sigma_Rbck,RWJets,sigma_RWJets,xi,sigma_xi,delta,sigma_delta,Afb,sigma_Afb,count_added,
			Loop(Rqqbar, Rbck, RWJets, xi, delta, Afb, 0));
	fclose(o);
	return 0;
}

//other loop function for making plots of the final fit compared to the data
//SEPARATED CASE
double angles_data::Loop(double Rqqbar_4jet, double sigma_Rqqbar_4jet, double Rbck_4jet, double sigma_Rbck_4jet, double RW4Jets, double sigma_RW4Jets,
						 double Rqqbar_5jet, double sigma_Rqqbar_5jet, double Rbck_5jet, double sigma_Rbck_5jet, double RW5Jets, double sigma_RW5Jets,
						 double xi_4jet, double sigma_xi_4jet, double xi_5jet, double sigma_xi_5jet,
						 double delta_4jet, double sigma_delta_4jet, double delta_5jet, double sigma_delta_5jet, 
						 double Afb_4jet, double sigma_Afb_4jet, double Afb_5jet, double sigma_Afb_5jet, char* rname) {
	//Set up loop
	if (fChain == 0) return 0.0;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	
	//make histograms to hold the data and the fit and their projections
	TH3D* data_hist = new TH3D("data_hist","Data Distribution; cos(#theta *); Feynman x (x_{F}); M_{t #bar{t}} (GeV)",
							   nbinsx,x_low,x_high,nbinsy,y_low,y_high,nbinsz,z_low,z_high);
	TH1D* data_x = new TH1D("data_x","Data x Projection",nbinsx,x_low,x_high);
	TH1D* data_y = new TH1D("data_y","Data y Projection",nbinsy,y_low,y_high);
	TH1D* data_z = new TH1D("data_z","Data z Projection",nbinsz,z_low,z_high);
	TH1D* gg_x = new TH1D("gg_x","gg x Projection",nbinsx,x_low,x_high);
	TH1D* gg_y = new TH1D("gg_y","gg y Projection",nbinsy,y_low,y_high);
	TH1D* gg_z = new TH1D("gg_z","gg z Projection",nbinsz,z_low,z_high);
	TH1D* qq_x = new TH1D("qq_x","qqbar x Projection",nbinsx,x_low,x_high);
	TH1D* qq_y = new TH1D("qq_y","qqbar y Projection",nbinsy,y_low,y_high);
	TH1D* qq_z = new TH1D("qq_z","qqbar z Projection",nbinsz,z_low,z_high);
	TH1D* bg_x = new TH1D("bg_x","background x Projection",nbinsx,x_low,x_high);
	TH1D* bg_y = new TH1D("bg_y","background y Projection",nbinsy,y_low,y_high);
	TH1D* bg_z = new TH1D("bg_z","background z Projection",nbinsz,z_low,z_high);
	TH1D* wj_x = new TH1D("wj_x","WJets x Projection",nbinsx,x_low,x_high);
	TH1D* wj_y = new TH1D("wj_y","WJets y Projection",nbinsy,y_low,y_high);
	TH1D* wj_z = new TH1D("wj_z","WJets z Projection",nbinsz,z_low,z_high);
	TH1D* event_numbers_data = new TH1D("event_numbers_data","lepton charge and jet multiplicity in data",6,0.,6.);
	TH1D* event_numbers_bck = new TH1D("event_numbers_bck","lepton charge and jet multiplicity in background",6,0.,6.);
	TH1D* event_numbers_WJets = new TH1D("event_numbers_WJets","lepton charge and jet multiplicity in WJets",6,0.,6.);
	TH1D* event_numbers_gg = new TH1D("event_numbers_gg","lepton charge and jet multiplicity in gg",6,0.,6.);
	TH1D* event_numbers_qq = new TH1D("event_numbers_qq","lepton charge and jet multiplicity in qq",6,0.,6.);
	event_numbers_bck->SetBarWidth(0.70);
	event_numbers_WJets->SetBarWidth(0.70);
	event_numbers_gg->SetBarWidth(0.70);
	event_numbers_qq->SetBarWidth(0.70);
	event_numbers_bck->SetBarOffset(0.15);
	event_numbers_WJets->SetBarOffset(0.15);
	event_numbers_gg->SetBarOffset(0.15);
	event_numbers_qq->SetBarOffset(0.15);

	TH3D* sideband = new TH3D("sideband","Sideband Distribution; cos(#theta *); Feynman x (x_{F}); M_{t #bar{t}} (GeV)",
							   nbinsx,x_low,x_high,nbinsy,y_low,y_high,nbinsz,z_low,z_high);

	int count_added_4jet = 0;
	int count_added_5jet = 0;
	//MAIN LOOP
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;

		//If Feynman x is negative, flip it.
		if (Feynman_x<0.0) 
			Feynman_x = -1.0*Feynman_x;
		
		//cut on whether the event is in range of the training set
		if (cos_theta_cs>=x_low && cos_theta_cs<=x_high && Feynman_x>=y_low && Feynman_x<=y_high && ttbar_mass>=z_low && ttbar_mass<=z_high && n_bTags >= num_b_tag_cuts && ln_L < lnL_cut) {
			//Fill the data histogram with the appropriate values
			double ev_fqqs_plus_4jet = fqqs_plus_4jet->GetBinContent(fqqs_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fqqs_minus_4jet = fqqs_minus_4jet->GetBinContent(fqqs_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fqqs_plus_5jet = fqqs_plus_5jet->GetBinContent(fqqs_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fqqs_minus_5jet = fqqs_minus_5jet->GetBinContent(fqqs_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fgg_plus_4jet  = fgg_plus_4jet->GetBinContent(fgg_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fgg_minus_4jet  = fgg_minus_4jet->GetBinContent(fgg_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fgg_plus_5jet  = fgg_plus_5jet->GetBinContent(fgg_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fgg_minus_5jet  = fgg_minus_5jet->GetBinContent(fgg_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fbck_plus_4jet = fbck_plus_4jet->GetBinContent(fbck_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fbck_minus_4jet = fbck_minus_4jet->GetBinContent(fbck_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fbck_plus_5jet = fbck_plus_5jet->GetBinContent(fbck_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_fbck_minus_5jet = fbck_minus_5jet->GetBinContent(fbck_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_WJets_plus_4jet = WJets_plus_4jet->GetBinContent(WJets_plus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_WJets_minus_4jet = WJets_minus_4jet->GetBinContent(WJets_minus_4jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_WJets_plus_5jet = WJets_plus_5jet->GetBinContent(WJets_plus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			double ev_WJets_minus_5jet = WJets_minus_5jet->GetBinContent(WJets_minus_5jet->FindFixBin(cos_theta_cs,Feynman_x,ttbar_mass));
			//double ev_gtt_4jet = gtt_4jet->GetBinContent(gtt_4jet->FindFixBin(ln_L));
			//double ev_gbk_4jet = gbk_4jet->GetBinContent(gbk_4jet->FindFixBin(ln_L));
			//double ev_gtt_5jet = gtt_5jet->GetBinContent(gtt_5jet->FindFixBin(ln_L));
			//double ev_gbk_5jet = gbk_5jet->GetBinContent(gbk_5jet->FindFixBin(ln_L));
			if ((Q_l>0 && n_valid_jets==4 && (ev_fgg_plus_4jet!=0 || ev_fqqs_plus_4jet!=0 || ev_fbck_plus_4jet!=0 || ev_WJets_plus_4jet)) ||
				(Q_l<0 && n_valid_jets==4 && (ev_fgg_minus_4jet!=0 || ev_fqqs_minus_4jet!=0 || ev_fbck_minus_4jet!=0 || ev_WJets_minus_4jet)) ||
			 	(Q_l>0 && n_valid_jets==5 && (ev_fgg_plus_5jet!=0 || ev_fqqs_plus_5jet!=0 || ev_fbck_plus_5jet!=0 || ev_WJets_plus_5jet)) ||
				(Q_l<0 && n_valid_jets==5 && (ev_fgg_minus_5jet!=0 || ev_fqqs_minus_5jet!=0 || ev_fbck_minus_5jet!=0 || ev_WJets_minus_5jet))) {
				data_hist->Fill(cos_theta_cs,Feynman_x,ttbar_mass);
				data_x->Fill(cos_theta_cs);
				data_y->Fill(Feynman_x);
				data_z->Fill(ttbar_mass);
				if (n_valid_jets==4 && Q_l>0) {
					++count_added_4jet;
					event_numbers_data->Fill(1.5);
				}
				else if (n_valid_jets==4 && Q_l<0) {
					++count_added_4jet;
					event_numbers_data->Fill(2.5);
				}
				else if (n_valid_jets==5 && Q_l>0) {
					++count_added_5jet;
					event_numbers_data->Fill(3.5);
				}
				else if (n_valid_jets==5 && Q_l<0) {
					++count_added_5jet;
					event_numbers_data->Fill(4.5);
				}
			}
			else 
				sideband->Fill(cos_theta_cs,Feynman_x,ttbar_mass);
		}//end cut
		nb = fChain->GetEntry(jentry);   nbytes += nb;
	}//end loop

	printf("BUILDING FIT COMPARISON PLOT\n");
	printf("# of events in data = %d\n", count_added_4jet+count_added_5jet);

	//build the total histograms by summing over the lepton charge and jet number and rescaling
	TH3D *fqqs = (TH3D*)fqqs_plus_4jet->Clone("fqqs");
	TH3D *fqqa = (TH3D*)fqqs_plus_4jet->Clone("fqqa");
	TH3D *fgg = (TH3D*)fqqs_plus_4jet->Clone("fgg");
	TH3D *fbck = (TH3D*)fqqs_plus_4jet->Clone("fbck");
	TH3D *fWJets = (TH3D*)fqqs_plus_4jet->Clone("fWJets");
	for (double x = x_low+(0.5*xbinwidth); x<x_high; x=x+xbinwidth) {
		for (double y = y_low+(0.5*ybinwidth); y<y_high; y=y+ybinwidth) {
			for (double z = z_low+(0.5*zbinwidth); z<z_high; z=z+zbinwidth) {
				double bck_4jet = Rbck_4jet*(fbck_plus_4jet->GetBinContent(fbck_plus_4jet->FindFixBin(x,y,z)) + fbck_minus_4jet->GetBinContent(fbck_minus_4jet->FindFixBin(x,y,z)));
				double bck_5jet = Rbck_5jet*(fbck_plus_5jet->GetBinContent(fbck_plus_5jet->FindFixBin(x,y,z)) + fbck_minus_5jet->GetBinContent(fbck_minus_4jet->FindFixBin(x,y,z)));
				double WJets_4jet = RW4Jets*(WJets_plus_4jet->GetBinContent(WJets_plus_4jet->FindFixBin(x,y,z)) + WJets_minus_4jet->GetBinContent(WJets_minus_4jet->FindFixBin(x,y,z)));
				double WJets_5jet = RW5Jets*(WJets_plus_5jet->GetBinContent(WJets_plus_5jet->FindFixBin(x,y,z)) + WJets_minus_5jet->GetBinContent(WJets_minus_4jet->FindFixBin(x,y,z)));
				double gg_4jet = (1.0-Rbck_4jet-RW4Jets)*(1.0-Rqqbar_4jet)*(fgg_plus_4jet->GetBinContent(fgg_plus_4jet->FindFixBin(x,y,z)) + fgg_minus_4jet->GetBinContent(fgg_minus_4jet->FindFixBin(x,y,z)));
				double gg_5jet = (1.0-Rbck_5jet-RW5Jets)*(1.0-Rqqbar_5jet)*(fgg_plus_5jet->GetBinContent(fgg_plus_5jet->FindFixBin(x,y,z)) + fgg_minus_5jet->GetBinContent(fgg_minus_5jet->FindFixBin(x,y,z)));
				double qqs_4jet = (1.0-Rbck_4jet-RW4Jets)*(Rqqbar_4jet/(1.0+xi_4jet*F_xi_4jet+delta_4jet*F_delta_4jet));
				double tmp = fqqs_plus_4jet->GetBinContent(fqqs_plus_4jet->FindFixBin(x,y,z)) + fqqs_minus_4jet->GetBinContent(fqqs_minus_4jet->FindFixBin(x,y,z));
				tmp = tmp + xi_4jet*(fqqs_xi_plus_4jet->GetBinContent(fqqs_xi_plus_4jet->FindFixBin(x,y,z)) + fqqs_xi_minus_4jet->GetBinContent(fqqs_xi_minus_4jet->FindFixBin(x,y,z)));
				tmp = tmp + delta_4jet*(fqqs_delta_plus_4jet->GetBinContent(fqqs_delta_plus_4jet->FindFixBin(x,y,z)) + fqqs_delta_minus_4jet->GetBinContent(fqqs_delta_minus_4jet->FindFixBin(x,y,z)));
				qqs_4jet = qqs_4jet*tmp;
				double qqs_5jet = (1.0-Rbck_5jet-RW5Jets)*(Rqqbar_5jet/(1.0+xi_5jet*F_xi_5jet+delta_5jet*F_delta_5jet));
				tmp = fqqs_plus_5jet->GetBinContent(fqqs_plus_5jet->FindFixBin(x,y,z)) + fqqs_minus_5jet->GetBinContent(fqqs_minus_5jet->FindFixBin(x,y,z));
				tmp = tmp + xi_5jet*(fqqs_xi_plus_5jet->GetBinContent(fqqs_xi_plus_5jet->FindFixBin(x,y,z)) + fqqs_xi_minus_5jet->GetBinContent(fqqs_xi_minus_5jet->FindFixBin(x,y,z)));
				tmp = tmp + delta_5jet*(fqqs_delta_plus_5jet->GetBinContent(fqqs_delta_plus_5jet->FindFixBin(x,y,z)) + fqqs_delta_minus_5jet->GetBinContent(fqqs_delta_minus_5jet->FindFixBin(x,y,z)));
				qqs_5jet = qqs_5jet*tmp;
				double qqa_4jet = Afb_4jet*(1.0-Rbck_4jet-RW4Jets)*(Rqqbar_4jet/(1.0+xi_4jet*F_xi_4jet+delta_4jet*F_delta_4jet));
				tmp = fqqa_plus_4jet->GetBinContent(fqqa_plus_4jet->FindFixBin(x,y,z)) + fqqa_minus_4jet->GetBinContent(fqqa_minus_4jet->FindFixBin(x,y,z));
				tmp = tmp + xi_4jet*(fqqa_xi_plus_4jet->GetBinContent(fqqa_xi_plus_4jet->FindFixBin(x,y,z)) + fqqa_xi_minus_4jet->GetBinContent(fqqa_xi_minus_4jet->FindFixBin(x,y,z)));
				tmp = tmp + delta_4jet*(fqqa_delta_plus_4jet->GetBinContent(fqqa_delta_plus_4jet->FindFixBin(x,y,z)) + fqqa_delta_minus_4jet->GetBinContent(fqqa_delta_minus_4jet->FindFixBin(x,y,z)));
				qqa_4jet = qqa_4jet*tmp;
				double qqa_5jet = Afb_5jet*(1.0-Rbck_5jet-RW5Jets)*(Rqqbar_5jet/(1.0+xi_5jet*F_xi_5jet+delta_5jet*F_delta_5jet));
				tmp = fqqa_plus_5jet->GetBinContent(fqqa_plus_5jet->FindFixBin(x,y,z)) + fqqa_minus_5jet->GetBinContent(fqqa_minus_5jet->FindFixBin(x,y,z));
				tmp = tmp + xi_5jet*(fqqa_xi_plus_5jet->GetBinContent(fqqa_xi_plus_5jet->FindFixBin(x,y,z)) + fqqa_xi_minus_5jet->GetBinContent(fqqa_xi_minus_5jet->FindFixBin(x,y,z)));
				tmp = tmp + delta_5jet*(fqqa_delta_plus_5jet->GetBinContent(fqqa_delta_plus_5jet->FindFixBin(x,y,z)) + fqqa_delta_minus_5jet->GetBinContent(fqqa_delta_minus_5jet->FindFixBin(x,y,z)));
				qqa_5jet = qqa_5jet*tmp;
				double qqs = count_added_4jet*qqs_4jet+count_added_5jet*qqs_5jet;
				double qqa = count_added_4jet*qqa_4jet+count_added_5jet*qqa_5jet;
				double gg  = count_added_4jet*gg_4jet+count_added_5jet*gg_5jet;
				double bck = count_added_4jet*bck_4jet+count_added_5jet*bck_5jet;
				double WJets = count_added_4jet*WJets_4jet+count_added_5jet*WJets_5jet;
				fqqs->SetBinContent(fqqs->FindFixBin(x,y,z),qqs);
				fqqa->SetBinContent(fqqa->FindFixBin(x,y,z),qqa);
				fgg->SetBinContent(fgg->FindFixBin(x,y,z),gg);
				fbck->SetBinContent(fbck->FindFixBin(x,y,z),bck);
				fWJets->SetBinContent(fWJets->FindFixBin(x,y,z),WJets);
			}
		}
	}//end loops over bins
	
	//second set of loops to make the fit histograms
	for (double x = x_low+(0.5*xbinwidth); x<x_high; x=x+xbinwidth) {
		for (double y = y_low+(0.5*ybinwidth); y<y_high; y=y+ybinwidth) {
			for (double z = z_low+(0.5*zbinwidth); z<z_high; z=z+zbinwidth) {
				double qqs = fqqs->GetBinContent(fqqs->FindFixBin(x,y,z));
				double qqa = fqqa->GetBinContent(fqqa->FindFixBin(x,y,z));
				double gg  = fgg->GetBinContent(fgg->FindFixBin(x,y,z));
				double bck = fbck->GetBinContent(fbck->FindFixBin(x,y,z));
				double WJets = fWJets->GetBinContent(fWJets->FindFixBin(x,y,z));
				gg_x->Fill(x,gg);
				gg_y->Fill(y,gg);
				gg_z->Fill(z,gg);
				qq_x->Fill(x,(qqs+qqa));
				qq_y->Fill(y,(qqs+qqa));
				qq_z->Fill(z,(qqs+qqa));
				bg_x->Fill(x,bck);
				bg_y->Fill(y,bck);
				bg_z->Fill(z,bck);
				wj_x->Fill(x,WJets);
				wj_y->Fill(y,WJets);
				wj_z->Fill(z,WJets);
			}
		}
	}//end loops over bins
	
	//and a few more loops over the histograms to set the errors on the projections for the data
	for (double x = x_low+(0.5*xbinwidth); x<x_high; x=x+xbinwidth) {
		double x_bin_var_data=0.0;
		for (double y = y_low+(0.5*ybinwidth); y<y_high; y=y+ybinwidth) {
			double y_bin_var_data=0.0;
			for (double z = z_low+(0.5*zbinwidth); z<z_high; z=z+zbinwidth) {
				y_bin_var_data=y_bin_var_data+data_hist->GetBinError(data_hist->FindFixBin(x,y,z))*data_hist->GetBinError(data_hist->FindFixBin(x,y,z));
			}
			x_bin_var_data=x_bin_var_data+y_bin_var_data;
		}
		data_x->SetBinError(data_x->FindFixBin(x),TMath::Sqrt(x_bin_var_data));
	}
	for (double z = z_low+(0.5*zbinwidth); z<z_high; z=z+zbinwidth) {
		double z_bin_var_data=0.0;
		for (double x = x_low+(0.5*xbinwidth); x<x_high; x=x+xbinwidth) {
			double x_bin_var_data=0.0;
			for (double y = y_low+(0.5*ybinwidth); y<y_high; y=y+ybinwidth) {
				x_bin_var_data=x_bin_var_data+data_hist->GetBinError(data_hist->FindFixBin(x,y,z))*data_hist->GetBinError(data_hist->FindFixBin(x,y,z));
			}
			z_bin_var_data=z_bin_var_data+x_bin_var_data;
		}
		data_z->SetBinError(data_z->FindFixBin(z),TMath::Sqrt(z_bin_var_data));
	}
	for (double y = y_low+(0.5*ybinwidth); y<y_high; y=y+ybinwidth) {
		double y_bin_var_data=0.0;
		for (double z = z_low+(0.5*zbinwidth); z<z_high; z=z+zbinwidth) {
			double z_bin_var_data=0.0;
			for (double x = x_low+(0.5*xbinwidth); x<x_high; x=x+xbinwidth) {
				z_bin_var_data=z_bin_var_data+data_hist->GetBinError(data_hist->FindFixBin(x,y,z))*data_hist->GetBinError(data_hist->FindFixBin(x,y,z));
			}
			y_bin_var_data=y_bin_var_data+z_bin_var_data;
		}
		data_y->SetBinError(data_y->FindFixBin(y),TMath::Sqrt(y_bin_var_data));
	}

	//Build the event type distributions for the background, gg, and qq
	double n_fit_bck_plus_4jet = (Rbck_4jet*fbck_plus_4jet->Integral()); 
	double n_fit_WJets_plus_4jet = (RW4Jets*WJets_plus_4jet->Integral()); 
	double n_fit_gg_plus_4jet  = (1.0-Rbck_4jet-RW4Jets)*((1.0-Rqqbar_4jet)*fgg_plus_4jet->Integral());
	double n_fit_qq_plus_4jet  = (1.0-Rbck_4jet-RW4Jets)*(Rqqbar_4jet/(1.0+xi_4jet*F_xi_4jet+delta_4jet*F_delta_4jet))*(fqqs_plus_4jet->Integral()+Afb_4jet*fqqa_plus_4jet->Integral());
		   n_fit_qq_plus_4jet += (1.0-Rbck_4jet-RW4Jets)*(Rqqbar_4jet/(1.0+xi_4jet*F_xi_4jet+delta_4jet*F_delta_4jet))*xi_4jet*(fqqs_xi_plus_4jet->Integral()+Afb_4jet*fqqa_xi_plus_4jet->Integral());
		   n_fit_qq_plus_4jet += (1.0-Rbck_4jet-RW4Jets)*(Rqqbar_4jet/(1.0+xi_4jet*F_xi_4jet+delta_4jet*F_delta_4jet))*delta_4jet*(fqqs_delta_plus_4jet->Integral()+Afb_4jet*fqqa_delta_plus_4jet->Integral());
	double n_fit_bck_minus_4jet = (Rbck_4jet*fbck_minus_4jet->Integral()); 
	double n_fit_WJets_minus_4jet = (RW4Jets*WJets_minus_4jet->Integral()); 
	double n_fit_gg_minus_4jet  = (1.0-Rbck_4jet-RW4Jets)*((1.0-Rqqbar_4jet)*fgg_minus_4jet->Integral());
	double n_fit_qq_minus_4jet  = (1.0-Rbck_4jet-RW4Jets)*(Rqqbar_4jet/(1.0+xi_4jet*F_xi_4jet+delta_4jet*F_delta_4jet))*(fqqs_minus_4jet->Integral()+Afb_4jet*fqqa_minus_4jet->Integral());
		   n_fit_qq_minus_4jet += (1.0-Rbck_4jet-RW4Jets)*(Rqqbar_4jet/(1.0+xi_4jet*F_xi_4jet+delta_4jet*F_delta_4jet))*xi_4jet*(fqqs_xi_minus_4jet->Integral()+Afb_4jet*fqqa_xi_minus_4jet->Integral());
		   n_fit_qq_minus_4jet += (1.0-Rbck_4jet-RW4Jets)*(Rqqbar_4jet/(1.0+xi_4jet*F_xi_4jet+delta_4jet*F_delta_4jet))*delta_4jet*(fqqs_delta_minus_4jet->Integral()+Afb_4jet*fqqa_delta_minus_4jet->Integral());
	double n_fit_bck_plus_5jet = (Rbck_5jet*fbck_plus_5jet->Integral()); 
	double n_fit_WJets_plus_5jet = (RW5Jets*WJets_plus_5jet->Integral()); 
	double n_fit_gg_plus_5jet  = (1.0-Rbck_5jet-RW5Jets)*((1.0-Rqqbar_5jet)*fgg_plus_5jet->Integral());
	double n_fit_qq_plus_5jet  = (1.0-Rbck_5jet-RW5Jets)*(Rqqbar_5jet/(1.0+xi_5jet*F_xi_5jet+delta_5jet*F_delta_5jet))*(fqqs_plus_5jet->Integral()+Afb_5jet*fqqa_plus_5jet->Integral());
		   n_fit_qq_plus_5jet += (1.0-Rbck_5jet-RW5Jets)*(Rqqbar_5jet/(1.0+xi_5jet*F_xi_5jet+delta_5jet*F_delta_5jet))*xi_5jet*(fqqs_xi_plus_5jet->Integral()+Afb_5jet*fqqa_xi_plus_5jet->Integral());
		   n_fit_qq_plus_5jet += (1.0-Rbck_5jet-RW5Jets)*(Rqqbar_5jet/(1.0+xi_5jet*F_xi_5jet+delta_5jet*F_delta_5jet))*delta_5jet*(fqqs_delta_plus_5jet->Integral()+Afb_5jet*fqqa_delta_plus_5jet->Integral());
	double n_fit_bck_minus_5jet = (Rbck_5jet*fbck_minus_5jet->Integral()); 
	double n_fit_WJets_minus_5jet = (RW5Jets*WJets_minus_5jet->Integral());
	double n_fit_gg_minus_5jet  = (1.0-Rbck_5jet-RW5Jets)*((1.0-Rqqbar_5jet)*fgg_minus_5jet->Integral());
	double n_fit_qq_minus_5jet  = (1.0-Rbck_5jet-RW5Jets)*(Rqqbar_5jet/(1.0+xi_5jet*F_xi_5jet+delta_5jet*F_delta_5jet))*(fqqs_minus_5jet->Integral()+Afb_5jet*fqqa_minus_5jet->Integral());
		   n_fit_qq_minus_5jet += (1.0-Rbck_5jet-RW5Jets)*(Rqqbar_5jet/(1.0+xi_5jet*F_xi_5jet+delta_5jet*F_delta_5jet))*xi_5jet*(fqqs_xi_minus_5jet->Integral()+Afb_5jet*fqqa_xi_minus_5jet->Integral());
		   n_fit_qq_minus_5jet += (1.0-Rbck_5jet-RW5Jets)*(Rqqbar_5jet/(1.0+xi_5jet*F_xi_5jet+delta_5jet*F_delta_5jet))*delta_5jet*(fqqs_delta_minus_5jet->Integral()+Afb_5jet*fqqa_delta_minus_5jet->Integral());
	event_numbers_bck->Fill(1.5,count_added_4jet*n_fit_bck_plus_4jet); event_numbers_bck->Fill(2.5,count_added_4jet*n_fit_bck_minus_4jet);
	event_numbers_bck->Fill(3.5,count_added_5jet*n_fit_bck_plus_5jet); event_numbers_bck->Fill(4.5,count_added_5jet*n_fit_bck_minus_5jet);
	event_numbers_WJets->Fill(1.5,count_added_4jet*n_fit_WJets_plus_4jet); event_numbers_WJets->Fill(2.5,count_added_4jet*n_fit_WJets_minus_4jet);
	event_numbers_WJets->Fill(3.5,count_added_5jet*n_fit_WJets_plus_5jet); event_numbers_WJets->Fill(4.5,count_added_5jet*n_fit_WJets_minus_5jet);
	event_numbers_gg->Fill(1.5,count_added_4jet*n_fit_gg_plus_4jet); event_numbers_gg->Fill(2.5,count_added_4jet*n_fit_gg_minus_4jet);
	event_numbers_gg->Fill(3.5,count_added_5jet*n_fit_gg_plus_5jet); event_numbers_gg->Fill(4.5,count_added_5jet*n_fit_gg_minus_5jet);
	event_numbers_qq->Fill(1.5,count_added_4jet*n_fit_qq_plus_4jet); event_numbers_qq->Fill(2.5,count_added_4jet*n_fit_qq_minus_4jet);
	event_numbers_qq->Fill(3.5,count_added_5jet*n_fit_qq_plus_5jet); event_numbers_qq->Fill(4.5,count_added_5jet*n_fit_qq_minus_5jet);
	
	//plot the new histograms on top of each other
	//build histogram stacks
	THStack *x_stack = new THStack("x_stack","cos(#theta *) Comparison; cos(#theta *)");
	THStack *y_stack = new THStack("y_stack","Feynman x Comparison; (x_{F})");
	THStack *z_stack = new THStack("z_stack","t #bar{t} Mass Comparison; M_{t #bar{t}}");
	THStack *event_numbers_stack = new THStack("event_numbers_stack","Lepton Charge and Jet Multiplicity Comparison");
	//add to stacks
	x_stack->Add(bg_x);
	y_stack->Add(bg_y);
	z_stack->Add(bg_z);
	x_stack->Add(wj_x);
	y_stack->Add(wj_y);
	z_stack->Add(wj_z);
	x_stack->Add(gg_x);
	y_stack->Add(gg_y);
	z_stack->Add(gg_z);
	x_stack->Add(qq_x);
	y_stack->Add(qq_y);
	z_stack->Add(qq_z);
	event_numbers_stack->Add(event_numbers_bck);
	event_numbers_stack->Add(event_numbers_WJets);
	event_numbers_stack->Add(event_numbers_gg);
	event_numbers_stack->Add(event_numbers_qq);
	//set drawing options for projections
	//glueglue
	gg_x->SetFillColor(38);
	gg_x->SetMarkerStyle(21);
	gg_x->SetMarkerColor(38);
	gg_y->SetFillColor(38);
	gg_y->SetMarkerStyle(21);
	gg_y->SetMarkerColor(38);
	gg_z->SetFillColor(38);
	gg_z->SetMarkerStyle(21);
	gg_z->SetMarkerColor(38);
	event_numbers_gg->SetFillColor(38);
	event_numbers_gg->SetMarkerStyle(21);
	event_numbers_gg->SetMarkerColor(38);
	//qqbar
	qq_x->SetFillColor(46);
	qq_x->SetMarkerStyle(21);
	qq_x->SetMarkerColor(46);
	qq_y->SetFillColor(46);
	qq_y->SetMarkerStyle(21);
	qq_y->SetMarkerColor(46);
	qq_z->SetFillColor(46);
	qq_z->SetMarkerStyle(21);
	qq_z->SetMarkerColor(46);
	event_numbers_qq->SetFillColor(46);
	event_numbers_qq->SetMarkerStyle(21);
	event_numbers_qq->SetMarkerColor(46);
	//background
	bg_x->SetFillColor(41);
	bg_x->SetMarkerStyle(21);
	bg_x->SetMarkerColor(41);
	bg_y->SetFillColor(41);
	bg_y->SetMarkerStyle(21);
	bg_y->SetMarkerColor(41);
	bg_z->SetFillColor(41);
	bg_z->SetMarkerStyle(21);
	bg_z->SetMarkerColor(41);
	event_numbers_bck->SetFillColor(41);
	event_numbers_bck->SetMarkerStyle(21);
	event_numbers_bck->SetMarkerColor(41);
	//WJets
	wj_x->SetFillColor(kGreen-3);
	wj_x->SetMarkerStyle(21);
	wj_x->SetMarkerColor(kGreen-3);
	wj_y->SetFillColor(kGreen-3);
	wj_y->SetMarkerStyle(21);
	wj_y->SetMarkerColor(kGreen-3);
	wj_z->SetFillColor(kGreen-3);
	wj_z->SetMarkerStyle(21);
	wj_z->SetMarkerColor(kGreen-3);
	event_numbers_WJets->SetFillColor(kGreen-3);
	event_numbers_WJets->SetMarkerStyle(21);
	event_numbers_WJets->SetMarkerColor(kGreen-3);
	//data
	data_x->SetMarkerColor(kBlack);
	data_y->SetMarkerColor(kBlack);
	data_z->SetMarkerColor(kBlack);
	data_x->SetLineWidth(2);
	data_y->SetLineWidth(2);
	data_z->SetLineWidth(2);
	event_numbers_data->SetMarkerColor(kBlack);
	event_numbers_data->SetLineWidth(2);
	//add a legend to the plot
	// TLegend *leg = new TLegend(0.56,0.73,0.88,0.88);
	TLegend *leg = new TLegend(0.14,0.73,0.46,0.88);
	leg->SetName("legend");
	leg->AddEntry(data_x,"Data","LPEX0");
	leg->AddEntry(qq_x,"q#bar{q} #rightarrow t#bar{t}(j)","F");
	leg->AddEntry(gg_x,"gg(qg) #rightarrow t#bar{t}(j)","F");
	leg->AddEntry(wj_x,"WJets Background","F");
	leg->AddEntry(bg_x,"Other Backgrounds","F");
	TCanvas* c = new TCanvas("c","Fit Comparison",1200,1200);
	c->Divide(2,2);
	c->cd(1);
	event_numbers_stack->SetMaximum(1.1*event_numbers_data->GetMaximum());
	event_numbers_stack->Draw("bar1");
	event_numbers_stack->GetXaxis()->SetBinLabel(event_numbers_stack->GetXaxis()->FindFixBin(1.5),"4jets, l+");
	event_numbers_stack->GetXaxis()->SetBinLabel(event_numbers_stack->GetXaxis()->FindFixBin(2.5),"4jets, l-");
	event_numbers_stack->GetXaxis()->SetBinLabel(event_numbers_stack->GetXaxis()->FindFixBin(3.5),"5jets, l+");
	event_numbers_stack->GetXaxis()->SetBinLabel(event_numbers_stack->GetXaxis()->FindFixBin(4.5),"5jets, l-");
	event_numbers_data->Draw("SAME PE1X0");
	leg->Draw();
	c->cd(2);
	y_stack->SetMaximum(1.1*data_y->GetMaximum());
	y_stack->Draw();
	data_y->Draw("SAME PE1X0");
	c->cd(3);
	x_stack->SetMaximum(1.1*data_x->GetMaximum());
	x_stack->Draw();
	data_x->Draw("SAME PE1X0");	
	c->cd(4);
	z_stack->SetMaximum(1.1*data_z->GetMaximum());
	z_stack->Draw();
	data_z->Draw("SAME PE1X0");
	//save the plots
	c->Print("fit_comparison.pdf","pdf");	

	// Save final stacks and data hist and legend into a root file for post processing
	printf("Writing final stack and data hists into a root file!\n");
	TFile* stack_file = new TFile("final_stack.root","Recreate");
	stack_file->cd();
	x_stack->Write();data_x->Write();
	y_stack->Write();data_y->Write();
	z_stack->Write();data_z->Write();
	event_numbers_stack->Write();event_numbers_data->Write();
	leg->Write();
	stack_file->Write();
	stack_file->Close();

	TCanvas* c1 = new TCanvas("c1","Sideband",900,900);
	c1->Divide(2,2,0.01,0.01,0);
	c1->cd(1);
	sideband->Draw("BOX");
	c1->cd(2);
	(sideband->ProjectionY())->Draw();
	c1->cd(3);
	(sideband->ProjectionX())->Draw();
	c1->cd(4);
	(sideband->ProjectionZ())->Draw();
	c1->Print("sideband.pdf","pdf");
	
	//finally, save the results of the fit to a text file.
	FILE* o = fopen("fit_results.txt","w");
	fprintf(o,"*************		FIT       RESULTS		*************\n");
	fprintf(o,"Rqqbar_4jet = %-.4f +/- %-.4f\n",Rqqbar_4jet,sigma_Rqqbar_4jet);
	fprintf(o,"Rbck_4jet = %-.4f +/- %-.4f\n",Rbck_4jet,sigma_Rbck_4jet);
	fprintf(o,"RW4Jets = %-.4f +/- %-.4f\n",RW4Jets,sigma_RW4Jets);
	fprintf(o,"Rqqbar_5jet = %-.4f +/- %-.4f\n",Rqqbar_5jet,sigma_Rqqbar_5jet);
	fprintf(o,"Rbck_5jet = %-.4f +/- %-.4f\n",Rbck_5jet,sigma_Rbck_5jet);
	fprintf(o,"RW5Jets = %-.4f +/- %-.4f\n",RW5Jets,sigma_RW5Jets);
	fprintf(o,"xi_4jet = %-.4f +/- %-.4f\n",xi_4jet,sigma_xi_4jet);
	fprintf(o,"xi_5jet = %-.4f +/- %-.4f\n",xi_5jet,sigma_xi_5jet);
	fprintf(o,"delta_4jet = %-.4f +/- %-.4f\n",delta_4jet,sigma_delta_4jet);
	fprintf(o,"delta_5jet = %-.4f +/- %-.4f\n",delta_5jet,sigma_delta_5jet);
	fprintf(o,"Afb_4jet    = %-.4f +/- %-.4f\n",Afb_4jet,sigma_Afb_4jet);
	fprintf(o,"Afb_5jet    = %-.4f +/- %-.4f\n",Afb_5jet,sigma_Afb_5jet);
	fprintf(o,"\n\n");
	fprintf(o,"Runname is %s\n",rname);
	fprintf(o,"\n");
	fprintf(o,"Binning and limits:\n");
	fprintf(o,"    %-.2f < costheta* < %-.2f in %d bins\n",x_low,x_high,nbinsx);
	fprintf(o,"    %-.2f < Feynman x < %-.2f in %d bins\n",y_low,y_high,nbinsy);
	fprintf(o,"    %-.2f < ttbarMass < %-.2f in %d bins\n",z_low,z_high,nbinsz);
	fclose(o);
	//and save an abridged version to the summary file
	o = fopen("summary_separated.txt","a");
	fprintf(o,"%s	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f\n",
			rname,Rqqbar_4jet,sigma_Rqqbar_4jet,Rbck_4jet,sigma_Rbck_4jet,RW4Jets,sigma_RW4Jets,Afb_4jet,sigma_Afb_4jet,
				  Rqqbar_5jet,sigma_Rqqbar_5jet,Rbck_5jet,sigma_Rbck_5jet,RW5Jets,sigma_RW5Jets,Afb_5jet,sigma_Afb_5jet);
	fprintf(o,"		%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%-.4f +/- %-.4f	%d	%.4f\n",
		xi_4jet,sigma_xi_4jet,delta_4jet,sigma_delta_4jet,xi_5jet,sigma_xi_5jet,delta_5jet,sigma_delta_5jet,(count_added_4jet+count_added_5jet),
		Loop(Rqqbar_4jet, Rbck_4jet, RW4Jets, Rqqbar_5jet, Rbck_5jet, RW5Jets, xi_4jet, xi_5jet, delta_4jet, delta_5jet, Afb_4jet, Afb_5jet));
	fclose(o);
	o = fopen("summary_separated.csv","a");
	fprintf(o,"%s,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%-.4f,%d,%.4f\n",
			rname,Rqqbar_4jet,sigma_Rqqbar_4jet,Rbck_4jet,sigma_Rbck_4jet,RW4Jets,sigma_RW4Jets,Afb_4jet,sigma_Afb_4jet,
				  Rqqbar_5jet,sigma_Rqqbar_5jet,Rbck_5jet,sigma_Rbck_5jet,RW5Jets,sigma_RW5Jets,Afb_5jet,sigma_Afb_5jet,
				  xi_4jet,sigma_xi_4jet,delta_4jet,sigma_delta_4jet,xi_5jet,sigma_xi_5jet,delta_5jet,sigma_delta_5jet,(count_added_4jet+count_added_5jet),
				  Loop(Rqqbar_4jet, Rbck_4jet, RW4Jets, Rqqbar_5jet, Rbck_5jet, RW5Jets, xi_4jet, xi_5jet, delta_4jet, delta_5jet, Afb_4jet, Afb_5jet));
	fclose(o);
	return 0;
}

//loads the histograms once for each object created
//COMBINED CASE
void angles_data::LoadHistogramsCombined() {
	//open the file
	TFile *f = new TFile(template_filename);
	//load histograms into object variables
	fqqs_plus_4jet = (TH3D*)f->Get("f_qqs_plus_4jet");
	fqqs_xi_plus_4jet = (TH3D*)f->Get("f_qqs_xi_plus_4jet");
	fqqs_delta_plus_4jet = (TH3D*)f->Get("f_qqs_delta_plus_4jet");
	fqqa_plus_4jet = (TH3D*)f->Get("f_qqa_plus_4jet");
	fqqa_xi_plus_4jet = (TH3D*)f->Get("f_qqa_xi_plus_4jet");
	fqqa_delta_plus_4jet = (TH3D*)f->Get("f_qqa_delta_plus_4jet");
	fqqs_minus_4jet = (TH3D*)f->Get("f_qqs_minus_4jet");
	fqqs_xi_minus_4jet = (TH3D*)f->Get("f_qqs_xi_minus_4jet");
	fqqs_delta_minus_4jet = (TH3D*)f->Get("f_qqs_delta_minus_4jet");
	fqqa_minus_4jet = (TH3D*)f->Get("f_qqa_minus_4jet");
	fqqa_xi_minus_4jet = (TH3D*)f->Get("f_qqa_xi_minus_4jet");
	fqqa_delta_minus_4jet = (TH3D*)f->Get("f_qqa_delta_minus_4jet");
	fqqs_plus_5jet = (TH3D*)f->Get("f_qqs_plus_5jet");
	fqqs_xi_plus_5jet = (TH3D*)f->Get("f_qqs_xi_plus_5jet");
	fqqs_delta_plus_5jet = (TH3D*)f->Get("f_qqs_delta_plus_5jet");
	fqqa_plus_5jet = (TH3D*)f->Get("f_qqa_plus_5jet");
	fqqa_xi_plus_5jet = (TH3D*)f->Get("f_qqa_xi_plus_5jet");
	fqqa_delta_plus_5jet = (TH3D*)f->Get("f_qqa_delta_plus_5jet");
	fqqs_minus_5jet = (TH3D*)f->Get("f_qqs_minus_5jet");
	fqqs_xi_minus_5jet = (TH3D*)f->Get("f_qqs_xi_minus_5jet");
	fqqs_delta_minus_5jet = (TH3D*)f->Get("f_qqs_delta_minus_5jet");
	fqqa_minus_5jet = (TH3D*)f->Get("f_qqa_minus_5jet");
	fqqa_xi_minus_5jet = (TH3D*)f->Get("f_qqa_xi_minus_5jet");
	fqqa_delta_minus_5jet = (TH3D*)f->Get("f_qqa_delta_minus_5jet");
	fgg_plus_4jet = (TH3D*)f->Get("f_gg_plus_4jet");
	fgg_minus_4jet = (TH3D*)f->Get("f_gg_minus_4jet");
	fgg_plus_5jet = (TH3D*)f->Get("f_gg_plus_5jet");
	fgg_minus_5jet = (TH3D*)f->Get("f_gg_minus_5jet");
	TFile *f2 = new TFile(histo_filename);
	gtt_4jet = (TH1D*)f2->Get("gtt_4jet_local");
	gbk_4jet = (TH1D*)f2->Get("gbk_4jet_local");
	gtt_5jet = (TH1D*)f2->Get("gtt_5jet_local");
	gbk_5jet = (TH1D*)f2->Get("gbk_5jet_local");
	//find limits of training set histograms
	nbins_lnL = gtt_4jet->GetNbinsX();
	lnL_low = gtt_4jet->GetXaxis()->GetXmin(); lnL_high = gtt_4jet->GetXaxis()->GetXmax();
	//printf("lnL_high = %.2f\n",lnL_high);
	nbinsx = fqqs_plus_4jet->GetNbinsX(); nbinsy = fqqs_plus_4jet->GetNbinsY(); nbinsz = fqqs_plus_4jet->GetNbinsZ();
	x_low = fqqs_plus_4jet->GetXaxis()->GetXmin(); x_high = fqqs_plus_4jet->GetXaxis()->GetXmax();
	y_low = fqqs_plus_4jet->GetYaxis()->GetXmin(); y_high = fqqs_plus_4jet->GetYaxis()->GetXmax();
	z_low = fqqs_plus_4jet->GetZaxis()->GetXmin(); z_high = fqqs_plus_4jet->GetZaxis()->GetXmax();
	xbinwidth = fqqs_plus_4jet->GetXaxis()->GetBinWidth(1);
	ybinwidth = fqqs_plus_4jet->GetYaxis()->GetBinWidth(1);
	zbinwidth = fqqs_plus_4jet->GetZaxis()->GetBinWidth(1);
	//load up the background distributions
	fbck_plus_4jet  = (TH3D*)f->Get("f_bck_plus_4jet");
	fbck_minus_4jet = (TH3D*)f->Get("f_bck_minus_4jet");
	fbck_plus_5jet  = (TH3D*)f->Get("f_bck_plus_5jet");
	fbck_minus_5jet = (TH3D*)f->Get("f_bck_minus_5jet");
	WJets_plus_4jet  = (TH3D*)f->Get("f_WJets_plus_4jet");
	WJets_minus_4jet = (TH3D*)f->Get("f_WJets_minus_4jet");
	WJets_plus_5jet  = (TH3D*)f->Get("f_WJets_plus_5jet");
	WJets_minus_5jet = (TH3D*)f->Get("f_WJets_minus_5jet");
	//Normalize the distributions
	printf("NORMALIZING DISTRIBUTIONS\n");
	double qq_rescale = 1.0/(fqqs_plus_4jet->Integral()+fqqs_minus_4jet->Integral()+fqqs_plus_5jet->Integral()+fqqs_minus_5jet->Integral());
	fqqs_plus_4jet->Scale(qq_rescale);
	fqqs_xi_plus_4jet->Scale(qq_rescale);
	fqqs_delta_plus_4jet->Scale(qq_rescale);
	fqqs_minus_4jet->Scale(qq_rescale);
	fqqs_xi_minus_4jet->Scale(qq_rescale);
	fqqs_delta_minus_4jet->Scale(qq_rescale);
	fqqa_plus_4jet->Scale(qq_rescale);
	fqqa_xi_plus_4jet->Scale(qq_rescale);
	fqqa_delta_plus_4jet->Scale(qq_rescale);
	fqqa_minus_4jet->Scale(qq_rescale);
	fqqa_xi_minus_4jet->Scale(qq_rescale);
	fqqa_delta_minus_4jet->Scale(qq_rescale);
	fqqs_plus_5jet->Scale(qq_rescale);
	fqqs_xi_plus_5jet->Scale(qq_rescale);
	fqqs_delta_plus_5jet->Scale(qq_rescale);
	fqqs_minus_5jet->Scale(qq_rescale);
	fqqs_xi_minus_5jet->Scale(qq_rescale);
	fqqs_delta_minus_5jet->Scale(qq_rescale);
	fqqa_plus_5jet->Scale(qq_rescale);
	fqqa_xi_plus_5jet->Scale(qq_rescale);
	fqqa_delta_plus_5jet->Scale(qq_rescale);
	fqqa_minus_5jet->Scale(qq_rescale);
	fqqa_xi_minus_5jet->Scale(qq_rescale);
	fqqa_delta_minus_5jet->Scale(qq_rescale);
	double gg_rescale = 1.0/(fgg_plus_4jet->Integral()+fgg_minus_4jet->Integral()+fgg_plus_5jet->Integral()+fgg_minus_5jet->Integral());
	fgg_plus_4jet->Scale(gg_rescale);
	fgg_minus_4jet->Scale(gg_rescale);
	fgg_plus_5jet->Scale(gg_rescale);
	fgg_minus_5jet->Scale(gg_rescale);
	double gtt_rescale = 1.0/(gtt_4jet->Integral()+gtt_5jet->Integral());
	gtt_4jet->Scale(gtt_rescale);
	gtt_5jet->Scale(gtt_rescale);
	double gbk_rescale = 1.0/(gbk_4jet->Integral()+gbk_5jet->Integral());
	gbk_4jet->Scale(gbk_rescale);
	gbk_5jet->Scale(gbk_rescale);
	double bck_rescale = 1.0/(fbck_plus_4jet->Integral()+fbck_minus_4jet->Integral()+fbck_plus_5jet->Integral()+fbck_minus_5jet->Integral());
	fbck_plus_4jet->Scale(bck_rescale);
	fbck_minus_4jet->Scale(bck_rescale);
	fbck_plus_5jet->Scale(bck_rescale);
	fbck_minus_5jet->Scale(bck_rescale);
	double WJets_rescale = 1.0/(WJets_plus_4jet->Integral()+WJets_minus_4jet->Integral()+WJets_plus_5jet->Integral()+WJets_minus_5jet->Integral());
	WJets_plus_4jet->Scale(WJets_rescale);
	WJets_minus_4jet->Scale(WJets_rescale);
	WJets_plus_5jet->Scale(WJets_rescale);
	WJets_minus_5jet->Scale(WJets_rescale);
	F_xi_comb = fqqs_xi_plus_4jet->Integral()+fqqs_xi_minus_4jet->Integral()+fqqs_xi_plus_5jet->Integral()+fqqs_xi_minus_5jet->Integral();
	F_delta_comb = fqqs_delta_plus_4jet->Integral()+fqqs_delta_minus_4jet->Integral()+fqqs_delta_plus_5jet->Integral()+fqqs_delta_minus_5jet->Integral();
	printf("4jet fraction in MC = %f, 5jet fraction in MC = %f\n", (gtt_4jet->Integral()+gbk_4jet->Integral())/2.0,(gtt_5jet->Integral()+gbk_5jet->Integral())/2.0);
}

//loads the histograms once for each object created
//SEPARATED CASE
void angles_data::LoadHistogramsSeparated() {
	//open the file
	TFile *f = new TFile(template_filename);
	//load histograms into object variables
	fqqs_plus_4jet = (TH3D*)f->Get("f_qqs_plus_4jet");
	fqqs_xi_plus_4jet = (TH3D*)f->Get("f_qqs_xi_plus_4jet");
	fqqs_delta_plus_4jet = (TH3D*)f->Get("f_qqs_delta_plus_4jet");
	fqqa_plus_4jet = (TH3D*)f->Get("f_qqa_plus_4jet");
	fqqa_xi_plus_4jet = (TH3D*)f->Get("f_qqa_xi_plus_4jet");
	fqqa_delta_plus_4jet = (TH3D*)f->Get("f_qqa_delta_plus_4jet");
	fqqs_minus_4jet = (TH3D*)f->Get("f_qqs_minus_4jet");
	fqqs_xi_minus_4jet = (TH3D*)f->Get("f_qqs_xi_minus_4jet");
	fqqs_delta_minus_4jet = (TH3D*)f->Get("f_qqs_delta_minus_4jet");
	fqqa_minus_4jet = (TH3D*)f->Get("f_qqa_minus_4jet");
	fqqa_xi_minus_4jet = (TH3D*)f->Get("f_qqa_xi_minus_4jet");
	fqqa_delta_minus_4jet = (TH3D*)f->Get("f_qqa_delta_minus_4jet");
	fqqs_plus_5jet = (TH3D*)f->Get("f_qqs_plus_5jet");
	fqqs_xi_plus_5jet = (TH3D*)f->Get("f_qqs_xi_plus_5jet");
	fqqs_delta_plus_5jet = (TH3D*)f->Get("f_qqs_delta_plus_5jet");
	fqqa_plus_5jet = (TH3D*)f->Get("f_qqa_plus_5jet");
	fqqa_xi_plus_5jet = (TH3D*)f->Get("f_qqa_xi_plus_5jet");
	fqqa_delta_plus_5jet = (TH3D*)f->Get("f_qqa_delta_plus_5jet");
	fqqs_minus_5jet = (TH3D*)f->Get("f_qqs_minus_5jet");
	fqqs_xi_minus_5jet = (TH3D*)f->Get("f_qqs_xi_minus_5jet");
	fqqs_delta_minus_5jet = (TH3D*)f->Get("f_qqs_delta_minus_5jet");
	fqqa_minus_5jet = (TH3D*)f->Get("f_qqa_minus_5jet");
	fqqa_xi_minus_5jet = (TH3D*)f->Get("f_qqa_xi_minus_5jet");
	fqqa_delta_minus_5jet = (TH3D*)f->Get("f_qqa_delta_minus_5jet");
	fgg_plus_4jet = (TH3D*)f->Get("f_gg_plus_4jet");
	fgg_minus_4jet = (TH3D*)f->Get("f_gg_minus_4jet");
	fgg_plus_5jet = (TH3D*)f->Get("f_gg_plus_5jet");
	fgg_minus_5jet = (TH3D*)f->Get("f_gg_minus_5jet");
	TFile *f2 = new TFile(histo_filename);
	gtt_4jet = (TH1D*)f2->Get("gtt_4jet_local");
	gbk_4jet = (TH1D*)f2->Get("gbk_4jet_local");
	gtt_5jet = (TH1D*)f2->Get("gtt_5jet_local");
	gbk_5jet = (TH1D*)f2->Get("gbk_5jet_local");
	//find limits of training set histograms
	nbins_lnL = gtt_4jet->GetNbinsX();
	lnL_low = gtt_4jet->GetXaxis()->GetXmin(); lnL_high = gtt_4jet->GetXaxis()->GetXmax();
	nbinsx = fqqs_plus_4jet->GetNbinsX(); nbinsy = fqqs_plus_4jet->GetNbinsY(); nbinsz = fqqs_plus_4jet->GetNbinsZ();
	x_low = fqqs_plus_4jet->GetXaxis()->GetXmin(); x_high = fqqs_plus_4jet->GetXaxis()->GetXmax();
	y_low = fqqs_plus_4jet->GetYaxis()->GetXmin(); y_high = fqqs_plus_4jet->GetYaxis()->GetXmax();
	z_low = fqqs_plus_4jet->GetZaxis()->GetXmin(); z_high = fqqs_plus_4jet->GetZaxis()->GetXmax();
	xbinwidth = fqqs_plus_4jet->GetXaxis()->GetBinWidth(1);
	ybinwidth = fqqs_plus_4jet->GetYaxis()->GetBinWidth(1);
	zbinwidth = fqqs_plus_4jet->GetZaxis()->GetBinWidth(1);
	//open background distributions
	fbck_plus_4jet  = (TH3D*)f->Get("f_bck_plus_4jet");
	fbck_minus_4jet = (TH3D*)f->Get("f_bck_minus_4jet");
	fbck_plus_5jet  = (TH3D*)f->Get("f_bck_plus_5jet");
	fbck_minus_5jet = (TH3D*)f->Get("f_bck_minus_5jet");
	WJets_plus_4jet  = (TH3D*)f->Get("f_WJets_plus_4jet");
	WJets_minus_4jet = (TH3D*)f->Get("f_WJets_minus_4jet");
	WJets_plus_5jet  = (TH3D*)f->Get("f_WJets_plus_5jet");
	WJets_minus_5jet = (TH3D*)f->Get("f_WJets_minus_5jet");
	//Normalize the distributions
	double qq_rescale_4jet = 1.0/(fqqs_plus_4jet->Integral()+fqqs_minus_4jet->Integral());
	// double qq_rescale_4jet = 0.0/(fqqs_plus_4jet->Integral()+fqqs_minus_4jet->Integral());
	fqqs_plus_4jet->Scale(qq_rescale_4jet);
	fqqs_xi_plus_4jet->Scale(qq_rescale_4jet);
	fqqs_delta_plus_4jet->Scale(qq_rescale_4jet);
	fqqs_minus_4jet->Scale(qq_rescale_4jet);
	fqqs_xi_minus_4jet->Scale(qq_rescale_4jet);
	fqqs_delta_minus_4jet->Scale(qq_rescale_4jet);
	fqqa_plus_4jet->Scale(qq_rescale_4jet);
	fqqa_xi_plus_4jet->Scale(qq_rescale_4jet);
	fqqa_delta_plus_4jet->Scale(qq_rescale_4jet);
	fqqa_minus_4jet->Scale(qq_rescale_4jet);
	fqqa_xi_minus_4jet->Scale(qq_rescale_4jet);
	fqqa_delta_minus_4jet->Scale(qq_rescale_4jet);
	double gg_rescale_4jet = 1.0/(fgg_plus_4jet->Integral()+fgg_minus_4jet->Integral());
	// double gg_rescale_4jet = 0.0/(fgg_plus_4jet->Integral()+fgg_minus_4jet->Integral());
	fgg_plus_4jet->Scale(gg_rescale_4jet);
	fgg_minus_4jet->Scale(gg_rescale_4jet);
	double bck_rescale_4jet = 1.0/(fbck_plus_4jet->Integral()+fbck_minus_4jet->Integral());
	// double bck_rescale_4jet = 0.0/(fbck_plus_4jet->Integral()+fbck_minus_4jet->Integral());
	fbck_plus_4jet->Scale(bck_rescale_4jet);
	fbck_minus_4jet->Scale(bck_rescale_4jet);
	double WJets_rescale_4jet = 1.0/(WJets_plus_4jet->Integral()+WJets_minus_4jet->Integral());
	// double WJets_rescale_4jet = 0.0/(WJets_plus_4jet->Integral()+WJets_minus_4jet->Integral());
	WJets_plus_4jet->Scale(WJets_rescale_4jet);
	WJets_minus_4jet->Scale(WJets_rescale_4jet);
	double qq_rescale_5jet = 1.0/(fqqs_plus_5jet->Integral()+fqqs_minus_5jet->Integral());
	// double qq_rescale_5jet = 0.0/(fqqs_plus_5jet->Integral()+fqqs_minus_5jet->Integral());
	fqqs_plus_5jet->Scale(qq_rescale_5jet);
	fqqs_xi_plus_5jet->Scale(qq_rescale_5jet);
	fqqs_delta_plus_5jet->Scale(qq_rescale_5jet);
	fqqs_minus_5jet->Scale(qq_rescale_5jet);
	fqqs_xi_minus_5jet->Scale(qq_rescale_5jet);
	fqqs_delta_minus_5jet->Scale(qq_rescale_5jet);
	fqqa_plus_5jet->Scale(qq_rescale_5jet);
	fqqa_xi_plus_5jet->Scale(qq_rescale_5jet);
	fqqa_delta_plus_5jet->Scale(qq_rescale_5jet);
	fqqa_minus_5jet->Scale(qq_rescale_5jet);
	fqqa_xi_minus_5jet->Scale(qq_rescale_5jet);
	fqqa_delta_minus_5jet->Scale(qq_rescale_5jet);
	double gg_rescale_5jet = 1.0/(fgg_plus_5jet->Integral()+fgg_minus_5jet->Integral());
	// double gg_rescale_5jet = 0.0/(fgg_plus_5jet->Integral()+fgg_minus_5jet->Integral());
	fgg_plus_5jet->Scale(gg_rescale_5jet);
	fgg_minus_5jet->Scale(gg_rescale_5jet);
	double bck_rescale_5jet = 1.0/(fbck_plus_5jet->Integral()+fbck_minus_5jet->Integral());
	// double bck_rescale_5jet = 0.0/(fbck_plus_5jet->Integral()+fbck_minus_5jet->Integral());
	fbck_plus_5jet->Scale(bck_rescale_5jet);
	fbck_minus_5jet->Scale(bck_rescale_5jet);
	double WJets_rescale_5jet = 1.0/(WJets_plus_5jet->Integral()+WJets_minus_5jet->Integral());
	// double WJets_rescale_5jet = 0.0/(WJets_plus_5jet->Integral()+WJets_minus_5jet->Integral());
	WJets_plus_5jet->Scale(WJets_rescale_5jet);
	WJets_minus_5jet->Scale(WJets_rescale_5jet);
	gtt_4jet->Scale(1.0/gtt_4jet->Integral());
	gbk_4jet->Scale(1.0/gbk_4jet->Integral());
	gtt_5jet->Scale(1.0/gtt_5jet->Integral());
	gbk_5jet->Scale(1.0/gbk_5jet->Integral());
	F_xi_4jet = fqqs_xi_plus_4jet->Integral()+fqqs_xi_minus_4jet->Integral();
	F_delta_4jet = fqqs_delta_plus_4jet->Integral()+fqqs_delta_minus_4jet->Integral();
	F_xi_5jet = fqqs_xi_plus_5jet->Integral()+fqqs_xi_minus_5jet->Integral();
	F_delta_5jet = fqqs_delta_plus_5jet->Integral()+fqqs_delta_minus_5jet->Integral();
}
