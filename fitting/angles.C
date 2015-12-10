#define angles_cxx
#include "angles.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>

///////////////////////////////////////////////////////////
// Given a root file list of ttbar_mass, cos_theta_cs,   //
// Feynman_x, and w_a, generates and saves a 3D histogram//
// of distribution function for events. If the given     //
// sample is qqbar and symmetric, an option exists to    //
// also generate the weighted asymmetric histogram       //
///////////////////////////////////////////////////////////

// If MG sample is used , the event weight should be modified using the cross section of mg signal sample and events of mg sample

/////////////////////////////////////////////////////////////////////////////////////////////
//Variable declarations that can be altered
bool genAsymmetricHisto = false; //Set to true if you want 
								//the asymmetric histogram
								//generated from given w_a
//limits on distributions of parameters
double min_ttbar_mass = 300.0;
double max_ttbar_mass = 1750.0;
double min_Feynman_x = 0.0;
double max_Feynman_x = 0.6;
const double min_cos_theta_cs = -1.0;
const double max_cos_theta_cs = 1.0;
const double min_lnL = -50.0;
const double max_lnL = 100.0;
//bin numbers for histograms
int nBins_cos_theta_cs = 20;
int nBins_Feynman_x  = 30;
int nBins_ttbar_mass = 40;
int nBins_lnL = 150;
//PDF member to be used
int pdf_index = -1;
// Number of tagged bs required
int b_tag_cuts = 2;
// kinfit chi2 cut
float kinfit_chi2_cut = 100;
// total number of templates
int num_of_temp = 5;
// template number
int temp_num = 5;
// signal events number scale.
float signal_scale = 1;

/////////////////////////////////////////////////////////////////////////////////////////////


//names and titles for the histograms
const char* sampleName = "testing";
const char* sampleName_formatted = "yeah testing";
const char* gtt_4jet_name = "gtt_4jet";
const char* gtt_4jet_title = "t#bar{t} likelihood distribution, 4jet";
const char* gbk_4jet_name = "gbk_4jet";
const char* gbk_4jet_title = "background likelihood distribution, 4jet";
const char* gtt_5jet_name = "gtt_5jet";
const char* gtt_5jet_title = "t#bar{t} likelihood distribution, 5jet";
const char* gbk_5jet_name = "gbk_5jet";
const char* gbk_5jet_title = "background likelihood distribution, 5jet";
//input filename
const char* infile = "angles.root";
//output filename
const char* outfile = "qqbar_s_and_a.root";
//treename to get (should be default except in the case of the data files)
const char* tree_name = "angles";
void angles::Loop(int is_for_dist) {
	Loop(is_for_dist,1,1.0);
}


//Working function containing main loop
void angles::Loop(int is_for_dist, int neventsgenerated, double crosssection)
{
	// signal_scale = 1/float(num_of_temp);
	// correct for the signal number of events generated. 
	if (is_for_dist == 1 || is_for_dist == 2)
	{
		neventsgenerated = neventsgenerated*signal_scale;	// the number of events need to be corrected only for signals
		cout << " Using "<<signal_scale<<" of signal templates:"<<endl;
	}
	//Print out which template is running for MC systematics
	// cout <<"Using template #"<<temp_num<<endl; 

	bool addTwice = (is_for_dist == 1 || is_for_dist == 2);
	//variables from input file
	float cs,xf,PT,mttbar,ln_L;
	float w_a,w_s_xi,w_a_xi,w_s_delta,w_a_delta;
	float w_a_opp,w_s_xi_opp,w_a_xi_opp,w_s_delta_opp,w_a_delta_opp;
	float btag_eff_reweight, tracking_reweight, lepID_reweight, lepIso_reweight, trigger_reweight, pileup_reweight, top_pT_reweight, GJR_reweight, CT10_reweight, cteq_reweight;
	float btag_eff_reweight_low, tracking_reweight_low, lepID_reweight_low, lepIso_reweight_low, trigger_reweight_low;
	float btag_eff_reweight_hi, tracking_reweight_hi, lepID_reweight_hi, lepIso_reweight_hi, trigger_reweight_hi;
	float pileup, pileup_real;
	float fitParValues[6];
	int lepton_charge,n_bTags,n_valid_jets;
	Int_t motherPIDs[2];
	//load input file
	TFile * in_f=new TFile(infile);
	TTree * tree=(TTree*)in_f->Get(tree_name);
	tree->SetBranchAddress("cos_theta_cs",&cs);
	tree->SetBranchAddress("Feynman_x",&xf);
	tree->SetBranchAddress("Qt",&PT);
	tree->SetBranchAddress("ttbar_mass",&mttbar);
	tree->SetBranchAddress("Q_l",&lepton_charge);
	tree->SetBranchAddress("n_valid_jets",&n_valid_jets);
	tree->SetBranchAddress("n_bTags",&n_bTags); 
	tree->SetBranchAddress("fitParValues",   &fitParValues);
	tree->SetBranchAddress("lnL",&ln_L);
	tree->SetBranchAddress("w_a",&w_a);
	tree->SetBranchAddress("w_a_opp",&w_a_opp);
	tree->SetBranchAddress("w_s_xi",&w_s_xi);
	tree->SetBranchAddress("w_s_xi_opp",&w_s_xi_opp);
	tree->SetBranchAddress("w_a_xi",&w_a_xi);
	tree->SetBranchAddress("w_a_xi_opp",&w_a_xi_opp);
	tree->SetBranchAddress("w_s_delta",&w_s_delta);
	tree->SetBranchAddress("w_s_delta_opp",&w_s_delta_opp);
	tree->SetBranchAddress("w_a_delta",&w_a_delta);
	tree->SetBranchAddress("w_a_delta_opp",&w_a_delta_opp);
	tree->SetBranchAddress("top_pT_reweight",&top_pT_reweight);
	tree->SetBranchAddress("motherPIDs",motherPIDs);
	tree->SetBranchAddress("pileup_reweight", &pileup_reweight);
	tree->SetBranchAddress("top_pT_reweight", &top_pT_reweight);
	tree->SetBranchAddress("GJR_reweight",    &GJR_reweight);
	tree->SetBranchAddress("CT10_reweight",   &CT10_reweight);
	tree->SetBranchAddress("cteq_reweight",   &cteq_reweight);
	tree->SetBranchAddress("btag_eff_reweight",     &btag_eff_reweight);
	tree->SetBranchAddress("tracking_reweight",     &tracking_reweight);
	tree->SetBranchAddress("lepID_reweight",        &lepID_reweight);
	tree->SetBranchAddress("lepIso_reweight",       &lepIso_reweight);
	tree->SetBranchAddress("trigger_reweight",      &trigger_reweight);
	tree->SetBranchAddress("btag_eff_reweight_low", &btag_eff_reweight_low);
	tree->SetBranchAddress("tracking_reweight_low", &tracking_reweight_low);
	tree->SetBranchAddress("lepID_reweight_low",    &lepID_reweight_low);
	tree->SetBranchAddress("lepIso_reweight_low",   &lepIso_reweight_low);
	tree->SetBranchAddress("trigger_reweight_low",  &trigger_reweight_low);
	tree->SetBranchAddress("btag_eff_reweight_hi",  &btag_eff_reweight_hi);
	tree->SetBranchAddress("tracking_reweight_hi",  &tracking_reweight_hi);
	tree->SetBranchAddress("lepID_reweight_hi",     &lepID_reweight_hi);
	tree->SetBranchAddress("lepIso_reweight_hi",    &lepIso_reweight_hi);
	tree->SetBranchAddress("trigger_reweight_hi",   &trigger_reweight_hi);	
	tree->SetBranchAddress("pileup",   &pileup);
	tree->SetBranchAddress("pileup_real",   &pileup_real);	
	
	//Set PDF weights address for signal templates
	int apply_PDF = 0 ; // By default no PDF weight will be applied, unless the template has a PDF weight branch
	double pdf_weights[100];
	float PDF_reweight = 1;
	if (pdf_index >= 0)
	{
		if(tree->GetListOfBranches()->FindObject("Pdf_weights")) // This will check if the branch "Pdf_weights" exists for this template root
		{
			tree->SetBranchAddress("Pdf_weights",&pdf_weights);
    		cout << "PDF member index is: " << pdf_index <<endl;
    		apply_PDF = 1;
    	}
	}
	//Set up histogram names and titles	
	char** names = (char**)malloc(25*sizeof(char*));
	char** titles = (char**)malloc(25*sizeof(char*));
	for (int i=0; i<25; i++) {
		names[i] = (char*)malloc(50*sizeof(char));
		titles[i] = (char*)malloc(250*sizeof(char));
		strcpy(names[i],sampleName); strcpy(titles[i],sampleName_formatted); 
	}
	strcat(names[0],"_s_p_4j"); 	   strcat(titles[0]," Unweighted Symmetric Distribution, + leptons, 4 jets");
	strcat(names[1],"_s_xi_p_4j"); 	   strcat(titles[1]," xi-weighted Symmetric Distribution, + leptons, 4 jets");
	strcat(names[2],"_s_delta_p_4j");  strcat(titles[2]," delta-weighted Symmetric Distribution, + leptons, 4 jets");
	strcat(names[3],"_a_p_4j"); 	   strcat(titles[3]," Weighted Asymmetric Distribution, + leptons, 4 jets");
	strcat(names[4],"_a_xi_p_4j"); 	   strcat(titles[4]," xi-weighted Asymmetric Distribution, + leptons, 4 jets");
	strcat(names[5],"_a_delta_p_4j");  strcat(titles[5]," delta-weighted Asymmetric Distribution, + leptons, 4 jets");
	strcat(names[6],"_s_m_4j"); 	   strcat(titles[6]," Unweighted Symmetric Distribution, - leptons, 4 jets");
	strcat(names[7],"_s_xi_m_4j"); 	   strcat(titles[7]," xi-weighted Symmetric Distribution, - leptons, 4 jets");
	strcat(names[8],"_s_delta_m_4j");  strcat(titles[8]," delta-weighted Symmetric Distribution, - leptons, 4 jets");
	strcat(names[9],"_a_m_4j"); 	   strcat(titles[9]," Weighted Asymmetric Distribution, - leptons, 4 jets");
	strcat(names[10],"_a_xi_m_4j");    strcat(titles[10]," xi-weighted Asymmetric Distribution, - leptons, 4 jets");
	strcat(names[11],"_a_delta_m_4j"); strcat(titles[11]," delta-weighted Asymmetric Distribution, - leptons, 4 jets");
	strcat(names[12],"_s_p_5j"); 	   strcat(titles[12]," Unweighted Symmetric Distribution, + leptons, 5 jets");
	strcat(names[13],"_s_xi_p_5j");    strcat(titles[13]," xi-weighted Symmetric Distribution, + leptons, 5 jets");
	strcat(names[14],"_s_delta_p_5j"); strcat(titles[14]," delta-weighted Symmetric Distribution, + leptons, 5 jets");
	strcat(names[15],"_a_p_5j");	   strcat(titles[15]," Weighted Asymmetric Distribution, + leptons, 5 jets");
	strcat(names[16],"_a_xi_p_5j");    strcat(titles[16]," xi-weighted Asymmetric Distribution, + leptons, 5 jets");
	strcat(names[17],"_a_delta_p_5j"); strcat(titles[17]," delta-weighted Asymmetric Distribution, + leptons, 5 jets");
	strcat(names[18],"_s_m_5j"); 	   strcat(titles[18]," Unweighted Symmetric Distribution, - leptons, 5 jets");
	strcat(names[19],"_s_xi_m_5j");    strcat(titles[19]," xi-weighted Symmetric Distribution, - leptons, 5 jets");
	strcat(names[20],"_s_delta_m_5j"); strcat(titles[20]," delta-weighted Symmetric Distribution, - leptons, 5 jets");
	strcat(names[21],"_a_m_5j"); 	   strcat(titles[21]," Weighted Asymmetric Distribution, - leptons, 5 jets");
	strcat(names[22],"_a_xi_m_5j");    strcat(titles[22]," xi-weighted Asymmetric Distribution, - leptons, 5 jets");
	strcat(names[23],"_a_delta_m_5j"); strcat(titles[23]," delta-weighted Asymmetric Distribution, - leptons, 5 jets");
	strcat(names[24],"_everything_no_pT_reweight"); strcat(titles[24]," everything distribution, no top pT reweighting");
	//Set up histograms themselves
	TH3D** histos = (TH3D**)malloc(25*sizeof(TH3D*));
	for (int i=0; i<25; i++) {
		histos[i] = new TH3D(names[i],titles[i],nBins_cos_theta_cs,min_cos_theta_cs,max_cos_theta_cs,
								nBins_Feynman_x,min_Feynman_x,max_Feynman_x,
								nBins_ttbar_mass,min_ttbar_mass,max_ttbar_mass);
		histos[i]->GetXaxis()->SetTitle("cos(#theta^{*})");
		histos[i]->GetYaxis()->SetTitle("Feynman x (x_{F})");
		histos[i]->GetZaxis()->SetTitle("M_{t #bar{t}} (GeV)");
	}
	TH1D *gtt_4jet = new TH1D(gtt_4jet_name,gtt_4jet_title,nBins_lnL,min_lnL,max_lnL);
	gtt_4jet->GetXaxis()->SetTitle("-2ln(L)");
	TH1D *gbk_4jet = new TH1D(gbk_4jet_name,gbk_4jet_title,nBins_lnL,min_lnL,max_lnL);
	gbk_4jet->GetXaxis()->SetTitle("-2ln(L)");
	TH1D *gtt_5jet = new TH1D(gtt_5jet_name,gtt_5jet_title,nBins_lnL,min_lnL,max_lnL);
	gtt_5jet->GetXaxis()->SetTitle("-2ln(L)");
	TH1D *gbk_5jet = new TH1D(gbk_5jet_name,gbk_5jet_title,nBins_lnL,min_lnL,max_lnL);
	gbk_5jet->GetXaxis()->SetTitle("-2ln(L)");

	//Set up TTree structure
	float costheta, xF, tt_M;
	float weight, wa, waxi, wadelta, wsxi, wsdelta;
	int nJets, Ql, dist_type;
	TFile * out_f=new TFile(outfile,"Recreate");
	TTree * outputTree=new TTree("output_tree","Recreate");
	outputTree->Branch("costheta",&costheta);
	outputTree->Branch("xF",&xF);
	outputTree->Branch("tt_M",&tt_M);
	outputTree->Branch("nJets",&nJets);
	outputTree->Branch("Ql",&Ql);
	outputTree->Branch("dist_type",&dist_type);
	outputTree->Branch("weight",&weight);
	outputTree->Branch("wa",&wa);
	outputTree->Branch("waxi",&waxi);
	outputTree->Branch("wadelta",&wadelta);
	outputTree->Branch("wsxi",&wsxi);
	outputTree->Branch("wsdelta",&wsdelta);

	//Set up loop
	Long64_t   nentries=tree->GetEntriesFast();

	//MAIN LOOP START
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		tree->GetEntry(jentry);

		// // for MC systematics
		// if (jentry%num_of_temp != temp_num && (is_for_dist == 1 || is_for_dist ==2) ) // create smaller templates only for signal which has much larger sample size
		// 	continue;

		// Test how size of signal templates change the fit results
		if (jentry > signal_scale*nentries && (is_for_dist == 1 || is_for_dist ==2) )
			continue;
		
		if (xf<0.0) {
			xf = xf*-1.0; //Use the absolute value of Feynman X
		}

		// Set PDF weights if the template root file has this branch
		if (apply_PDF == 1) {
			// Only use the PDF member if the weight is non-zero. The empty weight is set to be zero as default. 
			float tmp = 0 ;
	    	tmp = pdf_weights[pdf_index];
			// This is to avoid using the PDF member index out of the bound.
			if (tmp>0)
			{
				PDF_reweight = tmp ;
			}
			else {
				cout<<"The PDF index is out of bound for the PDF for this sample!"<<endl;
				cout<<"PDF weight = "<< tmp <<endl;
				cout<<"Will quit this run!"<<endl;
				break;
			}
		}		
		
		//MAKE CUTS FOR OUR ANALYSIS
		//cut that the event is within range of the plots
		if (xf<min_Feynman_x || xf>max_Feynman_x || cs<min_cos_theta_cs || cs>max_cos_theta_cs || mttbar<min_ttbar_mass || mttbar>max_ttbar_mass)
			continue;

		//cut for two bTags
		if (n_bTags<b_tag_cuts)
			continue;		

		//cut on chi^2
		if (ln_L > kinfit_chi2_cut)
			continue;

		//set the event weight
		// Comment: Here, the detector modeling sys can be evaluated
		float event_weight = pileup_reweight*btag_eff_reweight*tracking_reweight*lepID_reweight*lepIso_reweight*trigger_reweight*top_pT_reweight*PDF_reweight;
		float event_weight_no_top_pT_reweight = pileup_reweight*btag_eff_reweight*tracking_reweight*lepID_reweight*lepIso_reweight*trigger_reweight*PDF_reweight;

		//only add the q glu events once
		if (is_for_dist == 2 && (motherPIDs[0]!=21 || motherPIDs[1]!=21))
			addTwice = false;

		//half the event weight if we want to add it twice
		if (addTwice) {
			event_weight = 0.5*event_weight;
			event_weight_no_top_pT_reweight = 0.5*event_weight_no_top_pT_reweight;
		}

		//rescale the event weight by the cross section over the number of events
		
		// //For Powheg sample
		event_weight = event_weight * (crosssection/neventsgenerated) * (25523595./245.8)*signal_scale;
		event_weight_no_top_pT_reweight = event_weight_no_top_pT_reweight * (crosssection/neventsgenerated) * (25523595./245.8)*signal_scale;

		// For MG sample
		//event_weight = event_weight * (crosssection/neventsgenerated) * (25273288./108.0)*signal_scale;
		//event_weight_no_top_pT_reweight = event_weight_no_top_pT_reweight * (crosssection/neventsgenerated) * (25273288./108.0)*signal_scale;

		//add event to the symmetric histograms based on lepton charge
		if (lepton_charge > 0) {
			if (n_valid_jets==4) {
				histos[0]->Fill(cs,xf,mttbar,event_weight);
				histos[1]->Fill(cs,xf,mttbar,w_s_xi*event_weight);
				histos[2]->Fill(cs,xf,mttbar,w_s_delta*event_weight);
				histos[24]->Fill(cs,xf,mttbar,event_weight_no_top_pT_reweight);
				if (is_for_dist == 3 || is_for_dist == 4)
					gbk_4jet->Fill(ln_L,event_weight);
				else if (is_for_dist == 1 || is_for_dist == 2)
					gtt_4jet->Fill(ln_L,event_weight);
				if (addTwice) {
					histos[6]->Fill(-1.0*cs,xf,mttbar,event_weight);
					histos[7]->Fill(-1.0*cs,xf,mttbar,w_s_xi_opp*event_weight);
					histos[8]->Fill(-1.0*cs,xf,mttbar,w_s_delta_opp*event_weight);
					histos[24]->Fill(-1.0*cs,xf,mttbar,event_weight_no_top_pT_reweight);
					if (is_for_dist == 3 || is_for_dist == 4)
						gbk_4jet->Fill(ln_L,event_weight);
					else if (is_for_dist == 1 || is_for_dist == 2)
						gtt_4jet->Fill(ln_L,event_weight);
				}
			}
			else if (n_valid_jets==5) {
				histos[12]->Fill(cs,xf,mttbar,event_weight);
				histos[13]->Fill(cs,xf,mttbar,w_s_xi*event_weight);
				histos[14]->Fill(cs,xf,mttbar,w_s_delta*event_weight);
				histos[24]->Fill(cs,xf,mttbar,event_weight_no_top_pT_reweight);
				if (is_for_dist == 3 || is_for_dist == 4)
					gbk_5jet->Fill(ln_L,event_weight);
				else if (is_for_dist == 1 || is_for_dist == 2)
					gtt_5jet->Fill(ln_L,event_weight);
				if (addTwice) {
					histos[18]->Fill(-1.0*cs,xf,mttbar,event_weight);
					histos[19]->Fill(-1.0*cs,xf,mttbar,w_s_xi_opp*event_weight);
					histos[20]->Fill(-1.0*cs,xf,mttbar,w_s_delta_opp*event_weight);
					histos[24]->Fill(-1.0*cs,xf,mttbar,event_weight_no_top_pT_reweight);
					if (is_for_dist == 3 || is_for_dist == 4)
						gbk_5jet->Fill(ln_L,event_weight);
					else if (is_for_dist == 1 || is_for_dist == 2)
						gtt_5jet->Fill(ln_L,event_weight);
				}	
			}
		}
		else if (lepton_charge < 0) {
			if (n_valid_jets==4) {
				histos[6]->Fill(cs,xf,mttbar,event_weight);
				histos[7]->Fill(cs,xf,mttbar,w_s_xi*event_weight);
				histos[8]->Fill(cs,xf,mttbar,w_s_delta*event_weight);
				histos[24]->Fill(cs,xf,mttbar,event_weight_no_top_pT_reweight);
				if (is_for_dist == 3 || is_for_dist == 4)
					gbk_4jet->Fill(ln_L,event_weight);
				else if (is_for_dist == 1 || is_for_dist == 2)
					gtt_4jet->Fill(ln_L,event_weight);
				if (addTwice) {
					histos[0]->Fill(-1.0*cs,xf,mttbar,event_weight);
					histos[1]->Fill(-1.0*cs,xf,mttbar,w_s_xi_opp*event_weight);
					histos[2]->Fill(-1.0*cs,xf,mttbar,w_s_delta_opp*event_weight);
					histos[24]->Fill(-1.0*cs,xf,mttbar,event_weight_no_top_pT_reweight);
					if (is_for_dist == 3 || is_for_dist == 4)
						gbk_4jet->Fill(ln_L,event_weight);
					else if (is_for_dist == 1 || is_for_dist == 2)
						gtt_4jet->Fill(ln_L,event_weight);
				}
			}
			else if (n_valid_jets==5) {
				histos[18]->Fill(cs,xf,mttbar,event_weight);
				histos[19]->Fill(cs,xf,mttbar,w_s_xi*event_weight);
				histos[20]->Fill(cs,xf,mttbar,w_s_delta*event_weight);
				histos[24]->Fill(cs,xf,mttbar,event_weight_no_top_pT_reweight);
				if (is_for_dist == 3 || is_for_dist == 4)
					gbk_5jet->Fill(ln_L,event_weight);
				else if (is_for_dist == 1 || is_for_dist == 2)
					gtt_5jet->Fill(ln_L,event_weight);
				if (addTwice) {
					histos[12]->Fill(-1.0*cs,xf,mttbar,event_weight);
					histos[13]->Fill(-1.0*cs,xf,mttbar,w_s_xi_opp*event_weight);
					histos[14]->Fill(-1.0*cs,xf,mttbar,w_s_delta_opp*event_weight);
					histos[24]->Fill(-1.0*cs,xf,mttbar,event_weight_no_top_pT_reweight);
					if (is_for_dist == 3 || is_for_dist == 4)
						gbk_5jet->Fill(ln_L,event_weight);
					else if (is_for_dist == 1 || is_for_dist == 2)
						gtt_5jet->Fill(ln_L,event_weight);
				}
			}
		}
		//add events to the asymmetric histogram
		if (genAsymmetricHisto) {
			if (lepton_charge > 0) {
				if (n_valid_jets==4) {
					histos[3]->Fill(cs,xf,mttbar,w_a*event_weight);
					histos[4]->Fill(cs,xf,mttbar,w_a_xi*event_weight);
					histos[5]->Fill(cs,xf,mttbar,w_a_delta*event_weight);
					if (addTwice) {
						histos[9]->Fill(-1.0*cs,xf,mttbar,w_a_opp*event_weight);
						histos[10]->Fill(-1.0*cs,xf,mttbar,w_a_xi_opp*event_weight);
						histos[11]->Fill(-1.0*cs,xf,mttbar,w_a_delta_opp*event_weight);
					}
				}
				else if (n_valid_jets==5) {
					histos[15]->Fill(cs,xf,mttbar,w_a*event_weight);
					histos[16]->Fill(cs,xf,mttbar,w_a_xi*event_weight);
					histos[17]->Fill(cs,xf,mttbar,w_a_delta*event_weight);
					if (addTwice) {
						histos[21]->Fill(-1.0*cs,xf,mttbar,w_a_opp*event_weight);
						histos[22]->Fill(-1.0*cs,xf,mttbar,w_a_xi_opp*event_weight);
						histos[23]->Fill(-1.0*cs,xf,mttbar,w_a_delta_opp*event_weight);
					}
				}
			}
			else if (lepton_charge < 0) {
				if (n_valid_jets==4) {
					histos[9]->Fill(cs,xf,mttbar,w_a*event_weight);
					histos[10]->Fill(cs,xf,mttbar,w_a_xi*event_weight);
					histos[11]->Fill(cs,xf,mttbar,w_a_delta*event_weight);
					if (addTwice) {
						histos[3]->Fill(-1.0*cs,xf,mttbar,w_a_opp*event_weight);
						histos[4]->Fill(-1.0*cs,xf,mttbar,w_a_xi_opp*event_weight);
						histos[5]->Fill(-1.0*cs,xf,mttbar,w_a_delta_opp*event_weight);
					}
				}
				else if (n_valid_jets==5) {
					histos[21]->Fill(cs,xf,mttbar,w_a*event_weight);
					histos[22]->Fill(cs,xf,mttbar,w_a_xi*event_weight);
					histos[23]->Fill(cs,xf,mttbar,w_a_delta*event_weight);
					if (addTwice) {
						histos[15]->Fill(-1.0*cs,xf,mttbar,w_a_opp*event_weight);
						histos[16]->Fill(-1.0*cs,xf,mttbar,w_a_xi_opp*event_weight);
						histos[17]->Fill(-1.0*cs,xf,mttbar,w_a_delta_opp*event_weight);
					}
				}
			}
		}
		//fill the TTree
		costheta = cs;
		xF = xf;
		tt_M = mttbar;
		nJets = n_valid_jets;
		Ql = lepton_charge;
		dist_type = is_for_dist;
		weight = event_weight;
		wa = w_a;
		waxi = w_a_xi;
		wadelta = w_a_delta;
		wsxi = w_s_xi;
		wsdelta = w_s_delta;
		outputTree->Fill();
		if (addTwice) {
			costheta = -1.0*cs;
			xF = xf;
			tt_M = mttbar;
			nJets = n_valid_jets;
			Ql = -1.0*lepton_charge;
			dist_type = is_for_dist;
			weight = event_weight;
			wa = w_a_opp;
			waxi = w_a_xi_opp;
			wadelta = w_a_delta_opp;
			wsxi = w_s_xi_opp;
			wsdelta = w_s_delta_opp;
			outputTree->Fill();
		}
	}
	
	//Save completed histogram(s) to outputfile
	for (int i=0; i<25; i++) {
		histos[i]->Write();
	}
	gtt_4jet->Write();
	gbk_4jet->Write();
	gtt_5jet->Write();
	gbk_5jet->Write();
	//Save TTree to output file
	outputTree->Write();
	printf("		Done\n");
		
	//close the output file
	out_f->Close();
}

void angles::InitializeNames(const char* in_name, const char* out_name, const char* samplename, 
								const char* samplename_formatted, bool genasymhisto) {
	infile = in_name;
	outfile = out_name;
	sampleName = samplename;
	sampleName_formatted = samplename_formatted;
	genAsymmetricHisto = genasymhisto;
}

void angles::SetInputTreeName(const char* name) {
	tree_name = name;
}

void angles::SetBins(int xBins, int yBins, int zBins) {
	nBins_cos_theta_cs = xBins;
	nBins_Feynman_x = yBins;
	nBins_ttbar_mass = zBins;
}

void angles::SetMassLimits(double lowmass, double highmass) {
	min_ttbar_mass=lowmass;
	max_ttbar_mass=highmass;
}

void angles::SetxFLimits(double xfLow, double xfHigh) {
	min_Feynman_x = xfLow;
	max_Feynman_x = xfHigh;
}

void angles::SetPdfMembers(int pdf_member){
	pdf_index = pdf_member;
}
