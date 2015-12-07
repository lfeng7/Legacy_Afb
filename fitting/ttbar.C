#define ttbar_cxx
#include "ttbar.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"
#include <TLorentzVector.h>
#include <TLorentzRotation.h>

//Lorentz Declartions

const double sqrt_s=8000;
const double beam_energy=sqrt_s/2;

TLorentzRotation S;

//variable for ouput file name
const char* output_filename = "angles.root";
//variable for output tree name (should be default except for total data file)
const char* output_treename = "angles";

void ttbar::Loop(int data_or_mc)
{
	//Histograms for tests

	int kept=0;
	int count=0;
	//declarations of output variables
	float ttbar_mass,cos_theta_cs,x_f,Qt,cos_theta_mc,x_f_mc,ttbar_mass_mc,ln_L;
	float w_a,w_s_xi,w_a_xi,w_s_delta,w_a_delta;
	float w_a_opp,w_s_xi_opp,w_a_xi_opp,w_s_delta_opp,w_a_delta_opp;
	float btag_eff_reweight, tracking_reweight, lepID_reweight, lepIso_reweight, trigger_reweight, pileup_reweight, top_pT_reweight, GJR_reweight, CT10_reweight, cteq_reweight;
	float btag_eff_reweight_low, tracking_reweight_low, lepID_reweight_low, lepIso_reweight_low, trigger_reweight_low;
	float btag_eff_reweight_hi, tracking_reweight_hi, lepID_reweight_hi, lepIso_reweight_hi, trigger_reweight_hi;
	float pileup, pileup_real;
	float fitParValues[6];
	double Pdf_w[100];
	int Q_l, n_valid_jets, n_bTags;
	Int_t motherPIDs[2];
	//open output file and create branches
	TFile * file=new TFile(output_filename,"Recreate");
	TTree * output=new TTree(output_treename,"Recreate");
	output->Branch("ttbar_mass",&ttbar_mass);
	output->Branch("Qt",&Qt);
	output->Branch("cos_theta_cs",&cos_theta_cs);
	output->Branch("Feynman_x",&x_f);
	output->Branch("Q_l",&Q_l);
	output->Branch("cos_theta_mc",&cos_theta_mc);
	output->Branch("Feynman_x_mc",&x_f_mc);
	output->Branch("ttbar_mass_mc",&ttbar_mass_mc);
	output->Branch("lnL",&ln_L);
	output->Branch("n_valid_jets",&n_valid_jets);
	output->Branch("n_bTags",&n_bTags);
	output->Branch("fitParValues",fitParValues,"fitParValues[6]/F");
	output->Branch("w_a",          &w_a);
	output->Branch("w_a_opp",      &w_a_opp);
	output->Branch("w_s_xi",       &w_s_xi);
	output->Branch("w_s_xi_opp",   &w_s_xi_opp);
	output->Branch("w_a_xi",       &w_a_xi);
	output->Branch("w_a_xi_opp",   &w_a_xi_opp);
	output->Branch("w_s_delta",    &w_s_delta);
	output->Branch("w_s_delta_opp",&w_s_delta_opp);
	output->Branch("w_a_delta",    &w_a_delta);
	output->Branch("w_a_delta_opp",&w_a_delta_opp);
	output->Branch("motherPIDs",motherPIDs,"motherPIDs[2]/I");
	output->Branch("Pdf_weights",&Pdf_w,"Pdf_w[100]/D");
	output->Branch("pileup_reweight", &pileup_reweight);
	output->Branch("top_pT_reweight", &top_pT_reweight);
	output->Branch("GJR_reweight",    &GJR_reweight);
	output->Branch("CT10_reweight",   &CT10_reweight);
	output->Branch("cteq_reweight",   &cteq_reweight);
	output->Branch("btag_eff_reweight",     &btag_eff_reweight);
	output->Branch("tracking_reweight",     &tracking_reweight);
	output->Branch("lepID_reweight",        &lepID_reweight);
	output->Branch("lepIso_reweight",       &lepIso_reweight);
	output->Branch("trigger_reweight",      &trigger_reweight);
	output->Branch("btag_eff_reweight_low", &btag_eff_reweight_low);
	output->Branch("tracking_reweight_low", &tracking_reweight_low);
	output->Branch("lepID_reweight_low",    &lepID_reweight_low);
	output->Branch("lepIso_reweight_low",   &lepIso_reweight_low);
	output->Branch("trigger_reweight_low",  &trigger_reweight_low);
	output->Branch("btag_eff_reweight_hi",  &btag_eff_reweight_hi);
	output->Branch("tracking_reweight_hi",  &tracking_reweight_hi);
	output->Branch("lepID_reweight_hi",     &lepID_reweight_hi);
	output->Branch("lepIso_reweight_hi",    &lepIso_reweight_hi);
	output->Branch("trigger_reweight_hi",   &trigger_reweight_hi);	
	output->Branch("pileup",   &pileup);
	output->Branch("pileup_real",   &pileup_real);	
	
	//Setup loop
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	//MAIN LOOP
	for (Long64_t jentry=0; jentry< nentries ;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		count+=1;
		//intialize the rotations to the identity
		TLorentzRotation R_data;
		TLorentzRotation R_mc;
		
		//Supercool Load bar
		if (nentries > 100) {
			if ( count % (nentries/100) == 0 ){
				float ratio = count/(float)nentries;
				int   c     = ratio * 50;
				printf("%3d%% [", (int)(ratio*100) );
				for (int x=0; x<c; x++) {
					if (x+1>=c)
						printf(">");
					else
						printf("=");
				}
				for (int x=c; x<50; x++)
					printf(" ");
				printf("]\n\033[F\033[J");
			}//end supercool loadbar
		}		
		//Vectors of top and antitop
		TLorentzVector Top_MC;
		TLorentzVector ATop_MC;
		TLorentzVector Top_Data;
		TLorentzVector ATop_Data;
		TLorentzVector quark;
		TLorentzVector antiquark;
		
		int weight_is_valid = 1;
		//loop over particles in event and get the top and antitop fourvectors
		for (int j=6;j<8;j++) {
			if (PID[j] == 6) {
				Top_Data = *(new TLorentzVector(px[j],py[j],pz[j],E[j]));
				if (data_or_mc==1 && mc_E[j]!=0.0) 
					Top_MC = *(new TLorentzVector(mc_px[j],mc_py[j],mc_pz[j],mc_E[j]));
				else {
					Top_MC = *(new TLorentzVector(px[j],py[j],pz[j],E[j]));
					weight_is_valid = 0;
				}
			}
			if (PID[j] == -6) {
				ATop_Data = *(new TLorentzVector(px[j],py[j],pz[j],E[j]));
				if (data_or_mc==1 && mc_E[j]!=0.0) 
					ATop_MC = *(new TLorentzVector(mc_px[j],mc_py[j],mc_pz[j],mc_E[j]));
				else {
					ATop_MC = *(new TLorentzVector(px[j],py[j],pz[j],E[j]));
					weight_is_valid = 0;
				}
			}
		}//end of loop over particles
		//assign the quark and antiquark fourvectors from the last two particles in the list of mc_ particles
		if (data_or_mc==1 && mc_E[8]!=1.0 && mc_E[9]!=1.0) {
			quark = *(new TLorentzVector(mc_px[8],mc_py[8],mc_pz[8],mc_E[8]));
			antiquark = *(new TLorentzVector(mc_px[9],mc_py[9],mc_pz[9],mc_E[9]));
		}
		else {
			quark = *(new TLorentzVector(0.0,0.0,sqrt(beam_energy*beam_energy -1*1),beam_energy));
			antiquark = *(new TLorentzVector(0.0,0.0,-1*quark.Pz(),beam_energy));
			weight_is_valid = 0;
		}
		//set the lepton charge
		if (PID[0] > 0)
			Q_l = -1;
		else if (PID[0] < 0)
			Q_l = 1;
		else
			printf("Something went wrong in getting the lepton charge, not sure what.\n");
		//set the lnL
		ln_L = finalChi[0];
		//set the number of valid jets
		n_valid_jets = nValidJets;
		//set the number of btags
		n_bTags = nbTags;
		//set the PIDs of the mother particles
		motherPIDs[0] = motherParticles[0];
		motherPIDs[1] = motherParticles[1];
		
		//Make the 4-vector of the ttbar
		TLorentzVector Q_Data = Top_Data + ATop_Data;
		double ttbar_mass_data=Q_Data.Mag();
		TLorentzVector Q_MC = Top_MC + ATop_MC;
		ttbar_mass_mc=Q_MC.Mag();

		//defining the Px, Py,and Pz, and energies to boost into the ttbar rest frame
		double Bx_data = -1*Q_Data.Px()/Q_Data.E();  
		double By_data = -1*Q_Data.Py()/Q_Data.E();  
		double Bz_data = -1*Q_Data.Pz()/Q_Data.E();
		Qt = sqrt(Q_Data.Px()*Q_Data.Px()+Q_Data.Py()*Q_Data.Py());
		double Bx_mc = -1*Q_MC.Px()/Q_MC.E();  
		double By_mc = -1*Q_MC.Py()/Q_MC.E();  
		double Bz_mc = -1*Q_MC.Pz()/Q_MC.E();
		
		double beta_mc, beta_data;
		double M2_1_data = Top_Data.Mag2();
		double M2_2_data = ATop_Data.Mag2();
		double M2_1_mc   = Top_MC.Mag2();
		double M2_2_mc   = ATop_MC.Mag2();
		double num_data  = 1. - 2.*(M2_1_data+M2_2_data)/(ttbar_mass_data*ttbar_mass_data) + (M2_1_data-M2_2_data)*(M2_1_data-M2_2_data)/(ttbar_mass_data*ttbar_mass_data*ttbar_mass_data*ttbar_mass_data);
		double num_mc    = 1. - 2.*(M2_1_mc+M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc) + (M2_1_mc-M2_2_mc)*(M2_1_mc-M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc*ttbar_mass_mc*ttbar_mass_mc);
		double denom_1_data = (1. + (M2_1_data-M2_2_data)/(ttbar_mass_data*ttbar_mass_data))*(1. + (M2_1_data-M2_2_data)/(ttbar_mass_data*ttbar_mass_data));
		double denom_2_data = (1. + (M2_2_data-M2_1_data)/(ttbar_mass_data*ttbar_mass_data))*(1. + (M2_2_data-M2_1_data)/(ttbar_mass_data*ttbar_mass_data));
		double denom_1_mc   = (1. + (M2_1_mc-M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc))*(1. + (M2_1_mc-M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc));
		double denom_2_mc   = (1. + (M2_2_mc-M2_1_mc)/(ttbar_mass_mc*ttbar_mass_mc))*(1. + (M2_2_mc-M2_1_mc)/(ttbar_mass_mc*ttbar_mass_mc));
		beta_data = sqrt(sqrt((num_data*num_data)/(denom_1_data*denom_1_data) * (num_data*num_data)/(denom_2_data*denom_2_data)));
		beta_mc   = sqrt(sqrt((num_mc*num_mc)/(denom_1_mc*denom_1_mc) * (num_mc*num_mc)/(denom_2_mc*denom_2_mc)));
		
		//Need beta to be real 
		if (TMath::IsNaN(beta_data) || TMath::IsNaN(beta_mc))
			continue;

		//Feynman x
		x_f = 2*Q_Data.Pz()/sqrt_s;
		x_f_mc = 2*Q_MC.Pz()/sqrt_s;
		
		//Creating the Lorentz 4 vectors of the protons
		TLorentzVector Proton_data = *(new TLorentzVector(0,0,sqrt(beam_energy*beam_energy -1*1),beam_energy));
		TLorentzVector Proton2_data = *(new TLorentzVector(0,0,-1*Proton_data.Pz(),beam_energy));
		
		//Doing the boost
		R_data = R_data.Boost(Bx_data,By_data,Bz_data);
		Top_Data = R_data*Top_Data;
		ATop_Data = R_data*ATop_Data;
		Proton_data = R_data*Proton_data;
		Proton2_data = R_data*Proton2_data;
		R_mc = R_mc.Boost(Bx_mc,By_mc,Bz_mc);
		Top_MC = R_mc*Top_MC;
		ATop_MC = R_mc*ATop_MC;
		quark = R_mc*quark;
		antiquark = R_mc*antiquark;
		//Reset the boost
		R_data=S;
		R_mc=S;
		//Define three vectors for P,Pbar,top, quark and antiquark in ttbar c.m frame
		TVector3 top_data = *(new TVector3(Top_Data.Px(),Top_Data.Py(),Top_Data.Pz()));
		TVector3 proton_data = *(new TVector3(Proton_data.Px(),Proton_data.Py(),Proton_data.Pz()));
		TVector3 proton2_data = *(new TVector3(Proton2_data.Px(),Proton2_data.Py(),Proton2_data.Pz()));
		TVector3 top_mc = *(new TVector3(Top_MC.Px(),Top_MC.Py(),Top_MC.Pz()));
		TVector3 true_quark_direction = *(new TVector3(quark.Px(),quark.Py(),quark.Pz()));
		TVector3 true_antiquark_direction = *(new TVector3(antiquark.Px(),antiquark.Py(),antiquark.Pz()));
		
		//Flip the larger one between proton and proton2, and flip the antiquark direction
		if(proton_data.Mag()>proton2_data.Mag()){proton_data=-1.0*proton_data;}else{proton2_data=-1.0*proton2_data;}
		true_antiquark_direction = -1.0*true_antiquark_direction;
		
		//Normalize vectors
		top_data = top_data*(1.0/top_data.Mag());
		proton_data = proton_data*(1.0/proton_data.Mag());
		proton2_data = proton2_data*(1.0/proton2_data.Mag());
		top_mc = top_mc*(1.0/top_mc.Mag());
		true_quark_direction = true_quark_direction*(1.0/true_quark_direction.Mag());
		true_antiquark_direction = true_antiquark_direction*(1.0/true_antiquark_direction.Mag());
		//find the unit bisectors
		TVector3 bisector_data = (proton_data+proton2_data)*(1.0/(proton_data+proton2_data).Mag());
		TVector3 bisector_mc = (true_quark_direction+true_antiquark_direction)*(1.0/(true_quark_direction+true_antiquark_direction).Mag());
		//find the CS angle
		double cos_theta_cs_data=top_data*bisector_data;
		double cos_theta_cs_mc=top_mc*bisector_mc;

		//fill histograms
		cos_theta_mc = cos_theta_cs_mc;

		//printf("sintheta_top = %f, M_top = %f, pT_top = %f, SF_top = %f, sintheta_atop = %f, M_atop = %f, pT_atop = %f, SF_atop = %f, top_pT_weight = %f\n",sintheta_top, M_top, pT_top, SFtop, sintheta_atop, M_atop, pT_atop, SFatop, top_pT_weight);
		
		//weights (cf equation 5 in note)
		if (weight_is_valid != 1) {
			w_a = 0.0; w_a_opp = 0.0; w_s_xi = 0.0; w_s_xi_opp = 0.0; w_a_xi = 0.0; w_a_xi_opp = 0.0; 
			w_s_delta = 0.0; w_s_delta_opp = 0.0; w_a_delta = 0.0; w_a_delta_opp = 0.0;
		}
		else {
			double alpha = 0.0;
			////	THESE ALPHA VALUES MUST BE CHANGED FOR GENERATORS OTHER THAN MADGRAPH 5
			//if (nValidJets[0] == 4)
			//	alpha = -0.228;
			//else if (nValidJets[0] == 5)
			//	alpha = 0.010;
			//VALUES FOR CT10 POWHEG
			//if (nValidJets == 4)
			//	alpha = -0.256;
			//else if (nValidJets == 5)
			//	alpha = 0.143;
			//combined average alpha
			alpha = -0.129;
			double one_m_b2 = 1.0-beta_mc*beta_mc;
			double b2c2 = beta_mc*beta_mc*cos_theta_cs_mc*cos_theta_cs_mc;
			double otb2 = (1.0/3.0)*beta_mc*beta_mc;
			double denom = 1.0+b2c2+one_m_b2+alpha*(1.0-b2c2);
			w_a = 2.0 * ((1.0+otb2+one_m_b2+alpha*(1.0-otb2))/denom) * cos_theta_cs_mc; 
			w_s_xi = one_m_b2/denom;
			w_a_xi = 2.0*(one_m_b2/denom)*cos_theta_cs_mc;
			w_s_delta = (1.0-b2c2)/denom;
			w_a_delta = 2.0*((1.0-otb2)/denom)*cos_theta_cs_mc;
			w_a_opp = 2.0 * ((1.0+otb2+one_m_b2+alpha*(1.0-otb2))/denom) * (-1.0*cos_theta_cs_mc); 
			w_s_xi_opp = one_m_b2/denom;
			w_a_xi_opp = 2.0*(one_m_b2/denom)*(-1.0*cos_theta_cs_mc);
			w_s_delta_opp = (1.0-b2c2)/denom;
			w_a_delta_opp = 2.0*((1.0-otb2)/denom)*(-1.0*cos_theta_cs_mc);
		}
		//Set the other output variables to the correct values
		ttbar_mass=ttbar_mass_data;
		cos_theta_cs=cos_theta_cs_data;

		//copy over the reweighting factors, pileup, and kinematic fit results
		pileup_reweight       = weight_pileup;
		top_pT_reweight       = weight_top_pT;
		GJR_reweight          = weight_GJR_scale;
		CT10_reweight         = weight_CT10_scale;
		cteq_reweight         = weight_cteq_scale;
		btag_eff_reweight     = weight_btag_eff;
		tracking_reweight     = weight_tracking;
		lepID_reweight        = weight_lep_ID;
		lepIso_reweight       = weight_lep_iso;
		trigger_reweight      = weight_trig_eff;
		btag_eff_reweight_low = weight_btag_eff_err;
		tracking_reweight_low = weight_tracking_low;
		lepID_reweight_low    = weight_lep_ID_low;
		lepIso_reweight_low   = weight_lep_iso_low;
		trigger_reweight_low  = weight_trig_eff_low;
		btag_eff_reweight_hi  = weight_btag_eff_err;
		tracking_reweight_hi  = weight_tracking_hi;
		lepID_reweight_hi     = weight_lep_ID_hi;
		lepIso_reweight_hi    = weight_lep_iso_hi;
		trigger_reweight_hi   = weight_trig_eff_hi;
		pileup = pileup_events;
		pileup_real = mc_pileup_events;
		for (int i=0; i<6; i++) {
			fitParValues[i] = bestFitParValues[i];
		}

		//set PDF weights
		int pdf_member = sizeof(Pdf_weights)/sizeof(double);
		double pdf_n = Pdf_weights[0];
		for (int i = 0; i < 100; ++i) {
			if (i<pdf_member)
				Pdf_w[i] = Pdf_weights[i]/pdf_n;
			else
				Pdf_w[i] = 0;
		}

		//Write all data into leaves
		output->Fill();
		kept +=1;

	} //End of loop over entries

	cout << "		Done." << endl;
	
	//write the output file
	output->Write();
	
	//close the file
	file->Close();
	if (count != 0) 
		printf("# OF EVENTS CUT FROM SAMPLE = %d (%f) %%\n",count-kept,100.*(count-kept)/count);
}

void ttbar::SetOutputFilename(const char* outname) {
	output_filename = outname;
}

void ttbar::SetOutputTreeName(const char* treename) {
	output_treename = treename;
}
