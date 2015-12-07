#include "TH1.h"
#include <iostream>
#include <string>
#include <fstream>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFitter.h"
#include "TSystem.h"
#include "ttbar.C"
#include "angles.C"
#include "angles_data.C"
using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////

//total data title and callname
char mg_data_title[100] = {"Data Distribution"};
char mg_data_total_name[50] = {"all_data_angles.root"};

//for fitting
angles_data *ad;


//////////////////////////////////////////////////////////////////////////////////////////////////

//Shouldn't change anything below here between runs

//sample file struct declaration
struct sample_file {			//Holds all the necessary attributes of a sample file
	char* filepath;				//Path to the file
	char* sample_name;			//Callname of the file (used in naming files produced from this file)
	int is_for_distribution;	//integer indicating which distribution the file will be used in creating (1=qq,2=gg,3=bck)
	int nEventsGenerated;		//Number of events generated in the original MC sample
	double cross_section;		//Cross section for the type of sample generated
	char* sample_title;			//Title for histograms made from this file
	struct sample_file* next;
};

//data file struct declaration
struct data_file {			//Holds all the necessary attributes of a sample file
	char* filepath;				//Path to the file
	char* sample_name;			//Callname of the file (used in naming files produced from this file)
	struct data_file* next;	
};

//Bounds on x_F and mass
double xFLow = 0.0;
double xFHigh = 0.6;
double massLow = 350.0;
double massHigh = 1750.0;

//head and tail of sample file linked list
sample_file* MC_HEAD = NULL;
sample_file* MC_TAIL = NULL;
//number of sample files
int nSampleFiles = 0;

//head and tail of data file linked list
data_file* DATA_HEAD = NULL;
data_file* DATA_TAIL = NULL;
//number of sample files
int nDataFiles = 0;

//automatic filescope variables
int iters = 0;
char prepend[100];

//Global vars for running on the grid
//hard coded path
char run_env[400] = "/uscms_data/d3/lfeng7/Legacy_Afb/analysis/CMSSW_7_4_0/src/Legacy_Afb/templates/temp_angles";
int on_the_grid = 0 ;

//Function declarations
void handle_input(int which, char* filename);
void mergeHistoFiles();
void analyzeData();
void mergeAndPlotData();
void buildTemplates();
void fitCombined(char* runName);
void fitSeparated(char* runName);
void minuitfunccombined(int& nDim, double* gout, double& result, double* par, int flg);
void minuitfuncseparated(int& nDim, double* gout, double& result, double* par, int flg);
void move_stuff(char* runname);

//////////////////////////////////////////////////////////////////////////////////////////////////

//MAIN FUNCTION
int main_data(int analyze, int onGrid, int comb_or_sep, char* r, char* mc_input_filename, char* data_input_filename)
{
	char* runname;
	runname=r;
	//Set up the prepend for the directory of all the MC histogram files
	if (onGrid == 0) {
		strcpy(prepend,"./");
		//strcat(prepend,runname);
		strcat(prepend,"histo_files/");
	}
	else if (onGrid == 1) {
		strcpy(prepend,"tardir/");
		strcat(prepend,"histo_files/");
		on_the_grid = 1 ;
	}
	//Handle the input of all the MC histogram files
	handle_input(1,mc_input_filename);
	//Handle the input of all the data files given
	handle_input(2,data_input_filename);
	//analyze dataset to get angles file for the data, and plot histograms of the result
	if (analyze == 1) {
		printf("Analyzing Data Files\n");
		analyzeData();
	}
	//Make, merge, plot, and save the data histo files
	if (onGrid == 0) {
		//merge all the root histogram files to make one containing the four terms
		gSystem->cd(".");
		if (gSystem->Which(gSystem->pwd(),"templates.root") == 0) {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE NOT FOUND, MERGING MC HISTO FILES\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("Merging Histogram Files found at %s\n",prepend);
			mergeHistoFiles();
		}
		else {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE FOUND, NOT MERGING MC HISTO FILES\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		}
		printf("Making, merging, and plotting data _histo.root files\n");
		mergeAndPlotData(); // Use angles.C file to make several 3D distributions, the building block for likelihood
		gSystem->cd(".");
		if (gSystem->Which(gSystem->pwd(),"templates.root") == 0) {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE NOT FOUND, REMAKING SMOOTHED MC TEMPLATES\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("Building smoothed templates for fitting\n");
			buildTemplates();									// Smooth the MC templates
		}
		else {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE FOUND, USING LOADED PREMADE TEMPLATES!!!!\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		}
		//run the minimizer to find Rqqbar and Afb, plot results, and save to a text file.
		printf("Fitting data (get comfortable)\n");
		if (comb_or_sep==0) {
			printf("DOING COMBINED ANALYSIS\n");
			fitCombined(runname);
		}
		else if (comb_or_sep==1) {
			printf("DOING SEPARATED ANALYSIS\n");
			fitSeparated(runname);
		}
		printf("Done\n");
		//move stuff around so the home directory isn't so crowded with files and plots
		//DON'T DO THIS IF RUNNING ON CONDOR!!!
		printf("Moving Files, etc.\n");
		move_stuff(runname);
	}

if (onGrid == 1) {
		//merge all the root histogram files to make one containing the four terms
		gSystem->cd("."); // should i change it to tardir/ ?
		if (gSystem->Which(gSystem->pwd(),"templates.root") == 0) {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE NOT FOUND, MERGING MC HISTO FILES\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("Merging Histogram Files found at %s\n",prepend);
			mergeHistoFiles();
		}
		else {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE FOUND, NOT MERGING MC HISTO FILES\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		}
		printf("Making, merging, and plotting data _histo.root files\n");
		mergeAndPlotData();
		gSystem->cd(".");
		if (gSystem->Which(gSystem->pwd(),"templates.root") == 0) {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE NOT FOUND, REMAKING SMOOTHED MC TEMPLATES\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("Building smoothed templates for fitting\n");
			buildTemplates();
		}
		else {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE FOUND, USING LOADED PREMADE TEMPLATES!!!!\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		}
		//run the minimizer to find Rqqbar and Afb, plot results, and save to a text file.
		printf("Fitting data (get comfortable)\n");
		if (comb_or_sep==0) {
			printf("DOING COMBINED ANALYSIS\n");
			fitCombined(runname);
		}
		else if (comb_or_sep==1) {
			printf("DOING SEPARATED ANALYSIS\n");
			fitSeparated(runname);
		}
		printf("Done\n");
		//move stuff around so the home directory isn't so crowded with files and plots
		//DON'T DO THIS IF RUNNING ON CONDOR!!!
		// printf("Moving Files, etc.\n");
		// move_stuff(runname);
	}

	printf("All Done!\n");
	return 0;
	
}//end main

//make a linked list of all the MC sample_file structs or the data_file structs from the input file
void handle_input(int which, char* filename) {
	//open the file from the name
	FILE* ifp = fopen(filename,"r");
	if (ifp == NULL) {
		printf("ERROR: Inputted file %s cannot be opened. Check filename. Program will keep running (incorrectly)\n",filename);
		return;
	}
	char *fp = (char*)malloc(250*sizeof(char));
	char *sn = (char*)malloc(250*sizeof(char));
	char *ifd = (char*)malloc(10*sizeof(char));
	char *st = (char*)malloc(250*sizeof(char));
	int nEG;
	float c_s;
	//read in data from each line
	//if it's the MC input
	if (which==1) {
		while (fscanf(ifp,"%s %s %s %d %f %[^\n]s",fp,sn,ifd,&nEG,&c_s,st) != EOF) {
			//make a new sample file struct for each line
			sample_file* newFile = (sample_file*)malloc(sizeof(sample_file));
			//Insert it into the list
			newFile->next = NULL;
			if (MC_HEAD == NULL) {
				MC_HEAD = newFile;
				MC_TAIL = newFile;
			}
			else {
				MC_TAIL->next = newFile;
				MC_TAIL = newFile;
			}
			//Populate the new sample_file with the data from the line of the input file
			newFile->filepath = (char*)malloc(250*sizeof(char));
			strcpy(newFile->filepath,fp);
			newFile->sample_name = (char*)malloc(250*sizeof(char));
			strcpy(newFile->sample_name,sn);
			if (strcmp(ifd,"qq")==0)
				newFile->is_for_distribution = 1;
			else if (strcmp(ifd,"gg")==0)
				newFile->is_for_distribution = 2;
			else if (strcmp(ifd,"bck")==0 || strcmp(ifd,"bkg")==0)
				newFile->is_for_distribution = 3;
			else if (strcmp(ifd,"WJets")==0)
				newFile->is_for_distribution = 4;
			else if (strcmp(ifd,"sb")==0)
				newFile->is_for_distribution = -1;
			else 
				printf("ERROR: inputted distribution type %s is not recognized; Program will keep running (incorrectly)\n",ifd);
			newFile->nEventsGenerated = nEG;
			newFile->cross_section = c_s;
			newFile->sample_title = (char*)malloc(250*sizeof(char));
			strcpy(newFile->sample_title,st);
		}
		//get the number of sample files and print it out
		printf("Will merge histograms from the histo files produced from the following files:\n");
		for (sample_file *iter=MC_HEAD;iter!=NULL;iter=iter->next) {
			nSampleFiles = nSampleFiles+1;
			printf("	%s\n",iter->filepath);
		}
		printf("Total number of files is %d\n",nSampleFiles);
	}
	else if (which == 2) {
		while (fscanf(ifp,"%s %[^\n]s",fp,sn) != EOF) {
			//make a new sample file struct for each line
			data_file* newFile = (data_file*)malloc(sizeof(data_file));
			//Insert it into the list
			newFile->next = NULL;
			if (DATA_HEAD == NULL) {
				DATA_HEAD = newFile;
				DATA_TAIL = newFile;
			}
			else {
				DATA_TAIL->next = newFile;
				DATA_TAIL = newFile;
			}
			//Populate the new sample_file with the data from the line of the input file
			newFile->filepath = (char*)malloc(250*sizeof(char));
			strcpy(newFile->filepath,fp);
			newFile->sample_name = (char*)malloc(250*sizeof(char));
			strcpy(newFile->sample_name,sn);
		}
		//get the number of data files and print it out
		printf("Will analyze the following data files:\n");
		for (data_file *iter=DATA_HEAD;iter!=NULL;iter=iter->next) {
			nDataFiles = nDataFiles+1;
			printf("	%s\n",iter->filepath);
		}
		printf("Total number of data files is %d\n",nDataFiles);
	}
	fclose(ifp);
}//end handle_input

//Runs ttbar.C on the data files to get the angles information
void analyzeData() {
	//loop over data files
	int k=1;
	for (data_file *iter=DATA_HEAD;iter!=NULL;iter=iter->next) {
		printf("	Analyzing data file #%d (%s) to find angles data\n",k,iter->filepath);
		//Set up automatic file names
		char out[50];
		strcpy(out,"angles_data_");
		strcat(out,iter->sample_name);
		strcat(out,".root");
		//run ttbar.c to output angles files
		ttbar *t = new ttbar(iter->filepath);
		t->SetOutputFilename(out);
		t->SetOutputTreeName("angles_data");
		t->Loop(0);
		free(t);
		printf("	Done\n");
		k=k+1;
	}//end loop
}//end analyzeData


//Read all histogram files and reweight here
//merges all of the MC sample histogram files in the given directory 
void mergeHistoFiles() { 
	//Create the final file that will hold everything in the end
	TFile* out_file = new TFile("aggregated_distributions.root","Recreate");
	//Create the final TTree that will hold all of the MC events
	float costheta, xF, tt_M;
	float weight, wa, waxi, wadelta, wsxi, wsdelta;
	int nJets, Ql, dist_type;
	TTree *final_tree = new TTree("final_tree","Recreate");
	final_tree->Branch("costheta",&costheta);
	final_tree->Branch("xF",&xF);
	final_tree->Branch("tt_M",&tt_M);
	final_tree->Branch("nJets",&nJets);
	final_tree->Branch("Ql",&Ql);
	final_tree->Branch("dist_type",&dist_type);
	final_tree->Branch("weight",&weight);
	final_tree->Branch("wa",&wa);
	final_tree->Branch("waxi",&waxi);
	final_tree->Branch("wadelta",&wadelta);
	final_tree->Branch("wsxi",&wsxi);
	final_tree->Branch("wsdelta",&wsdelta);
	//get the first histogram to set bins and axis limits
	char file[250];
	strcpy(file,prepend);
	strcat(file,MC_HEAD->sample_name);
	strcat(file,"_histos.root");
	char h_s_plus_4jet[250];
	strcpy(h_s_plus_4jet,MC_HEAD->sample_name);
	strcat(h_s_plus_4jet,"_s_p_4j");
	TFile * f=new TFile(file);
	TH3D* histo = (TH3D*)f->Get(h_s_plus_4jet);
	int nbinsx_local, nbinsy_local, nbinsz_local, nbins_lnL_local;
	double x_low_local, x_high_local, y_low_local, y_high_local, z_low_local, z_high_local, lnL_low_local, lnL_high_local;
	nbinsx_local = histo->GetNbinsX(); nbinsy_local = histo->GetNbinsY(); nbinsz_local = histo->GetNbinsZ();
	x_low_local = histo->GetXaxis()->GetXmin(); x_high_local = histo->GetXaxis()->GetXmax();
	y_low_local = histo->GetYaxis()->GetXmin(); y_high_local = histo->GetYaxis()->GetXmax();
	z_low_local = histo->GetZaxis()->GetXmin(); z_high_local = histo->GetZaxis()->GetXmax();
	TH1D* histo2 = (TH1D*)f->Get("gtt_4jet");
	nbins_lnL_local = histo2->GetNbinsX(); lnL_low_local = histo2->GetXaxis()->GetXmin(); lnL_high_local = histo2->GetXaxis()->GetXmax();
	f->Close();
	//Set the bounds on Feynman x and mass
	xFLow = y_low_local;
	xFHigh = y_high_local;
	massLow = z_low_local;
	massHigh = z_high_local;
	//Set up new histogram names and titles
	char** new_names = (char**)malloc(36*sizeof(char*));
	char** new_titles = (char**)malloc(36*sizeof(char*));
	for (int i=0; i<36; i++) {
		new_names[i] = (char*)malloc(50*sizeof(char));
		new_titles[i] = (char*)malloc(250*sizeof(char));
		strcpy(new_names[i],"f_");
	}
	strcat(new_names[0],"qqs_plus_4jet");	strcat(new_names[1],"qqs_xi_plus_4jet");	strcat(new_names[2],"qqs_delta_plus_4jet");
	strcat(new_names[3],"qqa_plus_4jet");	strcat(new_names[4],"qqa_xi_plus_4jet");	strcat(new_names[5],"qqa_delta_plus_4jet");
	strcat(new_names[6],"qqs_minus_4jet");	strcat(new_names[7],"qqs_xi_minus_4jet");	strcat(new_names[8],"qqs_delta_minus_4jet");
	strcat(new_names[9],"qqa_minus_4jet");	strcat(new_names[10],"qqa_xi_minus_4jet");	strcat(new_names[11],"qqa_delta_minus_4jet");
	strcat(new_names[12],"qqs_plus_5jet");	strcat(new_names[13],"qqs_xi_plus_5jet");	strcat(new_names[14],"qqs_delta_plus_5jet");
	strcat(new_names[15],"qqa_plus_5jet");	strcat(new_names[16],"qqa_xi_plus_5jet");	strcat(new_names[17],"qqa_delta_plus_5jet");
	strcat(new_names[18],"qqs_minus_5jet");	strcat(new_names[19],"qqs_xi_minus_5jet");	strcat(new_names[20],"qqs_delta_minus_5jet");
	strcat(new_names[21],"qqa_minus_5jet");	strcat(new_names[22],"qqa_xi_minus_5jet");	strcat(new_names[23],"qqa_delta_minus_5jet");
	strcat(new_names[24],"gg_plus_4jet"); 	strcat(new_names[25],"bck_plus_4jet");		strcat(new_names[26],"WJets_plus_4jet");
	strcat(new_names[27],"gg_minus_4jet"); 	strcat(new_names[28],"bck_minus_4jet");		strcat(new_names[29],"WJets_minus_4jet");
	strcat(new_names[30],"gg_plus_5jet"); 	strcat(new_names[31],"bck_plus_5jet");		strcat(new_names[32],"WJets_plus_5jet");
	strcat(new_names[33],"gg_minus_5jet"); 	strcat(new_names[34],"bck_minus_5jet");		strcat(new_names[35],"WJets_minus_5jet");
	for (int i=0; i<36; i++) {
		strcpy(new_titles[i],new_names[i]);
		strcat(new_titles[i]," distribution; cos(#theta *); Feynman x (x_{F}); M_{t #bar{t}} (GeV)");
	}
	//Create new histograms
	TH3D** new_histos = (TH3D**)malloc(36*sizeof(TH3D*));
	for (int i=0; i<36; i++) {
		new_histos[i] = new TH3D(new_names[i],new_titles[i],nbinsx_local,x_low_local,x_high_local,
									nbinsy_local,y_low_local,y_high_local,nbinsz_local,z_low_local,z_high_local);
	}
	TH1D *gtt_4jet_local = new TH1D("gtt_4jet_local","t#bar{t} Likelihood distribution, 4jet; -2ln(L)",nbins_lnL_local,lnL_low_local,lnL_high_local);
	TH1D *gbk_4jet_local = new TH1D("gbk_4jet_local","Background Likelihood distribution, 4jet; -2ln(L)",nbins_lnL_local,lnL_low_local,lnL_high_local);
	TH1D *gtt_5jet_local = new TH1D("gtt_5jet_local","t#bar{t} Likelihood distribution, 5jet; -2ln(L)",nbins_lnL_local,lnL_low_local,lnL_high_local);
	TH1D *gbk_5jet_local = new TH1D("gbk_5jet_local","Background Likelihood distribution, 5jet; -2ln(L)",nbins_lnL_local,lnL_low_local,lnL_high_local);
	//add all the appropriate histograms to the final ones
	int k=1;
	for (sample_file *iter=MC_HEAD;iter!=NULL;iter=iter->next) {
		printf("	Adding histograms from sample file #%d (%s)\n",k,iter->filepath);
		//get filename and histogram names
		strcpy(file,prepend);
		strcat(file,iter->sample_name);
		strcat(file,"_histos.root");
		char** file_histo_names = (char**)malloc(25*sizeof(char*));
		for (int i=0; i<25; i++) {
			file_histo_names[i] = (char*)malloc(50*sizeof(char*));
			strcpy(file_histo_names[i],iter->sample_name);
		}
		strcat(file_histo_names[0],"_s_p_4j");			strcat(file_histo_names[1],"_s_xi_p_4j");
		strcat(file_histo_names[2],"_s_delta_p_4j");	strcat(file_histo_names[3],"_a_p_4j");
		strcat(file_histo_names[4],"_a_xi_p_4j");		strcat(file_histo_names[5],"_a_delta_p_4j");
		strcat(file_histo_names[6],"_s_m_4j");			strcat(file_histo_names[7],"_s_xi_m_4j");
		strcat(file_histo_names[8],"_s_delta_m_4j");	strcat(file_histo_names[9],"_a_m_4j");
		strcat(file_histo_names[10],"_a_xi_m_4j");		strcat(file_histo_names[11],"_a_delta_m_4j");
		strcat(file_histo_names[12],"_s_p_5j");			strcat(file_histo_names[13],"_s_xi_p_5j");
		strcat(file_histo_names[14],"_s_delta_p_5j");	strcat(file_histo_names[15],"_a_p_5j");
		strcat(file_histo_names[16],"_a_xi_p_5j");		strcat(file_histo_names[17],"_a_delta_p_5j");
		strcat(file_histo_names[18],"_s_m_5j");			strcat(file_histo_names[19],"_s_xi_m_5j");
		strcat(file_histo_names[20],"_s_delta_m_5j");	strcat(file_histo_names[21],"_a_m_5j");
		strcat(file_histo_names[22],"_a_xi_m_5j");		strcat(file_histo_names[23],"_a_delta_m_5j");
		strcat(file_histo_names[24],"_everything_no_pT_reweight");
		//Pull histograms from file and add to final histograms
		TH3D** file_histos = (TH3D**)malloc(25*sizeof(TH3D*));
		TFile * f2=new TFile(file);
		for (int i=0; i<25; i++)
			file_histos[i] = (TH3D*)f2->Get(file_histo_names[i]);
		double integral_corrected = file_histos[0]->Integral()+file_histos[6]->Integral()+file_histos[12]->Integral()+file_histos[18]->Integral();
		double integral_uncorrected = file_histos[24]->Integral();
		printf("		uncorrected_integral/corrected_integral = %.4f / %.4f = %.4f\n",integral_uncorrected,integral_corrected,integral_uncorrected/integral_corrected);
		for (int i=0; i<25; i++)
			file_histos[i]->Scale(1.0*(integral_uncorrected/integral_corrected));
		TH1D* histo_gtt_4jet = (TH1D*)f2->Get("gtt_4jet");
		TH1D* histo_gbk_4jet = (TH1D*)f2->Get("gbk_4jet");
		TH1D* histo_gtt_5jet = (TH1D*)f2->Get("gtt_5jet");
		TH1D* histo_gbk_5jet = (TH1D*)f2->Get("gbk_5jet");
		histo_gtt_4jet->Scale(1.0*(integral_uncorrected/integral_corrected));
		histo_gbk_4jet->Scale(1.0*(integral_uncorrected/integral_corrected));
		histo_gtt_5jet->Scale(1.0*(integral_uncorrected/integral_corrected));
		histo_gbk_5jet->Scale(1.0*(integral_uncorrected/integral_corrected));
		if (iter->is_for_distribution==1) {
			for (int i=0; i<24; i++)
				new_histos[i]->Add(file_histos[i]);
		}
		else if (iter->is_for_distribution==2) {
			new_histos[24]->Add(file_histos[0]); 
			new_histos[27]->Add(file_histos[6]); 
			new_histos[30]->Add(file_histos[12]); 
			new_histos[33]->Add(file_histos[18]); 
		}
		else if (iter->is_for_distribution==3) {
			new_histos[25]->Add(file_histos[0]); 
			new_histos[28]->Add(file_histos[6]); 
			new_histos[31]->Add(file_histos[12]); 
			new_histos[34]->Add(file_histos[18]); 
		}
		else if (iter->is_for_distribution==4) {
			new_histos[26]->Add(file_histos[0]); 
			new_histos[29]->Add(file_histos[6]); 
			new_histos[32]->Add(file_histos[12]); 
			new_histos[35]->Add(file_histos[18]); 
		}
		else if (iter->is_for_distribution==-1) {
			printf("lol, this is a sideband file\n");
		}
		else {
			printf("ERROR: distribution for sample file #%d (%s) not specified!\n",k,iter->filepath);
			return;
		}
		if (iter->is_for_distribution!=4) {
			gtt_4jet_local->Add(histo_gtt_4jet);
			gbk_4jet_local->Add(histo_gbk_4jet);
			gtt_5jet_local->Add(histo_gtt_5jet);
			gbk_5jet_local->Add(histo_gbk_5jet);
		}
		printf("	Done\n");
		//Open TTree from this file, read through it, rescale the event weights, and add it to the final tree
		printf("	Adding TTree from sample file #%d (%s)\n",k,iter->filepath);
		TTree * thistree=(TTree*)f2->Get("output_tree");
		float this_costheta, this_xF, this_tt_M;
		float this_weight, this_wa, this_waxi, this_wadelta, this_wsxi, this_wsdelta;
		int this_nJets, this_Ql, this_dist_type;
		thistree->SetBranchAddress("costheta",&this_costheta);
		thistree->SetBranchAddress("xF",&this_xF);
		thistree->SetBranchAddress("tt_M",&this_tt_M);
		thistree->SetBranchAddress("nJets",&this_nJets);
		thistree->SetBranchAddress("Ql",&this_Ql);
		thistree->SetBranchAddress("dist_type",&this_dist_type);
		thistree->SetBranchAddress("weight",&this_weight);
		thistree->SetBranchAddress("wa",&this_wa);
		thistree->SetBranchAddress("waxi",&this_waxi);
		thistree->SetBranchAddress("wadelta",&this_wadelta);
		thistree->SetBranchAddress("wsxi",&this_wsxi);
		thistree->SetBranchAddress("wsdelta",&this_wsdelta);
		Long64_t   nentries=thistree->GetEntriesFast();
		for (Long64_t jentry=0; jentry<nentries;jentry++) {
			thistree->GetEntry(jentry);
			costheta = this_costheta;
			xF = this_xF;
			tt_M = this_tt_M;
			nJets = this_nJets;
			Ql = this_Ql;
			dist_type = this_dist_type;
			weight = this_weight*(integral_uncorrected/integral_corrected);
			wa = this_wa;
			waxi = this_waxi;
			wadelta = this_wadelta;
			wsxi = this_wsxi;
			wsdelta = this_wsdelta;
			out_file->cd();
			final_tree->Fill();
		}

		f2->Close();
		printf("	Done\n");
		k=k+1;
	}//end loop over sample files
	//output the new histograms to a new file
	printf("	Aggregating histograms to aggregated_distributions.root\n");
	out_file->cd();
	for (int i=0; i<36; i++)
		new_histos[i]->Write();
	gtt_4jet_local->Write(); gbk_4jet_local->Write();
	gtt_5jet_local->Write(); gbk_5jet_local->Write();
	final_tree->Write();
	//close the file
	out_file->Close();
	//plot the histograms and save them with projections
	char** fn = (char**)malloc(36*sizeof(char*));
	for (int i=0; i<36; i++) {
		fn[i] = (char*)malloc(100*sizeof(char*));
		strcpy(fn[i],new_names[i]);
		strcat(fn[i],"_distribution.pdf");
	}
	for (int i=0; i<36; i++) {
		TCanvas* c = new TCanvas("c",new_names[i],900,900);
		c->Divide(2,2,0.01,0.01,0);
		c->cd(1);
		new_histos[i]->Draw("BOX");
		c->cd(2);
		new_histos[i]->ProjectionY()->Draw();
		c->cd(3);
		new_histos[i]->ProjectionX()->Draw();
		c->cd(4);
		new_histos[i]->ProjectionZ()->Draw();
		c->Print(fn[i],"pdf");
	}//end loop over histograms
	printf("	Done\n");
	printf("Done\n");
}//end mergeHistoFiles

//Runs angles.C on the data files, merges, and saves histograms and projections of the result
void mergeAndPlotData() {
	//Make the TChain to hold all the results and eventually for merging them
	TChain *ch = new TChain("angles_data");
	// TChain *ch = new TChain("angles");
	//loop over data files
	int k=1;
	for (data_file *iter=DATA_HEAD;iter!=NULL;iter=iter->next) {
		//Set up automatic file names
		char out[250];
		//add new angles file to chain
		char in[250];
		if (on_the_grid == 1) {
			sprintf(in,"%s/angles_data_%s.root",run_env,iter->sample_name); //grid					
		} 
		else {
		strcpy(in,"../angles_files/angles_data_");
		strcat(in,iter->sample_name);
		strcat(in,".root");
		}

		ch->Add(in);
		//use angles.C to make histogram files (skip plotting them though)
		printf("	Making histogram file from data file #%d (%s)\n",k,iter->filepath);
		//make filenames
		strcpy(out,iter->sample_name);
		strcat(out,"_data_histos.root");
		//run angles.c on the file to produce a histogram file
		angles *a = new angles();
		a->InitializeNames(in,out,iter->sample_name,iter->sample_name,false);
		a->SetInputTreeName("angles_data");
		// a->SetInputTreeName("angles");
		a->SetMassLimits(massLow,massHigh);
		a->SetxFLimits(xFLow,xFHigh);
		a->Loop(5);
		free(a);
		printf("	Done\n");
		k=k+1;
	}//end loop
	printf("	Merging All %d Data Files\n",nDataFiles);
	//merge all the data files together
	ch->Merge(mg_data_total_name);
	printf("	Done\n");
	printf("	Making histogram file from total data file\n");
	//run angles.c on the file to produce a histogram file
	angles *a = new angles();
	a->InitializeNames(mg_data_total_name,"all_data_histos.root","all_data","All Data",false);
	a->SetInputTreeName("angles_data");
	// a->SetInputTreeName("angles");
	a->SetMassLimits(massLow,massHigh);
	a->SetxFLimits(xFLow,xFHigh);
	a->Loop(5);
	free(a);
	printf("	Done\n");
	//plot and save total data histogram and projections
	TFile *f = new TFile("all_data_histos.root");
	TH3D* hist_plus_4jet = (TH3D*)f->Get("all_data_s_p_4j");
	TH1D* h_x_plus_4jet = (TH1D*)hist_plus_4jet->ProjectionX()->Clone("h_x_plus_4jet");
	TH1D* h_y_plus_4jet = (TH1D*)hist_plus_4jet->ProjectionY()->Clone("h_y_plus_4jet");
	TH1D* h_z_plus_4jet = (TH1D*)hist_plus_4jet->ProjectionZ()->Clone("h_z_plus_4jet");
	TH3D* hist_minus_4jet = (TH3D*)f->Get("all_data_s_m_4j");
	TH1D* h_x_minus_4jet = (TH1D*)hist_minus_4jet->ProjectionX()->Clone("h_x_minus_4jet");
	TH1D* h_y_minus_4jet = (TH1D*)hist_minus_4jet->ProjectionY()->Clone("h_y_minus_4jet");
	TH1D* h_z_minus_4jet = (TH1D*)hist_minus_4jet->ProjectionZ()->Clone("h_z_minus_4jet");
	TH3D* hist_plus_5jet = (TH3D*)f->Get("all_data_s_p_5j");
	TH1D* h_x_plus_5jet = (TH1D*)hist_plus_5jet->ProjectionX()->Clone("h_x_plus_5jet");
	TH1D* h_y_plus_5jet = (TH1D*)hist_plus_5jet->ProjectionY()->Clone("h_y_plus_5jet");
	TH1D* h_z_plus_5jet = (TH1D*)hist_plus_5jet->ProjectionZ()->Clone("h_z_plus_5jet");
	TH3D* hist_minus_5jet = (TH3D*)f->Get("all_data_s_m_5j");
	TH1D* h_x_minus_5jet = (TH1D*)hist_minus_5jet->ProjectionX()->Clone("h_x_minus_5jet");
	TH1D* h_y_minus_5jet = (TH1D*)hist_minus_5jet->ProjectionY()->Clone("h_y_minus_5jet");
	TH1D* h_z_minus_5jet = (TH1D*)hist_minus_5jet->ProjectionZ()->Clone("h_z_minus_5jet");
	//rename the projections
	char name[150];
	strcpy(name,"All Data (Positively Charged Leptons, 4 jets)");
	strcat(name, " X-Projection");
	h_x_plus_4jet->SetTitle(name);
	strcpy(name,"All Data (Positively Charged Leptons, 4 jets)");
	strcat(name, " Y-Projection");
	h_y_plus_4jet->SetTitle(name);
	strcpy(name,"All Data (Positively Charged Leptons, 4 jets)");
	strcat(name, " Z-Projection");
	h_z_plus_4jet->SetTitle(name);
	strcpy(name,"All Data (Negatively Charged Leptons, 4 jets)");
	strcat(name, " X-Projection");
	h_x_minus_4jet->SetTitle(name);
	strcpy(name,"All Data (Negatively Charged Leptons, 4 jets)");
	strcat(name, " Y-Projection");
	h_y_minus_4jet->SetTitle(name);
	strcpy(name,"All Data (Negatively Charged Leptons, 4 jets)");
	strcat(name, " Z-Projection");
	h_z_minus_4jet->SetTitle(name);
	strcpy(name,"All Data (Positively Charged Leptons, 4 jets)");
	strcat(name, " X-Projection");
	h_x_plus_5jet->SetTitle(name);
	strcpy(name,"All Data (Positively Charged Leptons, 5 jets)");
	strcat(name, " Y-Projection");
	h_y_plus_5jet->SetTitle(name);
	strcpy(name,"All Data (Positively Charged Leptons, 5 jets)");
	strcat(name, " Z-Projection");
	h_z_plus_5jet->SetTitle(name);
	strcpy(name,"All Data (Negatively Charged Leptons, 5 jets)");
	strcat(name, " X-Projection");
	h_x_minus_5jet->SetTitle(name);
	strcpy(name,"All Data (Negatively Charged Leptons, 5 jets)");
	strcat(name, " Y-Projection");
	h_y_minus_5jet->SetTitle(name);
	strcpy(name,"All Data (Negatively Charged Leptons, 5 jets)");
	strcat(name, " Z-Projection");
	h_z_minus_5jet->SetTitle(name);
	TCanvas *c1 = new TCanvas("c1","All Data Histogram (Positively Charged Leptons, 4 jets)",900,900);
	c1->Divide(2,2,0.01,0.01,0);
	c1->cd(1);
	hist_plus_4jet->Draw("BOX");
	c1->cd(2);
	h_y_plus_4jet->Draw();
	c1->cd(3);
	h_x_plus_4jet->Draw();
	c1->cd(4);
	h_z_plus_4jet->Draw();
	c1->Print("all_data_histos_plus_4jet.pdf","pdf");
	TCanvas *c2 = new TCanvas("c2","All Data Histogram (Negatively Charged Leptons, 4 jets)",900,900);
	c2->Divide(2,2,0.01,0.01,0);
	c2->cd(1);
	hist_minus_4jet->Draw("BOX");
	c2->cd(2);
	h_y_minus_4jet->Draw();
	c2->cd(3);
	h_x_minus_4jet->Draw();
	c2->cd(4);
	h_z_minus_4jet->Draw();
	c2->Print("all_data_histos_minus_4jet.pdf","pdf");
	TCanvas *c3 = new TCanvas("c3","All Data Histogram (Positively Charged Leptons, 5 jets)",900,900);
	c3->Divide(2,2,0.01,0.01,0);
	c3->cd(1);
	hist_plus_5jet->Draw("BOX");
	c3->cd(2);
	h_y_plus_5jet->Draw();
	c3->cd(3);
	h_x_plus_5jet->Draw();
	c3->cd(4);
	h_z_plus_5jet->Draw();
	c3->Print("all_data_histos_plus_5jet.pdf","pdf");
	TCanvas *c4 = new TCanvas("c4","All Data Histogram (Negatively Charged Leptons, 5 jets)",900,900);
	c4->Divide(2,2,0.01,0.01,0);
	c4->cd(1);
	hist_minus_5jet->Draw("BOX");
	c4->cd(2);
	h_y_minus_5jet->Draw();
	c4->cd(3);
	h_x_minus_5jet->Draw();
	c4->cd(4);
	h_z_minus_5jet->Draw();
	c4->Print("all_data_histos_minus_5jet.pdf","pdf");
	printf("Done\n");
}//end mergeAndPlotData

void buildTemplates() {
	//Build the JSON file that will be parsed by TemplateBuilder
	printf("	Building JSON file for TemplateBuilder. . .\n");
	//create and open the new JSON file; put the first few lines in
	ofstream json;
	json.open("templates.json");
	json << "//Configuration Options\n";
	json << "{\n";
	//input directory is the current directory
	gSystem->ChangeDirectory(".");
	json << "	\"inputDirectory\":\"" << gSystem->pwd() << "/\",\n";
	//outputfile is in the same directory, and is called templates.root
	json << "	\"outputFile\":\"" << gSystem->pwd() << "/templates.root\",\n";
	//start definining all templates
	json << "	//template definitions\n";
	json << "	\"templates\":[\n";
	//lists of template names, weighting factors, selections, and entriesperbin
	char** names = (char**)malloc(36*sizeof(char*));
	char** weights = (char**)malloc(36*sizeof(char*));
	char** selections = (char**)malloc(36*sizeof(char*));
	int entriesperbin[36];
	for (int i=0; i<36; i++) {
		names[i] = (char*)malloc(50*sizeof(char));
		weights[i] = (char*)malloc(100*sizeof(char));
		selections[i] = (char*)malloc(300*sizeof(char));
		strcpy(names[i],"f_");
		strcpy(weights[i],"weight");
	}
	strcat(names[0],"qqs_plus_4jet");			strcpy(selections[0],"dist_type == 1 && nJets==4 && Ql>0");		//no extra weighting factor
	strcat(names[1],"qqs_xi_plus_4jet");		strcpy(selections[1],"dist_type == 1 && nJets==4 && Ql>0");		strcat(weights[1],"*wsxi");
	strcat(names[2],"qqs_delta_plus_4jet");		strcpy(selections[2],"dist_type == 1 && nJets==4 && Ql>0");		strcat(weights[2],"*wsdelta");
	strcat(names[3],"qqa_plus_4jet");			strcpy(selections[3],"dist_type == 1 && nJets==4 && Ql>0");		strcat(weights[3],"*wa");
	strcat(names[4],"qqa_xi_plus_4jet");		strcpy(selections[4],"dist_type == 1 && nJets==4 && Ql>0");		strcat(weights[4],"*waxi");
	strcat(names[5],"qqa_delta_plus_4jet");		strcpy(selections[5],"dist_type == 1 && nJets==4 && Ql>0");		strcat(weights[5],"*wadelta");
	strcat(names[6],"qqs_minus_4jet");			strcpy(selections[6],"dist_type == 1 && nJets==4 && Ql<0");		//no extra weighting factor
	strcat(names[7],"qqs_xi_minus_4jet");		strcpy(selections[7],"dist_type == 1 && nJets==4 && Ql<0");		strcat(weights[7],"*wsxi");
	strcat(names[8],"qqs_delta_minus_4jet");	strcpy(selections[8],"dist_type == 1 && nJets==4 && Ql<0");		strcat(weights[8],"*wsdelta");
	strcat(names[9],"qqa_minus_4jet");			strcpy(selections[9],"dist_type == 1 && nJets==4 && Ql<0");		strcat(weights[9],"*wa");
	strcat(names[10],"qqa_xi_minus_4jet");		strcpy(selections[10],"dist_type == 1 && nJets==4 && Ql<0");	strcat(weights[10],"*waxi");
	strcat(names[11],"qqa_delta_minus_4jet");	strcpy(selections[11],"dist_type == 1 && nJets==4 && Ql<0");	strcat(weights[11],"*wadelta");
	strcat(names[12],"qqs_plus_5jet");			strcpy(selections[12],"dist_type == 1 && nJets==5 && Ql>0");	//no extra weighting factor
	strcat(names[13],"qqs_xi_plus_5jet");		strcpy(selections[13],"dist_type == 1 && nJets==5 && Ql>0");	strcat(weights[13],"*wsxi");
	strcat(names[14],"qqs_delta_plus_5jet");	strcpy(selections[14],"dist_type == 1 && nJets==5 && Ql>0");	strcat(weights[14],"*wsdelta");
	strcat(names[15],"qqa_plus_5jet");			strcpy(selections[15],"dist_type == 1 && nJets==5 && Ql>0");	strcat(weights[15],"*wa");
	strcat(names[16],"qqa_xi_plus_5jet");		strcpy(selections[16],"dist_type == 1 && nJets==5 && Ql>0");	strcat(weights[16],"*waxi");
	strcat(names[17],"qqa_delta_plus_5jet");	strcpy(selections[17],"dist_type == 1 && nJets==5 && Ql>0");	strcat(weights[17],"*wadelta");
	strcat(names[18],"qqs_minus_5jet");			strcpy(selections[18],"dist_type == 1 && nJets==5 && Ql<0");	//no extra weighting factor
	strcat(names[19],"qqs_xi_minus_5jet");		strcpy(selections[19],"dist_type == 1 && nJets==5 && Ql<0");	strcat(weights[19],"*wsxi");
	strcat(names[20],"qqs_delta_minus_5jet");	strcpy(selections[20],"dist_type == 1 && nJets==5 && Ql<0");	strcat(weights[20],"*wsdelta");
	strcat(names[21],"qqa_minus_5jet");			strcpy(selections[21],"dist_type == 1 && nJets==5 && Ql<0");	strcat(weights[21],"*wa");
	strcat(names[22],"qqa_xi_minus_5jet");		strcpy(selections[22],"dist_type == 1 && nJets==5 && Ql<0");	strcat(weights[22],"*waxi");
	strcat(names[23],"qqa_delta_minus_5jet");	strcpy(selections[23],"dist_type == 1 && nJets==5 && Ql<0");	strcat(weights[23],"*wadelta");
	strcat(names[24],"gg_plus_4jet");			strcpy(selections[24],"dist_type == 2 && nJets==4 && Ql>0");	//no extra weighting factor
	strcat(names[25],"bck_plus_4jet");			strcpy(selections[25],"dist_type == 3 && nJets==4 && Ql>0");	//no extra weighting factor
	strcat(names[26],"WJets_plus_4jet");		strcpy(selections[26],"dist_type == 4 && nJets==4 && Ql>0");	//no extra weighting factor
	strcat(names[27],"gg_minus_4jet");			strcpy(selections[27],"dist_type == 2 && nJets==4 && Ql<0");	//no extra weighting factor
	strcat(names[28],"bck_minus_4jet");			strcpy(selections[28],"dist_type == 3 && nJets==4 && Ql<0");	//no extra weighting factor
	strcat(names[29],"WJets_minus_4jet");		strcpy(selections[29],"dist_type == 4 && nJets==4 && Ql<0");	//no extra weighting factor
	strcat(names[30],"gg_plus_5jet");			strcpy(selections[30],"dist_type == 2 && nJets==5 && Ql>0");	//no extra weighting factor
	strcat(names[31],"bck_plus_5jet");			strcpy(selections[31],"dist_type == 3 && nJets==5 && Ql>0");	//no extra weighting factor
	strcat(names[32],"WJets_plus_5jet");		strcpy(selections[32],"dist_type == 4 && nJets==5 && Ql>0");	//no extra weighting factor
	strcat(names[33],"gg_minus_5jet");			strcpy(selections[33],"dist_type == 2 && nJets==5 && Ql<0");	//no extra weighting factor
	strcat(names[34],"bck_minus_5jet");			strcpy(selections[34],"dist_type == 3 && nJets==5 && Ql<0");	//no extra weighting factor
	strcat(names[35],"WJets_minus_5jet");		strcpy(selections[35],"dist_type == 4 && nJets==5 && Ql<0");	//no extra weighting factor
	//filling the entriesperbin is a bit annoying since we need to know the number of effective entries in the existing unsmoothed distributions
	TFile *f = new TFile("aggregated_distributions.root");
	// for (int i=0; i<6; i++)
		// entriesperbin[i] = (int)(((TH3F*)f->Get(names[0]))->GetEffectiveEntries()*0.0001)+1;
	// for (int i=6; i<12; i++)
		// entriesperbin[i] = (int)(((TH3F*)f->Get(names[6]))->GetEffectiveEntries()*0.0001)+1;
	// for (int i=12; i<18; i++)
		// entriesperbin[i] = (int)(((TH3F*)f->Get(names[12]))->GetEffectiveEntries()*0.0001)+1;
	// for (int i=18; i<24; i++)
		// entriesperbin[i] = (int)(((TH3F*)f->Get(names[18]))->GetEffectiveEntries()*0.0001)+1;
	for (int i=0; i<24; i++)
		entriesperbin[i] = 1;
	for (int i=24; i<36; i++)
		// entriesperbin[i] = (int)(((TH3F*)f->Get(names[i]))->GetEffectiveEntries()*0.000025)+1;
		entriesperbin[i] = 2;
	int nBinsXLocal = ((TH3F*)f->Get(names[0]))->GetNbinsX();
	int nBinsYLocal = ((TH3F*)f->Get(names[0]))->GetNbinsY();
	int nBinsZLocal = ((TH3F*)f->Get(names[0]))->GetNbinsZ();
	double costhetalow = ((TH3F*)f->Get(names[0]))->GetXaxis()->GetXmin();
	double costhetahigh = ((TH3F*)f->Get(names[0]))->GetXaxis()->GetXmax();
	f->Close();
	//build the "templates" section of the JSON file
	for (int i=0; i<36; i++) {
		if (i==3 || i==4 || i==5 || i==9 || i==10 || i==11 || i==15 || i==16 || i==17 || i==21 || i==22 || i==23)
			continue;
		//build just the histograms with exclusively positive weights
		json << "		//" << names[i] << "\n";
		json << "		{\n";
		json << "			\"name\":\"" << names[i] << "\",\n";
		json << "			\"files\":[\"aggregated_distributions.root\"],\n";
		json << "			\"tree\":\"final_tree\",\n";
		json << "			\"variables\":[\"costheta\",\"xF\",\"tt_M\"],\n";
		json << "			\"weight\":\"" << weights[i] << "\",\n";
		json << "			\"conserveSumOfWeights\":true,\n";
		json << "			\"selection\":\"" << selections[i] << "\",\n";
		json << "			\"assertion\":\"nJets>3 && nJets<6 && Ql != 0\",\n";
		json << "			\"binning\":{\n";
		json << "				\"type\":\"fixed\",\n";
		json << "				\"bins\":[" << nBinsXLocal << "," << costhetalow << "," << costhetahigh << "," << nBinsYLocal << "," << xFLow << "," << xFHigh << "," << nBinsZLocal << "," << massLow << "," << massHigh << "]\n";
		if (i<24 || i%3==0) {
		// if (i<24) {	
		// if (false) {
		// if (true) {
			json << "			}\n";
		}
		else {
			json << "			},\n";
			json << "			\"postprocessing\":[\n";
			json << "				{\"type\":\"smooth\", \"kernel\":\"adaptive\", \"entriesperbin\":" << entriesperbin[i] << "},\n";
			json << "				{\"type\":\"reweight\", \"axes\":[0,1,2]}\n";
			json << "			]\n";
			printf("MINIMUM NUMBER OF EVENTS PER BIN FOR DISTRIBUTION %s IS %d\n",names[i],entriesperbin[i]);
		}
		json << "		},\n";
	}
	for (int i=0; i<24; i++) {
		if (i!=3 && i!=4 && i!=5 && i!=9 && i!=10 && i!=11 && i!=15 && i!=16 && i!=17 && i!=21 && i!=22 && i!=23)
			continue;
		// build the histograms with negative weights: actually a linear combination of two templates because damn
		// positive weights
		json << "		//" << names[i] << "_positive_weights\n";
		json << "		{\n";
		json << "			\"name\":\"" << names[i] << "_positive_weights\",\n";
		json << "			\"files\":[\"aggregated_distributions.root\"],\n";
		json << "			\"tree\":\"final_tree\",\n";
		json << "			\"variables\":[\"costheta\",\"xF\",\"tt_M\"],\n";
		json << "			\"weight\":\"" << weights[i] << "\",\n";
		json << "			\"conserveSumOfWeights\":true,\n";
		json << "			\"selection\":\"" << selections[i] << " && (" << weights[i] << ")>0\",\n";
		json << "			\"assertion\":\"nJets>3 && nJets<6 && Ql != 0\",\n";
		json << "			\"binning\":{\n";
		json << "				\"type\":\"fixed\",\n";
		json << "				\"bins\":[" << nBinsXLocal << "," << costhetalow << "," << costhetahigh << "," << nBinsYLocal << "," << xFLow << "," << xFHigh << "," << nBinsZLocal << "," << massLow << "," << massHigh << "]\n";
		json << "			}\n"; //FIX THIS LINE WHEN YOU ADD SMOOTHING BACK IN
		// json << "			\"postprocessing\":[\n";
		// json << "				{\"type\":\"smooth\", \"kernel\":\"adaptive\", \"entriesperbin\":" << entriesperbin[i] << "}\n";
		// json << "			]\n";
		json << "		},\n";
		// negative weights
		json << "		//" << names[i] << "_negative_weights\n";
		json << "		{\n";
		json << "			\"name\":\"" << names[i] << "_negative_weights\",\n";
		json << "			\"files\":[\"aggregated_distributions.root\"],\n";
		json << "			\"tree\":\"final_tree\",\n";
		json << "			\"variables\":[\"costheta\",\"xF\",\"tt_M\"],\n";
		json << "			\"weight\":\"" << weights[i] << "*-1.0\",\n";
		json << "			\"conserveSumOfWeights\":true,\n";
		json << "			\"selection\":\"" << selections[i] << " && (" << weights[i] << ")<0\",\n";
		json << "			\"assertion\":\"nJets>3 && nJets<6 && Ql != 0\",\n";
		json << "			\"binning\":{\n";
		json << "				\"type\":\"fixed\",\n";
		json << "				\"bins\":[" << nBinsXLocal << "," << costhetalow << "," << costhetahigh << "," << nBinsYLocal << "," << xFLow << "," << xFHigh << "," << nBinsZLocal << "," << massLow << "," << massHigh << "]\n";
		json << "			}\n"; //FIX THIS LINE WHEN YOU ADD SMOOTHING BACK IN
		// json << "			\"postprocessing\":[\n";
		// json << "				{\"type\":\"smooth\", \"kernel\":\"adaptive\", \"entriesperbin\":" << entriesperbin[i] << "}\n";
		// json << "			]\n";
		json << "		},\n";
		// sum of the two
		json << "		//" << names[i] << "\n";
		json << "		{\n";
		json << "		\"name\":\"" << names[i] << "\",\n";
		json << "		\"templatesum\":[\n";
		json << "			{\"name\":\"" << names[i] <<"_positive_weights\",\"factor\":1.0},\n";
		json << "			{\"name\":\"" << names[i] <<"_negative_weights\",\"factor\":-1.0}\n";
		json << "		],\n";
		json << "		\"conserveSumOfWeights\":true\n";
		json << "		}";
		if (i<23)
			json << ",";
		json << "\n";
	}
	//finish the json file
	json << "	]\n";
	json << "}\n";
	json.close();
	printf("	Done\n");
	//call templateBuilder
	printf("	Building Smoothed Templates. . .\n");
	gSystem->Exec("$CMSSW_BASE/src/TemplateBuilder/buildTemplate.exe templates.json");
	printf("	Done.\n");
}//end buildTemplates

//fits the data with the expected distribution
void fitCombined(char* runName) {
	//Create the angles_data object and load in histograms
	ad = new angles_data(mg_data_total_name);
	ad->LoadHistogramsCombined();
	//Create the minimizer
	TFitter* minimizer = new TFitter();
	//Set the function to be minimized
	minimizer->SetFCN(minuitfunccombined);
	//Set the parameters
	//arg1 - parameter number
	//arg2 - parameter name
	//arg3 - first guess at the parameter value
	//arg4 - estimated distance to minimum
	//arg5 - lower value
	//arg6 - upper value 
	minimizer->SetParameter(0,"Rqqbar", 0.06, 0.03,  0.0, 1.0);

	minimizer->SetParameter(1,"Rbck",   0.18, 0.03,  0.0, 1.0); //floating
	// minimizer->SetParameter(1,"Rbck",   0.171, 0.03,  0.0, 1.0); //fixed
	// minimizer->SetParameter(1,"Rbck",   0.3633, 0.03,  0.0, 1.0); //fixed, 1.1*Rbck
	// minimizer->SetParameter(1,"Rbck",   0.2837, 0.03,  0.0, 1.0); //fixed, restricted eta
	// minimizer->SetParameter(1,"Rbck",   0.3121, 0.03,  0.0, 1.0); //fixed, 1.1*Rbck, restricted eta
	// minimizer->SetParameter(1,"Rbck",   0.2871, 0.03,  0.0, 1.0); //fixed, harder jet cuts
	// minimizer->SetParameter(1,"Rbck",   0.3158, 0.03,  0.0, 1.0); //fixed, 1.1*Rbck, harder jet cuts
	// minimizer->SetParameter(1,"Rbck",   0.1738, 0.03,  0.0, 1.0); //fixed, restricted costheta

	minimizer->SetParameter(2,"RWjets", 0.15, 0.05,0.0,1.0); //floating
	// minimizer->SetParameter(2,"RWjets", 0.037, 0.05,0.0,1.0); //fixed
	// minimizer->SetParameter(2,"RWjets", 0.1298, 0.05,0.0,1.0); //fixed, restricted costheta

	minimizer->SetParameter(3,"xi", 	0.00, 0.50, -1.0, 1.0);
	minimizer->SetParameter(4,"delta",  0.00, 0.50, -1.0, 1.0);
	minimizer->SetParameter(5,"Afb",    0.05, 0.10, -1.0, 1.0);
	//minimizer->FixParameter(0);
	// minimizer->FixParameter(1);	//Rbck
	// minimizer->FixParameter(2);	//RWjets
	minimizer->FixParameter(3);
	minimizer->FixParameter(4);
	//minimizer->FixParameter(5);
	//set up arguments and minimize
	double arg[2];
	arg[0] = 5000; //maximum number of iterations
	arg[1] = 1000;  //function tolerance*1000
	minimizer->ExecuteCommand("MIGRAD",arg,2);
	double arg2[1]; 
	arg2[0] = 5000;
	minimizer->ExecuteCommand("HESSE",arg2,1);
	//Get the resulting parameters and errors
	double fit_Rqqbar     = minimizer->GetParameter(0);
	double fit_Rqqbar_err = minimizer->GetParError(0);
	double fit_Rbck = minimizer->GetParameter(1);
	double fit_Rbck_err = minimizer->GetParError(1);
	double fit_RWJets  = minimizer->GetParameter(2);
	double fit_RWJets_err = minimizer->GetParError(2);
	double fit_xi = minimizer->GetParameter(3);
	double fit_xi_err = minimizer->GetParError(3);
	double fit_delta = minimizer->GetParameter(4);
	double fit_delta_err = minimizer->GetParError(4);
	double fit_Afb     = minimizer->GetParameter(5);
	double fit_Afb_err = minimizer->GetParError(5);
	//make a new object and call the other Loop function to make plots of the results
	free(ad);
	angles_data *ad2 = new angles_data(mg_data_total_name);
	ad2->LoadHistogramsCombined();
	ad2->Loop(fit_Rqqbar,fit_Rqqbar_err,fit_Rbck,fit_Rbck_err,fit_RWJets,fit_RWJets_err,
			  fit_xi,fit_xi_err,fit_delta,fit_delta_err,fit_Afb,fit_Afb_err,runName);
	free(ad2);
}

//fits the data with the expected distribution and changes the values of the global parameters Rqqbar and Afb
void fitSeparated(char* runName) {
	//Create the angles_data object and load in histograms
	ad = new angles_data(mg_data_total_name);
	ad->LoadHistogramsSeparated();
	//Create the minimizer
	TFitter* minimizer = new TFitter();
	//Set the function to be minimized
	minimizer->SetFCN(minuitfuncseparated);
	//Set the parameters
	//arg1 - parameter number
	//arg2 - parameter name
	//arg3 - first guess at the parameter value
	//arg4 - estimated distance to minimum
	//arg5 - lower value
	//arg6 - upper value 
	minimizer->SetParameter(0,"Rqqbar_4jet", 0.10, 0.05,  0.0, 1.0);

	minimizer->SetParameter(1,"Rbck_4jet",   0.20, 0.05,  0.0, 1.0);   //floating
	// minimizer->SetParameter(1,"Rbck_4jet",   0.193, 0.05,  0.0, 1.0); //fixed
	// minimizer->SetParameter(1,"Rbck_4jet",   0.2170, 0.05,  0.0, 1.0);  //fixed, 1.1*Rbck

	minimizer->SetParameter(2,"R_W4Jets",    0.18, 0.05,  0.0, 1.0); //floating
	// minimizer->SetParameter(2,"R_W4Jets",    0.043, 0.05,  0.0, 1.0); //fixed

	minimizer->SetParameter(3,"Rqqbar_5jet", 0.05, 0.05,  0.0, 1.0);

	minimizer->SetParameter(4,"Rbck_5jet",   0.15, 0.05,  0.0, 1.0);   //floating
	// minimizer->SetParameter(4,"Rbck_5jet",   0.140, 0.05,  0.0, 1.0); //fixed
	// minimizer->SetParameter(4,"Rbck_5jet",   0.2823, 0.05,  0.0, 1.0);   //fixed, 1.1*Rbck
	
	minimizer->SetParameter(5,"R_W5Jets",    0.11, 0.05,  0.0, 1.0); //floating
	// minimizer->SetParameter(5,"R_W5Jets",    0.028, 0.05,  0.0, 1.0); //fixed

	minimizer->SetParameter(6,"xi_4jet", 	 0.00, 0.50, -1.0, 1.0);
	minimizer->SetParameter(7,"xi_5jet", 	 0.00, 0.50, -1.0, 1.0);
	minimizer->SetParameter(8,"delta_4jet",  0.00, 0.50, -1.0, 1.0);
	minimizer->SetParameter(9,"delta_5jet",  0.00, 0.50, -1.0, 1.0);
	minimizer->SetParameter(10,"Afb_4jet",    0.05, 0.10, -1.0, 1.0);
	minimizer->SetParameter(11,"Afb_5jet",    0.05, 0.10, -1.0, 1.0);
	// minimizer->FixParameter(0);
	// minimizer->FixParameter(1);			//Rbck_4jet
	// minimizer->FixParameter(2);		//R_W4Jets
	// minimizer->FixParameter(3);		
	// minimizer->FixParameter(4);			//Rbck_5jet
	// minimizer->FixParameter(5);		//R_W5Jets
	minimizer->FixParameter(6);
	minimizer->FixParameter(7);
	minimizer->FixParameter(8);
	minimizer->FixParameter(9);
	// minimizer->FixParameter(10);
	// minimizer->FixParameter(11);
	//set up arguments and minimize
	double arg[2];
	arg[0] = 5000; //maximum number of iterations
	arg[1] = 1000;  //function tolerance*1000
	minimizer->ExecuteCommand("MIGRAD",arg,2);
	double arg2[1]; 
	arg2[0] = 5000;
	minimizer->ExecuteCommand("HESSE",arg2,1);
	//Get the resulting parameters and errors
	double fit_Rqqbar_4jet     = minimizer->GetParameter(0);
	double fit_Rqqbar_4jet_err = minimizer->GetParError(0);
	double fit_Rbck_4jet 	   = minimizer->GetParameter(1);
	double fit_Rbck_4jet_err   = minimizer->GetParError(1);
	double fit_RW4Jets		   = minimizer->GetParameter(2);
	double fit_RW4Jets_err	   = minimizer->GetParError(2);
	double fit_Rqqbar_5jet     = minimizer->GetParameter(3);
	double fit_Rqqbar_5jet_err = minimizer->GetParError(3);
	double fit_Rbck_5jet 	   = minimizer->GetParameter(4);
	double fit_Rbck_5jet_err   = minimizer->GetParError(4);
	double fit_RW5Jets		   = minimizer->GetParameter(5);
	double fit_RW5Jets_err	   = minimizer->GetParError(5);
	double fit_xi_4jet 		   = minimizer->GetParameter(6);
	double fit_xi_4jet_err 	   = minimizer->GetParError(6);
	double fit_xi_5jet 		   = minimizer->GetParameter(7);
	double fit_xi_5jet_err 	   = minimizer->GetParError(7);
	double fit_delta_4jet 	   = minimizer->GetParameter(8);
	double fit_delta_4jet_err  = minimizer->GetParError(8);
	double fit_delta_5jet 	   = minimizer->GetParameter(9);
	double fit_delta_5jet_err  = minimizer->GetParError(9);
	double fit_Afb_4jet        = minimizer->GetParameter(10);
	double fit_Afb_4jet_err    = minimizer->GetParError(10);
	double fit_Afb_5jet        = minimizer->GetParameter(11);
	double fit_Afb_5jet_err    = minimizer->GetParError(11);
	//make a new object and call the other Loop function to make plots of the results
	free(ad);
	angles_data *ad2 = new angles_data(mg_data_total_name);
	ad2->LoadHistogramsSeparated();
	ad2->Loop(fit_Rqqbar_4jet,fit_Rqqbar_4jet_err,fit_Rbck_4jet,fit_Rbck_4jet_err,fit_RW4Jets,fit_RW4Jets_err,
			  fit_Rqqbar_5jet,fit_Rqqbar_5jet_err,fit_Rbck_5jet,fit_Rbck_5jet_err,fit_RW5Jets,fit_RW5Jets_err,
			  fit_xi_4jet,fit_xi_4jet_err,fit_xi_5jet,fit_xi_5jet_err,
			  fit_delta_4jet,fit_delta_4jet_err,fit_delta_5jet,fit_delta_5jet_err,
			  fit_Afb_4jet,fit_Afb_4jet_err,fit_Afb_5jet,fit_Afb_5jet_err,runName);
	free(ad2);
}

//The actual Likelihood fitting function: returns -2ln(L) for a given set of parameters
//COMBINED CASE
double myfunc(double R_qqbar, double R_bck, double R_WJets, double xi, double delta, double A_fb){
	iters = iters+1;
	
	//holds -2ln(L) value to be incremented
	double logL = 0;
	
	//get the log likelihood from the angles_data object's loop
	logL = ad->Loop(R_qqbar,R_bck,R_WJets,xi,delta,A_fb,0);
	
	if (iters%5==0) {
		printf("-2ln(L) at iteration %d = %.10f\n",iters,logL);
		printf("----------------------------------------------------------\n");
	}
	
	//return the -2ln(L) value for this iteration of parameters
	return logL;
}

//SEPARATED CASE
double myfunc(double R_qqbar_4jet, double R_bck_4jet, double R_W4Jets, double R_qqbar_5jet, double R_bck_5jet, double R_W5Jets,
				double xi_4jet, double xi_5jet, double delta_4jet, double delta_5jet, double A_fb_4jet, double A_fb_5jet){
	iters = iters+1;
	
	//holds -2ln(L) value to be incremented
	double logL = 0;
	
	//get the log likelihood from the angles_data object's loop
	logL = ad->Loop(R_qqbar_4jet,R_bck_4jet,R_W4Jets,R_qqbar_5jet,R_bck_5jet,R_W5Jets,xi_4jet,xi_5jet,delta_4jet,delta_5jet,A_fb_4jet,A_fb_5jet);
	
	if (iters%5==0) {
		printf("-2ln(L) at iteration %d = %.10f\n",iters,logL);
		printf("----------------------------------------------------------\n");
	}
	
	//return the -2ln(L) value for this iteration of parameters
	return logL;
}

//runs a few shell commands to organize the output of the run.
void move_stuff(char* runname) {
	char dirname[150];
	strcpy(dirname,runname);
	strcat(dirname,"_output");
	char cmd[150];
	strcpy(cmd,"mkdir ");
	strcat(cmd,dirname);
	//gSystem->Exec(cmd);
	//gSystem->Exec("mkdir data_angles_files");
	gSystem->Exec("mkdir data_histo_files");
	gSystem->Exec("mkdir final_fit_stuff");
	//gSystem->Exec("mv angles_*.root data_angles_files");
	gSystem->Exec("mv *_data_histos.root data_histo_files");
	gSystem->Exec("mv f_*_distribution.pdf final_fit_stuff");
	gSystem->Exec("mv all_data_histos*.* final_fit_stuff");
	gSystem->Exec("mv all_data_angles.* final_fit_stuff");
	// gSystem->Exec("mv aggregated_distributions.root final_fit_stuff");
	gSystem->Exec("mv fit_results.txt final_fit_stuff");
	gSystem->Exec("mv fit_comparison.pdf final_fit_stuff");
	gSystem->Exec("mv sideband.pdf final_fit_stuff");
	gSystem->Exec("mv ev_lnL_*.pdf final_fit_stuff");
	//strcpy(cmd,"rsync -avh data_angles_files ");
	//strcat(cmd,dirname);
	//gSystem->Exec(cmd);
	//gSystem->Exec("rm -rf data_angles_files");
	strcpy(cmd,"rsync -aq data_histo_files ");
	strcat(cmd,dirname);
	gSystem->Exec(cmd);
	gSystem->Exec("rm -rf data_histo_files");
	strcpy(cmd,"rsync -aq final_fit_stuff ");
	strcat(cmd,dirname);
	gSystem->Exec(cmd);
	gSystem->Exec("rm -rf final_fit_stuff");
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	        Accessory Functions for the Minimizer; not important and shouldn't be changed.                 //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Function should be defined like this but we'll only use the result and the par arguments
void minuitfunccombined(int& nDim, double* gout, double& result, double* par, int flg){
	//kill some bloody compiler warnings
	nDim=nDim;
	gout=gout;
	par=par;
	flg=flg;
	//calls our function to get the chi^2 for current set of parameters
	result= myfunc(par[0],par[1],par[2],par[3],par[4],par[5]);
}

//function for separated case
void minuitfuncseparated(int& nDim, double* gout, double& result, double* par, int flg){
	//kill some bloody compiler warnings
	nDim=nDim;
	gout=gout;
	par=par;
	flg=flg;
	//calls our function to get the chi^2 for current set of parameters
	result= myfunc(par[0],par[1],par[2],par[3],par[4],par[5],par[6],par[7],par[8],par[9],par[10],par[11]);
}
