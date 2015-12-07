#include "TH1.h"
#include <iostream>
#include <string>
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFitter.h"
#include "TSystem.h"
#include "ttbar.C"
#include "angles.C"
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////


// Tone of detector reweighting
const char* detector_tune = "The tune of detector reweighting is : nominal ";

//binning, mass range, and xf range defaults
int XBINS = 20;
int YBINS = 30;
int ZBINS = 40;
double LOWMASS = 350.0;
double HIGHMASS = 1750.0;
double XFLOW = 0.0;
double XFHIGH = 0.6;
int pdf_member = -1; // pdf weight: Default pdf member

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

//head and tail of sample file linked list
sample_file* HEAD = NULL;
sample_file* TAIL = NULL;
//number of sample files
int nSampleFiles = 0;
//wether on the grid
int on_the_grid = 0;
char dir_angles[400] = "/Users/ray/HEP/analysis/Legacy_Afb/templates/temp_angles";

//Function declarations
//Add pdf_index as a parameter
int main_MC(int analyze, int onGrid, char* r, char* input_filename, int index);
int main_MC(int analyze, int onGrid, char* r, char* input_filename, int index, int xbins, int ybins, int zbins, double lowmass, double highmass, double xflow, double xfhigh);
void handle_input(char* filename);
void makeAnglesFiles();
void makeHistosFiles(int xBins,int yBins,int zBins,double mass_lower,double mass_upper, double xflow, double xfhigh);
void move_stuff(char* runname);

//////////////////////////////////////////////////////////////////////////////////////////////////


//MAIN FUNCTION
int main_MC(int analyze, int onGrid, char* r, char* input_filename, int index)
{
	return main_MC(analyze,onGrid,r,input_filename,index,XBINS,YBINS,ZBINS,LOWMASS,HIGHMASS,XFLOW,XFHIGH);
}
int main_MC(int analyze, int onGrid, char* r, char* input_filename, int index, int xbins, int ybins, int zbins, double lowmass, double highmass, double xflow, double xfhigh)
{

	// For systematics evaluation only. To remind that we are using the modified reweighting factor indeed.
	cout<<"\n"<<detector_tune<<"\n"<<endl;

	char* runname;
	runname=r;

	//pdf weight: take input of pdf_index and asign it to the global var pdf_member. The pdf_member will be used later on 
	//in angles::SetPdfMembers(int pdf_member)
	pdf_member = index;

	//take in the input file and make the linked list of sample file structs
	printf("Handling input file %s :\n",input_filename);
	handle_input(input_filename);
	//analyze the mg5 samples to make "angles" files
	if (analyze == 1) {
		printf("Analyzing raw sample root files\n");
		makeAnglesFiles();
	}
	if (onGrid == 0) {
		//make symmetric and weighted asymmetric histogram files
		printf("Making individual histograms for training set files\n");
		gSystem->cd(".");
		if (gSystem->Which(gSystem->pwd(),"templates.root") == 0) {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE NOT FOUND, REMAKING MC HISTO FILES\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			makeHistosFiles(xbins,ybins,zbins,lowmass,highmass,xflow,xfhigh);
		}
		else {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE FOUND, MC HISTO FILES WILL NOT BE REMADE\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		}
		//move the angles and histogram files into the appropriate locations
		//DO NOT DO THIS IF RUNNING ON CONDOR: THE OUTPUT WILL NOT COME BACK
		printf("Moving files into place for later use. . .\n");
		move_stuff(runname);

	}

	if (onGrid == 1) {
		//make symmetric and weighted asymmetric histogram files
		printf("Making individual histograms for training set files\n");
		on_the_grid = 1;
		cout << " Code will run on the grid!" << endl ;		
		gSystem->cd(".");
		if (gSystem->Which(gSystem->pwd(),"templates.root") == 0) {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE NOT FOUND, REMAKING MC HISTO FILES\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			makeHistosFiles(xbins,ybins,zbins,lowmass,highmass,xflow,xfhigh);
		}
		else {
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			printf("TEMPLATE FILE FOUND, MC HISTO FILES WILL NOT BE REMADE\n");
			printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		}


	}


	printf("All done!\n");
	return 0;
}//end main

//take the input file, parse the lines, and build up the linked list of sample file structs
void handle_input(char* filename) {
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
	while (fscanf(ifp,"%s %s %s %d %f %[^\n]s",fp,sn,ifd,&nEG,&c_s,st) != EOF) {
		//make a new sample file struct for each line
		sample_file* newFile = (sample_file*)malloc(sizeof(sample_file));
		//Insert it into the list
		newFile->next = NULL;
		if (HEAD == NULL) {
			HEAD = newFile;
			TAIL = newFile;
		}
		else {
			TAIL->next = newFile;
			TAIL = newFile;
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
		newFile->cross_section = 1.0*c_s;
		newFile->sample_title = (char*)malloc(250*sizeof(char));
		strcpy(newFile->sample_title,st);
	}
	//get the number of sample files and print it out
	printf("Will analyze the following files:\n");
	for (sample_file *iter=HEAD;iter!=NULL;iter=iter->next) {
		nSampleFiles = nSampleFiles+1;
		printf("	%s\n",iter->filepath);
	}
	printf("Total number of files is %d\n",nSampleFiles);
	fclose(ifp);
	
}//end handle_input

//loops over mg5 sample files, runs ttbar.C to make angles files
void makeAnglesFiles() {
	//loop over sample files
	int i = 1;
	for (sample_file *iter=HEAD;iter!=NULL;iter=iter->next) {
		printf("	Analyzing template sample file #%d (%s)\n",i,iter->filepath);
		//Set up automatic file names
		char out[50];
		strcpy(out,"angles_");
		strcat(out,iter->sample_name);
		strcat(out,".root");
		//run ttbar.c to output angles files
		ttbar *t = new ttbar(iter->filepath);
		t->SetOutputFilename(out);
		t->Loop(1);
		free(t);
		i = i+1;
	}//end loop
	printf("Done\n");
}//end makeAnglesFiles

//loops over angles files, runs angles.C to make histogram files
void makeHistosFiles(int xBins, int yBins, int zBins, double lowmass, double highmass, double xflow, double xfhigh) {
	int i = 1;
	//loop over mg5 sample files
	for (sample_file *iter=HEAD;iter!=NULL;iter=iter->next) {
		printf("	Making histograms from sample file #%d (%s)\n",i,iter->filepath);
		//make filenames
		char in[200];
		if (on_the_grid == 1)
		{
			sprintf(in,"%s/angles_%s.root",dir_angles,iter->sample_name); //grid
		}
		else {
		strcpy(in,"../angles_files/angles_");
		strcat(in,iter->sample_name);
		strcat(in,".root");
		}
		char out[200];
		strcpy(out,iter->sample_name);
		strcat(out,"_histos.root");
		//run angles.c on all the files to produce histogram files
		angles *a = new angles();
		if (iter->is_for_distribution==1)
			a->InitializeNames(in,out,iter->sample_name,iter->sample_title,true);
		else
			a->InitializeNames(in,out,iter->sample_name,iter->sample_title,false);
		a->SetBins(xBins,yBins,zBins);
		a->SetMassLimits(lowmass,highmass);
		a->SetxFLimits(xflow,xfhigh);
		// Set PDF member
		a->SetPdfMembers(pdf_member);
		a->Loop(iter->is_for_distribution,iter->nEventsGenerated,iter->cross_section);
		free(a);
		i=i+1;
	}//end loop
	printf("Done\n");
}//end makeHistosFiles

//runs a few shell commands to organize the output of the run.
void move_stuff(char* runname) {
	char dirname[200];
	strcpy(dirname,runname);
	strcat(dirname,"_output");
	char cmd[200];
	strcpy(cmd,"mkdir ");
	strcat(cmd,dirname);
	gSystem->Exec(cmd);
	//gSystem->Exec("mkdir angles_files");
	gSystem->Exec("mkdir histo_files");
	//gSystem->Exec("mv angles_*.root angles_files");
	gSystem->Exec("mv *_histos.root histo_files");
	//strcpy(cmd,"rsync -avh angles_files ");
	//strcat(cmd,dirname);
	//gSystem->Exec(cmd);
	//gSystem->Exec("rm -rf angles_files");
	strcpy(cmd,"rsync -aq histo_files ");
	strcat(cmd,dirname);
	gSystem->Exec(cmd);
	gSystem->Exec("rm -rf histo_files");
}//end move_stuff
