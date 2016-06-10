//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 23 12:33:01 2013 by ROOT version 5.30/06
// from TTree angles/Recreate
// found on file: angles.root
//////////////////////////////////////////////////////////

#ifndef angles_h
#define angles_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class angles {
	public :
	
	virtual void     Loop(int is_for_dist, int neventsgenerated, double crosssection);
	virtual void	 Loop(int is_for_dist);
	virtual void     InitializeNames(const char* in_name, const char* out_name, const char* samplename, 
										const char* samplename_formatted, bool genasymhisto);
	void SetInputTreeName(const char* name);
	void SetBins(int xBins, int yBins, int zBins);
	void SetMassLimits(double lowmass, double highmass);
	void SetxFLimits(double xfLow, double xfHigh);
	void SetPdfMembers(int pdf_member);
	void Set_SF_sys(const char* SF_sys_type);

};

#endif