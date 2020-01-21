#include "TRint.h"
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TROOT.h>
#include <TChain.h>  
#include <TFile.h>
#include "TObject.h"
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <TVector3.h>
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TH1I.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveStats.h"
#include "TVector3.h"
#include <stdlib.h>

using namespace std;
 
int main(int __argc,char *__argv[]){

	char *outFileName = "Event_Analysis.root";
	extern int optind;   
	TROOT troot();  
	TFile outFile(outFileName,"recreate");  

	// Declares parameters.
	Int_t nentries, j;

       Float_t     p1px, p1py, p1pz, p1pe;
       Float_t     p2px, p2py, p2pz, p2pE;
       Float_t     pbpx, pbpy, pbpz, pbpE;
       Float_t     pIpx, pIpy, pIpz, pIpE;

	Float_t     boost;
	Float_t     p1bcos, p2bcos, pbbcos;
	Float_t     p1bphi, p2bphi, pbbphi;

	Float_t xmm, xmm2;

	TLorentzVector P4p1, P4p2, P4pb, P4miss, P4target, P4pho, P4inc, P4all, P4prots, P4prMM;

	const Float_t xmpr = 0.938272; // Proton mass.

	TH1F *Hegam = new TH1F("Hegam","Photon energy",200,0.0,6.0);
	TH1F *HMMsq = new TH1F("MMsq","MM^{2}(p p)",500,-2.0,2.0);
	TH1F *HMMsqall = new TH1F("MMsqa","MM^{2}(p p pbar)",500,-0.25,0.25);
	TH1F *HMM = new TH1F("2MM","MM(p p pbar)",500,-2.0,2.0);

	TH1F *p1gcos = new TH1F("P1g_CosTheta","CosTheta(p)",20,-1.2,1.2);
	TH1F *p2gcos = new TH1F("P2g_CosTheta","CosTheta(p)",20,-1.2,1.2);
	TH1F *pbgcos = new TH1F("Pbg_CosTheta","CosTheta(pbar)",20,-1.2,1.2);

	TH1F *p1gphi = new TH1F("P1g_Phi","Phi(p)",20,-3.2,3.2);
	TH1F *p2gphi = new TH1F("P2g_Phi","Phi(p)",20,-3.2,3.2);
	TH1F *pbgphi = new TH1F("Pbg_Phi","Phi(pbar)",20,-3.2,3.2);

	TH2F *p1vsegam = new TH2F("P1_EgamvsCosTheta","Egam(GeV) vs. CosTheta",20,1,-1,20,4.5,5.5);
	TH2F *p2vsegam = new TH2F("P2_EgamvsCosTheta","Egam(GeV) vs. CosTheta",20,1,-1,20,4.5,5.5);
	TH2F *pbvsegam = new TH2F("Pb_EgamvsCosTheta","Egam(GeV) vs. CosTheta",20,1,-1,20,4.5,5.5);


/* ***************************MAIN PROGRAM*************************** */

	/* Loops over data files. */
	for(int n_arg = optind; n_arg < __argc; n_arg++){
   
	TFile inFile(__argv[n_arg]);
		
	/* Opens each input file. */        
	if(TTree *pg = (TTree*)inFile.Get("T")){
			
	pg->SetBranchAddress("p1px",&p1px);
	pg->SetBranchAddress("p1py",&p1py);
	pg->SetBranchAddress("p1pz",&p1pz);
	pg->SetBranchAddress("p1pe",&p1pe);
	pg->SetBranchAddress("p2px",&p2px);
	pg->SetBranchAddress("p2py",&p2py);
	pg->SetBranchAddress("p2pz",&p2pz);
	pg->SetBranchAddress("p2pE",&p2pE);
	pg->SetBranchAddress("pbpx",&pbpx);
	pg->SetBranchAddress("pbpy",&pbpy);
	pg->SetBranchAddress("pbpz",&pbpz);
	pg->SetBranchAddress("pbpE",&pbpE);
	pg->SetBranchAddress("pIpx",&pIpx);
	pg->SetBranchAddress("pIpy",&pIpy);
	pg->SetBranchAddress("pIpz",&pIpz);
	pg->SetBranchAddress("pIpE",&pIpE);

	nentries = (Int_t)pg->GetEntries();
			
	/* Loops over events. */
	for (j=0;j<nentries;j++){ 

	pg->GetEntry(j);

	boost = -pIpE/(pIpE + xmpr);
		
	// Set the target and photon 4-Vectors.
	P4target.SetPxPyPzE(0.0,0.0,0.0,xmpr);		
	P4pho.SetPxPyPzE(pIpx,pIpy,pIpz,pIpE);
	P4inc =  P4target + P4pho;
       
	// Sets the 4-Vectors of the particles. 
	P4p1.SetPxPyPzE(p1px,p1py,p1pz,sqrt(p1px*p1px+p1py*p1py+p1pz*p1pz+xmpr*xmpr));
	P4p2.SetPxPyPzE(p2px,p2py,p2pz,sqrt(p2px*p2px+p2py*p2py+p2pz*p2pz+xmpr*xmpr));
	P4pb.SetPxPyPzE(pbpx,pbpy,pbpz,sqrt(pbpx*pbpx+pbpy*pbpy+pbpz*pbpz+xmpr*xmpr));


	// Calculate various 4-vectors.
	P4all = P4p1 + P4p2 + P4pb;
	P4prots = P4p1 + P4p2;
	P4miss = P4inc - P4all;
	P4prMM = P4inc - P4prots;


	//Boost the vectors
	P4p1.Boost(0,0,boost);
	P4p2.Boost(0,0,boost);
	P4pb.Boost(0,0,boost);

	p1bcos = P4p1.Vect().CosTheta();
	p2bcos = P4p2.Vect().CosTheta();
	pbbcos = P4pb.Vect().CosTheta();
	
	p1bphi = P4p1.Vect().Phi();
	p2bphi = P4p2.Vect().Phi();
	pbbphi = P4pb.Vect().Phi();

	// Get masses.
	xmm  = P4miss.M();
	xmm2 = P4prMM.M();


	//Fill Histograms
	Hegam->Fill(pIpE);

	p1gcos->Fill(p1bcos);
	p2gcos->Fill(p2bcos);
	pbgcos->Fill(pbbcos);

	p1gphi->Fill(p1bphi);
	p2gphi->Fill(p2bphi);
	pbgphi->Fill(pbbphi);

	p1vsegam->Fill(p1bcos,pIpE);
	p2vsegam->Fill(p2bcos,pIpE);
	pbvsegam->Fill(pbbcos,pIpE);


	    HMM->Fill(xmm);
	    HMMsq->Fill(P4prMM.M2());
	    HMMsqall->Fill(P4miss.M2());


       

	
	}// entries
	}//if TTree *nbar = (TTree*)
		
	else {
	  cout << "File has no TTree named pipr or does not exist!!" << endl;
	  cout << __argv[n_arg] << endl;
	}  
	
	}// files

 outFile.Write(); // Write to the output file.
 outFile.Close(); // Close the output file.

}
