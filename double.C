#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TPad.h>
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
 
using namespace std;
void dbt(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

  TH1F* h_delR=new TH1F("","",30,0,6);

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);

    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
    Float_t*  fatjetSDmass = data.GetPtrFloat("FATjetSDmass");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    Int_t*   nSubSoftDropJet = data.GetPtrInt("FATnSubSDJet");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    Int_t*   FATnSubSDJet   = data.GetPtrInt("FATnSubSDJet");
    vector<float>* FATsubjetSDCSV       = data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);

    
    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    
    vector<int> fatjet;
    for(int ij=0; ij<nFJets; ij++)
      {
	
	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
	if(thisJet->Pt()<300)continue;
	if(fabs(thisJet->Eta())>2.5)continue;
	if(!passFatJetLooseID[ij])continue;
	if( FATnSubSDJet[ij] != 2 ) continue;
        if( FATsubjetSDCSV[ij][0] < 0.605 || FATsubjetSDCSV[ij][1] < 0.605 ) continue;  	

	fatjet.push_back(ij);

      }
    if(fatjet.size()<2)continue;
    TLorentzVector *leadfatjet = (TLorentzVector*)fatjetP4->At(fatjet[0]);
    
    vector<int> addjet;
    Float_t del_R;
    for(int ad=0; ad<nAJets; ad++)
      {
	TLorentzVector* theseJet = (TLorentzVector*)addjetP4->At(ad);
	del_R = theseJet->DeltaR(*leadfatjet);
      }
   
    h_delR->Fill(del_R);
  }

  h_delR->Draw();
}
