#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TCanvas.h"
#include <TLegend.h>
#include <TPad.h>
#include "TF1.h"
#include "TF2.h"
#include <TStyle.h>
#include <iostream>
#include <cstring>
#include <string>
#include <Riostream.h>
#include "TFile.h"
#include <TChain.h>
#include <TClonesArray.h>
#include <Selection.h>
#include "TApplication.h"
#include <tuple>
#include <set>



//include C++ library classes
#include <sstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

//include other parts of the code
#include "plotCode_new.h"



using namespace std;




ClassImp(Selection)

//_______________________________________________________default constructor_____
Selection::Selection():TObject()

{
}
//_______________________________________________________ constructor_____
Selection::Selection(TString FileNameTree_in):TObject()

{
}
//________________________________________________________________distruttore_____
Selection::~Selection()	 {
    // destructor
}


//==================================================================
void Selection::selezione(   Int_t _gen_nL,
		   Int_t           _nL,
		    Double_t        _lPt[7], 
		    Double_t        _lEta[7],   
		    Double_t        _lE[7],   
		    Bool_t          _lPOGMedium[7],   
		    Bool_t          _lPOGLoose[7],   
		    UInt_t                 _lFlavor[20],   
		    Int_t        _lCharge[7],   
		    Double_t        _relIso[7],
		    Double_t        _ipPV[7],   
		    Double_t        _ipZPV[7],   
		    Double_t        _3dIP[7],   
		    Double_t        _3dIPsig[7], 
		      
		    Bool_t           _lIsPrompt[7],
		    double index[3],
		    double number_tight[1],
		    double number_prompt[1],
		    int skip_event[1],
			      
		    int effsam,
		  
		    double faxtore[1],
		  
		     Bool_t          _lpassConversionVeto[20],   
			     Double_t        _muNumberInnerHits[20],  
		    Double_t        _eleNumberInnerHitsMissing[20],
		     Double_t        _sumChargedHadronPt03[20]

			    
			     


			      ){
    
     
   
    
    
    
    
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PARAMETERS AND CUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
    const double MVA_cuts_pt15[3] = {0.77, 0.56, 0.48};
    const double MVA_cuts_pt25[3] = {0.52, 0.11, -0.01};
    const double newMVALooseFR[3]= {-0.2, -0.7, -0.8};   
    //const double newMVALooseFR[3]= {-0.02, -0.52, -0.52};   
    //const double newMVALooseFR[3]= {-0.99, -0.99, -0.99};


    const double isolation_loose=0.3;
    const double isolation_tight=0.1;  
    const int number_veto_leptons=3;



    number_tight[0] = 0;
    number_prompt[0] = 0;
    
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< PARAMETERS AND CUTS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  
    for(Int_t et=0;et< 1;et++){
      
      //  cout<<"----------------------- in selction"<<endl;
      skip_event[0] = -1000;
       
      //number of leptons in tree
      if(_nL < number_veto_leptons) skip_event[0] = -100;
      if(_nL < number_veto_leptons) continue;
      // cout<<"3 leptons"<<endl;
      index[0] =  0;  
      index[1] =  1;   
      index[2] =  2;   



      number_prompt[0] = 0;
     
   


  


      faxtore[0] = 1;
      
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VECTORS AND VARIABLES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      unsigned*         ind = new unsigned[_nL];	//new indices of good leptons
      unsigned          lCount = 0;	//Count number of FO leptons that are not taus
      unsigned*         _isFO= new unsigned[_nL];

      unsigned          lCount1 = 0;	//Count number of FO leptons that are not taus
      unsigned*         _isFO1= new unsigned[_nL];
            
      Bool_t            _passedMVA90[_nL];
      double*           conePt = new double[_nL];
            
            
      unsigned*          ordind = new unsigned[_nL];	//Order FO leptons by Pt
      std::set<unsigned> usedLep;
            
      unsigned           tightC = 0;//Count number of T leptons
      unsigned*          _isT= new unsigned[_nL];
            
      unsigned           promptC = 0;
            
      double            low_mass_pt_base[1];


     
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< VECTORS AND VARIABLES <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
     
       skip_event[0] = -90;
      //FO
      for(unsigned l = 0; l < _nL; ++l){
	_isFO[l] = false;
	_isFO1[l] = false;
	//if (_lFlavor[l] == 0 && _lPt[l] < 5 ) continue;
	//if (_lFlavor[l] == 1 && _lPt[l] < 3.5 ) continue;
	if (_lFlavor[l]== 1 && _lPOGLoose[l] )       _isFO1[l] = true; // tracker or global muon --> loose definition from POG
	if (_lFlavor[l]== 0 && _lpassConversionVeto[l] && _eleNumberInnerHitsMissing[l]<=1 && _lPOGLoose[l])    _isFO1[l] = true; //no conversion and number missing hits
	if (_lFlavor[l]== 1 && _lPOGLoose[l]  && _sumChargedHadronPt03[l] < 10 )       _isFO[l] = true; // tracker or global muon --> loose definition from POG
	if (_lFlavor[l]== 0 && _lpassConversionVeto[l] && _eleNumberInnerHitsMissing[l]<=1  && _sumChargedHadronPt03[l] < 10 && _lPOGLoose[l])    _isFO[l] = true; //no conversion and number missing hits
	if(_isFO1[l] && _lFlavor[l] != 2){
	  ++lCount1;
	}
	if(_isFO[l] && _lFlavor[l] != 2){
	  ind[lCount] = l;
	  ++lCount;
	}
      } 
      if(lCount1 < number_veto_leptons) continue;
      skip_event[0] = -80;
      //if(lCount1 < number_veto_leptons) skip_event[0] = -80;
      //if (lCount != number_veto_leptons ) skip_event[0] = -1;
      if(lCount < number_veto_leptons) continue;
      skip_event[0] = -75;


  
      // if (lCount != number_veto_leptons ) continue;

      
      for(unsigned l = 0; l < lCount; ++l){
	conePt[ind[l]] =0.;
	conePt[ind[l]] = _lPt[ind[l]];
      }

    
      
      for(unsigned k =0; k < lCount; ++k){
	//unsigned maxI = 999;
	double maxConePt = 0;
	for(unsigned l = 0; l < lCount; ++l){
	  if(usedLep.find(ind[l]) == usedLep.end()){
	    if(conePt[ind[l]] > maxConePt){
	      maxConePt = conePt[ind[l]];
	      ordind[k] = ind[l];
	    }
	  }
	}
	usedLep.insert(ordind[k]);
      }
      for(unsigned i = 0; i < lCount; ++i){
	ind[i] = ordind[i];
      }

      int pt_treshold = 100;
      for(unsigned l = 0; l < lCount; ++l){
	if((_lFlavor[ind[l]] == 0 && conePt[ind[l]] < 5)  ) pt_treshold = -1;
	if((_lFlavor[ind[l]] == 1 && conePt[ind[l]] < 3.5)  ) pt_treshold = -1;
      }
      if(pt_treshold == -1)       continue;
      skip_event[0] = -70;
  
            
      

      _isT[ind[0]]=true;
      _isT[ind[1]]=true;
      _isT[ind[2]]=true;
      

      if (_lPt[ind[0]] < 20 || _relIso[ind[0]] > 0.1 || TMath::Abs(_ipPV[ind[0]]) > 0.05 || TMath::Abs(_ipZPV[ind[0]]) > 0.1 || fabs(_3dIPsig [ind[0]]) > 4  ) continue;
      skip_event[0] = -65; 

          
      


      /*  if ((TMath::Abs(_ipPV[ind[1]]) < 0.05 && TMath::Abs(_ipZPV[ind[1]]) < 0.1 && TMath::Abs(_3dIPsig[ind[1]] < 4) ) || (TMath::Abs(_ipPV[ind[2]]) < 0.05 && TMath::Abs(_ipZPV[ind[2]]) < 0.1 && TMath::Abs(_3dIPsig[ind[2]] < 4) ) )  skip_event[0] = -1; 
       if ((TMath::Abs(_ipPV[ind[1]]) < 0.05 && TMath::Abs(_ipZPV[ind[1]]) < 0.1 && TMath::Abs(_3dIPsig[ind[1]] < 4))  || (TMath::Abs(_ipPV[ind[2]]) < 0.05 && TMath::Abs(_ipZPV[ind[2]]) < 0.1 && TMath::Abs(_3dIPsig[ind[2]] < 4) ) )  continue;

      */

      if ( (TMath::Abs(_ipPV[ind[1]]) < 0.05 && TMath::Abs(_ipZPV[ind[1]]) < 0.1 && _3dIPsig[ind[1]] < 4)  )  continue;
      skip_event[0] = -60; 
      if ( (TMath::Abs(_ipPV[ind[2]]) < 0.05 && TMath::Abs(_ipZPV[ind[2]]) < 0.1 && _3dIPsig[ind[2]] < 4)  )  continue;
      skip_event[0] = -55; 


      if( ( !_lPOGLoose[ind[1]] || !_lPOGMedium[ind[1]] || _relIso[ind[1]] > 0.6)  || ( !_lPOGLoose[ind[2]] || !_lPOGMedium[ind[2]] || _relIso[ind[2]] > 0.6)) continue;
	 skip_event[0] = 20;




      if (_isT[ind[0]]  && _isT[ind[1]]  && _isT[ind[2]]  ) tightC=3;
      number_tight[0] = tightC;
      if ( tightC==3 ) number_tight[0]=3;
   

      // only prompt for MC 
      // cout<<"prima 3 prompts"<<endl;

      if (_lIsPrompt[ind[0]]  && _lIsPrompt[ind[1]]  && _lIsPrompt[ind[2]]  ) promptC=3;

      /*
      if (promptC != 3){
	cout<<"lCount:   "<<lCount<<endl;
	cout<<_lIsPrompt[ind[0]]  <<"  -   "<< _lIsPrompt[ind[1]] <<"  -   "<<_lIsPrompt[ind[2]] <<endl;
	cout<<"pt:  "<<_lPt[ind[0]]  <<"  -   "<< _lPt[ind[1]] <<"  -   "<<_lPt[ind[2]] <<endl;
	cout<<"eta:  "<<_lEta[ind[0]]  <<"  -   "<< _lEta[ind[1]] <<"  -   "<<_lEta[ind[2]] <<endl;
	cout<<"FO:  "<<_lHNFO[ind[0]]<<"  -  "<<_lHNFO[ind[1]]<<"  -  "<<_lHNFO[ind[2]]<<endl;
	 for(unsigned l = 0; l < _gen_nL; ++l){
	   cout<<l<<")   pt_gen: "<<_gen_lPt[l]<<"   eta_gen: "<<_gen_lEta[l]<<"    flov: "<< _gen_lFlavor[l]<<"   prompt: "<< _gen_lIsPrompt[l]<<endl;
	 }
     }
      */


      // if ( promptC != 3) skip_event[0] = -1;			   
      // if ( promptC != 3) continue;
      //	cout<<"3 prompts"<<endl;


 
   
            

    
      index[0] =  ind[0];  
      index[1] =  ind[1];   
      index[2] =  ind[2];   

      if (promptC == 0)number_prompt[0] = 0;
      if (promptC == 1)number_prompt[0] = 1;
      if (promptC == 2)number_prompt[0] = 2;
      if (promptC == 3)number_prompt[0] = 3;
      
   


  


      faxtore[0] = 1;
      
      /*faxtore[0] = 1;
      if(tightFail ){
	if (effsam  == 0 )faxtore[0] *= -1;
	if (effsam  != 0 )faxtore[0] = 1 ; 
	for(unsigned l = 0; l < lCount; ++l){	
	 
	  if(!_isT[ind[l]] && _isFO[ind[l]]){
	    double fr = FR_factor (*&fakeRate_mu, *&fakeRate_e ,_lEta[ind[l]], _flavors[ind[l]], conePt[ind[l]]);
	    faxtore[0] *= -fr/(1-fr);
	  }
	}   
      }
      if (!tightFail) faxtore[0] = 1;*/

        
    }    
    
    
}// end analisi






//==================================================================
double Selection::maximum(double a, double b){
    double massimo=0;
    massimo = a;
    if (massimo < b) massimo = b;
    return massimo+1;
}
//==================================================================
double Selection::find_eta(double b){
    double etaa=-1;
    if(b < 0.8 ) etaa = 0;
    else if(b < 1.479 ) etaa = 1;
    else etaa = 2;
    return etaa;
}
//==================================================================
void Selection::check_low_mass_pt(double check[1], double pt_cone_1,double pt_cone_2,double pt_cone_3,double flavor_3){ 
  check[0] = 1;
  if(pt_cone_1 < 15 || pt_cone_2 < 10 || (pt_cone_3 < 5 && flavor_3 == 1) || (pt_cone_3 < 10 && flavor_3 == 0)) check[0] = 0;
}   

//==================================================================
void Selection::from_TGraph_to_TH1D (TGraphAsymmErrors *graph, TH1D *histo, int number_point){

  const int numero= number_point;

  double x_graph[numero];
  double y_graph[numero];  
  for (int i =0; i <number_point; i ++){
    x_graph[i]=0;
    y_graph[i]=0;
  }
  for (int i =0; i <number_point; i ++){
    graph -> GetPoint(i, x_graph[i], y_graph[i]);
    histo->SetBinContent (i+1, x_graph[i],  y_graph[i]);

    //cout<<i<<") "<<y_graph[i]<<"  "<< histo->GetBinContent (i+1)<<endl;

  }
}


//==================================================================
double Selection::FR_factor(TGraphAsymmErrors *fakeRate_mu[3],
                              TGraphAsymmErrors *fakeRate_e[3],
                              double eta,
                              double flavors,
                              double lptcone
                              ){
    
    
    eta = fabs(eta);
    
    
    TH1D *fakeRate_mu_histo[3];
    TH1D *fakeRate_e_histo[3];
    Double_t newBins_mu1[7] = {5,10, 15, 25, 35, 50, 70};
    Double_t newBins_e1[6] = {10, 15, 25, 35, 50, 70};
    fakeRate_mu_histo[0]= new TH1D("fake_rate_mu_histo_eta1","",6,newBins_mu1);
    fakeRate_mu_histo[1]= new TH1D("fake_rate_mu_histo_eta2","",6,newBins_mu1);
    fakeRate_mu_histo[2]= new TH1D("fake_rate_mu_histo_eta3","",6,newBins_mu1);
    fakeRate_e_histo[0]= new TH1D("fake_rate_e_histo_eta1","",5,newBins_e1);
    fakeRate_e_histo[1]= new TH1D("fake_rate_e_histo_eta2","",5,newBins_e1);
    fakeRate_e_histo[2]= new TH1D("fake_rate_e_histo_eta3","",5,newBins_e1);
    
    for (int i=0; i< 3; i++){
        from_TGraph_to_TH1D(*&fakeRate_mu[i],*&fakeRate_mu_histo[i],6);
        from_TGraph_to_TH1D(*&fakeRate_e[i],*&fakeRate_e_histo[i],5);
    }
    
    //double momentum = part.Pt() * maximum( 1, iso - 0.1);
    double momentum = lptcone;
    
    double factor=0;
    double factore=0;
    if (flavors == 0){
        if (momentum < 49){
            if (eta < 0.8){
                factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(momentum));
            }//eta1
            else if (eta > 0.8 && eta<1.479){
                factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(momentum));
            }//eta1
            else {
                factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(momentum));
            }//eta1
        }// <70
        else {
            if (eta < 0.8){
                factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(45));
            }//eta1
            else if (eta > 0.8 && eta<1.479){
                factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(45));
            }//eta1
            else {
                factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(45));
            }//eta1
        }
    }//e
    
    if (flavors == 1){
        if (momentum < 49){
            if (eta < 0.8){
                factore = fakeRate_mu_histo[0]->GetBinContent(fakeRate_mu_histo[0]->FindBin(momentum));
            }//eta1
            else if (eta > 0.8 && eta<1.479){
                factore = fakeRate_mu_histo[1]->GetBinContent(fakeRate_mu_histo[1]->FindBin(momentum));
            }//eta1
            else {
                factore = fakeRate_mu_histo[2]->GetBinContent(fakeRate_mu_histo[2]->FindBin(momentum));
            }//eta1
        }// <70
        else {
            if (eta < 0.8){
                factore = fakeRate_mu_histo[0]->GetBinContent(fakeRate_mu_histo[0]->FindBin(45));
            }//eta1
            else if (eta > 0.8 && eta<1.479){
                factore = fakeRate_mu_histo[1]->GetBinContent(fakeRate_mu_histo[1]->FindBin(45));
            }//eta1
            else {
                factore = fakeRate_mu_histo[2]->GetBinContent(fakeRate_mu_histo[2]->FindBin(45));
            }//eta1
        }
    }//e
    
    delete fakeRate_mu_histo[0];
    delete fakeRate_mu_histo[1];
    delete fakeRate_mu_histo[2];
    delete fakeRate_e_histo[0];
    delete fakeRate_e_histo[1];
    delete fakeRate_e_histo[2];
    
    
    
    return factore;
    
}
