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
#include "TLine.h"
#include "TProfile.h"
#include "TAttFill.h"
#include <iostream>
#include <cstring>
#include <string>
#include <TGraphErrors.h>
#include <Riostream.h>
#include "TFile.h"
#include <TChain.h>
#include <TClonesArray.h>
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>
#include <THStack.h>
#include <TPaveText.h>
#include <Analysis_mc.h>
#include "TApplication.h"
#include "TColor.h"
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
#include "tdrstyle.h"
#include "plotCode_new.h"
//#include "Selection.h"




using namespace std;




static Double_t pigreco= TMath::ACos(-1);
ClassImp(Analysis_mc)

//_______________________________________________________default constructor_____
Analysis_mc::Analysis_mc():TObject()

{
}
//_______________________________________________________ constructor_____
Analysis_mc::Analysis_mc(int selezione, string FileNameTree_in):TObject()

{
}
//________________________________________________________________distruttore_____
Analysis_mc::~Analysis_mc()	 {
  // destructor
}

//_______________________________________________________ constructor_____
void Analysis_mc::printProgress(double progress){
  const unsigned barWidth = 100;
  std::cout << "[";
  unsigned pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << unsigned(progress * 100.0) << " %\r" << std::flush;
}





//==================================================================
void Analysis_mc::analisi(int selezione, int num_histo_kin
                          ){
    
  cout<<"in analisi"<<endl;
  cout<<"---------------------------"<<endl;   
  setTDRStyle();
    
    
  
  const double glugluToZZkFactor = 1;
  const double WZSF = 1;
  const double ZZSF = 1;
  const double XgammaSF = 1;
  const double low_coupling= 0.0001;
  const double low_coupling_2= 0.00001;
  const double high_coupling = 0.001;
  const double coupling_factor_low= low_coupling / 0.0001;
  const double coupling_factor_low_2= low_coupling_2 / 0.0001;
 
  const int nSamples= 50;
  const int nSamples_eff = 17;
  const int nSamples_signal=10;
    

  const TString fileList[nSamples] = {  "ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",

					"gev1.root", "gev2.root", "gev3.root", "gev4.root", "gev5.root", "gev5_5.root", "gev6.root", "gev7.root", "gev8.root", "gev9.root",

					"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",

					"TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", "TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", "TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root ",	"ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root" ,"ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root" ,"ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root" ,"ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root" ,"ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",
					
					"WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root",


					"ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root", "VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root",

					"WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root","ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root","WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "WZG_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root","WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "WGGJets_TuneCUETP8M1_13TeV_madgraphMLM_pythia8.root",
					
					"GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root","WWTo2L2Nu_13TeV-powheg.root", "WWTo2L2Nu_DoubleScattering_13TeV-pythia8.root","ZZTo4L_13TeV_powheg_pythia8.root ", "WZTo3LNu_mllmin01_13TeV-powheg-pythia8_ext1.root",


					"WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root","ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root", 

					"TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8.root", "TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8.root", "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root "    , "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8___1.root", "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","tZq_ll_4f_13TeV-amcatnlo-pythia8.root"


				 
  };

       
  const TString names[nSamples]   =  {  "total",      
					"trilMaj1", "trilMaj2", "trilMaj3","trilMaj4","trilMaj5", "trilMaj5_5", "trilMaj6","trilMaj7", "trilMaj8","trilMaj9",    
					"DY",    "DY",   
					"TTbar","TTbar","TTbar","TTbar", "TTbar","TTbar", "TTbar","TTbar",
					"WJets",
					"multiboson","multiboson",
					"multiboson", "multiboson", "multiboson","multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson",
					"X+#gamma","X+#gamma",      
					"TT/T + X","TT/T + X","TT/T + X", "TT/T + X","TT/T + X","TT/T + X","TT/T + X"
				
					
  };
  const TString eff_names[nSamples_eff +1 ] = { "total",      
						"trilMaj1", "trilMaj2", "trilMaj3","trilMaj4","trilMaj5", "trilMaj5_5", "trilMaj6","trilMaj7", "trilMaj8","trilMaj9",    
						"DY",  
						"t#bar{t}",
						"WJets",
						"multiboson", 
						"X+#gamma",    
						"TT/T + X",		
						"no prompt"};

  // 9 0.05173

  const double xSections[nSamples]= {0,        
				     0.5201 * coupling_factor_low,0.5273 * coupling_factor_low, 0.5217 * coupling_factor_low,0.5216 * coupling_factor_low,0.05190 * coupling_factor_low_2, 0.05220* coupling_factor_low_2, 0.05226* coupling_factor_low_2, 0.05226* coupling_factor_low_2, 0.05193* coupling_factor_low_2, 0.05173* coupling_factor_low_2,
        
				     18610, 1921.8*3,
				     87.315, 182.175, 182.75,3.36 , 26.38 ,  44.33 , 35.85, 35.85 , 
				     61526.7,

				     0.752,0.001034,

				     0.2086, 0.1651,  0.01398,0.05565,0.04123 , 0.2147 , 1.711, 0.00319,0.00319,0.00319,0.00159,0.00159,0.00159,12.178, 0.1729, 1.256, 58.59*0.652,


				     405.271,123.9,
				     
				     2.967, 0.01731, 3.697, 0.2043, 0.2529,0.4719,0.0758
				     
				    

				    
  };
    
  double luminosity = 35.867;
    
  TFile *hfile[nSamples];
  TTree *inputTree[nSamples];
  //TApplication* rootapp = new TApplication("example",&argc, argv);
    
    
    
  // Declaration of leaf types
  ULong64_t       _runNb;
  ULong64_t       _lumiBlock;
  ULong64_t       _eventNb;
  UChar_t         _nVertex;
  Double_t        _weight;
  Double_t        _lheHTIncoming;
  Double_t        _ctauHN;
  UChar_t         _nLheWeights;
  _nLheWeights = 110;
  Double_t        _lheWeight[110];   //[_nLheWeights]
  Float_t         _nTrueInt;
  UChar_t         _ttgEventType;
  UChar_t         _zgEventType;
  Double_t        _gen_met;
  Double_t        _gen_metPhi;
  UChar_t         _gen_nPh;
  _gen_nPh = 7;
  Double_t        _gen_phPt[7];   //[_gen_nPh]
  Double_t        _gen_phEta[7];   //[_gen_nPh]
  Double_t        _gen_phPhi[7];   //[_gen_nPh]
  Double_t        _gen_phE[7];   //[_gen_nPh]
  Int_t           _gen_phMomPdg[7];   //[_gen_nPh]
  Bool_t          _gen_phIsPrompt[7];   //[_gen_nPh]
  Double_t        _gen_phMinDeltaR[7];   //[_gen_nPh]
  Bool_t          _gen_phPassParentage[7];   //[_gen_nPh]
  UChar_t         _gen_nL;
  _gen_nL = 12;
  Double_t        _gen_lPt[12];   //[_gen_nL]
  Double_t        _gen_lEta[12];   //[_gen_nL]
  Double_t        _gen_lPhi[12];   //[_gen_nL]
  Double_t        _gen_lE[12];   //[_gen_nL]
  //UInt_t          _gen_lFlavor[12];   //[_gen_nL]
  Int_t           _gen_lCharge[12];   //[_gen_nL]
  Int_t           _gen_lMomPdg[12];   //[_gen_nL]
  Bool_t          _gen_lIsPrompt[12];   //[_gen_nL]
  Bool_t          _passHN_1l;
  Bool_t          _HLT_Ele27_WPTight_Gsf;
  Int_t           _HLT_Ele27_WPTight_Gsf_prescale;
  Bool_t          _HLT_IsoMu24;
  Int_t           _HLT_IsoMu24_prescale;
  Bool_t          _HLT_IsoTkMu24;
  Int_t           _HLT_IsoTkMu24_prescale;
  Bool_t          _passHN_eee;
  Bool_t          _HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;
  Int_t           _HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale;
  Bool_t          _HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Int_t           _HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
  Bool_t          _passHN_eem;
  Bool_t          _HLT_Mu8_DiEle12_CaloIdL_TrackIdL;
  Int_t           _HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale;
  Bool_t          _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
  Int_t           _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
  Bool_t          _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;
  Int_t           _HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
  Bool_t          _HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _passHN_emm;
  Bool_t          _HLT_DiMu9_Ele9_CaloIdL_TrackIdL;
  Int_t           _HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  Int_t           _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  Int_t           _HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale;
  Bool_t          _HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  Int_t           _HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  Int_t           _HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  Int_t           _HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale;
  Bool_t          _HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  Int_t           _HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale;
  Bool_t          _passHN_mmm;
  Bool_t          _HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx;
  Int_t           _HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_prescale;
  Bool_t          _HLT_TripleMu_12_10_5;
  Int_t           _HLT_TripleMu_12_10_5_prescale;
  Bool_t          _passMET;
  Bool_t          _HLT_MET300;
  Int_t           _HLT_MET300_prescale;
  Bool_t          _HLT_HT350_MET100;
  Int_t           _HLT_HT350_MET100_prescale;
  Bool_t          _HLT_AllMET300;
  Int_t           _HLT_AllMET300_prescale;
  Bool_t          _HLT_AllMET170;
  Int_t           _HLT_AllMET170_prescale;
  Bool_t          _HLT_jet;
  Int_t           _HLT_jet_prescale;
  Bool_t          _HLT_dijet;
  Int_t           _HLT_dijet_prescale;
  Bool_t          _HLT_MET170_BeamHaloCleaned;
  Int_t           _HLT_MET170_BeamHaloCleaned_prescale;
  Bool_t          _HLT_MET170_NotCleaned;
  Int_t           _HLT_MET170_NotCleaned_prescale;
  Bool_t          _HLT_HT800;
  Int_t           _HLT_HT800_prescale;
  Bool_t          _HLT_HT900;
  Int_t           _HLT_HT900_prescale;
  Bool_t          _HLT_dijet55met110;
  Int_t           _HLT_dijet55met110_prescale;
  Bool_t          _HLT_dijet70met120;
  Int_t           _HLT_dijet70met120_prescale;
  Bool_t          _HLT_HT600;
  Int_t           _HLT_HT600_prescale;
  Bool_t          _HLT_HT475;
  Int_t           _HLT_HT475_prescale;
  Bool_t          _HLT_HT350;
  Int_t           _HLT_HT350_prescale;
  Bool_t          _passMETFilters;
  Bool_t          _Flag_HBHENoiseFilter;
  Bool_t          _Flag_HBHENoiseIsoFilter;
  Bool_t          _Flag_EcalDeadCellTriggerPrimitiveFilter;
  Bool_t          _Flag_goodVertices;
  Bool_t          _Flag_eeBadScFilter;
  Bool_t          _Flag_globalTightHalo2016Filter;
  Bool_t          _flag_badPFMuonFilter;
  Bool_t          _flag_badChCandFilter;
  Bool_t          _passTTG_e;
  Bool_t          _HLT_Ele105_CaloIdVT_GsfTrkIdT;
  Int_t           _HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale;
  Bool_t          _HLT_Ele115_CaloIdVT_GsfTrkIdT;
  Int_t           _HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale;
  Bool_t          _passTTG_ee;
  Bool_t          _HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  Int_t           _HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;
  Bool_t          _HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
  Int_t           _HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale;
  Bool_t          _HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW;
  Int_t           _HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale;
  Bool_t          _passTTG_em;
  Bool_t          _HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;
  Int_t           _HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale;
  Bool_t          _HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL;
  Int_t           _HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale;
  Bool_t          _passTTG_m;
  Bool_t          _HLT_Mu50;
  Int_t           _HLT_Mu50_prescale;
  Bool_t          _HLT_TkMu50;
  Int_t           _HLT_TkMu50_prescale;
  Bool_t          _HLT_Mu45_eta2p1;
  Int_t           _HLT_Mu45_eta2p1_prescale;
  Bool_t          _passTTG_mm;
  Bool_t          _HLT_Mu30_TkMu11;
  Int_t           _HLT_Mu30_TkMu11_prescale;
  UChar_t         _nL;
  _nL = 20;
  UChar_t         _nMu;
  _nMu = 20;
  UChar_t         _nEle;
  _nEle = 20;
  UChar_t         _nLight;
  _nLight = 20;
  UChar_t         _nTau;
  UChar_t         _nVFit;
  _nVFit  = 25;
  UChar_t         _nGoodLeading;
  Double_t        _lIndex[20];   //[_nL]
  Double_t        _vertices[25][12];
  Double_t        _lPt[20];   //[_nL]
  Double_t        _lEta[20];   //[_nL]
  Double_t        _lEtaSC[20];   //[_nLight]
  Double_t        _lPhi[20];   //[_nL]
  Double_t        _lE[20];   //[_nL]
  UInt_t          _lFlavor[20];   //[_nL]
  Int_t           _lCharge[20];   //[_nL]
  Double_t        _dxy[20];   //[_nL]
  Double_t        _dz[20];   //[_nL]
  Double_t        _3dIP[20];   //[_nL]
  Double_t        _3dIPSig[20];   //[_nL]
  Double_t        _2dIP[20];   //[_nL]
  Double_t        _2dIPSig[20];   //[_nL]
  Float_t         _lElectronMva[20];   //[_nLight]
  Bool_t          _lElectronPassEmu[20];   //[_nLight]
  Bool_t          _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[20];   //[_nL]
  Bool_t          _lPOGVeto[20];   //[_nL]
  Bool_t          _lPOGLoose[20];   //[_nL]
  Bool_t          _lPOGMedium[20];   //[_nL]
  Bool_t          _lPOGTight[20];   //[_nL]
  Bool_t          _lpassConversionVeto[20];   //[_nL]
  Double_t        _muNumberInnerHits[20];   //[_nL]
  Double_t        _eleNumberInnerHitsMissing[20];   //[_nL]
  Double_t        _relIso[20];   //[_nLight]
  Double_t        _puCorr[20];   //[_nL]
  Double_t        _absIso03[20];   //[_nL]
  Double_t        _absIso04[20];   //[_nL]
  Double_t        _sumNeutralHadronEt04[20];   //[_nL]
  Double_t        _sumChargedHadronPt04[20];   //[_nL]
  Double_t        _sumPhotonEt04[20];   //[_nL]
  Double_t        _sumNeutralHadronEt03[20];   //[_nL]
  Double_t        _sumChargedHadronPt03[20];   //[_nL]
  Double_t        _sumPhotonEt03[20];   //[_nL]
  Double_t        _trackIso[20];   //[_nL]
  Double_t        _ecalIso[20];   //[_nL]
  Double_t        _hcalIso[20];   //[_nL]
  Double_t        _deltaBIso[20];   //[_nL]
  Double_t        _ecalPFClusterIso[20];   //[_nL]
  Double_t        _hcalPFClusterIso[20];   //[_nL]
  Double_t        _ptRel[20];   //[_nLight]
  Double_t        _ptRatio[20];   //[_nLight]
  Double_t        _closestJetCsv[20];   //[_nLight]
  UInt_t          _selectedTrackMult[20];   //[_nLight]
  Double_t        _muonSegComp[10];   //[_nMu]
  Bool_t          _tauMuonVeto[20];   //[_nL]
  Bool_t          _tauEleVeto[20];   //[_nL]
  Bool_t          _decayModeFindingNew[20];   //[_nL]
  Bool_t          _tauVLooseMvaNew[20];   //[_nL]
  Bool_t          _tauLooseMvaNew[20];   //[_nL]
  Bool_t          _tauMediumMvaNew[20];   //[_nL]
  Bool_t          _tauTightMvaNew[20];   //[_nL]
  Bool_t          _tauVTightMvaNew[20];   //[_nL]
  Bool_t          _tauVTightMvaOld[20];   //[_nL]
  Bool_t          _lIsPrompt[20];   //[_nL]
  Int_t           _lMatchPdgId[20];   //[_nL]
  UChar_t         _nPh;
  _nPh = 20;
  Double_t        _phPt[10];   //[_nPh]
  Double_t        _phEta[10];   //[_nPh]
  Double_t        _phEtaSC[10];   //[_nPh]
  Double_t        _phPhi[10];   //[_nPh]
  Double_t        _phE[10];   //[_nPh]
  Bool_t          _phCutBasedLoose[10];   //[_nPh]
  Bool_t          _phCutBasedMedium[10];   //[_nPh]
  Bool_t          _phCutBasedTight[10];   //[_nPh]
  Double_t        _phMva[10];   //[_nPh]
  Double_t        _phRandomConeChargedIsolation[10];   //[_nPh]
  Double_t        _phChargedIsolation[10];   //[_nPh]
  Double_t        _phNeutralHadronIsolation[10];   //[_nPh]
  Double_t        _phPhotonIsolation[10];   //[_nPh]
  Double_t        _phSigmaIetaIeta[10];   //[_nPh]
  Double_t        _phSigmaIetaIphi[10];   //[_nPh]
  Double_t        _phHadronicOverEm[10];   //[_nPh]
  Bool_t          _phPassElectronVeto[10];   //[_nPh]
  Bool_t          _phHasPixelSeed[10];   //[_nPh]
  Bool_t          _phIsPrompt[10];   //[_nPh]
  Int_t           _phMatchMCPhotonAN15165[10];   //[_nPh]
  Int_t           _phMatchMCLeptonAN15165[10];   //[_nPh]
  Int_t           _phMatchPdgId[10];   //[_nPh]
  UChar_t         _nJets;
  _nJets = 20;
  Double_t        _jetPt[20];   //[_nJets]
  Double_t        _jetPt_JECUp[20];   //[_nJets]
  Double_t        _jetPt_JECDown[20];   //[_nJets]
  Double_t        _jetPt_JERUp[20];   //[_nJets]
  Double_t        _jetPt_JERDown[20];   //[_nJets]
  Double_t        _jetEta[20];   //[_nJets]
  Double_t        _jetPhi[20];   //[_nJets]
  Double_t        _jetE[20];   //[_nJets]
  Double_t        _jetCsvV2[20];   //[_nJets]
  Double_t        _jetDeepCsv_udsg[20];   //[_nJets]
  Double_t        _jetDeepCsv_b[20];   //[_nJets]
  Double_t        _jetDeepCsv_c[20];   //[_nJets]
  Double_t        _jetDeepCsv_bb[20];   //[_nJets]
  UInt_t          _jetHadronFlavor[20];   //[_nJets]
  UInt_t          _jetId[20];   //[_nJets]
  Double_t        _met;
  Double_t        _metJECDown;
  Double_t        _metJECUp;
  Double_t        _metUnclDown;
  Double_t        _metUnclUp;
  Double_t        _metPhi;
  Double_t        _metPhiJECDown;
  Double_t        _metPhiJECUp;
  Double_t        _metPhiUnclDown;
  Double_t        _metPhiUnclUp;
    
  // List of branches
  TBranch        *b__runNb;   //!
  TBranch        *b__lumiBlock;   //!
  TBranch        *b__eventNb;   //!
  TBranch        *b__nVertex;   //!
  TBranch        *b__weight;   //!
  TBranch        *b__lheHTIncoming;   //!
  TBranch        *b__ctauHN;   //!
  TBranch        *b__nLheWeights;   //!
  TBranch        *b__lheWeight;   //!
  TBranch        *b__nTrueInt;   //!
  TBranch        *b__ttgEventType;   //!
  TBranch        *b__zgEventType;   //!
  TBranch        *b__gen_met;   //!
  TBranch        *b__gen_metPhi;   //!
  TBranch        *b__gen_nPh;   //!
  TBranch        *b__gen_phPt;   //!
  TBranch        *b__gen_phEta;   //!
  TBranch        *b__gen_phPhi;   //!
  TBranch        *b__gen_phE;   //!
  TBranch        *b__gen_phMomPdg;   //!
  TBranch        *b__gen_phIsPrompt;   //!
  TBranch        *b__gen_phMinDeltaR;   //!
  TBranch        *b__gen_phPassParentage;   //!
  TBranch        *b__gen_nL;   //!
  TBranch        *b__gen_lPt;   //!
  TBranch        *b__gen_lEta;   //!
  TBranch        *b__gen_lPhi;   //!
  TBranch        *b__gen_lE;   //!
  //TBranch        *b__gen_lFlavor;   //!
  TBranch        *b__gen_lCharge;   //!
  TBranch        *b__gen_lMomPdg;   //!
  TBranch        *b__gen_lIsPrompt;   //!
  TBranch        *b__passHN_1l;   //!
  TBranch        *b__HLT_Ele27_WPTight_Gsf;   //!
  TBranch        *b__HLT_Ele27_WPTight_Gsf_prescale;   //!
  TBranch        *b__HLT_IsoMu24;   //!
  TBranch        *b__HLT_IsoMu24_prescale;   //!
  TBranch        *b__HLT_IsoTkMu24;   //!
  TBranch        *b__HLT_IsoTkMu24_prescale;   //!
  TBranch        *b__passHN_eee;   //!
  TBranch        *b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;   //!
  TBranch        *b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale;   //!
  TBranch        *b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   //!
  TBranch        *b__passHN_eem;   //!
  TBranch        *b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL;   //!
  TBranch        *b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__passHN_emm;   //!
  TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL;   //!
  TBranch        *b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale;   //!
  TBranch        *b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
  TBranch        *b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale;   //!
  TBranch        *b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL;   //!
  TBranch        *b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale;   //!
  TBranch        *b__passHN_mmm;   //!
  TBranch        *b__HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx;   //!
  TBranch        *b__HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_prescale;   //!
  TBranch        *b__HLT_TripleMu_12_10_5;   //!
  TBranch        *b__HLT_TripleMu_12_10_5_prescale;   //!
  TBranch        *b__passMET;   //!
  TBranch        *b__HLT_MET300;   //!
  TBranch        *b__HLT_MET300_prescale;   //!
  TBranch        *b__HLT_HT350_MET100;   //!
  TBranch        *b__HLT_HT350_MET100_prescale;   //!
  TBranch        *b__HLT_AllMET300;   //!
  TBranch        *b__HLT_AllMET300_prescale;   //!
  TBranch        *b__HLT_AllMET170;   //!
  TBranch        *b__HLT_AllMET170_prescale;   //!
  TBranch        *b__HLT_jet;   //!
  TBranch        *b__HLT_jet_prescale;   //!
  TBranch        *b__HLT_dijet;   //!
  TBranch        *b__HLT_dijet_prescale;   //!
  TBranch        *b__HLT_MET170_BeamHaloCleaned;   //!
  TBranch        *b__HLT_MET170_BeamHaloCleaned_prescale;   //!
  TBranch        *b__HLT_MET170_NotCleaned;   //!
  TBranch        *b__HLT_MET170_NotCleaned_prescale;   //!
  TBranch        *b__HLT_HT800;   //!
  TBranch        *b__HLT_HT800_prescale;   //!
  TBranch        *b__HLT_HT900;   //!
  TBranch        *b__HLT_HT900_prescale;   //!
  TBranch        *b__HLT_dijet55met110;   //!
  TBranch        *b__HLT_dijet55met110_prescale;   //!
  TBranch        *b__HLT_dijet70met120;   //!
  TBranch        *b__HLT_dijet70met120_prescale;   //!
  TBranch        *b__HLT_HT600;   //!
  TBranch        *b__HLT_HT600_prescale;   //!
  TBranch        *b__HLT_HT475;   //!
  TBranch        *b__HLT_HT475_prescale;   //!
  TBranch        *b__HLT_HT350;   //!
  TBranch        *b__HLT_HT350_prescale;   //!
  TBranch        *b__passMETFilters;   //!
  TBranch        *b__Flag_HBHENoiseFilter;   //!
  TBranch        *b__Flag_HBHENoiseIsoFilter;   //!
  TBranch        *b__Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
  TBranch        *b__Flag_goodVertices;   //!
  TBranch        *b__Flag_eeBadScFilter;   //!
  TBranch        *b__Flag_globalTightHalo2016Filter;   //!
  TBranch        *b__flag_badPFMuonFilter;   //!
  TBranch        *b__flag_badChCandFilter;   //!
  TBranch        *b__passTTG_e;   //!
  TBranch        *b__HLT_Ele105_CaloIdVT_GsfTrkIdT;   //!
  TBranch        *b__HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale;   //!
  TBranch        *b__HLT_Ele115_CaloIdVT_GsfTrkIdT;   //!
  TBranch        *b__HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale;   //!
  TBranch        *b__passTTG_ee;   //!
  TBranch        *b__HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
  TBranch        *b__HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale;   //!
  TBranch        *b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale;   //!
  TBranch        *b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW;   //!
  TBranch        *b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale;   //!
  TBranch        *b__passTTG_em;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;   //!
  TBranch        *b__HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale;   //!
  TBranch        *b__HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL;   //!
  TBranch        *b__HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale;   //!
  TBranch        *b__passTTG_m;   //!
  TBranch        *b__HLT_Mu50;   //!
  TBranch        *b__HLT_Mu50_prescale;   //!
  TBranch        *b__HLT_TkMu50;   //!
  TBranch        *b__HLT_TkMu50_prescale;   //!
  TBranch        *b__HLT_Mu45_eta2p1;   //!
  TBranch        *b__HLT_Mu45_eta2p1_prescale;   //!
  TBranch        *b__passTTG_mm;   //!
  TBranch        *b__HLT_Mu30_TkMu11;   //!
  TBranch        *b__HLT_Mu30_TkMu11_prescale;   //!
  TBranch        *b__nL;   //!
  TBranch        *b__nMu;   //!
  TBranch        *b__nEle;   //!
  TBranch        *b__nLight;   //!
  TBranch        *b__nTau;   //!
  TBranch        *b__nVFit;   //!
  TBranch        *b__nGoodLeading;   //!
  TBranch        *b__lIndex;   //!
  TBranch        *b__vertices;   //!
  TBranch        *b__lPt;   //!
  TBranch        *b__lEta;   //!
  TBranch        *b__lEtaSC;   //!
  TBranch        *b__lPhi;   //!
  TBranch        *b__lE;   //!
  TBranch        *b__lFlavor;   //!
  TBranch        *b__lCharge;   //!
  TBranch        *b__dxy;   //!
  TBranch        *b__dz;   //!
  TBranch        *b__3dIP;   //!
  TBranch        *b__3dIPSig;   //!
  TBranch        *b__2dIP;   //!
  TBranch        *b__2dIPSig;   //!
  TBranch        *b__lElectronMva;   //!
  TBranch        *b__lElectronPassEmu;   //!
  TBranch        *b__lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto;   //!
  TBranch        *b__lPOGVeto;   //!
  TBranch        *b__lPOGLoose;   //!
  TBranch        *b__lPOGMedium;   //!
  TBranch        *b__lPOGTight;   //!
  TBranch        *b__lpassConversionVeto;   //!
  TBranch        *b__muNumberInnerHits;   //!
  TBranch        *b__eleNumberInnerHitsMissing;   //!
  TBranch        *b__relIso;   //!
  TBranch        *b__puCorr;   //!
  TBranch        *b__absIso03;   //!
  TBranch        *b__absIso04;   //!
  TBranch        *b__sumNeutralHadronEt04;   //!
  TBranch        *b__sumChargedHadronPt04;   //!
  TBranch        *b__sumPhotonEt04;   //!
  TBranch        *b__sumNeutralHadronEt03;   //!
  TBranch        *b__sumChargedHadronPt03;   //!
  TBranch        *b__sumPhotonEt03;   //!
  TBranch        *b__trackIso;   //!
  TBranch        *b__ecalIso;   //!
  TBranch        *b__hcalIso;   //!
  TBranch        *b__deltaBIso;   //!
  TBranch        *b__ecalPFClusterIso;   //!
  TBranch        *b__hcalPFClusterIso;   //!
  TBranch        *b__ptRel;   //!
  TBranch        *b__ptRatio;   //!
  TBranch        *b__closestJetCsv;   //!
  TBranch        *b__selectedTrackMult;   //!
  TBranch        *b__muonSegComp;   //!
  TBranch        *b__tauMuonVeto;   //!
  TBranch        *b__tauEleVeto;   //!
  TBranch        *b__decayModeFindingNew;   //!
  TBranch        *b__tauVLooseMvaNew;   //!
  TBranch        *b__tauLooseMvaNew;   //!
  TBranch        *b__tauMediumMvaNew;   //!
  TBranch        *b__tauTightMvaNew;   //!
  TBranch        *b__tauVTightMvaNew;   //!
  TBranch        *b__tauVTightMvaOld;   //!
  TBranch        *b__lIsPrompt;   //!
  TBranch        *b__lMatchPdgId;   //!
  TBranch        *b__nPh;   //!
  TBranch        *b__phPt;   //!
  TBranch        *b__phEta;   //!
  TBranch        *b__phEtaSC;   //!
  TBranch        *b__phPhi;   //!
  TBranch        *b__phE;   //!
  TBranch        *b__phCutBasedLoose;   //!
  TBranch        *b__phCutBasedMedium;   //!
  TBranch        *b__phCutBasedTight;   //!
  TBranch        *b__phMva;   //!
  TBranch        *b__phRandomConeChargedIsolation;   //!
  TBranch        *b__phChargedIsolation;   //!
  TBranch        *b__phNeutralHadronIsolation;   //!
  TBranch        *b__phPhotonIsolation;   //!
  TBranch        *b__phSigmaIetaIeta;   //!
  TBranch        *b__phSigmaIetaIphi;   //!
  TBranch        *b__phHadronicOverEm;   //!
  TBranch        *b__phPassElectronVeto;   //!
  TBranch        *b__phHasPixelSeed;   //!
  TBranch        *b__phIsPrompt;   //!
  TBranch        *b__phMatchMCPhotonAN15165;   //!
  TBranch        *b__phMatchMCLeptonAN15165;   //!
  TBranch        *b__phMatchPdgId;   //!
  TBranch        *b__nJets;   //!
  TBranch        *b__jetPt;   //!
  TBranch        *b__jetPt_JECUp;   //!
  TBranch        *b__jetPt_JECDown;   //!
  TBranch        *b__jetPt_JERUp;   //!
  TBranch        *b__jetPt_JERDown;   //!
  TBranch        *b__jetEta;   //!
  TBranch        *b__jetPhi;   //!
  TBranch        *b__jetE;   //!
  TBranch        *b__jetCsvV2;   //!
  TBranch        *b__jetDeepCsv_udsg;   //!
  TBranch        *b__jetDeepCsv_b;   //!
  TBranch        *b__jetDeepCsv_c;   //!
  TBranch        *b__jetDeepCsv_bb;   //!
  TBranch        *b__jetHadronFlavor;   //!
  TBranch        *b__jetId;   //!
  TBranch        *b__met;   //!
  TBranch        *b__metJECDown;   //!
  TBranch        *b__metJECUp;   //!
  TBranch        *b__metUnclDown;   //!
  TBranch        *b__metUnclUp;   //!
  TBranch        *b__metPhi;   //!
  TBranch        *b__metPhiJECDown;   //!
  TBranch        *b__metPhiJECUp;   //!
  TBranch        *b__metPhiUnclDown;   //!
  TBranch        *b__metPhiUnclUp;   //!
    
    
  double hcounter[nSamples];
    
  for(int sam = 0; sam < nSamples; ++sam){
        
    cout<<"================================= 1"<<endl;
    cout<<"--------> "<< "name " << names[sam] << endl;
    cout<<"--------> "<< "fileList[sam] " << fileList[sam] << endl;
    hfile[sam] = new TFile("/Users/Martina/Desktop/file/new_displaced/"+fileList[sam],"read");
    hfile[sam]->cd("blackJackAndHookers");
    inputTree[sam] = static_cast<TTree*>(hfile[sam]->Get("blackJackAndHookers/blackJackAndHookersTree"));
    inputTree[sam]->SetBranchAddress("_runNb", &_runNb, &b__runNb);
    inputTree[sam]->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
    inputTree[sam]->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
    inputTree[sam]->SetBranchAddress("_nVertex", &_nVertex, &b__nVertex);
    inputTree[sam]->SetBranchAddress("_weight", &_weight, &b__weight);
    inputTree[sam]->SetBranchAddress("_lheHTIncoming", &_lheHTIncoming, &b__lheHTIncoming);
    inputTree[sam]->SetBranchAddress("_ctauHN", &_ctauHN, &b__ctauHN);
    inputTree[sam]->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
    inputTree[sam]->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
    inputTree[sam]->SetBranchAddress("_nTrueInt", &_nTrueInt, &b__nTrueInt);
    inputTree[sam]->SetBranchAddress("_ttgEventType", &_ttgEventType, &b__ttgEventType);
    inputTree[sam]->SetBranchAddress("_zgEventType", &_zgEventType, &b__zgEventType);
    inputTree[sam]->SetBranchAddress("_gen_met", &_gen_met, &b__gen_met);
    inputTree[sam]->SetBranchAddress("_gen_metPhi", &_gen_metPhi, &b__gen_metPhi);
    inputTree[sam]->SetBranchAddress("_gen_nPh", &_gen_nPh, &b__gen_nPh);
    inputTree[sam]->SetBranchAddress("_gen_phPt", _gen_phPt, &b__gen_phPt);
    inputTree[sam]->SetBranchAddress("_gen_phEta", _gen_phEta, &b__gen_phEta);
    inputTree[sam]->SetBranchAddress("_gen_phPhi", _gen_phPhi, &b__gen_phPhi);
    inputTree[sam]->SetBranchAddress("_gen_phE", _gen_phE, &b__gen_phE);
    inputTree[sam]->SetBranchAddress("_gen_phMomPdg", _gen_phMomPdg, &b__gen_phMomPdg);
    inputTree[sam]->SetBranchAddress("_gen_phIsPrompt", _gen_phIsPrompt, &b__gen_phIsPrompt);
    inputTree[sam]->SetBranchAddress("_gen_phMinDeltaR", _gen_phMinDeltaR, &b__gen_phMinDeltaR);
    inputTree[sam]->SetBranchAddress("_gen_phPassParentage", _gen_phPassParentage, &b__gen_phPassParentage);
    inputTree[sam]->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
    inputTree[sam]->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
    inputTree[sam]->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
    inputTree[sam]->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
    inputTree[sam]->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
    //inputTree[sam]->SetBranchAddress("_gen_lFlavor", _gen_lFlavor, &b__gen_lFlavor);
    inputTree[sam]->SetBranchAddress("_gen_lCharge", _gen_lCharge, &b__gen_lCharge);
    inputTree[sam]->SetBranchAddress("_gen_lMomPdg", _gen_lMomPdg, &b__gen_lMomPdg);
    inputTree[sam]->SetBranchAddress("_gen_lIsPrompt", _gen_lIsPrompt, &b__gen_lIsPrompt);
    inputTree[sam]->SetBranchAddress("_passHN_1l", &_passHN_1l, &b__passHN_1l);
    inputTree[sam]->SetBranchAddress("_HLT_Ele27_WPTight_Gsf", &_HLT_Ele27_WPTight_Gsf, &b__HLT_Ele27_WPTight_Gsf);
    inputTree[sam]->SetBranchAddress("_HLT_Ele27_WPTight_Gsf_prescale", &_HLT_Ele27_WPTight_Gsf_prescale, &b__HLT_Ele27_WPTight_Gsf_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_IsoMu24", &_HLT_IsoMu24, &b__HLT_IsoMu24);
    inputTree[sam]->SetBranchAddress("_HLT_IsoMu24_prescale", &_HLT_IsoMu24_prescale, &b__HLT_IsoMu24_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_IsoTkMu24", &_HLT_IsoTkMu24, &b__HLT_IsoTkMu24);
    inputTree[sam]->SetBranchAddress("_HLT_IsoTkMu24_prescale", &_HLT_IsoTkMu24_prescale, &b__HLT_IsoTkMu24_prescale);
    inputTree[sam]->SetBranchAddress("_passHN_eee", &_passHN_eee, &b__passHN_eee);
    inputTree[sam]->SetBranchAddress("_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
    inputTree[sam]->SetBranchAddress("_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale", &_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale, &b__HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale, &b__HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_passHN_eem", &_passHN_eem, &b__passHN_eem);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &_HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale", &_HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale, &b__HLT_Mu8_DiEle12_CaloIdL_TrackIdL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale, &b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ, &b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale, &b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_passHN_emm", &_passHN_emm, &b__passHN_emm);
    inputTree[sam]->SetBranchAddress("_HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &_HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
    inputTree[sam]->SetBranchAddress("_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale", &_HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale, &b__HLT_DiMu9_Ele9_CaloIdL_TrackIdL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale, &b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale", &_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale, &b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale", &_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale, &b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale", &_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale, &b__HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale", &_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale, &b__HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale", &_HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale, &b__HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_prescale);
    inputTree[sam]->SetBranchAddress("_passHN_mmm", &_passHN_mmm, &b__passHN_mmm);
    inputTree[sam]->SetBranchAddress("_HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx", &_HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx, &b__HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx);
    inputTree[sam]->SetBranchAddress("_HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_prescale", &_HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_prescale, &b__HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_TripleMu_12_10_5", &_HLT_TripleMu_12_10_5, &b__HLT_TripleMu_12_10_5);
    inputTree[sam]->SetBranchAddress("_HLT_TripleMu_12_10_5_prescale", &_HLT_TripleMu_12_10_5_prescale, &b__HLT_TripleMu_12_10_5_prescale);
    inputTree[sam]->SetBranchAddress("_passMET", &_passMET, &b__passMET);
    inputTree[sam]->SetBranchAddress("_HLT_MET300", &_HLT_MET300, &b__HLT_MET300);
    inputTree[sam]->SetBranchAddress("_HLT_MET300_prescale", &_HLT_MET300_prescale, &b__HLT_MET300_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_HT350_MET100", &_HLT_HT350_MET100, &b__HLT_HT350_MET100);
    inputTree[sam]->SetBranchAddress("_HLT_HT350_MET100_prescale", &_HLT_HT350_MET100_prescale, &b__HLT_HT350_MET100_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_AllMET300", &_HLT_AllMET300, &b__HLT_AllMET300);
    inputTree[sam]->SetBranchAddress("_HLT_AllMET300_prescale", &_HLT_AllMET300_prescale, &b__HLT_AllMET300_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_AllMET170", &_HLT_AllMET170, &b__HLT_AllMET170);
    inputTree[sam]->SetBranchAddress("_HLT_AllMET170_prescale", &_HLT_AllMET170_prescale, &b__HLT_AllMET170_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_jet", &_HLT_jet, &b__HLT_jet);
    inputTree[sam]->SetBranchAddress("_HLT_jet_prescale", &_HLT_jet_prescale, &b__HLT_jet_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_dijet", &_HLT_dijet, &b__HLT_dijet);
    inputTree[sam]->SetBranchAddress("_HLT_dijet_prescale", &_HLT_dijet_prescale, &b__HLT_dijet_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_MET170_BeamHaloCleaned", &_HLT_MET170_BeamHaloCleaned, &b__HLT_MET170_BeamHaloCleaned);
    inputTree[sam]->SetBranchAddress("_HLT_MET170_BeamHaloCleaned_prescale", &_HLT_MET170_BeamHaloCleaned_prescale, &b__HLT_MET170_BeamHaloCleaned_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_MET170_NotCleaned", &_HLT_MET170_NotCleaned, &b__HLT_MET170_NotCleaned);
    inputTree[sam]->SetBranchAddress("_HLT_MET170_NotCleaned_prescale", &_HLT_MET170_NotCleaned_prescale, &b__HLT_MET170_NotCleaned_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_HT800", &_HLT_HT800, &b__HLT_HT800);
    inputTree[sam]->SetBranchAddress("_HLT_HT800_prescale", &_HLT_HT800_prescale, &b__HLT_HT800_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_HT900", &_HLT_HT900, &b__HLT_HT900);
    inputTree[sam]->SetBranchAddress("_HLT_HT900_prescale", &_HLT_HT900_prescale, &b__HLT_HT900_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_dijet55met110", &_HLT_dijet55met110, &b__HLT_dijet55met110);
    inputTree[sam]->SetBranchAddress("_HLT_dijet55met110_prescale", &_HLT_dijet55met110_prescale, &b__HLT_dijet55met110_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_dijet70met120", &_HLT_dijet70met120, &b__HLT_dijet70met120);
    inputTree[sam]->SetBranchAddress("_HLT_dijet70met120_prescale", &_HLT_dijet70met120_prescale, &b__HLT_dijet70met120_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_HT600", &_HLT_HT600, &b__HLT_HT600);
    inputTree[sam]->SetBranchAddress("_HLT_HT600_prescale", &_HLT_HT600_prescale, &b__HLT_HT600_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_HT475", &_HLT_HT475, &b__HLT_HT475);
    inputTree[sam]->SetBranchAddress("_HLT_HT475_prescale", &_HLT_HT475_prescale, &b__HLT_HT475_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_HT350", &_HLT_HT350, &b__HLT_HT350);
    inputTree[sam]->SetBranchAddress("_HLT_HT350_prescale", &_HLT_HT350_prescale, &b__HLT_HT350_prescale);
    inputTree[sam]->SetBranchAddress("_passMETFilters", &_passMETFilters, &b__passMETFilters);
    inputTree[sam]->SetBranchAddress("_Flag_HBHENoiseFilter", &_Flag_HBHENoiseFilter, &b__Flag_HBHENoiseFilter);
    inputTree[sam]->SetBranchAddress("_Flag_HBHENoiseIsoFilter", &_Flag_HBHENoiseIsoFilter, &b__Flag_HBHENoiseIsoFilter);
    inputTree[sam]->SetBranchAddress("_Flag_EcalDeadCellTriggerPrimitiveFilter", &_Flag_EcalDeadCellTriggerPrimitiveFilter, &b__Flag_EcalDeadCellTriggerPrimitiveFilter);
    inputTree[sam]->SetBranchAddress("_Flag_goodVertices", &_Flag_goodVertices, &b__Flag_goodVertices);
    inputTree[sam]->SetBranchAddress("_Flag_eeBadScFilter", &_Flag_eeBadScFilter, &b__Flag_eeBadScFilter);
    inputTree[sam]->SetBranchAddress("_Flag_globalTightHalo2016Filter", &_Flag_globalTightHalo2016Filter, &b__Flag_globalTightHalo2016Filter);
    inputTree[sam]->SetBranchAddress("_flag_badPFMuonFilter", &_flag_badPFMuonFilter, &b__flag_badPFMuonFilter);
    inputTree[sam]->SetBranchAddress("_flag_badChCandFilter", &_flag_badChCandFilter, &b__flag_badChCandFilter);
    inputTree[sam]->SetBranchAddress("_passTTG_e", &_passTTG_e, &b__passTTG_e);
    inputTree[sam]->SetBranchAddress("_HLT_Ele105_CaloIdVT_GsfTrkIdT", &_HLT_Ele105_CaloIdVT_GsfTrkIdT, &b__HLT_Ele105_CaloIdVT_GsfTrkIdT);
    inputTree[sam]->SetBranchAddress("_HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale", &_HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale, &b__HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Ele115_CaloIdVT_GsfTrkIdT", &_HLT_Ele115_CaloIdVT_GsfTrkIdT, &b__HLT_Ele115_CaloIdVT_GsfTrkIdT);
    inputTree[sam]->SetBranchAddress("_HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale", &_HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale, &b__HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale);
    inputTree[sam]->SetBranchAddress("_passTTG_ee", &_passTTG_ee, &b__passTTG_ee);
    inputTree[sam]->SetBranchAddress("_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b__HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
    inputTree[sam]->SetBranchAddress("_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale", &_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale, &b__HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL, &b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
    inputTree[sam]->SetBranchAddress("_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale", &_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale, &b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW", &_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW, &b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW);
    inputTree[sam]->SetBranchAddress("_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale", &_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale, &b__HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_prescale);
    inputTree[sam]->SetBranchAddress("_passTTG_em", &_passTTG_em, &b__passTTG_em);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", &_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL, &b__HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale", &_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale, &b__HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL", &_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL, &b__HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL);
    inputTree[sam]->SetBranchAddress("_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale", &_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale, &b__HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_prescale);
    inputTree[sam]->SetBranchAddress("_passTTG_m", &_passTTG_m, &b__passTTG_m);
    inputTree[sam]->SetBranchAddress("_HLT_Mu50", &_HLT_Mu50, &b__HLT_Mu50);
    inputTree[sam]->SetBranchAddress("_HLT_Mu50_prescale", &_HLT_Mu50_prescale, &b__HLT_Mu50_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu50", &_HLT_TkMu50, &b__HLT_TkMu50);
    inputTree[sam]->SetBranchAddress("_HLT_TkMu50_prescale", &_HLT_TkMu50_prescale, &b__HLT_TkMu50_prescale);
    inputTree[sam]->SetBranchAddress("_HLT_Mu45_eta2p1", &_HLT_Mu45_eta2p1, &b__HLT_Mu45_eta2p1);
    inputTree[sam]->SetBranchAddress("_HLT_Mu45_eta2p1_prescale", &_HLT_Mu45_eta2p1_prescale, &b__HLT_Mu45_eta2p1_prescale);
    inputTree[sam]->SetBranchAddress("_passTTG_mm", &_passTTG_mm, &b__passTTG_mm);
    inputTree[sam]->SetBranchAddress("_HLT_Mu30_TkMu11", &_HLT_Mu30_TkMu11, &b__HLT_Mu30_TkMu11);
    inputTree[sam]->SetBranchAddress("_HLT_Mu30_TkMu11_prescale", &_HLT_Mu30_TkMu11_prescale, &b__HLT_Mu30_TkMu11_prescale);
    inputTree[sam]->SetBranchAddress("_nL", &_nL, &b__nL);
    inputTree[sam]->SetBranchAddress("_nMu", &_nMu, &b__nMu);
    inputTree[sam]->SetBranchAddress("_nEle", &_nEle, &b__nEle);
    inputTree[sam]->SetBranchAddress("_nLight", &_nLight, &b__nLight);
    inputTree[sam]->SetBranchAddress("_nTau", &_nTau, &b__nTau);
    inputTree[sam]->SetBranchAddress("_nVFit", &_nVFit, &b__nVFit);
    inputTree[sam]->SetBranchAddress("_nGoodLeading", &_nGoodLeading, &b__nGoodLeading);
    inputTree[sam]->SetBranchAddress("_lIndex", _lIndex, &b__lIndex);
    inputTree[sam]->SetBranchAddress("_vertices", _vertices, &b__vertices);
    inputTree[sam]->SetBranchAddress("_lPt", _lPt, &b__lPt);
    inputTree[sam]->SetBranchAddress("_lEta", _lEta, &b__lEta);
    inputTree[sam]->SetBranchAddress("_lEtaSC", _lEtaSC, &b__lEtaSC);
    inputTree[sam]->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
    inputTree[sam]->SetBranchAddress("_lE", _lE, &b__lE);
    inputTree[sam]->SetBranchAddress("_lFlavor", _lFlavor, &b__lFlavor);
    inputTree[sam]->SetBranchAddress("_lCharge", _lCharge, &b__lCharge);
    inputTree[sam]->SetBranchAddress("_dxy", _dxy, &b__dxy);
    inputTree[sam]->SetBranchAddress("_dz", _dz, &b__dz);
    inputTree[sam]->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
    inputTree[sam]->SetBranchAddress("_3dIPSig", _3dIPSig, &b__3dIPSig);
    inputTree[sam]->SetBranchAddress("_2dIP", _2dIP, &b__2dIP);
    inputTree[sam]->SetBranchAddress("_2dIPSig", _2dIPSig, &b__2dIPSig);
    inputTree[sam]->SetBranchAddress("_lElectronMva", _lElectronMva, &b__lElectronMva);
    inputTree[sam]->SetBranchAddress("_lElectronPassEmu", _lElectronPassEmu, &b__lElectronPassEmu);
    inputTree[sam]->SetBranchAddress("_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto", _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, &b__lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto);
    inputTree[sam]->SetBranchAddress("_lPOGVeto", _lPOGVeto, &b__lPOGVeto);
    inputTree[sam]->SetBranchAddress("_lPOGLoose", _lPOGLoose, &b__lPOGLoose);
    inputTree[sam]->SetBranchAddress("_lPOGMedium", _lPOGMedium, &b__lPOGMedium);
    inputTree[sam]->SetBranchAddress("_lPOGTight", _lPOGTight, &b__lPOGTight);
    inputTree[sam]->SetBranchAddress("_lpassConversionVeto", _lpassConversionVeto, &b__lpassConversionVeto);
    inputTree[sam]->SetBranchAddress("_muNumberInnerHits", _muNumberInnerHits, &b__muNumberInnerHits);
    inputTree[sam]->SetBranchAddress("_eleNumberInnerHitsMissing", _eleNumberInnerHitsMissing, &b__eleNumberInnerHitsMissing);
    inputTree[sam]->SetBranchAddress("_relIso", _relIso, &b__relIso);
    inputTree[sam]->SetBranchAddress("_puCorr", _puCorr, &b__puCorr);
    inputTree[sam]->SetBranchAddress("_absIso03", _absIso03, &b__absIso03);
    inputTree[sam]->SetBranchAddress("_absIso04", _absIso04, &b__absIso04);
    inputTree[sam]->SetBranchAddress("_sumNeutralHadronEt04", _sumNeutralHadronEt04, &b__sumNeutralHadronEt04);
    inputTree[sam]->SetBranchAddress("_sumChargedHadronPt04", _sumChargedHadronPt04, &b__sumChargedHadronPt04);
    inputTree[sam]->SetBranchAddress("_sumPhotonEt04", _sumPhotonEt04, &b__sumPhotonEt04);
    inputTree[sam]->SetBranchAddress("_sumNeutralHadronEt03", _sumNeutralHadronEt03, &b__sumNeutralHadronEt03);
    inputTree[sam]->SetBranchAddress("_sumChargedHadronPt03", _sumChargedHadronPt03, &b__sumChargedHadronPt03);
    inputTree[sam]->SetBranchAddress("_sumPhotonEt03", _sumPhotonEt03, &b__sumPhotonEt03);
    inputTree[sam]->SetBranchAddress("_trackIso", _trackIso, &b__trackIso);
    inputTree[sam]->SetBranchAddress("_ecalIso", _ecalIso, &b__ecalIso);
    inputTree[sam]->SetBranchAddress("_hcalIso", _hcalIso, &b__hcalIso);
    inputTree[sam]->SetBranchAddress("_deltaBIso", _deltaBIso, &b__deltaBIso);
    inputTree[sam]->SetBranchAddress("_ecalPFClusterIso", _ecalPFClusterIso, &b__ecalPFClusterIso);
    inputTree[sam]->SetBranchAddress("_hcalPFClusterIso", _hcalPFClusterIso, &b__hcalPFClusterIso);
    inputTree[sam]->SetBranchAddress("_ptRel", _ptRel, &b__ptRel);
    inputTree[sam]->SetBranchAddress("_ptRatio", _ptRatio, &b__ptRatio);
    inputTree[sam]->SetBranchAddress("_closestJetCsv", _closestJetCsv, &b__closestJetCsv);
    inputTree[sam]->SetBranchAddress("_selectedTrackMult", _selectedTrackMult, &b__selectedTrackMult);
    inputTree[sam]->SetBranchAddress("_muonSegComp", _muonSegComp, &b__muonSegComp);
    inputTree[sam]->SetBranchAddress("_tauMuonVeto", _tauMuonVeto, &b__tauMuonVeto);
    inputTree[sam]->SetBranchAddress("_tauEleVeto", _tauEleVeto, &b__tauEleVeto);
    inputTree[sam]->SetBranchAddress("_decayModeFindingNew", _decayModeFindingNew, &b__decayModeFindingNew);
    inputTree[sam]->SetBranchAddress("_tauVLooseMvaNew", _tauVLooseMvaNew, &b__tauVLooseMvaNew);
    inputTree[sam]->SetBranchAddress("_tauLooseMvaNew", _tauLooseMvaNew, &b__tauLooseMvaNew);
    inputTree[sam]->SetBranchAddress("_tauMediumMvaNew", _tauMediumMvaNew, &b__tauMediumMvaNew);
    inputTree[sam]->SetBranchAddress("_tauTightMvaNew", _tauTightMvaNew, &b__tauTightMvaNew);
    inputTree[sam]->SetBranchAddress("_tauVTightMvaNew", _tauVTightMvaNew, &b__tauVTightMvaNew);
    inputTree[sam]->SetBranchAddress("_tauVTightMvaOld", _tauVTightMvaOld, &b__tauVTightMvaOld);
    inputTree[sam]->SetBranchAddress("_lIsPrompt", _lIsPrompt, &b__lIsPrompt);
    inputTree[sam]->SetBranchAddress("_lMatchPdgId", _lMatchPdgId, &b__lMatchPdgId);
    inputTree[sam]->SetBranchAddress("_nPh", &_nPh, &b__nPh);
    inputTree[sam]->SetBranchAddress("_phPt", _phPt, &b__phPt);
    inputTree[sam]->SetBranchAddress("_phEta", _phEta, &b__phEta);
    inputTree[sam]->SetBranchAddress("_phEtaSC", _phEtaSC, &b__phEtaSC);
    inputTree[sam]->SetBranchAddress("_phPhi", _phPhi, &b__phPhi);
    inputTree[sam]->SetBranchAddress("_phE", _phE, &b__phE);
    inputTree[sam]->SetBranchAddress("_phCutBasedLoose", _phCutBasedLoose, &b__phCutBasedLoose);
    inputTree[sam]->SetBranchAddress("_phCutBasedMedium", _phCutBasedMedium, &b__phCutBasedMedium);
    inputTree[sam]->SetBranchAddress("_phCutBasedTight", _phCutBasedTight, &b__phCutBasedTight);
    inputTree[sam]->SetBranchAddress("_phMva", _phMva, &b__phMva);
    inputTree[sam]->SetBranchAddress("_phRandomConeChargedIsolation", _phRandomConeChargedIsolation, &b__phRandomConeChargedIsolation);
    inputTree[sam]->SetBranchAddress("_phChargedIsolation", _phChargedIsolation, &b__phChargedIsolation);
    inputTree[sam]->SetBranchAddress("_phNeutralHadronIsolation", _phNeutralHadronIsolation, &b__phNeutralHadronIsolation);
    inputTree[sam]->SetBranchAddress("_phPhotonIsolation", _phPhotonIsolation, &b__phPhotonIsolation);
    inputTree[sam]->SetBranchAddress("_phSigmaIetaIeta", _phSigmaIetaIeta, &b__phSigmaIetaIeta);
    inputTree[sam]->SetBranchAddress("_phSigmaIetaIphi", _phSigmaIetaIphi, &b__phSigmaIetaIphi);
    inputTree[sam]->SetBranchAddress("_phHadronicOverEm", _phHadronicOverEm, &b__phHadronicOverEm);
    inputTree[sam]->SetBranchAddress("_phPassElectronVeto", _phPassElectronVeto, &b__phPassElectronVeto);
    inputTree[sam]->SetBranchAddress("_phHasPixelSeed", _phHasPixelSeed, &b__phHasPixelSeed);
    inputTree[sam]->SetBranchAddress("_phIsPrompt", _phIsPrompt, &b__phIsPrompt);
    inputTree[sam]->SetBranchAddress("_phMatchMCPhotonAN15165", _phMatchMCPhotonAN15165, &b__phMatchMCPhotonAN15165);
    inputTree[sam]->SetBranchAddress("_phMatchMCLeptonAN15165", _phMatchMCLeptonAN15165, &b__phMatchMCLeptonAN15165);
    inputTree[sam]->SetBranchAddress("_phMatchPdgId", _phMatchPdgId, &b__phMatchPdgId);
    inputTree[sam]->SetBranchAddress("_nJets", &_nJets, &b__nJets);
    inputTree[sam]->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
    inputTree[sam]->SetBranchAddress("_jetPt_JECUp", _jetPt_JECUp, &b__jetPt_JECUp);
    inputTree[sam]->SetBranchAddress("_jetPt_JECDown", _jetPt_JECDown, &b__jetPt_JECDown);
    inputTree[sam]->SetBranchAddress("_jetPt_JERUp", _jetPt_JERUp, &b__jetPt_JERUp);
    inputTree[sam]->SetBranchAddress("_jetPt_JERDown", _jetPt_JERDown, &b__jetPt_JERDown);
    inputTree[sam]->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
    inputTree[sam]->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
    inputTree[sam]->SetBranchAddress("_jetE", _jetE, &b__jetE);
    inputTree[sam]->SetBranchAddress("_jetCsvV2", _jetCsvV2, &b__jetCsvV2);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_udsg", _jetDeepCsv_udsg, &b__jetDeepCsv_udsg);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_b", _jetDeepCsv_b, &b__jetDeepCsv_b);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_c", _jetDeepCsv_c, &b__jetDeepCsv_c);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_bb", _jetDeepCsv_bb, &b__jetDeepCsv_bb);
    //inputTree[sam]->SetBranchAddress("_jetHadronFlavor", _jetHadronFlavor, &b__jetHadronFlavor);
    inputTree[sam]->SetBranchAddress("_jetId", _jetId, &b__jetId);
    inputTree[sam]->SetBranchAddress("_met", &_met, &b__met);
    inputTree[sam]->SetBranchAddress("_metJECDown", &_metJECDown, &b__metJECDown);
    inputTree[sam]->SetBranchAddress("_metJECUp", &_metJECUp, &b__metJECUp);
    inputTree[sam]->SetBranchAddress("_metUnclDown", &_metUnclDown, &b__metUnclDown);
    inputTree[sam]->SetBranchAddress("_metUnclUp", &_metUnclUp, &b__metUnclUp);
    inputTree[sam]->SetBranchAddress("_metPhi", &_metPhi, &b__metPhi);
    inputTree[sam]->SetBranchAddress("_metPhiJECDown", &_metPhiJECDown, &b__metPhiJECDown);
    inputTree[sam]->SetBranchAddress("_metPhiJECUp", &_metPhiJECUp, &b__metPhiJECUp);
    inputTree[sam]->SetBranchAddress("_metPhiUnclDown", &_metPhiUnclDown, &b__metPhiUnclDown);
    inputTree[sam]->SetBranchAddress("_metPhiUnclUp", &_metPhiUnclUp, &b__metPhiUnclUp);
    TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    _hCounter->Read("hCounter");
    hcounter[sam] = _hCounter->GetBinContent(1);
  }//end for on tree
    
  //******************* HISTO **********************
  const int nCat=42;
  const int nDist = 120;  //Number of distributions to plo
  const int FinalState=5;

  TH1D* Histos[nDist][nCat][nSamples_eff +1];

  const TString statusNames[FinalState] = {"_ALL", "_EEE", "_EMM", "_EEM", "_MMM"};
  const TString catNames[nCat]= {"_0", "_1", "_2", "_3", "_4", "_5", "_final", "_0_ossf", "_1_ossf", "_2_ossf", "_3_ossf", "_4_ossf", "_5_ossf", "_final_ossf", "_0_no_ossf", "_1_no_ossf", "_2_no_ossf", "_3_no_ossf", "_4_no_ossf", "_5_no_ossf", "_final_no_ossf", 
				 "__prompt__0", "__prompt__1", "__prompt__2", "__prompt__3", "__prompt__4", "__prompt__5", "__prompt__final", "__prompt__0_ossf", "__prompt__1_ossf", "__prompt__2_ossf", "__prompt__3_ossf", "__prompt__4_ossf", "__prompt__5_ossf", "__prompt__final_ossf", "__prompt__0_no_ossf", "__prompt__1_no_ossf", "__prompt__2_no_ossf", "__prompt__3_no_ossf", "__prompt__4_no_ossf", "__prompt__5_no_ossf", "__prompt__final_no_ossf"};
  const TString Histnames_ossf[nDist] = {"categories","categories_talk", "prompt",
					 "LeptonPt_le","LeptonPt_subl", "LeptonPt_tr", "Sum3Pt","Sum2Pt_lt","Sum2Pt_st","Sum2Pt_ls",
					 "Mlll","Mll_min","Mll_min_os", "Mll_l2l3", "Mll_pair", "MT_pair", "MT_3body", "MT_t", "MET", "MET_phi", "NJets", "NbJets","HT",
					 "dxy_l","dz_l","3dIP_l","2dIP_l", "3dIPSig_l", "2dIPSig_l",
					 "dxy_s","dz_s","3dIP_s","2dIP_s", "3dIPSig_s", "2dIPSig_s",
					 "dxy_t","dz_t","3dIP_t","2dIP_t", "3dIPSig_t", "2dIPSig_t",
					 "relIso03_l","absIso03_l","relIso04_l","absIso04_l","trackIso_l", "deltaBIso_l", "sumChargedHadronPt03_l",
					 "relIso03_s","absIso03_s","relIso04_s","absIso04_s","trackIso_s", "deltaBIso_s", "sumChargedHadronPt03_s",
					 "relIso03_t","absIso03_t","relIso04_t","absIso04_t","trackIso_t", "deltaBIso_t", "sumChargedHadronPt03_t",
					 "DeltaR_pair","DeltaR_lt","DeltaR_st", "DeltaPhi_pair","DeltaPhi_lt","DeltaPhi_st",
					 "vertex_X_ls", "vertex_Y_ls", "vertex_Z_ls", "vertex_R_ls","vertex_Rsign_ls","vertex_R2D_ls","vertex_R2Dsign_ls", "vertex_chi2_ls", "vertex_normalized_chi2_ls","vertex_chi2_ls_talk", "vertex_normalized_chi2_ls_talk","chi2_Y_ls", "R_X_ls",
					 "vertex_X_lt", "vertex_Y_lt", "vertex_Z_lt", "vertex_R_lt","vertex_Rsign_lt","vertex_R2D_lt","vertex_R2Dsign_lt", "vertex_chi2_lt", "vertex_normalized_chi2_lt","vertex_chi2_lt_talk", "vertex_normalized_chi2_lt_talk","chi2_Y_lt", "R_X_lt",
					 "vertex_X_st", "vertex_Y_st", "vertex_Z_st", "vertex_R_st","vertex_Rsign_st","vertex_R2D_st","vertex_R2Dsign_st", "vertex_chi2_st", "vertex_normalized_chi2_st","vertex_chi2_st_talk", "vertex_normalized_chi2_st_talk","chi2_Y_st", "R_X_st",
					 "minDeltaPhi", "ratioDeltaPhi1", "ratioDeltaPhi2","ratioDeltaPhi_withmin1", "ratioDeltaPhi_withmin2","RatioDeltaRPhi1","RatioDeltaRPhi2","mlll_met", "trasverseMass_met",
					 "SR_2DSIP", "SR_vtxR", "SR_vtxR2D", "SR_vtxX"};


  const TString Xaxes[nDist] = {"categories"," ","#prompt",
				"P_{T} #left(l_{1} #right) (GeV)","P_{T} #left(l_{2} #right) (GeV)", "P_{T} #left(l_{3} #right) (GeV)", "SumP_{T}#left(3leptons#right) (GeV)","SumP_{T}#left(l+l_{3} #right) (GeV)","SumP_{T}#left(l_{2}+l_{3} #right) (GeV)","SumP_{T}#left(l+l_{2} #right) (GeV)",
				"M_{lll} (GeV)","M_{ll}min (GeV)","M_{ll}minOS (GeV)","M_{ll} #left(l_{2}+l_{3} #right) (GeV)", "M_{ll} #left(M_{Z} pair#right) (GeV)", "M_{T} #left(no M_{Z} pair#right) (GeV)","M_{T} #left(l_{3} #right) (GeV)","M_{T} #left(3l+MET#right) (GeV)","MET (GeV)", "MET_phi (GeV)", "number of jets", "number of b-jets","HT",
				"#left|#it{d}_{xy}#right| #left(l_{1} #right) (cm)","#left|#it{d}_{z}#right| #left(l_{1} #right) (cm)","3D IP #left(l_{1} #right) (cm)","2D IP #left(l_{1} #right) (cm)", "3D SIP #left(l_{1} #right)", "2D SIP #left(l_{1} #right)",
				"#left|#it{d}_{xy}#right| #left(l_{2} #right) (cm)","#left|#it{d}_{z}#right| #left(l_{2} #right) (cm)","3D IP #left(l_{2} #right) (cm)","2D IP #left(l_{2} #right) (cm)", "3D SIP #left(l_{2} #right)", "2D SIP #left(l_{2} #right)",
				"#left|#it{d}_{xy}#right| #left(l_{3} #right) (cm)","#left|#it{d}_{z}#right| #left(l_{3} #right) (cm)","3D IP #left(l_{3} #right) (cm)","2D IP #left(l_{3} #right) (cm)", "3D SIP #left(l_{3} #right)", "2D SIP #left(l_{3} #right)",
				"relative Iso_{03} #left(l_{1} #right)","absolute Iso_{03} #left(l_{1} #right)", "relative Iso_{04} #left(l_{1} #right)","absolute Iso_{04} #left(l_{1} #right)","trackIso #left(l_{1} #right)", "delta#beta Iso_{03} #left(l_{1} #right)", "sumChargedHadronPt03 #left(l_{1} #right)",
				"relative Iso_{03} #left(l_{2} #right)","absolute Iso_{03} #left(l_{2} #right)","relative Iso_{04} #left(l_{2} #right)","absolute Iso_{04} #left(l_{2} #right)","trackIso #left(l_{2} #right)", "delta#beta Iso_{03} #left(l_{2} #right)", "sumChargedHadronPt03 #left(l_{2} #right)",
				"relative Iso_{03} #left(l_{3} #right)","absolute Iso_{03} #left(l_{3} #right)","relative Iso_{04} #left(l_{3} #right)","absolute Iso_{04} #left(l_{3} #right)","trackIso #left(l_{3} #right)", "delta#beta Iso_{03} #left(l_{3} #right)", "sumChargedHadronPt03 #left(l_{3} #right)",
				"#Delta#it{R} #left(M_{Z} pair#right)","#Delta#it{R} #left(l_{1}-l_{3} #right)","#Delta#it{R} #left(l_{2}-l_{3} #right)", "#Delta#phi #left(M_{Z} pair#right)","#Delta#phi #left(l_{1}-l_{3} #right)","#Delta#phi #left(l_{2}-l_{3} #right)",
				"vertex X (l_{1}-l_{2} ) (cm)", "vertex Y (l_{1}-l_{2} ) (cm)", "vertex Z (l_{1}-l_{2} ) (cm)", "vertex R (l_{1}-l_{2} ) (cm)","vertex R significance (l_{1}-l_{2} )", "vertex R2D (l_{1}-l_{2} ) (cm)","vertex R2D significance (l_{1}-l_{2} )", "vertex #chi ^{2} (l_{1}-l_{2} )", "vertex normalize#it{d}_#chi ^{2} (l_{1}-l_{2} )","vertex #chi ^{2} (l_{1}-l_{2} )", "vertex normalize#it{d}_#chi ^{2} (l_{1}-l_{2} )","#chi^{2}/vertex Y (l_{1}-l_{2} )", "R/vertex X  (l_{1}-l_{2} )",
				"vertex X (l_{1}-l_{3} ) (cm)", "vertex Y (l_{1}-l_{3} ) (cm)", "vertex Z (l_{1}-l_{3} ) (cm)", "vertex R (l_{1}-l_{3} ) (cm)","vertex R significance (l_{1}-l_{3} )", "vertex R2D (l_{1}-l_{3} ) (cm)","vertex R2D significance (l_{1}-l_{3} )", "vertex #chi ^{2} (l_{1}-l_{3} )", "vertex normalize#it{d}_#chi ^{2} (l_{1}-l_{3} )","vertex #chi ^{2} (l_{1}-l_{3} )", "vertex normalize#it{d}_#chi ^{2} (l_{1}-l_{3} )","#chi^{2}/vertex Y (l_{1}-l_{3} )", "R/vertex X  (l_{1}-l_{3} )",
				"vertex X (l_{2}-l_{3} ) (cm)", "vertex Y (l_{2}-l_{3} ) (cm)", "vertex Z (l_{2}-l_{3} ) (cm)", "vertex R (l_{2}-l_{3} ) (cm)","vertex R significance (l_{2}-l_{3} )", "vertex R2D (l_{2}-l_{3} ) (cm)","vertex R2D significance (l_{2}-l_{3} )", "vertex #chi ^{2} (l_{2}-l_{3} )", "vertex normalize#it{d}_#chi ^{2} (l_{2}-l_{3} )","vertex #chi ^{2} (l_{2}-l_{3} )", "vertex normalize#it{d}_#chi ^{2} (l_{2}-l_{3} )", "#chi^{2}/vertex Y (l_{2}-l_{3} )", "R/vertex X  (l_{2}-l_{3} )",
				"min #Delta#phi #left(l_{1}-other#right)",
				"(#Delta#phi(l_{2}-l_{3}))/(#Delta#phi(l_{1}-l_{2}))",
				"(#Delta#phi(l_{1}-l_{2}))/(#Delta#phi(l_{2}-l_{3}))",
				"(#Delta#phi(l_{2}-l_{3}))/(min #Delta#phi (l_{1}-other))",
				"(min #Delta#phi (l_{1}-other))/(#Delta#phi(l_{2}-l_{3}))",
				"(#Delta#it{R}(l_{2}-l_{3}))/(min #Delta#phi (l_{1}-other))",
				"(min #Delta#phi (l_{1}-other))/(#Delta#it{R}(l_{2}-l_{3}))",




				"M_{lll}+MET (GeV)", "MT+MET(GeV)",


				"SR_2DSIP", "SR_vtxR", "SR_vtxR2D", "SR_vtxX"


  };


  const TString Units[nDist] = {" "," "," ",
				"GeV", "GeV", "GeV", "GeV", "GeV","GeV","GeV", 
				"GeV","GeV","GeV", "GeV", "GeV", "GeV", "GeV","GeV","GeV","GeV", "","","GeV",
				"cm","cm","cm","cm","","",
				"cm","cm","cm","cm","","",
				"cm","cm","cm","cm","","",
				"", "GeV", "", "GeV", "GeV","", "GeV",
				"", "GeV", "", "GeV", "GeV","", "GeV",
				"", "GeV", "", "GeV", "GeV","", "GeV",
				"","","","","","",
				"cm","cm","cm", "cm","","cm","", "","","","","cm^{-1}", "",
				"cm","cm","cm", "cm","","cm","", "","","","","cm^{-1}", "",
				"cm","cm","cm", "cm","","cm","", "","","","","cm^{-1}", "",
				"", "", "","","","","","GeV", "GeV",
				"","","",""};
  const double HistMin[nDist] = {-0.5,-0.5,-0.5,
				 0, 0, 0, 0, 0, 0, 0,
				 20, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 40,0,
				 0, 0, 0, 0, 0, 0,
				 0, 0, 0, 0, 0, 0,
				 0, 0, 0, 0, 0, 0,
				 0, 0, 0, 0, 0, 0,0,
				 0, 0, 0, 0, 0, 0,0,
				 0, 0, 0, 0, 0, 0,0,
				 0,0,0,0,0,0,
				 -50, -50, -200, 0,0,0,0, -20,-20,0,0,0,0,
				 -50, -50, -200, 0,0,0,0, -20,-20,0,0,0,0,
				 -50, -50, -200, 0,0,0,0, -20,-20,0,0,0,0,
				 0,0,0,0,0,0,0,0,0,
				 -0.5,-0.5,-0.5,-0.5};
  const double HistMax[nDist] = { 41.5,11.5, 3.5,
				  200, 100, 100, 100, 100, 100, 100,
				  140,100,100, 100, 100, 100, 100,  100, 150,3.5, 10,10, 100,
				  0.5, 1, 20, 2, 100, 100,
				  0.5, 1, 20, 2, 100, 100,
				  0.5, 1, 20, 2, 100, 100,
				  1, 20, 1, 30, 40, 20,20,
				  1, 20, 1, 30, 40, 20,20,
				  1, 20, 1, 30, 40, 20,20,
				  4,4,4,3.3,3.3,3.3,
				  50, 50, 200, 50,1000,30,1000, 200,50,60,100,100,100,
				  50, 50, 200, 50,1000,30,1000, 200,50,60,100,100,100,
				  50, 50, 200, 50,1000,30,1000, 200,50,60,100,100,100,
				  3.15,0.3, 300,0.3,300,0.4,400,200, 200,
				  5.5,8.5,8.5,6.5};
  const int nBins[nDist] =      { 42,12,4,
				  25, 25, 25, 25, 25, 25, 25,
				  21,25,25, 25, 25, 25, 25,  25, 25,20, 10,10, 100,
				  30, 50, 30, 50, 50, 40,
				  30, 50, 30, 50, 50, 40,
				  30, 50, 30, 50, 50, 40,
				  50, 100, 50, 100, 100, 50,100,
				  50, 100, 50, 100, 100, 50,100,
				  50, 100, 50, 100, 100, 50,100,
				  20,20,20,20,20,20,
				  25, 25, 100, 50,100,30,100, 50,100,20,10, 100,100,
				  25, 25, 100, 50,100,30,100, 50,100,20,10, 100,100,
				  25, 25, 100, 50,100,30,100, 50,100 ,20,10,100,100,
				  20,30,60,30,60,40,80,50,50,
				  6,9,9,7};    
  cout<<"------ 1"<<endl;
  for(int i = 0; i < nDist; ++i){
    float BinWidth = (HistMax[i] - HistMin[i])/nBins[i];
    std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
    for(int effsam = 0; effsam < nSamples_eff + 1; ++effsam){
      for(int cat = 0; cat < nCat; ++cat){               
	Histos[i][cat][effsam] = new TH1D(eff_names[effsam] + catNames[cat] + Histnames_ossf[i] , eff_names[effsam] + catNames[cat] + Histnames_ossf[i]  + ";" + Xaxes[i] + "; events /" + Yaxis + Units[i], nBins[i], HistMin[i], HistMax[i]);
	Histos[i][cat][effsam]->Sumw2();
      }      }
  }
}
//Calculate the center of the maximum bin of each histogram
double maxBinC[nDist];
for(int i = 0; i < nDist; ++i){
  maxBinC[i] = Histos[i][0][0]->GetBinCenter(Histos[i][0][0]->GetNbinsX());
 }
    
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PARAMETERS AND CUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
const double met_cuts =80;
const double mlll_cuts = 85;
const int number_veto_leptons=3;
const double b_jets_wp= 0.5426;
const double b_jets_pt= 25;

const double MVA_cuts_pt15[3] = {0.77, 0.56, 0.48};
const double MVA_cuts_pt25[3] = {0.52, 0.11, -0.01};
const double newMVALooseFR[3]= {-0.02, -0.52, -0.52}; 

const double isolation_loose=0.6;
const double isolation_tight=0.1;

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
//Loop over all samples
Double_t Counter[nSamples];
for(int i = 0; i  < nSamples; ++i){
  Counter[i] = 0;
 }
Double_t scale[nSamples];
    
    
//-------------->  LOOP ON ALL SAMPLES
for(int sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
  Long64_t nEntries = inputTree[sam]->GetEntries();
  cout<<"----------------"<<endl;
  if(sam != 0){
    cout<<"=========  effsam: "<<effsam<<" "<<names[sam]<<"   "<<sam<<endl;
    if(names[sam] == names[sam -1]) --effsam;
    cout<<"+++++++++  effsam: "<<effsam<<" "<<names[sam]<<"   "<<sam<<endl;

  }
  if(sam >= 0){ // 1 data set
    scale[sam] = xSections[sam]*luminosity*1000/(hcounter[sam]);
  }
  //if (effsam == 0) continue;
  std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
  cout << effsam << endl;
        
  //if (effsam !=1) continue;
  //  if (sam > 7) continue;
        
  double counters_cut[20] ;
  double counters_cut_ossf[20] ;
  double counters_cut_no_ossf[20] ;


  // if (sam <= 12 ){

  for (int i =0; i< 20; i++){
    counters_cut[i] =0.;
    counters_cut_ossf[i] =0.;
    counters_cut_no_ossf[i] =0.;
  }
  // }
  double progress = 0; 	//For printing progress bar
        
        
  //--------------> LOOOP ON THE TREE"S ENTRIES
  for (Long64_t it = 0; it < nEntries; ++it){
    inputTree[sam]->GetEntry(it);
      
    //if (it%10000 == 0) cout<<'.'<<flush;
            
    if(it%100 == 0 && it != 0){
      progress += (double) (100./nEntries);
      printProgress(progress);
    } else if(it == nEntries -1){
      progress = 1.;
      printProgress(progress);
    }
            
            
    double scal = 0;
    scal = scale[sam]*_weight;
     
            
    if (effsam == 0) continue;
    if (effsam == 17) continue;   
            
            
            
    // if(_nL < number_veto_leptons) continue;
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VECTORS AND VARIABLES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    int               nBjets = 0;
    unsigned*         ind = new unsigned[_nL];	//new indices of good leptons
    //double*           conePt = new double[_nL];
    double            faxtore_FR=0;            
    double            prov_index[3]={-1,-1,-1};
    double            prov_number_tight[1]= {-1};
    double            prov_number_prompt[1]= {-1};
    double            prov_fattore[1]= {-1};
    int               skip_event[1]= {-1};
    double            faxtore[1]= {-100};            
    TLorentzVector    lepton_reco[3];
    int               flavors_3l[3];
    int               charge_3l[3];
    TLorentzVector    sum_3l_rec;	//M_3l
    TLorentzVector    pair [3];
    int               kind[1] ={-1}; // 0 no-ossf
    TLorentzVector    sum_2l_rec_pair; 	//M_2l best Z candidate
    int               event_clas[1]={-1}; // 1* = eee   2* = emm   3* = eem  4* = mmm
    int               check_mt= -1;
    int               check_deltaR= -1;
    int               check_mt_os= -1;
    Double_t          delta_R_max=-1;
    Double_t          delta_R_min=-1;
    Double_t          _mll_min=50000;
    Double_t          _mll_min_os=50000;
    TLorentzVector    lepton_transv[3];
    TLorentzVector    METvec;
    double            m_T=0;
    double            MET=0;
    double            MET_PHI=0;           
    int               nominal=-1;
    int               jec_check=-1;
    int               unclustered_met=-1;
    int               up=-1;
    int               down=-1;
    double            new_met[1]= {0};
    double            new_met_phi[1]= {0};
    int               new_number_bjets[1]= {0};            
    int               new_pt_checks= -1;
    bool              trigger_fired = false;
    bool              low_pt_event=false;
    bool              high_pt_event = false;            
    unsigned*         _SR_low= new unsigned[8];
    double             search_region_fill[1]={-1};
    bool              data_control_region=false;



    unsigned          lCount = 0;	//Count number of FO leptons that are not taus
    unsigned*         _isFO= new unsigned[_nL];
    Bool_t            _passedMVA90[_nL];   
      
     
    unsigned           promptC = 0;
    double             low_mass_pt_base[1];
      
    unsigned*          _isWithTrack= new unsigned[_nL];
    unsigned           wTrack=0;

    double            iV_ls=0;
    double            iV_lt=0;
    double            iV_st=0;

    double            _vertex_X[3];
    double            _vertex_Y[3];
    double            _vertex_Z[3];
    double            _vertex_R2D[3];
    double            _vertex_sR2D[3];
    double            _vertex_R[3];
    double            _vertex_sR[3];
    double            _vertex_chi2[3];
    double            _vertex_normchi2[3];

    unsigned         ind_new_leading=0;
    unsigned         ind_new_p=0;
    unsigned         ind_new_pp=0;

    bool             isAll = false;
    bool             isAll_met_mlll = false;

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<            
    //------------------------------------------------------------ selection class
    //FO
    for(unsigned l = 0; l < _nL; ++l){
      _isFO[l] = false;
      if (_lFlavor[l] == 0 && _lPt[l] < 10 ) continue;
      if (_lFlavor[l] == 1 && _lPt[l] < 5 ) continue;
      if ( _relIso[l] > 0.2) continue;
      if (_lFlavor[l] == 0 && !_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[l]) continue;
      if (_lFlavor[l] == 1 && !_lPOGLoose[l]) continue; 
      _isFO[l] = true;
      if(_isFO[l] && _lFlavor[l] != 2){
	ind[lCount] = l;
	++lCount;
      }
    } 

  
    double*  conePt = new double[_nL];
    for(unsigned l = 0; l < lCount; ++l){
      conePt[ind[l]] =0.;
      conePt[ind[l]] = _lPt[ind[l]];
    }
      

    if(lCount < 3) continue;
    //Order FO leptons by Pt
    unsigned* ordind = new unsigned[lCount];
    std::set<unsigned> usedLep;      
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
    //======================= new Selection!
    ind_new_leading = -1;
    int counter_leading=0;
    for(unsigned l = 0; l < lCount; ++l){
      if (counter_leading == 0){
	if ((_lFlavor[ind[l]]== 1 && _lPOGLoose[ind[l]] && _lPOGMedium[ind[l]] && _relIso[ind[l]] < 0.1 && TMath::Abs(_dxy[ind[l]]) < 0.05 && TMath::Abs(_dz[ind[l]]) < 0.1 && fabs(_3dIPSig [ind[l]]) < 4 && _lPt[ind[l]]> 25 ) || (_lFlavor[ind[l]]== 0 && _lPOGLoose[ind[l]] && _lPOGMedium[ind[l]] && _relIso[ind[l]] < 0.1 && TMath::Abs(_dxy[ind[l]]) < 0.05 && TMath::Abs(_dz[ind[l]]) < 0.1 && fabs(_3dIPSig [ind[l]]) < 4 && _lPt[ind[l]]> 30)   ) {
	  ++counter_leading;
	  ind_new_leading = ind[l];
	}
      }
    }
    if (counter_leading <1)continue;
    unsigned displacedC = 0;
    unsigned* _isDisplaced= new unsigned[_nL];
    for(unsigned l = 0; l < lCount; ++l){
      _isDisplaced[ind[l]] = false;
      //if ((TMath::Abs(_dxy[ind[l]]) > 0.05 || TMath::Abs(_dz[ind[l]]) > 0.1 || fabs(_3dIPSig[ind[l]]) > 4))	_isDisplaced[ind[l]] = true;
      //if ( fabs(_3dIPSig[ind[l]]) > 4)	_isDisplaced[ind[l]] = true;
      if (TMath::Abs(_dxy[ind[l]]) > 0.05 )	_isDisplaced[ind[l]] = true;


    }     
    for(unsigned l = 0; l < lCount; ++l){
      if(_isDisplaced[ind[l]]) ++displacedC;	
    }
    if (displacedC < 2) continue;
    //if (displacedC == 2) continue;
    // ===================== 3 leptons selected ======================
     
    // leading lepton
    for (int i =0; i < lCount; i ++){
      if (ind[i] == ind_new_leading) {
	lepton_reco[0].SetPtEtaPhiE( _lPt[ind[i]],  _lEta[ind[i]], _lPhi[ind[i]], _lE[ind[i]]);
	flavors_3l[0]=_lFlavor[ind[i]];
	charge_3l[0]=_lCharge[ind[i]];
	conePt[0] =  _lPt[ind[i]];
      }
    }
    // easy case == just 2 leptons
    TLorentzVector   lepton_tobeselected[20];
    int              index_displaced[20];
    int index_l[2];


    if (displacedC == 2) {
      int check_disp= 0;
      for (int i =0; i < lCount; i ++){
	if (_isDisplaced[ind[i]] && check_disp==0) {
	  lepton_reco[1].SetPtEtaPhiE( _lPt[ind[i]],  _lEta[ind[i]], _lPhi[ind[i]], _lE[ind[i]]);
	  flavors_3l[1]=_lFlavor[ind[i]];
	  charge_3l[1]=_lCharge[ind[i]];
	  conePt[1] =  _lPt[ind[i]];
	  check_disp = 1;
	  index_l[0] =ind[i]; 
	  index_displaced[0] = ind[i];
	  lepton_tobeselected[0].SetPtEtaPhiE( _lPt[ind[i]],  _lEta[ind[i]], _lPhi[ind[i]], _lE[ind[i]]);
	}
	if (_isDisplaced[ind[i]] && check_disp!=0) {
	  lepton_reco[2].SetPtEtaPhiE( _lPt[ind[i]],  _lEta[ind[i]], _lPhi[ind[i]], _lE[ind[i]]);
	  flavors_3l[2]=_lFlavor[ind[i]];
	  charge_3l[2]=_lCharge[ind[i]];
	  conePt[2] =  _lPt[ind[i]];
	  index_l[1] =ind[i]; 
	  index_displaced[1] = ind[i];
	  lepton_tobeselected[1].SetPtEtaPhiE( _lPt[ind[i]],  _lEta[ind[i]], _lPhi[ind[i]], _lE[ind[i]]);

	}
      }
    }
     
    // more than 2 leptons!!! 
    if (displacedC > 2){
	
      for (int i =0; i < 20; i ++){
	index_displaced[i] = -50;
	lepton_tobeselected[i].SetPtEtaPhiE( 0,  0, 0, 0);
      }
      int displacedC_check=0;
      for (int i =0; i < lCount; i ++){
	if (_isDisplaced[ind[i]]){	    
	  lepton_tobeselected[displacedC_check].SetPtEtaPhiE( _lPt[ind[i]],  _lEta[ind[i]], _lPhi[ind[i]], _lE[ind[i]]);
	  index_displaced[displacedC_check] = ind[i];
	  displacedC_check++;
	}
      }
      find_leptons(selezione,displacedC, lepton_tobeselected,index_displaced, index_l);
      lepton_reco[1].SetPtEtaPhiE( _lPt[index_l[0]],  _lEta[index_l[0]], _lPhi[index_l[0]], _lE[index_l[0]]);
      flavors_3l[1]=_lFlavor[index_l[0]];
      charge_3l[1]=_lCharge[index_l[0]];
      conePt[1] =  _lPt[index_l[0]];
      lepton_reco[2].SetPtEtaPhiE( _lPt[index_l[1]],  _lEta[index_l[1]], _lPhi[index_l[1]], _lE[index_l[1]]);
      flavors_3l[2]=_lFlavor[index_l[1]];
      charge_3l[2]=_lCharge[index_l[1]];
      conePt[2] =  _lPt[index_l[1]];
    }// end >2

      

    /////////////////////////////// --------------- ////////////////////////
    /////////////////////////////// --------------- ////////////////////////

    if (_lIsPrompt[ind_new_leading]) promptC++;
    if (_lIsPrompt[index_l[1]])promptC++;
    if (_lIsPrompt[index_l[0]])promptC++;

       
    //calculate the index for the vertex
    iV_ls= _lIndex[ind_new_leading] * 100 +  _lIndex[index_l[0]];
    iV_st= _lIndex[index_l[0]] * 100 +  _lIndex[index_l[1]];
    iV_lt= _lIndex[ind_new_leading] * 100 +  _lIndex[index_l[1]];

    //[0] == ls
    //[1] == st
    //[2] == lt

    _vertex_chi2[0]     = -10;
    _vertex_normchi2[0] = -10;
    _vertex_chi2[1]     = -10;
    _vertex_normchi2[1] = -10;
    _vertex_chi2[2]     = -10;
    _vertex_normchi2[2] = -10;

      
    //M_3l
    sum_3l_rec.SetPtEtaPhiE(0,0,0,0);
    sum_3l_rec= (lepton_reco[0]+ lepton_reco[1]+lepton_reco[2] );


    if (charge_3l[0] == charge_3l[1] ){
      _vertex_chi2[0]     = -12;
      _vertex_normchi2[0] = -12;	
    }
    if (charge_3l[1] == charge_3l[2] ){
      _vertex_chi2[1]     = -14;
      _vertex_normchi2[1] = -14;
    }
    if (charge_3l[0] == charge_3l[2] ){
      _vertex_chi2[2]     = -16;
      _vertex_normchi2[2] = -16;	
    }

    if (_vertices[0][0] == 102 && _vertices[1][0] == 305){
      _vertices[0][0] = 102;
      _vertices[1][0] = 302; 
    }
    if (_vertices[0][0] == 103 && _vertices[1][0] == 206){
      _vertices[0][0] = 103;
      _vertices[1][0] = 203; 
    }
    if (_vertices[0][0] == 103 && _vertices[1][0] == 204){
      _vertices[0][0] = 103;
      _vertices[1][0] = 203; 
    }


    for(unsigned v = 0; v < _nVFit; ++v){
      if ((_vertices[v][0] == (_lIndex[ind_new_leading] * 100 +  _lIndex[index_l[0]])) || (_vertices[v][0] == (_lIndex[ind_new_leading]  +  _lIndex[index_l[0]]*100))) {
	_vertex_X[0]        = _vertices[v][1];
	_vertex_Y[0]        = _vertices[v][2];
	_vertex_Z[0]        = _vertices[v][3];
	_vertex_R2D[0]      = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2]);
	_vertex_sR2D[0]     = TMath::Sqrt(derivate2_with_sigmaR2D(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2])  + derivate2_with_sigmaR2D(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2])   + 2*_vertices[v][7]*_vertices[v][7] *derivateR2D(_vertices[v][1],_vertices[v][1], _vertices[v][2])*derivateR2D(_vertices[v][2],_vertices[v][1], _vertices[v][2]) );
	_vertex_R[0]        = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2] + _vertices[v][3]*_vertices[v][3]);
	_vertex_sR[0]= TMath::Sqrt(derivate2_with_sigmaR(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   derivate2_with_sigmaR(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   derivate2_with_sigmaR(_vertices[v][3], _vertices[v][6],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   2*_vertices[v][7]*_vertices[v][7]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   2*_vertices[v][9]*_vertices[v][9]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   2*_vertices[v][8]*_vertices[v][8]*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3])  );
	_vertex_chi2[0]     = _vertices[v][11];
	_vertex_normchi2[0] = _vertices[v][11]/_vertices[v][10];
      }
      if ((_vertices[v][0] == (_lIndex[index_l[0]] * 100 +  _lIndex[index_l[1]] )) ||(_vertices[v][0] == (_lIndex[index_l[0]] +  _lIndex[index_l[1]] *100) )) {
	_vertex_X[1]        = _vertices[v][1];
	_vertex_Y[1]        = _vertices[v][2];
	_vertex_Z[1]        = _vertices[v][3];
	_vertex_R2D[1]      = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2]);
	_vertex_sR2D[1]     = TMath::Sqrt(derivate2_with_sigmaR2D(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2])  + derivate2_with_sigmaR2D(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2])   + 2*_vertices[v][7]*_vertices[v][7] *derivateR2D(_vertices[v][1],_vertices[v][1], _vertices[v][2])*derivateR2D(_vertices[v][2],_vertices[v][1], _vertices[v][2]) );
	_vertex_R[1]        = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2] + _vertices[v][3]*_vertices[v][3]);
	_vertex_sR[1]= TMath::Sqrt(derivate2_with_sigmaR(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   derivate2_with_sigmaR(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   derivate2_with_sigmaR(_vertices[v][3], _vertices[v][6],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   2*_vertices[v][7]*_vertices[v][7]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   2*_vertices[v][9]*_vertices[v][9]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   2*_vertices[v][8]*_vertices[v][8]*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3])  );
	_vertex_chi2[1]     = _vertices[v][11];
	_vertex_normchi2[1] = _vertices[v][11]/_vertices[v][10];
      }
      if ((_vertices[v][0] == (_lIndex[ind_new_leading] * 100 +  _lIndex[index_l[1]])) || (_vertices[v][0] == (_lIndex[ind_new_leading] +  _lIndex[index_l[1]]*100))) {
	_vertex_X[2]        = _vertices[v][1];
	_vertex_Y[2]        = _vertices[v][2];
	_vertex_Z[2]        = _vertices[v][3];
	_vertex_R2D[2]      = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2]);
	_vertex_sR2D[2]     = TMath::Sqrt(derivate2_with_sigmaR2D(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2])  + derivate2_with_sigmaR2D(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2])   + 2*_vertices[v][7]*_vertices[v][7] *derivateR2D(_vertices[v][1],_vertices[v][1], _vertices[v][2])*derivateR2D(_vertices[v][2],_vertices[v][1], _vertices[v][2]) );
	_vertex_R[2]        = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2] + _vertices[v][3]*_vertices[v][3]);
	_vertex_sR[2]= TMath::Sqrt(derivate2_with_sigmaR(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   derivate2_with_sigmaR(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   derivate2_with_sigmaR(_vertices[v][3], _vertices[v][6],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   2*_vertices[v][7]*_vertices[v][7]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   2*_vertices[v][9]*_vertices[v][9]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				   2*_vertices[v][8]*_vertices[v][8]*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3])  );
	_vertex_chi2[2]     = _vertices[v][11];
	_vertex_normchi2[2] = _vertices[v][11]/_vertices[v][10];
      }	           
    }// end loop vertices




    //--------------------------------------------------------------------------------   
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   

            
    //================== event classification ========================
    // ---------------- > OSSF or NO_OSSF
    ossf_no_ossf( kind, pair,lepton_reco[0], lepton_reco[1], lepton_reco[2], flavors_3l, charge_3l);
    if (kind[0]  == -1) continue;

    //---------------- > M_2l best Z candidate
    sum_2l_rec_pair.SetPtEtaPhiE(0,0,0,0);
    sum_2l_rec_pair= (pair[0]+pair[1] );
    bool ossf_event= false;
    if (kind[0] == 1) ossf_event = true;

      
    METvec.SetPtEtaPhiE(_met, 0, _metPhi,_met);    
    lepton_transv[0].SetPtEtaPhiE(lepton_reco[0].Pt(),0, lepton_reco[0].Phi(), lepton_reco[0].Pt());
    lepton_transv[1].SetPtEtaPhiE(lepton_reco[1].Pt(),0, lepton_reco[1].Phi(), lepton_reco[1].Pt());
    lepton_transv[2].SetPtEtaPhiE(lepton_reco[2].Pt(),0, lepton_reco[2].Phi(), lepton_reco[2].Pt());	

    double mT_met = 0;
    mT_met = (lepton_transv[0] +lepton_transv[1] +lepton_transv[2] + METvec).Mag();
    double mlll_met = 0;
    mlll_met=  (lepton_reco[0] +lepton_reco[1] +lepton_reco[2]).Mag() +_met;
  
    


    // ---------------- > CHANNELS
    class_os( event_clas,  flavors_3l, charge_3l);
    if (event_clas[0] == -1) continue;
    // 1* = eee
    // 2* = emm
    // 3* = eem
    // 4* = mmm
    isAll = true;

    // ---------------- > M_min
    check_mt=-1;    
    _mll_min = (lepton_reco[0]+lepton_reco[1]).M();
    check_mt= 2;    
    if ( (lepton_reco[0]+lepton_reco[2]).M() < _mll_min){
      _mll_min = (lepton_reco[0]+lepton_reco[2]).M();
      check_mt= 1;
    }
    if ( (lepton_reco[1]+lepton_reco[2]).M() < _mll_min) {
      _mll_min = (lepton_reco[1]+lepton_reco[2]).M();
      check_mt= 0;
    }   
    // 2 == ls
    // 1 == lt
    // 0 == st
    // ---------------- > M_min OS
    check_mt_os=-1;
    if ((charge_3l[0] != charge_3l[1] )) {
      _mll_min_os = (lepton_reco[0]+lepton_reco[1]).M();
      check_mt_os= 2;
    }
    if ((charge_3l[0] != charge_3l[2] ) && (lepton_reco[0]+lepton_reco[2]).M() < _mll_min_os){
      _mll_min_os = (lepton_reco[0]+lepton_reco[2]).M();
      check_mt_os= 1;
    }
    if ((charge_3l[1] != charge_3l[2] ) && (lepton_reco[1]+lepton_reco[2]).M() < _mll_min_os) {
      _mll_min_os = (lepton_reco[1]+lepton_reco[2]).M();
      check_mt_os= 0;
    }   
    check_deltaR=-1;
    delta_R_min = lepton_reco[0].DeltaR(lepton_reco[1]);
    check_deltaR=2;
    if (lepton_reco[0].DeltaR(lepton_reco[2]) < delta_R_min) {
      delta_R_min = lepton_reco[0].DeltaR(lepton_reco[2]);
      check_deltaR=1;
    }
    if (lepton_reco[1].DeltaR(lepton_reco[2]) < delta_R_min) {
      delta_R_min = lepton_reco[1].DeltaR(lepton_reco[2]);
      check_deltaR=0;
    }



    double min_delta_phi = 0;
    min_delta_phi = fabs(lepton_reco[0].DeltaPhi(lepton_reco[1]));
    if (fabs(lepton_reco[0].DeltaPhi(lepton_reco[2])) < min_delta_phi)  min_delta_phi = fabs(lepton_reco[0].DeltaPhi(lepton_reco[2]));
    double ration1_deltaphi = fabs(lepton_reco[1].DeltaPhi(lepton_reco[2]))/fabs(lepton_reco[0].DeltaPhi(lepton_reco[1]));
    double ration2_deltaphi = fabs(lepton_reco[0].DeltaPhi(lepton_reco[1]))/fabs(lepton_reco[2].DeltaPhi(lepton_reco[1]));

    double ration1_deltaphimin = fabs(lepton_reco[1].DeltaPhi(lepton_reco[2]))/min_delta_phi;
    double ration2_deltaphimin = min_delta_phi/fabs(lepton_reco[2].DeltaPhi(lepton_reco[1]));

    double ration1_deltaphiRmin = fabs(lepton_reco[1].DeltaR(lepton_reco[2]))/min_delta_phi;
    double ration2_deltaphiRmin = min_delta_phi/fabs(lepton_reco[2].DeltaR(lepton_reco[1]));


    //=============================================================
      
      
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      HISTOGRAMS     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    double pt_cone_leading=          lepton_reco[0].Pt() ;
    double pt_cone_sub_leading=      lepton_reco[1].Pt();
    double pt_cone_trailing=         lepton_reco[2].Pt();
  
    double trackIso_subleading =_trackIso[index_l[0]];
    double trackIso_trailing= _trackIso[index_l[1]];

    if (lepton_reco[1].DeltaR(lepton_reco[2]) < 0.3){
      _trackIso[index_l[1]] =  _trackIso[index_l[1]] - pt_cone_sub_leading;
      _trackIso[index_l[0]] =  _trackIso[index_l[0]] - pt_cone_trailing;
    }
            

    bool _ss_st = false;
    bool _ss_mll = false;
    bool _ss_delta = false;
    bool _os_st = false;
    bool _os_mll = false;
    bool _os_delta = false;
    bool _no_prompt = false;
    bool _prompt = false;
    bool _delta = false;
    bool _chi2 = false;

    bool _delta2=false;
    bool _mlll_cuts = false;
    bool _mll_cut = false;
    bool _minDeltaphi=false;
      
    if (effsam < 11 && promptC < 3) _no_prompt = true;
    if ((effsam < 11 && promptC == 3) ||(effsam >=11) ) _prompt = true;

    bool selection_0=false;
    bool selection_1=false;
    bool selection_2=false;
    bool selection_3=false;
    bool selection_4=false;
    bool selection_5=false;
    bool selection_final=false;


    // cout<<min_delta_phi<<endl;


    if (charge_3l[2] != charge_3l[1])                                        selection_0 = true;
    if ( selection_0 && lepton_reco[1].DeltaR(lepton_reco[2]) < 1)           selection_1 = true;
    if ( selection_1 &&  min_delta_phi > 2)                                  selection_2 = true;
    if ( selection_2 &&  sum_3l_rec.M() > 45 && sum_3l_rec.M() < 85)         selection_3 = true;
    if ( selection_3 &&  _met < 85)                                          selection_4 = true;
    if ( selection_4 &&  _vertex_chi2[1] < 50)                               selection_5 = true;
    if ( selection_5 )                                                       selection_final = true;
     
    //************************   SR    ********************************
    bool _sr_bool_2dsip[6] ;
    bool _sr_bool_vtxR[9];
    bool _sr_bool_vtx2R[9];
    bool _sr_bool_vtxX[7];
      
    for (int i = 0; i < 6 ; i++){
      _sr_bool_2dsip[i]=false;
    }
    for (int i = 0; i < 9 ; i++){
      _sr_bool_vtxR[i]=false;
      _sr_bool_vtx2R[i]=false;
    }
    for (int i = 0; i < 7 ; i++){
      _sr_bool_vtxX[i]=false;
    }
    if (fabs(_2dIPSig[index_l[0]]) < 10)                                       _sr_bool_2dsip[0] = true;
    if (fabs(_2dIPSig[index_l[0]]) > 10 &&  fabs(_2dIPSig[index_l[0]]) < 20)   _sr_bool_2dsip[1] = true;
    if (fabs(_2dIPSig[index_l[0]]) > 20 &&  fabs(_2dIPSig[index_l[0]]) < 40)   _sr_bool_2dsip[2] = true;
    if (fabs(_2dIPSig[index_l[0]]) > 40 &&  fabs(_2dIPSig[index_l[0]]) < 60)   _sr_bool_2dsip[3] = true;
    if (fabs(_2dIPSig[index_l[0]]) > 60 &&  fabs(_2dIPSig[index_l[0]]) < 80)   _sr_bool_2dsip[4] = true;
    if (fabs(_2dIPSig[index_l[0]]) > 80)                                       _sr_bool_2dsip[5] = true;

    if (_vertex_R[1] < 2 )                                                     _sr_bool_vtxR[0] = true;
    if (_vertex_R[1] > 2  && _vertex_R[1] < 4)                                 _sr_bool_vtxR[1] = true;
    if (_vertex_R[1] > 4  && _vertex_R[1] < 6)                                 _sr_bool_vtxR[2] = true;
    if (_vertex_R[1] > 6  && _vertex_R[1] < 8)                                 _sr_bool_vtxR[3] = true;
    if (_vertex_R[1] > 8  && _vertex_R[1] < 10)                                _sr_bool_vtxR[4] = true;
    if (_vertex_R[1] > 10  && _vertex_R[1] < 20)                               _sr_bool_vtxR[5] = true;
    if (_vertex_R[1] > 20  && _vertex_R[1] < 30)                               _sr_bool_vtxR[6] = true;
    if (_vertex_R[1] > 30  && _vertex_R[1] < 50)                               _sr_bool_vtxR[7] = true;
    if (_vertex_R[1] > 50  )                                                   _sr_bool_vtxR[8] = true;

    if (_vertex_R2D[1] < 2 )                                                    _sr_bool_vtx2R[0] = true;
    if (_vertex_R2D[1] > 2  && _vertex_R2D[1] < 4)                              _sr_bool_vtx2R[1] = true;
    if (_vertex_R2D[1] > 4  && _vertex_R2D[1] < 6)                              _sr_bool_vtx2R[2] = true;
    if (_vertex_R2D[1] > 6  && _vertex_R2D[1] < 8)                              _sr_bool_vtx2R[3] = true;
    if (_vertex_R2D[1] > 8  && _vertex_R2D[1] < 10)                             _sr_bool_vtx2R[4] = true;
    if (_vertex_R2D[1] > 10  && _vertex_R2D[1] < 15)                            _sr_bool_vtx2R[5] = true;
    if (_vertex_R2D[1] > 15  && _vertex_R2D[1] < 20)                            _sr_bool_vtx2R[6] = true;
    if (_vertex_R2D[1] > 20  && _vertex_R2D[1] < 30)                            _sr_bool_vtx2R[7] = true;
    if (_vertex_R2D[1] > 30  )                                                  _sr_bool_vtx2R[8] = true;

    if (_vertex_X[1] < 2 )                                                      _sr_bool_vtxX[0] = true;
    if (_vertex_X[1] > 2  && _vertex_X[1] < 5)                                  _sr_bool_vtxX[1] = true;
    if (_vertex_X[1] > 5  && _vertex_X[1] < 10)                                 _sr_bool_vtxX[2] = true;
    if (_vertex_X[1] > 10  && _vertex_X[1] < 20)                                _sr_bool_vtxX[3] = true;
    if (_vertex_X[1] > 20  && _vertex_X[1] < 30)                                _sr_bool_vtxX[4] = true;
    if (_vertex_X[1] > 30  && _vertex_X[1] < 40)                                _sr_bool_vtxX[5] = true;
    if (_vertex_X[1] > 40)                                                      _sr_bool_vtxX[6] = true;
     

    //****************************************************************
    bool _final_state_bool[4] ;
    for (int i =0; i < 4; i++){
      _final_state_bool[i]= false;
    }     
    bool _eee = false;
    bool _emm = false;	    
    bool _eem = false;
    bool _mmm = false;
    if (event_clas[0] == 1 || event_clas[0] == 10) _eee= true;
    if (event_clas[0] == 2 || event_clas[0] == 20) _emm= true;
    if (event_clas[0] == 3 || event_clas[0] == 30) _eem= true;
    if (event_clas[0] == 4 || event_clas[0] == 40) _mmm= true;

    if (_eee) _final_state_bool[0] = true;
    if (_emm) _final_state_bool[1] = true;
    if (_eem) _final_state_bool[2] = true;
    if (_mmm) _final_state_bool[3] = true;

            
    double values[nDist] ={static_cast<double>(0) ,static_cast<double>(0) , static_cast<double>(promptC),
			   
			   pt_cone_leading, pt_cone_sub_leading, pt_cone_trailing,
			   pt_cone_leading+ pt_cone_sub_leading+ pt_cone_trailing,
			   pt_cone_leading+ pt_cone_trailing,
			   pt_cone_sub_leading+ pt_cone_trailing,
			   pt_cone_leading+ pt_cone_sub_leading,
			   sum_3l_rec.M(),_mll_min,_mll_min_os,(lepton_reco[1]+lepton_reco[2]).M(), sum_2l_rec_pair.M(),0.,0.,0., _met,_metPhi, static_cast<double>(_nJets), static_cast<double>(nBjets),0.,
			   fabs(_dxy[ind_new_leading]),fabs(_dz[ind_new_leading]), fabs(_3dIP[ind_new_leading]), fabs(_2dIP[ind_new_leading]), fabs(_3dIPSig[ind_new_leading]), fabs(_2dIPSig[ind_new_leading]),
			   fabs(_dxy[index_l[0]]),fabs(_dz[index_l[0]]), fabs(_3dIP[index_l[0]]), fabs(_2dIP[index_l[0]]), fabs(_3dIPSig[index_l[0]]), fabs(_2dIPSig[index_l[0]]),
			   fabs(_dxy[index_l[1]]),fabs(_dz[index_l[1]]), fabs(_3dIP[index_l[1]]), fabs(_2dIP[index_l[1]]), fabs(_3dIPSig[index_l[1]]), fabs(_2dIPSig[index_l[1]]),
			   _relIso[ind_new_leading], _absIso03[ind_new_leading], _absIso04[ind_new_leading]/lepton_reco[0].Pt(), _absIso04[ind_new_leading], _trackIso[ind_new_leading], _deltaBIso[ind_new_leading],_sumChargedHadronPt03[ind_new_leading],
			   _relIso[index_l[0]], _absIso03[index_l[0]], _absIso04[index_l[0]]/lepton_reco[1].Pt(), _absIso04[index_l[0]], _trackIso[index_l[0]], _deltaBIso[index_l[0]],_sumChargedHadronPt03[index_l[0]],
			   _relIso[index_l[1]], _absIso03[index_l[1]], _absIso04[index_l[1]]/lepton_reco[2].Pt(), _absIso04[index_l[1]], _trackIso[index_l[1]], _deltaBIso[index_l[1]],_sumChargedHadronPt03[index_l[1]],
			   pair[0].DeltaR(pair[1]),lepton_reco[0].DeltaR(lepton_reco[2]),lepton_reco[1].DeltaR(lepton_reco[2]), TMath::Abs(pair[0].DeltaPhi(pair[1])),TMath::Abs(lepton_reco[0].DeltaPhi(lepton_reco[2])),TMath::Abs(lepton_reco[1].DeltaPhi(lepton_reco[2])),
			   _vertex_X[0], _vertex_Y[0],_vertex_Z[0], _vertex_R[0], _vertex_R[0]/ _vertex_sR[0], _vertex_R2D[0], _vertex_R2D[0]/ _vertex_sR2D[0],_vertex_chi2[0],_vertex_normchi2[0],_vertex_chi2[0],_vertex_normchi2[0],_vertex_chi2[0]/fabs(_vertex_Y[0]),_vertex_R[0]/fabs(_vertex_X[0]), 
			   _vertex_X[2], _vertex_Y[2],_vertex_Z[2], _vertex_R[2], _vertex_R[2]/ _vertex_sR[2], _vertex_R2D[2], _vertex_R2D[2]/ _vertex_sR2D[2],_vertex_chi2[2],_vertex_normchi2[2],_vertex_chi2[2],_vertex_normchi2[2],_vertex_chi2[2]/fabs(_vertex_Y[2]),_vertex_R[2]/fabs(_vertex_X[2]), 
			   _vertex_X[1], _vertex_Y[1],_vertex_Z[1], _vertex_R[1], _vertex_R[1]/ _vertex_sR[1], _vertex_R2D[1], _vertex_R2D[1]/ _vertex_sR2D[1],_vertex_chi2[1],_vertex_normchi2[1] ,_vertex_chi2[1],_vertex_normchi2[1] , _vertex_chi2[1]/fabs(_vertex_Y[1]), _vertex_R[1]/fabs(_vertex_X[1]),		
			   min_delta_phi, ration1_deltaphi, ration2_deltaphi, ration1_deltaphimin, ration2_deltaphimin, ration1_deltaphiRmin, ration2_deltaphiRmin, mlll_met,mT_met,
			   static_cast<double>(0) ,static_cast<double>(0) ,static_cast<double>(0) ,static_cast<double>(0) 
    };
            



            
    //      double values_sr[nDist_sr] = {search_region_fill[0], search_region_fill[0]};
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
 
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  FILLING  HISTOGRAMS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    unsigned fill = effsam;

    // ------------------- Histo kinematics
    for(int numero_histo = 0; numero_histo < nDist; ++numero_histo){
      if (numero_histo == 0  || numero_histo == 1 ||numero_histo == 119 ||numero_histo == 118 ||numero_histo == 117 ||numero_histo == 116 ) continue;
      if (selection_0)      Histos[numero_histo][0][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_1)      Histos[numero_histo][1][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_2)      Histos[numero_histo][2][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_3)      Histos[numero_histo][3][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_4)      Histos[numero_histo][4][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_5)      Histos[numero_histo][5][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_final)  Histos[numero_histo][6][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (ossf_event ){
	if (selection_0)      Histos[numero_histo][7][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_1)      Histos[numero_histo][8][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_2)      Histos[numero_histo][9][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_3)      Histos[numero_histo][10][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_4)      Histos[numero_histo][11][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_5)      Histos[numero_histo][12][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_final)  Histos[numero_histo][13][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      }//end ossf
      if (!ossf_event ){
	if (selection_0)      Histos[numero_histo][14][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_1)      Histos[numero_histo][15][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_2)      Histos[numero_histo][16][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_3)      Histos[numero_histo][17][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_4)      Histos[numero_histo][18][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_5)      Histos[numero_histo][19][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_final)  Histos[numero_histo][20][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      }//end ossf

      if (_prompt){
	if (selection_0)      Histos[numero_histo][21][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_1)      Histos[numero_histo][22][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_2)      Histos[numero_histo][23][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_3)      Histos[numero_histo][24][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_4)      Histos[numero_histo][25][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_5)      Histos[numero_histo][26][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_final)  Histos[numero_histo][27][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (ossf_event ){
	  if (selection_0)      Histos[numero_histo][28][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_1)      Histos[numero_histo][29][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_2)      Histos[numero_histo][30][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_3)      Histos[numero_histo][31][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_4)      Histos[numero_histo][32][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_5)      Histos[numero_histo][33][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_final)  Histos[numero_histo][34][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	}//end ossf
	if (!ossf_event ){
	  if (selection_0)      Histos[numero_histo][35][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_1)      Histos[numero_histo][36][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_2)      Histos[numero_histo][37][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_3)      Histos[numero_histo][38][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_4)      Histos[numero_histo][39][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_5)      Histos[numero_histo][40][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_final)  Histos[numero_histo][41][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	}//end ossf
      }//end prompt
    }//end histo

    for(int numero_histo = 0; numero_histo < nDist; ++numero_histo){
      if (numero_histo == 0  || numero_histo == 1 ||numero_histo == 119 ||numero_histo == 118 ||numero_histo == 117 ||numero_histo == 116 ) continue;
      if (selection_0)      Histos[numero_histo][0][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_1)      Histos[numero_histo][1][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_2)      Histos[numero_histo][2][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_3)      Histos[numero_histo][3][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_4)      Histos[numero_histo][4][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_5)      Histos[numero_histo][5][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (selection_final)  Histos[numero_histo][6][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      if (ossf_event ){
	if (selection_0)      Histos[numero_histo][7][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_1)      Histos[numero_histo][8][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_2)      Histos[numero_histo][9][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_3)      Histos[numero_histo][10][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_4)      Histos[numero_histo][11][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_5)      Histos[numero_histo][12][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_final)  Histos[numero_histo][13][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      }//end ossf
      if (!ossf_event ){
        if (selection_0)      Histos[numero_histo][14][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_1)      Histos[numero_histo][15][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_2)      Histos[numero_histo][16][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_3)      Histos[numero_histo][17][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_4)      Histos[numero_histo][18][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_5)      Histos[numero_histo][19][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_final)  Histos[numero_histo][20][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
      }//end ossf
    
      if (_prompt){
        if (selection_0)      Histos[numero_histo][21][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_1)      Histos[numero_histo][22][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_2)      Histos[numero_histo][23][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_3)      Histos[numero_histo][24][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_4)      Histos[numero_histo][25][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_5)      Histos[numero_histo][26][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (selection_final)  Histos[numero_histo][27][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        if (ossf_event ){
	  if (selection_0)      Histos[numero_histo][28][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_1)      Histos[numero_histo][29][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_2)      Histos[numero_histo][30][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_3)      Histos[numero_histo][31][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_4)      Histos[numero_histo][32][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_5)      Histos[numero_histo][33][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_final)  Histos[numero_histo][34][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        }//end ossf
        if (!ossf_event ){
	  if (selection_0)      Histos[numero_histo][35][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_1)      Histos[numero_histo][36][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_2)      Histos[numero_histo][37][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_3)      Histos[numero_histo][38][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_4)      Histos[numero_histo][39][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_5)      Histos[numero_histo][40][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_final)  Histos[numero_histo][41][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
        }//end ossf
      }//end prompt
    }//end histo





    if (selection_0)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(0), maxBinC[0]), scal);
    if (selection_1)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(1), maxBinC[0]), scal);
    if (selection_2)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(2), maxBinC[0]), scal);
    if (selection_3)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(3), maxBinC[0]), scal);
    if (selection_4)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(4), maxBinC[0]), scal);
    if (selection_5)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(5), maxBinC[0]), scal);
    if (selection_final)  Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(6), maxBinC[0]), scal);
    if (ossf_event ){
      if (selection_0)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(7), maxBinC[0]), scal);
      if (selection_1)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(8), maxBinC[0]), scal);
      if (selection_2)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(9), maxBinC[0]), scal);
      if (selection_3)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(10), maxBinC[0]), scal);
      if (selection_4)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(11), maxBinC[0]), scal);
      if (selection_5)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(12), maxBinC[0]), scal);
      if (selection_final)  Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(13), maxBinC[0]), scal);
    }//end ossf
    if (!ossf_event ){
      if (selection_0)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(14), maxBinC[0]), scal);
      if (selection_1)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(15), maxBinC[0]), scal);
      if (selection_2)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(16), maxBinC[0]), scal);
      if (selection_3)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(17), maxBinC[0]), scal);
      if (selection_4)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(18), maxBinC[0]), scal);
      if (selection_5)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(19), maxBinC[0]), scal);
      if (selection_final)  Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(20), maxBinC[0]), scal);
    }//end ossf

    if (_prompt){
      if (selection_0)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(21), maxBinC[0]), scal);
      if (selection_1)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(22), maxBinC[0]), scal);
      if (selection_2)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(23), maxBinC[0]), scal);
      if (selection_3)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(24), maxBinC[0]), scal);
      if (selection_4)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(25), maxBinC[0]), scal);
      if (selection_5)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(26), maxBinC[0]), scal);
      if (selection_final)  Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(27), maxBinC[0]), scal);
      if (ossf_event ){
	if (selection_0)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(28), maxBinC[0]), scal);
	if (selection_1)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(29), maxBinC[0]), scal);
	if (selection_2)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(30), maxBinC[0]), scal);
	if (selection_3)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(31), maxBinC[0]), scal);
	if (selection_4)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(32), maxBinC[0]), scal);
	if (selection_5)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(33), maxBinC[0]), scal);
	if (selection_final)  Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(34), maxBinC[0]), scal);
      }//end ossf
      if (!ossf_event ){
        if (selection_0)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(35), maxBinC[0]), scal);
        if (selection_1)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(36), maxBinC[0]), scal);
        if (selection_2)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(37), maxBinC[0]), scal);
        if (selection_3)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(38), maxBinC[0]), scal);
        if (selection_4)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(39), maxBinC[0]), scal);
        if (selection_5)      Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(40), maxBinC[0]), scal);
        if (selection_final)  Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(41), maxBinC[0]), scal);
      }//end ossf
    }//end prompt


    if (_prompt){
      if (ossf_event ){
        if (selection_0)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(0), maxBinC[0]), scal);
        if (selection_1)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(1), maxBinC[0]), scal);
        if (selection_2)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(2), maxBinC[0]), scal);
        if (selection_3)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(3), maxBinC[0]), scal);
        if (selection_4)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(4), maxBinC[0]), scal);
        if (selection_5)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(5), maxBinC[0]), scal);
      }//end ossf
      if (!ossf_event ){
        if (selection_0)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(6), maxBinC[0]), scal);
        if (selection_1)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(7), maxBinC[0]), scal);
        if (selection_2)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(8), maxBinC[0]), scal);
        if (selection_3)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(9), maxBinC[0]), scal);
        if (selection_4)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(10), maxBinC[0]), scal);
        if (selection_5)      Histos[1][0][fill]->Fill(TMath::Min(static_cast<double>(11), maxBinC[0]), scal);
      }//end ossf
    }//end prompt








    if (selection_final){
      if (_prompt){
        for (int i =0; i < 6; i++){
	  if (_sr_bool_2dsip[i])                                         Histos[116][27][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
	  if (ossf_event && _sr_bool_2dsip[i])                           Histos[116][34][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
	  if (!ossf_event && _sr_bool_2dsip[i])                          Histos[116][41][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
        }//end 2dsip
        for (int i =0; i < 9; i++){
	  if (_sr_bool_vtxR[i])                                         Histos[117][27][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
	  if (ossf_event && _sr_bool_vtxR[i])                           Histos[117][34][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
	  if (!ossf_event && _sr_bool_vtxR[i])                          Histos[117][41][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
            
	  if (_sr_bool_vtx2R[i])                                         Histos[118][27][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
	  if (ossf_event && _sr_bool_vtx2R[i])                           Histos[118][34][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
	  if (!ossf_event && _sr_bool_vtx2R[i])                          Histos[118][41][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
        }//end 2dsip
        for (int i =0; i < 7; i++){
	  if (_sr_bool_vtxX[i])                                         Histos[119][27][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
	  if (ossf_event && _sr_bool_vtxX[i])                           Histos[119][34][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
	  if (!ossf_event && _sr_bool_vtxX[i])                          Histos[119][41][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
        }//end 2dsip
      }//ed prompt
    
      for (int i =0; i < 6; i++){
        if (_sr_bool_2dsip[i])                                         Histos[116][6][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
      }//end 2dsip
      for (int i =0; i < 9; i++){
        if (_sr_bool_vtxR[i])                                         Histos[117][6][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
        if (_sr_bool_vtx2R[i])                                         Histos[118][6][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
      }//end 2dsip
      for (int i =0; i < 7; i++){
        if (_sr_bool_vtxX[i])                                         Histos[119][6][fill]->Fill(TMath::Min(static_cast<double>(i), maxBinC[0]), scal);
      }//end 2dsip
    }//end final
     
   
  }
    
 } 
//Split data and MC histograms for plotting and propagating uncertainties
TH1D* dataYields[nDist][nCat];
for(unsigned dist = 0; dist < nDist; ++dist){
  for(unsigned cat = 0; cat < nCat; ++cat){
    dataYields[dist][cat] = (TH1D*) Histos[dist][cat][11]->Clone();
  }
 }
    
// cout<< "ok 1"<<endl;
    
TH1D* bkgYields[nDist][nCat][nSamples_eff -10]; //change to nSamples_eff if sig is removed
for(unsigned dist = 0; dist < nDist; ++dist){
  for(unsigned cat = 0; cat < nCat; ++cat){
    for(unsigned effsam1 = 11; effsam1 < nSamples_eff +1 ; ++effsam1){
                
      //	cout<< effsam1<<"   "<<nSamples_eff<<endl;
      bkgYields[dist][cat][effsam1 -11] = (TH1D*) Histos[dist][cat][effsam1]->Clone();
      if(effsam1 > 11 && effsam1 < 17){
	dataYields[dist][cat]->Add(bkgYields[dist][cat][effsam1 -11]);
      }
    }
  }
 }
    
    
    
const TString sigNames[ nSamples_signal] = {"m_{N} = 1 ", "m_{N} = 2 ", "m_{N} = 3 ", "m_{N} = 4 ","m_{N} = 5 |V|^{2} = 10^{-5}", "m_{N} = 5.5 |V|^{2} = 10^{-5}", "m_{N} = 6 ", "m_{N} = 7 |V|^{2} = 10^{-5}", "m_{N} = 8 |V|^{2} = 10^{-5}", "m_{N} = 9 "};
TH1D* signals[nSamples_signal];
//Plot the yields as a function of the search region
for(unsigned dist = 0; dist < nDist; ++dist){
  for(unsigned cat = 0; cat < nCat; ++cat){
     
    signals[0] = (TH1D*) Histos[dist][cat][1]->Clone() ;
    signals[1] = (TH1D*) Histos[dist][cat][2]->Clone() ;
    signals[2] = (TH1D*) Histos[dist][cat][3]->Clone() ;
    signals[3] = (TH1D*) Histos[dist][cat][4]->Clone() ;
    signals[4] = (TH1D*) Histos[dist][cat][5]->Clone() ;
    signals[5] = (TH1D*) Histos[dist][cat][6]->Clone() ;
    signals[6] = (TH1D*) Histos[dist][cat][7]->Clone() ;
    signals[7] = (TH1D*) Histos[dist][cat][8]->Clone() ;
    signals[8] = (TH1D*) Histos[dist][cat][9]->Clone() ;
    signals[9] = (TH1D*) Histos[dist][cat][10]->Clone() ;
    // if (dist != 0 && cat !=0){      
    plotDataVSMC(cat,dist,dataYields[dist][cat], bkgYields[dist][cat], eff_names, nSamples_eff -  nSamples_signal - 1, Histnames_ossf[dist] + "_" +  catNames[cat], catNames[cat], true, 2, true, signals,  sigNames , nSamples_signal, false);
    //   cout<< "ok 4"<<endl;
    //}
    /*if (dist == 0 && cat ==0){
      cout<<"------------------"<<endl;
      plotDataVSMC(dataYields[dist][cat], bkgYields[dist][cat], eff_names, nSamples_eff -  nSamples_signal - 1, Histnames_ossf[dist] + "_" +  catNames[cat], true, 0, true, signals,  sigNames , nSamples_signal, false);
      for (int i = 0; i < 22; i ++){
	  
      cout<<"===="<<endl;
      cout<<"bin# "<<i+1<<") signal: "<<signals[0] -> GetBinContent (i+1)<<" "<<signals[1] -> GetBinContent (i+1)<<" "<<signals[2] -> GetBinContent (i+1)<<" "<<signals[3] -> GetBinContent (i+1)<<" "<<signals[4] -> GetBinContent (i+1)<<" "<<signals[5] -> GetBinContent (i+1)<<" "<<signals[6] -> GetBinContent (i+1)<<" "<<signals[7] -> GetBinContent (i+1)<<" "<<signals[8] -> GetBinContent (i+1)<<" "<<signals[9] -> GetBinContent (i+1)<<endl;
	  
      }    
      }*/


  }
 }
 
    
    
    
    
    
    
}// end analisi

//==================================================================

void Analysis_mc::find_leptons(int selezione,  unsigned displacedC, TLorentzVector lepton_tobeselected[20], int index_displaced[20], int index_s[2]){
  // ordine in pT
  int index_1=-5;
  int index_2= -5;
  if (selezione == 0 ){
    index_1 = index_displaced[0];
    index_2 = index_displaced[1];	
  }

  Double_t          delta_R_min=-1;
  Double_t          _mll_min=50000;
 
  //min mass
  if (selezione == 1 ){
    for (int h = 0; h <displacedC; h ++ ){
      for (int g = 0; g < displacedC; g ++){
	if (h != g) {
	  if (h == 0 && g ==1) {
	    _mll_min = (lepton_tobeselected[h]+lepton_tobeselected[g]).M();
	    index_1 = index_displaced[g];
	    index_2 = index_displaced[h];

	  }
	  if((lepton_tobeselected[h]+lepton_tobeselected[g]).M() < _mll_min) {
	    _mll_min = (lepton_tobeselected[h]+lepton_tobeselected[g]).M();
	    index_1 = index_displaced[g];
	    index_2 = index_displaced[h];

	  }
	}
      }// loop1
    }//loop 2
  }//selezione1

  //min delta
  if (selezione == 2 ){
    for (int h = 0; h <displacedC; h ++ ){
      for (int g = 0; g < displacedC; g ++){
	if (h != g) {
	  if (h == 0 && g ==1) {
	    delta_R_min=  lepton_tobeselected[h].DeltaR(lepton_tobeselected[g]);
	    index_1 = index_displaced[g];
	    index_2 = index_displaced[h];
	  }
	  if(lepton_tobeselected[h].DeltaR(lepton_tobeselected[g]) < delta_R_min) {
	    delta_R_min=  lepton_tobeselected[h].DeltaR(lepton_tobeselected[g]);
	    index_1 = index_displaced[g];
	    index_2 = index_displaced[h];
	  }
	}	
      }// loop1
    }//loop 2
  }
  index_s [0] = index_1;
  index_s[1]  = index_2; 

}//end funciton





//==================================================================
double Analysis_mc::maximum(double a, double b){
  double massimo=0;
  massimo = a;
  if (massimo < b) massimo = b;
  return massimo+1;
}
//___________________________________________________________________
double Analysis_mc::derivateR(double a, double x, double y, double z){
  double result =0;
  result = a /((x*x + y*y + z*z)*(x*x + y*y + z*z));
  return result;
}
//___________________________________________________________________
double Analysis_mc::derivate2_with_sigmaR(double a, double sa, double x, double y, double z){
  double result =0;
  result = (a*a)*(sa*sa) /(x*x + y*y + z*z);
  return result;
}
//___________________________________________________________________
double Analysis_mc::derivateR2D(double a, double x, double y){
  double result =0;
  result = a /((x*x + y*y )*(x*x + y*y ));
  return result;
}
//___________________________________________________________________
double Analysis_mc::derivate2_with_sigmaR2D(double a, double sa, double x, double y){
  double result =0;
  result = (a*a)*(sa*sa) /(x*x + y*y );
  return result;
}


//___________________________________________________________________

void Analysis_mc::from_TGraph_to_TH1D (TGraphAsymmErrors *graph, TH1D *histo, int number_point){
    
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
/*
//==================================================================
double Analysis_mc::fakeWeight(const unsigned ind, const int flavors, const double conePt, const double eta, const bool tight, TH2D* frMap, const unsigned lCount){
unsigned nFO = 0;
double fr[4];
for(unsigned l = 0; l < lCount; ++l){
if(!tight[ind[l]]){
fr[nFO] = frMap[flavors[ind[l]]]->GetBinContent(frMap[flavors[ind[l]]]->FindBin(TMath::Min(conePt[l], 99.), fabs(eta[ind[l]])));
++nFO;
}
}
double weight = 1;
for(unsigned f = 0; f < nFO; ++f){
weight *= fr[f]/(1-fr[f]);
}
if(nFO == 2) weight*= -1;
return weight;
}
*/






//==================================================================
double Analysis_mc::FR_factor(TGraphAsymmErrors *fakeRate_mu[3],
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





//___________________________________________________________________
void Analysis_mc::class_os(int event_clas[1], int  flavors_3l[3], int  charge_3l[3]){
    
  int ch_lepton1=charge_3l[0];
  int ch_lepton2=charge_3l[1];
  int ch_lepton3=charge_3l[2];
  int fl_lepton1=flavors_3l[0];
  int fl_lepton2=flavors_3l[1];
  int fl_lepton3=flavors_3l[2];
    
    
  // 1* = eee
  // 2* = emm
  // 3* = eem
  // 4* = mmm
    
    
    
    
    
  event_clas[0]=-1;
    
    
  if (ch_lepton1 == ch_lepton2 && ch_lepton1 == ch_lepton3 && ch_lepton3 == ch_lepton2)   event_clas[0]=-1;
    
    
  if (fl_lepton2 == 0 ||  fl_lepton3  == 0 ||  fl_lepton1== 0) {
    if ((fl_lepton2 + fl_lepton3 + fl_lepton1) == 0 ) event_clas[0] = 10; //e e e
    if ((fl_lepton2 + fl_lepton3 + fl_lepton1) == 1 ) event_clas[0] = 30; //e e mu
    if ((fl_lepton2 + fl_lepton3 + fl_lepton1) == 2 ) {
      if (fl_lepton2 == 1 ||  fl_lepton3  == 1 ||  fl_lepton1== 1) event_clas[0]=20; // e mu mu
    }
  }// at least an electron
  else if (fl_lepton2 == 1 &&  fl_lepton3  == 1 &&  fl_lepton1 == 1) {
    event_clas[0] = 40;
  }
  else {
    event_clas[0] =-1;
  }
    
  if (event_clas[0]  == 30){
    if ((fl_lepton1 == 0 && fl_lepton2 == 0 && fl_lepton3 == 1) ) { // e e mu
      if ((ch_lepton1 + ch_lepton2) == 0) event_clas[0] = 3;
    }
    if (fl_lepton1 == 0 && fl_lepton2 == 1 && fl_lepton3 == 0) { // e mu e
      if ((ch_lepton1 + ch_lepton3) == 0) event_clas[0] = 3;
    }
    if (fl_lepton1 == 1 && fl_lepton2 == 0 && fl_lepton3 == 0 ) { // mu e e
      if ((ch_lepton2 + ch_lepton3) == 0) event_clas[0] = 3;
    }
  }
    
    
  if (event_clas[0]  == 20){
    if ((fl_lepton1 == 1 && fl_lepton2 == 1 && fl_lepton3 == 0) ) { // e e mu
      if ((ch_lepton1 + ch_lepton2) == 0) event_clas[0] = 2;
    }
    if (fl_lepton1 == 1 && fl_lepton2 == 0 && fl_lepton3 == 1 ) { // e mu e
      if ((ch_lepton1 + ch_lepton3) == 0) event_clas[0] = 2;
    }
    if (fl_lepton1 == 0 && fl_lepton2 == 1 && fl_lepton3 == 1 ) { // mu e e
      if  ((ch_lepton2 + ch_lepton3) == 0)event_clas[0] = 2;
    }
  }
    
    
  if (event_clas[0]  == 10 ) {
    if ((ch_lepton1 + ch_lepton2 + ch_lepton3) == 1 || (ch_lepton1 + ch_lepton2 + ch_lepton3) == -1)  event_clas[0] =1;
  }
    
  if (event_clas[0]  == 40 ) {
    if ((ch_lepton1 + ch_lepton2 + ch_lepton3) == 1 || (ch_lepton1 + ch_lepton2 + ch_lepton3) == -1)  event_clas[0] =4;
  }
    
    
}



//___________________________________________________________________
void Analysis_mc::ossf_no_ossf(int kind[1],TLorentzVector pair[3],TLorentzVector leep1, TLorentzVector leep2,TLorentzVector leep3, int  flavors_3l[3], int  charge_3l[3]){
    
  int ch_lepton1=charge_3l[0];
  int ch_lepton2=charge_3l[1];
  int ch_lepton3=charge_3l[2];
  int fl_lepton1=flavors_3l[0];
  int fl_lepton2=flavors_3l[1];
  int fl_lepton3=flavors_3l[2];
    
    
  kind[0] = -1;
    
  if (ch_lepton1 == ch_lepton2 && ch_lepton1 == ch_lepton3 && ch_lepton3 == ch_lepton2)   kind[0] = -1;
    
    
  // OSSF
  if (     ((ch_lepton1 != ch_lepton2)    && (fl_lepton1 == fl_lepton2))  || ((ch_lepton1 != ch_lepton3)   && (fl_lepton1 == fl_lepton3)) || ((ch_lepton2 != ch_lepton3)  && (fl_lepton3 == fl_lepton2)) ){ // ossf
    //cout<<"in function where kind is 1: "<<kind[0]<<endl;
        
        
    kind[0] = 1;
    double i_m[3]={33333,33333,33333};
    double mass_inv=0;
    int index_inv=100;
    double min_mass=999;
    if ((ch_lepton1 != ch_lepton2)  && (fl_lepton1 == fl_lepton2)) i_m[0]= TMath:: Abs((leep1 + leep2).Mag() - 91.1876);
    if ((ch_lepton1 != ch_lepton3)  && (fl_lepton1 == fl_lepton3)) i_m[1]= TMath:: Abs((leep1 + leep3).Mag() - 91.1876);
    if ((ch_lepton2 != ch_lepton3)  && (fl_lepton3 == fl_lepton2)) i_m[2]= TMath:: Abs((leep2 + leep3).Mag() - 91.1876);
    for (int i =0; i < 3; i++){
      if (i_m[i] == 33333) continue;
      mass_inv = i_m[i];
      if (min_mass > mass_inv ){
	min_mass = mass_inv;
	index_inv = i;
      }
    }
    if (index_inv == 0) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[2].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
    }
    if (index_inv == 1) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
            
    }
    if (index_inv == 2) {
      pair[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
            
    }
  }// end ossf
  // No_OSSF
  else if (   ((ch_lepton1 + ch_lepton2) == 0  )  || ((ch_lepton1 + ch_lepton3) == 0   ) || ((ch_lepton3 + ch_lepton2) == 0   )   ){
    //cout<<"in function where kind is 0: "<<kind[0]<<endl;
    kind[0] = 0;
    double i_m[3]={33333,33333,33333};
    double mass_inv=0;
    int index_inv=100;
    double min_mass=999;
    if ((ch_lepton1 != ch_lepton2) ) i_m[0]= TMath:: Abs((leep1 + leep2).Mag() - 91.1876);
    if ((ch_lepton1 != ch_lepton3)  ) i_m[1]= TMath:: Abs((leep1 + leep3).Mag() - 91.1876);
    if ((ch_lepton2 != ch_lepton3) ) i_m[2]= TMath:: Abs((leep2 + leep3).Mag() - 91.1876);
    for (int i =0; i < 3; i++){
      if (i_m[i] == 33333) continue;
      mass_inv = i_m[i];
      if (min_mass > mass_inv ){
	min_mass = mass_inv;
	index_inv = i;
      }
    }
    if (index_inv == 0) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[2].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
    }
    if (index_inv == 1) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
            
    }
    if (index_inv == 2) {
      pair[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
            
    }
        
        
  }//end no-ossf
    
  /*
    cout<< ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   "<<kind<<endl;
    cout<<"1: "<<ch_lepton1<<"  "<<fl_lepton1<< "  "<< leep1.Pt()<< endl;
    cout<<"2: "<<ch_lepton2<<"  "<<fl_lepton2<< "  "<< leep2.Pt()<< endl;
    cout<<"3: "<<ch_lepton3<<"  "<<fl_lepton3<< "  "<< leep3.Pt()<< endl;
    cout<<"pair 1: "<<pair1.Pt()<< endl;
    cout<<"pair 2: "<<pair2.Pt()<< endl;
  */
  //cout<<"in function "<<kind[0]<<endl;
}





//___________________________________________________________________
void Analysis_mc::fr_selection(int number, TLorentzVector lepton_fake_order[3],TLorentzVector leep1, TLorentzVector leep2,TLorentzVector leep3, int index_leptons[3],  int flavor_leptons[3], int origin_leptons[3],int index_3l[3],  int flavor_3l[3], int origin_3l[3]){
    
  lepton_fake_order[0].SetPtEtaPhiE(0,0,0,0);
  lepton_fake_order[1].SetPtEtaPhiE(0,0,0,0);
  lepton_fake_order[2].SetPtEtaPhiE(0,0,0,0);
  for(int i =0; i< 3; i++){
    index_leptons[i]=  -5;
    flavor_leptons[i]= -5;
    origin_leptons[i]= -5;
  }
    
  if (number == 3) {
    lepton_fake_order[0].SetPtEtaPhiE(leep1.Pt(),leep1.Eta(),leep1.Phi(),leep1.E());
    lepton_fake_order[1].SetPtEtaPhiE(leep2.Pt(),leep2.Eta(),leep2.Phi(),leep2.E());
    lepton_fake_order[2].SetPtEtaPhiE(leep3.Pt(),leep3.Eta(),leep3.Phi(),leep3.E());
    index_leptons[0]=index_3l[0];
    index_leptons[1]=index_3l[1];
    index_leptons[2]=index_3l[2];
    flavor_leptons[0]=flavor_3l[0];
    flavor_leptons[1]=flavor_3l[1];
    flavor_leptons[2]=flavor_3l[2];
    origin_leptons[0]=origin_3l[0];
    origin_leptons[1]=origin_3l[1];
    origin_leptons[2]=origin_3l[2];
  }
  if (number == 0) {
    lepton_fake_order[0].SetPtEtaPhiE(leep1.Pt(),leep1.Eta(),leep1.Phi(),leep1.E());
    lepton_fake_order[1].SetPtEtaPhiE(leep2.Pt(),leep2.Eta(),leep2.Phi(),leep2.E());
    lepton_fake_order[2].SetPtEtaPhiE(leep3.Pt(),leep3.Eta(),leep3.Phi(),leep3.E());
    index_leptons[0]=index_3l[0];
    index_leptons[1]=index_3l[1];
    index_leptons[2]=index_3l[2];
    flavor_leptons[0]=flavor_3l[0];
    flavor_leptons[1]=flavor_3l[1];
    flavor_leptons[2]=flavor_3l[2];
    origin_leptons[0]=origin_3l[0];
    origin_leptons[1]=origin_3l[1];
    origin_leptons[2]=origin_3l[2];
  }
    
    
    
}//end fr
