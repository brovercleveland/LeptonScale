#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TLorentzVector.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <vector>                   // STL vector class
#include <utility>                  // For STL pair class
#include <map>                      // STL map class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "Math/LorentzVector.h"     // 4-vector class

#include "CPlot.hh"                 // helper class for plots
#include "MitStyleRemix.hh"         // style settings for drawing
#endif

// RooFit headers
#include "RooRealVar.h"
#include "RooLinearVar.h"
#include "RooFormulaVar.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooFitResult.h"


//typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LorentzVector;


//=== MAIN MACRO ================================================================================================= 

void EleScale(const Int_t opt=3) {
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  TString outputDir = "";
  TString pufname = "";
  UInt_t runNumberMin=0;
  UInt_t runNumberMax=0;
  
  if(opt==0) {
    outputDir = "Zee_Full2011";
    pufname   = "PileupReweighting.Fall11_To_Full2011.root";
    runNumberMin = 0;
    runNumberMax = 999999;
  
  } else if(opt==1) {
    outputDir = "Zee_Run2011A";
    pufname   = "PileupReweighting.Fall11_To_Run2011A.root";
    runNumberMin = 0;
    runNumberMax = 173692;
  
  } else if(opt==2) {
    outputDir = "Zee_Run2011B";
    pufname   = "PileupReweighting.Fall11_To_Run2011B.root";
    runNumberMin = 173693;
    runNumberMax = 999999;

  } else if(opt==3) {
    outputDir = "Zee_Run2012";
    pufname   = "PileupReweighting.Fall11_To_Run2011B.root";
    runNumberMin = 173693;
    runNumberMax = 999999;
  
  
  } else {
    cout << "Invalid option! Aborting..." << endl;
    return;
  }  
  
  vector<TString> infilenamev;
  infilenamev.push_back("data_select.raw.root");  // data
  infilenamev.push_back("zee_select.root");       // MC

  vector<TString> treenamev;
  treenamev.push_back("eleTree_DATA");  // data
  treenamev.push_back("eleTree_DYJets");       // MC

  TString infilename = "/tthome/bpollack/CMSSW_5_3_11_patch6/src/HZG_Analyzer/HiggsZGAnalyzer/eleFiles/eleFile_EE2012ABCD_07-14-14_ele.root";
  
  const Double_t MASS_LOW  = 60;
  const Double_t MASS_HIGH = 120;
  const Double_t PT_CUT    = 10;
  const Double_t ETA_CUT   = 2.5;
  const Double_t ELE_MASS  = 0.000511;  
  
  vector<pair<Double_t,Double_t> > scEta_limits;
  scEta_limits.push_back(make_pair(0.0,0.4));
  scEta_limits.push_back(make_pair(0.4,0.8));
  scEta_limits.push_back(make_pair(0.8,1.2));
  scEta_limits.push_back(make_pair(1.2,1.4442));
  //scEta_limits.push_back(make_pair(1.566,2.1));
  scEta_limits.push_back(make_pair(1.4442,2.1));
  scEta_limits.push_back(make_pair(2.1,2.5));

  vector<pair<Double_t,Double_t> > pT_limits;
  pT_limits.push_back(make_pair(10.0,20.0));
  pT_limits.push_back(make_pair(20.0,30.0));
  pT_limits.push_back(make_pair(30.0,100000.0));

  CPlot::sOutDir = outputDir;
  
  const TString format("png");
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
   
  enum { eData=0, eMC };
  
  TFile *pufile = new TFile(pufname); assert(pufile);
  TH1D  *puWeights = (TH1D*)pufile->Get("puWeights");
  
  char hname[100];

  vector<TH1D*> hMCv, hDatav;  
  int nMap[6][6][3][3];

  int n =0;
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=0; jbin<scEta_limits.size(); jbin++) {
      for(UInt_t mbin=0; mbin<pT_limits.size(); mbin++) {
        for(UInt_t nbin=mbin; nbin<pT_limits.size(); nbin++) {
          sprintf(hname,"mc_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          hMCv.push_back(new TH1D(hname,"",100,MASS_LOW,MASS_HIGH));
          hMCv.back()->Sumw2();
          
          sprintf(hname,"data_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          hDatav.push_back(new TH1D(hname,"",100,MASS_LOW,MASS_HIGH));
          hDatav.back()->Sumw2();
          nMap[ibin][jbin][mbin][nbin] = n;
          n++;
        }
      }
    }
  }
  
  //
  // Declare ntuple variables
  //
  UInt_t runNumber, lumiSection;
  ULong64_t  eventNumber;
  UInt_t  npv, npu;
  Int_t   q1, q2;
  TLorentzVector *dilep=0, *ele1=0, *ele2=0;
  Float_t SCetaEl1, SCetaEl2, eventWeight;


  cout << "Processing " << infilename << "..." << endl;
  TFile *infile = new TFile(infilename); assert(infile);
  
  for(UInt_t itree=0; itree<treenamev.size(); itree++) {
    cout << "Processing " << treenamev[itree] << "..." << endl;
    TTree *intree = (TTree*)infile->Get(treenamev[itree]); assert(intree);
  
    intree->SetBranchAddress("runNumber", &runNumber);   // event run number
    intree->SetBranchAddress("lumiSection",&lumiSection);  // event lumi section
    intree->SetBranchAddress("eventNumber", &eventNumber);   // event number
    intree->SetBranchAddress("ele1",   &ele1);     // lead lepton 4-vector
    intree->SetBranchAddress("ele2",   &ele2);     // trail lepton 4-vector
    intree->SetBranchAddress("SCetaEl1",    &SCetaEl1);	   // lead Supercluster eta 
    intree->SetBranchAddress("SCetaEl2",    &SCetaEl2);	   // trail Supercluster eta 
    intree->SetBranchAddress("eventWeight",    &eventWeight);	   // PU weight 
  
    int nEvents = 0;
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      nEvents++;
      if (nEvents%1000000 == 0) cout<<"Events passed: "<<nEvents<<endl; 
      
        
      //if(q1 == q2) continue;
      //if(dilep->M()	  < MASS_LOW)  continue;
      //if(dilep->M()	  > MASS_HIGH) continue;
      if(ele1->Pt()	  < PT_CUT)    continue;
      if(ele2->Pt()	  < PT_CUT)    continue;
      if(fabs(SCetaEl1) > ETA_CUT)   continue;      
      if(fabs(SCetaEl2) > ETA_CUT)   continue;
    
      TLorentzVector vele1 = (*ele1);
      TLorentzVector vele2 = (*ele2);
      TLorentzVector vDilep = vele1 + vele2;
    
      Int_t bin1=-1, bin2=-1;
      for(UInt_t i=0; i<scEta_limits.size(); i++) {
        Double_t etalow  = scEta_limits.at(i).first;
        Double_t etahigh = scEta_limits.at(i).second;
        if(fabs(SCetaEl1)>=etalow && fabs(SCetaEl1)<=etahigh) bin1=i;
        if(fabs(SCetaEl2)>=etalow && fabs(SCetaEl2)<=etahigh) bin2=i;
      }
      
      assert(bin1>=0);
      assert(bin2>=0);
      Int_t ibin= bin1; 
      Int_t jbin= bin2;

      Int_t bin3=-1, bin4=-1;
      for(UInt_t i=0; i<pT_limits.size(); i++) {
        Double_t ptlow  = pT_limits.at(i).first;
        Double_t pthigh = pT_limits.at(i).second;
        if(fabs(ele1->Pt())>=ptlow && fabs(ele1->Pt())<=pthigh) bin3=i;
        if(fabs(ele2->Pt())>=ptlow && fabs(ele2->Pt())<=pthigh) bin4=i;
      }

      assert(bin3>=0);
      assert(bin4>=0);

      Int_t mbin= (bin3<=bin4) ? bin3 : bin4;
      Int_t nbin= (bin3<=bin4) ? bin4 : bin3;
      
      /*
      UInt_t n=jbin-ibin;
      for(Int_t k=0; k<ibin; k++)
        n+=(scEta_limits.size()-k);
      
      if(itree==eData) hDatav[n]->Fill(vDilep.M(),eventWeight);
      if(itree==eMC)   hMCv[n]->Fill(vDilep.M(),eventWeight);
      */
      if(itree==eData) hDatav[nMap[ibin][jbin][mbin][nbin]]->Fill(vDilep.M(),eventWeight);
      if(itree==eMC)   hMCv[nMap[ibin][jbin][mbin][nbin]]->Fill(vDilep.M(),eventWeight);
    }  
    intree=0;
  }
  infile=0;
  delete infile;
  
  //
  // Fit for energy scale and resolution corrections
  //
  char vname[100];  // buffer for RooFit object names
  
  char pname[100];
  char str1[100];
  char str2[100];
  char str3[100];
  char str4[100];
  TCanvas *c = MakeCanvas("c","c",800,600);
  
  // Dummy histograms for TLegend (I can't figure out how to properly pass RooFit objects...)
  TH1D *hDummyData = new TH1D("hDummyData","",0,0,10);
  hDummyData->SetMarkerStyle(kFullCircle);
  hDummyData->SetMarkerSize(0.9);
  TH1D *hDummyMC = new TH1D("hDummyMC","",0,0,10);
  hDummyMC->SetLineColor(kBlue);
  hDummyMC->SetFillColor(kBlue);
  hDummyMC->SetFillStyle(3002);  
  TH1D *hDummyFit = new TH1D("hDummyFit","",0,0,10);
  hDummyFit->SetLineColor(kGreen+2);

  RooRealVar mass("mass","M_{ee}",60.0,120.0,"GeV") ;
  mass.setBins(1600,"cache");
  
  RooRealVar massmc("massmc","massmc",0.0,150.0,"GeV");  // mass variable for building MC template

  RooCategory zscEta_cat("zscEta_cat","zscEta_cat");  
  RooSimultaneous combscalefit("combscalefit","combscalefit",zscEta_cat);
  
  map<string,TH1*> hmap;  // Mapping of category labels and data histograms
  
  //RooArgList scalebins;   // List of RooRealVars storing per bin energy scale corrections
  //RooArgList sigmabins;   // List of RooRealVars storing per bin energy resolution corrections
  RooArgList scalebins[3];   // List of RooRealVars storing per bin energy scale corrections
  RooArgList sigmabins[3];   // List of RooRealVars storing per bin energy resolution corrections
  Int_t intOrder = 1;     // Interpolation order for       
  
  for(UInt_t mbin=0; mbin<pT_limits.size(); mbin++) {
    for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
      sprintf(vname,"scale_%i_%i",ibin,mbin);
      RooRealVar *scalebinned = new RooRealVar(vname,vname,1.0,0.8,1.2);      
      scalebins[mbin].add(*scalebinned);
      
      sprintf(vname,"sigma_%i_%i",ibin,mbin);
      RooRealVar *sigmabinned = new RooRealVar(vname,vname,0.8,0.0,5.0);      
      sigmabins[mbin].add(*sigmabinned);
    }
  }
    
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=0; jbin<scEta_limits.size(); jbin++) {
      for(UInt_t mbin=0; mbin<pT_limits.size(); mbin++) {
        for(UInt_t nbin=mbin; nbin<pT_limits.size(); nbin++) {
      //UInt_t n=jbin-ibin;
      //for(UInt_t k=0; k<ibin; k++)
      //  n+=(scEta_limits.size()-k);
      //
      
          cout<<"bins "<<ibin<<jbin<<mbin<<nbin<<endl;
          sprintf(vname,"masslinearshifted_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          RooFormulaVar *masslinearshifted = new RooFormulaVar(vname,vname,"sqrt(@0*@1)",RooArgList(*scalebins[mbin].at(ibin),*scalebins[nbin].at(jbin)));

          sprintf(vname,"massshiftedscEta_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          RooLinearVar *massshiftedscEta = new RooLinearVar(vname,vname,mass,*masslinearshifted,RooConst(0.0));        

          // MC-based template
          sprintf(vname,"zmassmcscEta_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          RooDataHist *zmassmcscEta = new RooDataHist(vname,vname,RooArgList(massmc),hMCv[nMap[ibin][jbin][mbin][nbin]]);      
          sprintf(vname,"masstemplatescEta_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          RooHistPdf *masstemplatescEta = new RooHistPdf(vname,vname,RooArgList(*massshiftedscEta),RooArgList(massmc),*zmassmcscEta,intOrder);             

          // Gaussian smearing function
          sprintf(vname,"sigmascEta_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          RooFormulaVar *sigmascEta = new RooFormulaVar(vname,vname,"TMath::Max(0.01,sqrt(@0*@0+@1*@1))",RooArgList(*sigmabins[mbin].at(ibin),*sigmabins[nbin].at(jbin)));        
          sprintf(vname,"resscEta_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          RooGaussian *resscEta = new RooGaussian(vname,vname,mass,RooConst(0.),*sigmascEta);

          // Fit model: MC-template convoluted with Gaussian
          sprintf(vname,"fftscEta_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          RooFFTConvPdf *fftscEta = new RooFFTConvPdf(vname,vname,mass,*masstemplatescEta,*resscEta);
          fftscEta->setBufferStrategy(RooFFTConvPdf::Flat);

          // Add bin as a category
          char zscEta_catname[100];
          sprintf(zscEta_catname,"zscEta_cat_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          zscEta_cat.defineType(zscEta_catname); 
          zscEta_cat.setLabel(zscEta_catname);
          hmap.insert(pair<string,TH1*>(zscEta_catname,hDatav[nMap[ibin][jbin][mbin][nbin]]));      
          combscalefit.addPdf(*fftscEta,zscEta_catname);

        }
      }
    }
  }
  
  // perform fit
  RooDataHist zdatascEta_comb("zdatascEta_comb","zdatascEta_comb",RooArgList(mass),zscEta_cat,hmap,1.0/(hDatav.size()));
  combscalefit.fitTo(zdatascEta_comb,PrintEvalErrors(kFALSE),Minos(kFALSE),Strategy(0),Minimizer("Minuit2",""));

  Double_t xval[pT_limits.size()][scEta_limits.size()];
  Double_t xerr[pT_limits.size()][scEta_limits.size()];
  Double_t scaleDatatoMC[pT_limits.size()][scEta_limits.size()];
  Double_t scaleDatatoMCerr[pT_limits.size()][scEta_limits.size()];
  Double_t scaleMCtoData[pT_limits.size()][scEta_limits.size()];
  Double_t scaleMCtoDataerr[pT_limits.size()][scEta_limits.size()];
  Double_t sigmaMCtoData[pT_limits.size()][scEta_limits.size()];
  Double_t sigmaMCtoDataerr[pT_limits.size()][scEta_limits.size()];
  
  for(UInt_t mbin=0; mbin<pT_limits.size(); mbin++) {
    for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
        Double_t etalow  = scEta_limits.at(ibin).first;
        Double_t etahigh = scEta_limits.at(ibin).second;
        
        xval[mbin][ibin] = 0.5*(etahigh+etalow);
        xerr[mbin][ibin] = 0.5*(etahigh-etalow);
        
        scaleDatatoMC[mbin][ibin]    = ((RooRealVar*)scalebins[mbin].at(ibin))->getVal();
        scaleDatatoMCerr[mbin][ibin] = ((RooRealVar*)scalebins[mbin].at(ibin))->getError();

        scaleMCtoData[mbin][ibin]    = 1.0/scaleDatatoMC[ibin][mbin];
        scaleMCtoDataerr[mbin][ibin] = scaleDatatoMCerr[ibin][mbin]/scaleDatatoMC[ibin][mbin]/scaleDatatoMC[ibin][mbin];
        
        sigmaMCtoData[mbin][ibin]    = ((RooRealVar*)sigmabins[mbin].at(ibin))->getVal();
        sigmaMCtoDataerr[mbin][ibin] = ((RooRealVar*)sigmabins[mbin].at(ibin))->getError();
    }
    TGraphErrors *grScaleDatatoMC = new TGraphErrors(scEta_limits.size(),xval[mbin],scaleDatatoMC[mbin],xerr[mbin],scaleDatatoMCerr[mbin]);
    TGraphErrors *grScaleMCtoData = new TGraphErrors(scEta_limits.size(),xval[mbin],scaleMCtoData[mbin],xerr[mbin],scaleMCtoDataerr[mbin]);
    TGraphErrors *grSigmaMCtoData = new TGraphErrors(scEta_limits.size(),xval[mbin],sigmaMCtoData[mbin],xerr[mbin],sigmaMCtoDataerr[mbin]);
    

    char ptrange[100];
    sprintf(ptrange,"ele_scale_datatomc_pT_%i",mbin);
    CPlot plotScale1(ptrange,"","Supercluster |#eta|","Data scale correction");
    plotScale1.AddGraph(grScaleDatatoMC,"",kBlue);
    plotScale1.SetYRange(0.99,1.01);
    plotScale1.AddLine(0,1,2.75,1,kBlack,7);
    plotScale1.Draw(c,kTRUE,format);
    
    sprintf(ptrange,"ele_scale_mctodata_%i",mbin);
    CPlot plotScale2(ptrange,"","Supercluster |#eta|","MC#rightarrowData scale correction");
    plotScale2.AddGraph(grScaleMCtoData,"",kBlue);
    plotScale2.SetYRange(0.99,1.01);
    plotScale2.AddLine(0,1,2.75,1,kBlack,7);
    plotScale2.Draw(c,kTRUE,format);

    sprintf(ptrange,"ele_res_mctodata_%i",mbin);
    CPlot plotRes(ptrange,"","Supercluster |#eta|","MC#rightarrowData additional smear [GeV]");
    plotRes.AddGraph(grSigmaMCtoData,"",kBlue);
    plotRes.SetYRange(0,1.6);
    plotRes.Draw(c,kTRUE,format);
  }

  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    for(UInt_t jbin=0; jbin<scEta_limits.size(); jbin++) {
      for(UInt_t mbin=0; mbin<pT_limits.size(); mbin++) {
        for(UInt_t nbin=mbin; nbin<pT_limits.size(); nbin++) {
          //UInt_t n=jbin-ibin;
          //for(UInt_t k=0; k<ibin; k++)
          // n+=(scEta_limits.size()-k);

          // Post-fit plot
          RooPlot *frame = mass.frame();
          char catname[100]; sprintf(catname,"zscEta_cat_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          char cutstr[100];  sprintf(cutstr,"zscEta_cat==zscEta_cat::%s",catname); 
          RooDataHist zmc(catname,catname,RooArgList(mass),hMCv[nMap[ibin][jbin][mbin][nbin]]);
          RooHistPdf mctemplate(catname,catname,RooArgList(mass),zmc,intOrder);          
          mctemplate.plotOn(frame,LineColor(kBlue),LineWidth(1),Normalization(hDatav[nMap[ibin][jbin][mbin][nbin]]->GetEntries()));            
          mctemplate.plotOn(frame,LineColor(kBlue),FillColor(kBlue),FillStyle(3002),DrawOption("F"),Normalization(hDatav[nMap[ibin][jbin][mbin][nbin]]->GetEntries()));      
          zdatascEta_comb.plotOn(frame,Cut(cutstr),MarkerStyle(kFullCircle),MarkerSize(1.0),DrawOption("ZP"));
          combscalefit.plotOn(frame,Slice(zscEta_cat,catname),ProjWData(RooArgSet(mass,catname),zdatascEta_comb),
              LineColor(kGreen+2));           
          sprintf(pname,"postfit_%i_%i_%i_%i",ibin,jbin,mbin,nbin);
          sprintf(str1,"ele1Eta: [%.1f, %.1f]",scEta_limits.at(ibin).first,scEta_limits.at(ibin).second);
          sprintf(str2,"ele2Eta: [%.1f, %.1f]",scEta_limits.at(jbin).first,scEta_limits.at(jbin).second);
          sprintf(str3,"ele1Pt: [%.1f, %.1f]",pT_limits.at(mbin).first,pT_limits.at(mbin).second);
          sprintf(str4,"ele2Pt: [%.1f, %.1f]",pT_limits.at(nbin).first,pT_limits.at(nbin).second);
          CPlot plot(pname,frame,"","m(e^{+}e^{-}) [GeV/c^{2}]","Events / 0.6 GeV/c^{2}");
          plot.AddTextBox(str1,0.21,0.80,0.45,0.87,0,kBlack,-1);
          plot.AddTextBox(str2,0.21,0.73,0.45,0.80,0,kBlack,-1);
          plot.AddTextBox(str3,0.21,0.66,0.45,0.80,0,kBlack,-1);
          plot.AddTextBox(str4,0.21,0.51,0.45,0.80,0,kBlack,-1);
          plot.SetLegend(0.75,0.64,0.93,0.88);
          plot.GetLegend()->AddEntry(hDummyData,"Data","PL");
          plot.GetLegend()->AddEntry(hDummyMC,"Sim","FL");
          plot.GetLegend()->AddEntry(hDummyFit,"Fit","L");
          plot.Draw(c,kTRUE,format);
        }
      }
    }
  }

      
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
   
  ofstream txtfile;
  char txtfname[100];    
  sprintf(txtfname,"%s/summary.txt",outputDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  txtfile << "  MC->Data scale correction" << endl;
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    Double_t etalow  = scEta_limits.at(ibin).first;
    Double_t etahigh = scEta_limits.at(ibin).second;
    txtfile << "$" << etalow << " < |\\eta| < " << etahigh << "$ & ";
    //txtfile << "$" << ((RooRealVar*)scalebins.at(ibin))->getVal() << "$ \\pm $" << ((RooRealVar*)scalebins.at(ibin))->getError() << "$ \\\\" << endl;
    txtfile << "$" << 1.0/((RooRealVar*)scalebins[3].at(ibin))->getVal() << "$ \\pm $" << "-1" << "$ \\\\" << endl;
  }
  txtfile << endl;
  txtfile << "  MC->Data resolution correction [GeV]" << endl;
  for(UInt_t ibin=0; ibin<scEta_limits.size(); ibin++) {
    Double_t etalow  = scEta_limits.at(ibin).first;
    Double_t etahigh = scEta_limits.at(ibin).second;
    txtfile << etalow << " < |\\eta| < " << etahigh << " & ";
    txtfile << "$" << ((RooRealVar*)sigmabins[3].at(ibin))->getVal() << "$ \\pm $" << ((RooRealVar*)sigmabins[3].at(ibin))->getError() << "$ \\\\" << endl;
  }
  txtfile.close();
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
}
