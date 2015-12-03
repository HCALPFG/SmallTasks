#include <algorithm> 
#include <ctime>
#include <fstream>
#include <iomanip> // for setw()
#include <iostream>
#include <sstream>
#include <vector>

#include "TBranch.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TError.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "HFTiming.h"

TString fc_thresh = "test";
TString which_bx = "first"; //Choose "iso", "first", or "last" (in a train)
TString outdir ="/afs/cern.ch/user/r/rbhandar/www/hcal/hftiming/"+which_bx+"/fc"+fc_thresh+"/";

TChain* ch = new TChain("hcalTupleTree/tree");

vector<unsigned int> run  = { //254790,254852,254879,254906,254907,254914,                                           // 2015C
			      //			      256630,256675,256676,256677,256801,256843,256867,256868,256926,256941,257461,257531, // 2015D
			      //			      257599,257613,257645,257682,257722,257723,257735,257751,257805,257816,257819,257968,
			      //			      257969,258157,258158,258159,
			      258177,258656,258694,258702,258703,258705,258706,258712,
			      258713,258742,258745,258749,258750
			      
			      //			      259152,259157,259158,259159,259161,259162,259163,259164,259167,259199,259200,259201, // Totem
			      //			      259207,259208,259236,259237,259351,259352,259384,259385,259388,259399,259429,259431

			      // 259626,259681,259683,259685,259686,259721,259809,259810,259811,259818,259820,259821, // 2015D
			      //			      259822,259862,259890,259891,260373,260427,260431,260532,260533,260534,260536,260538,
			      //			      260541,260575,260576,260577,260593,260627
                             }; //End of 2015 pp collisions

const unsigned int Nrun=run.size();

//int Ethres[4]={fc_thresh.Atoi(),100000,100001,100002};
int Ethres[4]={50,100000,100001,100002};
//int Ethres[4]={50,100,150,200};
//if(fcthresh=="binned"){  Ethres[0]=0; Ethresh[1]=50; Ethresh[2]=100; Ethresh[3]=150;}
//if(fcthresh=="test"){  Ethres[0]=50; Ethresh[1]=100000; Ethresh[2]=100001; Ethresh[3]=100002;}

int HistColor[4]={kBlack,kRed,kBlue,kGreen};
bool DoNorm=false;

// Macro
void HFTiming(){ 
  //Set Up
  gROOT->SetBatch(true);  //Don't show canvases
  gErrorIgnoreLevel=1001; //Don't print message about canvas being saved

  //  ch->Add("/afs/cern.ch/user/j/jaehyeok/public/JetHT_Run2015C-v1_RAW/*.root");
  //  ch->Add("/afs/cern.ch/user/j/jaehyeok/public/JetHT_Run2015D-v1_RAW/*.root");
  ch->Add("/afs/cern.ch/user/j/jaehyeok/public/JetHT_Run2015D-v1_RAW_258177_258750/*.root");
  //  ch->Add("/afs/cern.ch/work/r/rbhandar/public/hcaltuples/ZeroBias_Run2015D-v1_RAW_259152_259431/*.root"); // Totem
  //  ch->Add("/afs/cern.ch/work/r/rbhandar/public/hcaltuples/JetHT_Run2015D-v1_RAW_259626_259891/*.root");
  //  ch->Add("/afs/cern.ch/work/r/rbhandar/public/hcaltuples/JetHT_Run2015D-v1_RAW_260373_260426/*.root");
  //  ch->Add("/afs/cern.ch/work/r/rbhandar/public/hcaltuples/JetHT_Run2015D-v1_RAW_260427_260627/*.root");

  int dir = gSystem->mkdir(outdir,true);
  if(dir==0){
    cout<<"\n[HF Timing] Configuration: Created directory: "<<outdir<<"\n"<<endl;
    gSystem->CopyFile("/afs/cern.ch/user/r/rbhandar/www/index.php",outdir+"index.php");
  }
  else if(dir==-1) 
    cout<<"\n[HF Timing] Configuration: Saving to "<<outdir<<"\n"<<endl;

  //Analysis
  vector< vector<int> > fillScheme = findFillScheme(ch, which_bx, true);

  HFTimingOne(fillScheme, 2, 1, 41, 3, 2);
  HFTimingOne(fillScheme, 2, 1, -41, 3, 2);
}

// Main Looper
void HFTimingOne(vector< vector<int> > goodbxs, int TStoCheck = 2, int TSadjacent = 1, int IETA=999, int IPHI=999, int DEPTH=999){ 
  time_t begtime, endtime;
  time(&begtime);

  init_chain(ch);

  // histograms
  TH1F *h1[Nrun][4], *h2[Nrun][4], *h12[Nrun][4], *h2over12[Nrun][4], *ht[Nrun][4], *h1over2[Nrun][4], *havgtime[Nrun][4]; 
  TProfile *h2profile[Nrun][4];

  for(unsigned int irun=0; irun<Nrun; irun++){
    for(int i=0; i<4; i++){
      int Ethreslow=Ethres[i];
      int Ethreshigh;
      if(i==3) Ethreshigh=999;
      else Ethreshigh=Ethres[i+1];

      h1[irun][i]        = new TH1F(Form("h1_run%i_E%iTo%i", run[irun],         Ethreslow,Ethreshigh),  Form("TS%i",TSadjacent),                                 50, 0, 300);
      h2[irun][i]        = new TH1F(Form("h2_run%i_E%iTo%i", run[irun],         Ethreslow,Ethreshigh),  Form("TS%i",TStoCheck),                                  50, 0, 300);
      h12[irun][i]       = new TH1F(Form("h12_run%i_E%iTo%i", run[irun],        Ethreslow,Ethreshigh),  Form("TS%i+TS%i",TSadjacent,TStoCheck),                  50, 0, 600);
      h2over12[irun][i]  = new TH1F(Form("h2over12_run%i_E%iTo%i", run[irun],    Ethreslow,Ethreshigh),  Form("TS%i/(TS%i+TS%i)",TStoCheck,TSadjacent,TStoCheck), 20, 0, 1);
      h1over2[irun][i]   = new TH1F(Form("h1over2_run%i_E%iTo%i", run[irun],    Ethreslow,Ethreshigh),  Form("TS%i/TS%i",TSadjacent,TStoCheck),                  20, 0, 5);
      havgtime[irun][i]  = new TH1F(Form("havgtime_run%i_E%iTo%i", run[irun],   Ethreslow,Ethreshigh),  "Energy-avg timinig (in the unit of TS)",                30, 0, 3);
      h2profile[irun][i]  = new TProfile(Form("h2profile_run%i_E%iTo%i", run[irun],   Ethreslow,Ethreshigh),  "h2profile",        1000,0,1000,0,1);
        
      h1[irun][i]->Sumw2();
      h2[irun][i]->Sumw2();
      h12[irun][i]->Sumw2();
      h2over12[irun][i]->Sumw2();
      h1over2[irun][i]->Sumw2();
      havgtime[irun][i]->Sumw2();

      h1cosmetic(h1[irun][i],       HistColor[i], 2, 0); 
      h1cosmetic(h2[irun][i],       HistColor[i], 2, 0); 
      h1cosmetic(h12[irun][i],      HistColor[i], 2, 0); 
      h1cosmetic(h2over12[irun][i], HistColor[i], 2, 0); 
      h1cosmetic(h1over2[irun][i],  HistColor[i], 2, 0); 
      h1cosmetic(havgtime[irun][i], HistColor[i], 2, 0); 
    } 
  } 
   
  // main event loop
  unsigned int skiprun = 9999;
  unsigned int nentries = (Int_t)ch->GetEntries();
  for(unsigned int ievent = 0; ievent<nentries; ievent++){
    ch->GetEntry(ievent); 
    
    // Status
    if((ievent<100&&ievent%10==0) || (ievent<1000&&ievent%100==0) || (ievent<10000&&ievent%1000==0) || 
       (ievent<100000&&ievent%10000==0) || (ievent%100000==0) || (ievent+1==nentries)){
      if (isatty(1)) {
	printf("\r[HF Timing] Processing Events: %i / %i (%i%%)",ievent,nentries,static_cast<int>((ievent+1)*100./nentries));
	fflush(stdout);
	if(ievent+1==nentries) printf("\n");
      }
    }        
	
    // Choose bunch crossings
    if(!isGoodBX(run_, bx_, goodbxs)) continue;

    // Skip List
    if((run_>=254231 && run_<=254323)) continue;           // Pre-calibration
    if((run_>=256673 && run_<=256674)) continue;           // No Entries
    if(run_==256842) continue;                             // No Entries
    if(run_==256866) continue;                             // Low stats
    if(run_==256869) continue;                             // Low stats
    if(run_==257614) continue;                             // Low stats
    if(run_==257804) continue;                             // Low stats
    if(run_==258129) continue;                             // Low stats
    if(run_==258136) continue;                             // Low stats
    if((run_>=258211 && run_<=258448)) continue;           // No iso bx
    if(run_==258655) continue;                             // No Entries
    if(run_==258714) continue;                             // Low stats
    if(run_==258741) continue;                             // Low stats
    if(run_==259158) continue;                             // Low stats
    if(run_==259637) continue;                             // Low stats
    if(run_==259813) continue;                             // Low stats
    if(run_==259817) continue;                             // Low stats
    if(run_==260425) continue;                             // No iso bx
    if(run_==260426) continue;                             // No iso bx

    // loop over channels
    for(unsigned int ich=0; ich<HFDigiSubdet_->size(); ich++){ 
      // Selected only interesting channel
      if( IETA!=999 && !(HFDigiIEta_->at(ich)==IETA && HFDigiIPhi_->at(ich)==IPHI && HFDigiDepth_->at(ich)==DEPTH)) continue; 

      // Fill histograms 
      for(int i=0; i<4; i++){  // for different E or Q cuts
	if(i<3 && ((HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))<Ethres[i] || 
		   (HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))>Ethres[i+1]) ) continue;
	if(i==3 && ((HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))<Ethres[i]) ) continue;

	if(0){ // DEBUG
	  cout << Form("[HF Timing] Subdet=%i IEta=%i IPhi=%i Depth=%i RecEnergy=%.3f\n",
		       HFDigiSubdet_->at(ich), 
		       HFDigiIEta_->at(ich),    
		       HFDigiIPhi_->at(ich),    
		       HFDigiDepth_->at(ich),   
		       HFDigiRecEnergy_->at(ich) ); 
	} 
                
	//Make sure this run is part of the run vector above
	int ithisrun=-1;
	for(unsigned int irun=0; irun<Nrun; irun++) 
	  if(run_==run[irun]) ithisrun = irun;

	if(ithisrun==-1) {
	  if(run_!=skiprun){
	    cout<<"\n\e[31m[HF Timing]\e[0m WARNING: Run = "<<run_<<" is being skipped!"<<endl;
	    skiprun = run_;
	    continue;
	  }
	  else		    
	    continue;
	}

	//Filling Histograms
	h1[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TSadjacent),299.999));
	h2[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TStoCheck),299.999)); 
	h12[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck),599.999));
	h1over2[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TSadjacent)/HFDigiFC_->at(ich).at(TStoCheck),4.999));
	h2over12[ithisrun][i]->Fill(Min(HFDigiFC_->at(ich).at(TStoCheck)/(HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck)),9.999));
	h2profile[ithisrun][i]->Fill(ls_, HFDigiFC_->at(ich).at(TStoCheck)/(HFDigiFC_->at(ich).at(TSadjacent)+HFDigiFC_->at(ich).at(TStoCheck))); 
                

	// energy average timing (considering only TS1 and TS2) 
	float avgtime=1; 
	float num =0;
	float deno =0;
	for(int icap=1; icap<3; icap++){ 
	  num = num + icap*(HFDigiFC_->at(ich).at(icap));
	  deno = deno + HFDigiFC_->at(ich).at(icap);
	}
	avgtime = num/deno;
	havgtime[ithisrun][i]->Fill(avgtime);  
      }

    } //for(unsigned int ich=0; ich<HFDigiSubdet_->size(); ich++)
  } //for(unsigned int ievent = 0; ievent<nentries; ievent++)

   
    // Draw 
  if(DoNorm){
    for(unsigned int irun=0; irun<Nrun; irun++){ 
      for(int i=0; i<4; i++){ 
	h1[irun][i]->Scale(1./h1[irun][i]->Integral());
	h2[irun][i]->Scale(1./h2[irun][i]->Integral());
	h12[irun][i]->Scale(1./h12[irun][i]->Integral());
	h1over2[irun][i]->Scale(1./h1over2[irun][i]->Integral());
	h2over12[irun][i]->Scale(1./h2over12[irun][i]->Integral());
	havgtime[irun][i]->Scale(1./havgtime[irun][i]->Integral());
      }
    }
  }

  TLatex *tex_RUN = new TLatex(0.4,0.85,Form("Run number = %i", run_));
  tex_RUN->SetNDC();
  tex_RUN->SetTextSize(0.04);
  tex_RUN->SetLineWidth(2);
  tex_RUN->SetTextAlign(12);
    
  TLatex *tex_ch = new TLatex(0.4,0.8,Form("ieta=%i iphi=%i depth=%i",IETA,IPHI,DEPTH));
  tex_ch->SetNDC();
  tex_ch->SetTextSize(0.04);
  tex_ch->SetLineWidth(2);
  tex_ch->SetTextAlign(12);

  for(unsigned int irun=0; irun<Nrun; irun++){  

    // Legend 
    TLegend *l1 = new TLegend(0.4, 0.5, 0.9, 0.75);
    l1->SetBorderSize(0);
    l1->SetFillColor(0);
    l1->SetFillStyle(0);
    l1->SetFillColor(kWhite);
    l1->SetLineColor(kWhite);
    l1->SetShadowColor(kWhite);
    l1->AddEntry(h1[irun][0],        Form("%i < E_{rechit} < %i GeV ",Ethres[0],Ethres[1]),    "lp");
    l1->AddEntry(h1[irun][1],        Form("%i < E_{rechit} < %i GeV ",Ethres[1],Ethres[2]),    "lp");
    l1->AddEntry(h1[irun][2],        Form("%i < E_{rechit} < %i GeV ",Ethres[2],Ethres[3]),    "lp");
    l1->AddEntry(h1[irun][3],        Form("E_{rechit} > %i GeV",Ethres[3]),    "lp");
        
    // Canvas 
    TCanvas *c = new TCanvas("c", "c", 1200, 800);
    c->Divide(3,2);
    c->cd(1);
    h1[irun][0]->SetMaximum(h1[irun][0]->GetMaximum()*2);
    for(int i=0; i<4; i++) h1[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));  
    tex_RUN->Draw("SAME");
    tex_ch->Draw("SAME");
    l1->Draw("same");
    c->cd(2);
    h2[irun][0]->SetMaximum(h2[irun][0]->GetMaximum()*2);
    for(int i=0; i<4; i++) h2[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
    tex_RUN->Draw("SAME");
    tex_ch->Draw("SAME");
    c->cd(3);
    h12[irun][0]->SetMaximum(h12[irun][0]->GetMaximum()*2);
    for(int i=0; i<4; i++) h12[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
    tex_RUN->Draw("SAME");
    tex_ch->Draw("SAME");
    c->cd(4);
    c->cd(4)->SetLogy(1);
    h1over2[irun][0]->SetMaximum(h1over2[irun][0]->GetMaximum()*50);
    for(int i=0; i<4; i++) h1over2[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
    tex_RUN->Draw("SAME");
    tex_ch->Draw("SAME");
    c->cd(5);
    c->cd(5)->SetLogy(0);
    h2over12[irun][0]->SetMaximum(h2over12[irun][0]->GetMaximum()*2);
    for(int i=0; i<4; i++) h2over12[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
    tex_RUN->Draw("SAME");
    tex_ch->Draw("SAME");
    c->cd(6);
    havgtime[irun][0]->SetMaximum(havgtime[irun][0]->GetMaximum()*2);
    for(int i=0; i<4; i++) havgtime[irun][i]->Draw(Form("E %s",i==0?"":"SAME"));
    tex_RUN->Draw("SAME");
    tex_ch->Draw("SAME");  
        

    //cout << Form("RUN=%i IETA=%i IPHI=%i DEPTH=%i",run[irun],IETA,IPHI,DEPTH) << endl;
    //cout << "Peak : " << havgtime[irun][0]->GetMean() << " +/- " << havgtime[irun][0]->GetRMS()<< endl;

    c->Print(Form(outdir+"/RUN%i_IETA%i_IPHI%i_DEPTH%i.pdf",run[irun],IETA,IPHI,DEPTH));
    delete c;
        
    TCanvas *cprofile = new TCanvas("cprofile", "cprofile", 600, 400);
    cprofile->cd(1);
    h2profile[irun][0]->Draw();
    if(h2profile[irun][0]->GetEntries()<=11) cout<<"\e[31m[HF Timing]\e[0m WARNING: Low stats run = "<<run[irun]<<endl;
    cprofile->Print(Form(outdir+"/Profile_RUN%i_IETA%i_IPHI%i_DEPTH%i.pdf",run[irun],IETA,IPHI,DEPTH));
	
    delete cprofile;
  }

  //
  // Make a summary plot 
  //
  TH1F *hsummary = new TH1F("hsummary",  "hsummary", Nrun, 0, Nrun); 
  for(unsigned int irun=0; irun<Nrun; irun++){   
    hsummary->SetBinContent(irun+1,h2over12[irun][0]->GetMean());
    hsummary->SetBinError(irun+1,h2over12[irun][0]->GetMeanError());  
    //hsummary->SetBinError(irun+1,0);  
    hsummary->GetXaxis()->SetBinLabel(irun+1,Form("%i",(run[irun]-250000)));  
    hsummary->GetXaxis()->SetLabelSize(0.06);  
    hsummary->GetXaxis()->LabelsOption("v");
    hsummary->GetYaxis()->SetTitle("<Q2/(Q1+Q2)>");
    hsummary->GetYaxis()->SetTitleSize(0.06);
    hsummary->GetYaxis()->SetTitleOffset(0.8);
    hsummary->GetYaxis()->SetLabelSize(0.06);  
    hsummary->SetStats(0);
    hsummary->SetTitle(Form("iEta=%i, iPhi=%i, Depth=%i",IETA,IPHI,DEPTH));

  } 
  TCanvas *csum = new TCanvas("csum", "csum", 1200, 400);
  csum->cd(1);
  csum->SetGridy(1);
  hsummary->SetMinimum(0);
  hsummary->SetMaximum(1);
  hsummary->SetMarkerStyle(20);
  hsummary->SetMarkerSize(1);
  hsummary->Draw("ep");
  csum->Print(Form(outdir+"/Summary_IETA%i_IPHI%i_DEPTH%i.pdf",IETA,IPHI,DEPTH));
  csum->Print(Form(outdir+"/Summary_IETA%i_IPHI%i_DEPTH%i.root",IETA,IPHI,DEPTH));

  delete csum;

  makeTextFile(h2over12,IETA);
    
  // clean 
  for(unsigned int irun=0; irun<Nrun; irun++) {
    for(int i=0; i<4; i++){
      delete h1[irun][i]; 
      delete h2[irun][i]; 
      delete h12[irun][i]; 
      delete h1over2[irun][i]; 
      delete h2over12[irun][i]; 
      delete havgtime[irun][i]; 
      delete h2profile[irun][i]; 
    } 
  }
  delete hsummary; 
  
  time(&endtime);
  int seconds = difftime(endtime, begtime);
  float hertz = nentries; hertz /= seconds;
  cout<<"[HF Timing] Processed channel "<<IETA<<":"<<IPHI<<":"<<DEPTH<<" in "<<seconds<<" seconds ("<<hoursMinSec(seconds)<<") for "<<nentries
      <<" events -> "<<roundNumber(hertz,1,1000)<<" kHz, "<<roundNumber(1000,2,hertz)<<" ms per event"<<endl<<endl;
}


vector< vector<int> > findFillScheme(TChain *ch, TString bxtype="iso", bool printbx=false){
  time_t begtime, endtime;
  time(&begtime);

  if(bxtype!="iso" && bxtype!="first" && bxtype!="last")
    cout<<"Fill scheme type invalid. bxtype must be \"iso\", \"first\", or \"last\""<<endl;  

  unsigned int   run_ = 0;
  ch->SetBranchAddress("run", &run_);
  unsigned int   bx_ = 0;
  ch->SetBranchAddress("bx", &bx_);
  vector<vector<float> >   *HFDigiFC_ = 0;
  ch->SetBranchAddress("HFDigiFC", &HFDigiFC_);

  //Initialize vector of TH2Ds the size of Nrun
  vector<TH2D*> h_bxs;
  for(unsigned int irun=0; irun<Nrun; irun++)
    h_bxs.push_back(new TH2D(Form("run[%i]",irun),"",3563,0,3563,1,0,4000));

  //Loop over entries to fill h_bxs
  unsigned int currentrun = 666;
  unsigned int runidx = 99999;
  unsigned int nentries = ch->GetEntries();
  for(unsigned int iloop=0; iloop<nentries; iloop++){
    ch->GetEntry(iloop);   

    if((iloop<100&&iloop%10==0) || (iloop<1000&&iloop%100==0) || (iloop<10000&&iloop%1000==0) || 
       (iloop<100000&&iloop%10000==0) || (iloop%100000==0) || (iloop+1==nentries)){
      if (isatty(1)) {
	printf("\r[HF Timing] Analyzing fill scheme: %i / %i (%i%%)",iloop,nentries,static_cast<int>((iloop)*100./nentries));
	fflush(stdout);
	if(iloop+1==nentries){
	  printf("\r[HF Timing] Analyzing fill scheme: %i / %i (%i%%)",iloop+1,nentries,static_cast<int>((iloop+1.)*100./nentries));
	  fflush(stdout);
	  printf("\n");
	} 
      }
    }        

    //Get total charge in channels
    double totalFC12=0;
    for(unsigned int ich=0; ich<HFDigiFC_->size(); ich++)
      totalFC12 += HFDigiFC_->at(ich).at(1) + HFDigiFC_->at(ich).at(2);      

    //Get index of current run only if it is a new run
    if(run_ != currentrun){
      runidx = find(run.begin(),run.end(),run_)-run.begin();
      currentrun = run_;
    }
    if(runidx!=Nrun)
      if(totalFC12>100)
	h_bxs[runidx]->Fill(bx_,totalFC12);
  }

  // Loop over histograms and get bxs
  vector< vector<int> > selectedbxs(Nrun, vector<int>());
  for(unsigned int ihist=0; ihist<h_bxs.size(); ihist++){

    if((ihist%1)==0){
      if (isatty(1)) {
	printf("\r[HF Timing] Extracting fill scheme: %i / %i (%i%%)",ihist,(int)h_bxs.size(),static_cast<int>((ihist)*100./h_bxs.size()));
	fflush(stdout);
	if(ihist+1>=h_bxs.size()){
	  printf("\r[HF Timing] Extracting fill scheme: %i / %i (%i%%)",ihist+1,(int)h_bxs.size(),static_cast<int>((ihist+1)*100./h_bxs.size()));
	  fflush(stdout);
	  printf("\n");
	}
      }
    }
    
    for(int ibx=1; ibx<=3563; ibx++){   // account for th2 binning starting at 1
      
      bool isfirst=false, islast=false;
      
      double intbx = h_bxs[ihist]->Integral(ibx,ibx,1,3563);
      double intbx_prev1 = h_bxs[ihist]->Integral(ibx-1,ibx-1,1,3563);
      double intbx_prev2 = h_bxs[ihist]->Integral(ibx-2,ibx-2,1,3563);
      double intbx_next1 = h_bxs[ihist]->Integral(ibx+1,ibx+1,1,3563);
      double intbx_next2 = h_bxs[ihist]->Integral(ibx+2,ibx+2,1,3563);

      // Three methods to try
      // 1. Individual comparison: Charge[i]/Charge[i+1]<0.1 && Charge[i]/Charge[i+2]<0.1 && ... 
      // 2. Summed comparison: Charge[i]/(Charge[i+1]+Charge[i+2]+...)<0.1
      // 3. Averaged comparison: (n*Charge[i])/(Charge[i+1]+Charge[i+2]+...Charge[i+n])<0.1

      if(intbx_prev1/intbx<0.01 && intbx_prev2/intbx<0.01) //if prev two bunches empty
      	isfirst=true;
      if(intbx_next1/intbx<0.01 && intbx_next2/intbx<0.01) //if next two bunches empty
      	islast=true;

      if(isfirst&&islast && bxtype=="iso")
	selectedbxs[ihist].push_back(ibx-1);
      else if(isfirst && !islast && bxtype=="first")
	selectedbxs[ihist].push_back(ibx-1);
      else if(islast && !isfirst && bxtype=="last")
	selectedbxs[ihist].push_back(ibx-1);
    }
  }
    
  if(printbx){
    ofstream file(outdir+"fillscheme.txt");
    for(unsigned int i=0;i<selectedbxs.size();i++){
      file<<"Run: "<<run[i]<<", Nbxs: "<<selectedbxs[i].size()<<endl;
      for(unsigned int j=0;j<selectedbxs[i].size();j++){
	file<<"bx["<<i<<"]["<<j<<"] = "<<selectedbxs[i][j]<<endl;
      }
    }
    file.close();
    cout<<"[HF Timing] Writing bxs: fillscheme.txt"<<endl;      
  }
  
  time(&endtime);
  int seconds = difftime(endtime, begtime);
  float hertz = nentries; hertz /= seconds;

  cout<<"[HF Timing] Found fill scheme in "<<seconds<<" seconds ("<<hoursMinSec(seconds)<<") for "<<addCommas(nentries)
      <<" events -> "<<roundNumber(hertz,1,1000)<<" kHz, "<<roundNumber(1000,2,hertz)<<" ms per event"<<endl<<endl;

  return selectedbxs;
}

bool isGoodBX(int run_num, int bx, vector< vector<int> > fillscheme){
  
  bool isgood=false;
  unsigned int runidx = find(run.begin(),run.end(),run_num)-run.begin();

  if(runidx != run.size()){
    for(unsigned int ibx=0; ibx<fillscheme[runidx].size(); ibx++){
      if(bx==fillscheme[runidx][ibx])
	isgood=true;
    }
  }
  
  return isgood;
}


//Print all the values for the summary plot
void makeTextFile(TH1F* timing[][4], int ieta){

  TString name = "HFp_Timing";
  if(ieta<0) name = "HFm_Timing";
  name += ".txt";

  ofstream file(outdir+name);
  for(unsigned int i=0; i<Nrun; i++){
    //run declared up top
    file<<run[i]<<", "<<timing[i][0]->GetMean()<<", "<<timing[i][0]->GetMeanError()<<endl;
  }
  file.close();

  cout<<"[HF Timing] Writing values: "<<name<<endl;
}
