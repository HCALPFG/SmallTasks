// Code to produce 50/50 curves for specified channels.
#include <vector>
#include <iostream>
#include "TH1D.h"
#include "TChain.h"
#include "TTree.h"
#include <TError.h>
#include "TString.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TSystem.h"

namespace {
  TString fc_thresh("50");
  TString outdir("/afs/cern.ch/user/r/rbhandar/www/hcal/ffcurves/fc"+fc_thresh+"/");

  const int noffset = 6;
  double offset[2][noffset] = { {-1,0,1,2,3,4},    // x axis
				{ 0,0,0,0,0,0} };  // error
  int runlist[noffset] = {250902,250897,250897,250897,250897,250902};
  int lumilist[2][noffset] = { {73 ,11,72 ,140,195,1 },   // First lumi
			       {103,64,132,187,230,40} }; // Last lumi

  TChain ff("hcalTupleTree/tree");
}

void ffcurveone(TString ieta, TString iphi, TString depth);
void getMeanErr(TString histname, TString var, TString cut, double mean_err[][noffset], int idx);


void ffcurve(){
  gROOT->SetBatch(true);  //Don't show canvases
  gErrorIgnoreLevel=1001; //Don't print message about canvas being saved

  ff.Add("~jaehyeok/public/PFGtuple5050Curve/*.root");
  const int nentries = ff.GetEntries();

  int dir = gSystem->mkdir(outdir,true);
  if(dir==0) cout<<"Created directory: "<<outdir<<endl;

  ffcurveone("41","3","2");
  ffcurveone("-41","3","2");
}

void ffcurveone(TString ieta, TString iphi, TString depth){

  double mean_charge212[2][noffset]; //0th array for mean. 1st array for error.

  for(int ioffset=0; ioffset<noffset; ioffset++){
    
    TString hist("h_"+ieta+"_"+iphi+"_"+depth+"_"+to_string((int)offset[0][ioffset]));

    TString sign("1");
    if(ieta.Contains("-")) sign = "0";

    TString pass_runlist = "(run=="+to_string(runlist[ioffset])+"&&ls>="+to_string(lumilist[0][ioffset])+"&&ls<="+to_string(lumilist[1][ioffset])+")";
    TString pass_channel = "(HFDigiIEta=="+ieta+"&&HFDigiIPhi=="+iphi+"&&HFDigiDepth=="+depth+")";
    TString pass_fcthresh = "((HFDigiFC["+sign+"][1]+HFDigiFC["+sign+"][2])>"+fc_thresh+")";

    TString variable = "HFDigiFC["+sign+"][2]/(HFDigiFC["+sign+"][1]+HFDigiFC["+sign+"][2])";

    // Get mean and associated errors for the offset
    getMeanErr(hist, variable, pass_runlist+"&&"+pass_channel+"&&"+pass_fcthresh, mean_charge212, ioffset);
  }
  
  TCanvas can;
  TString curvename("curve5050p");
  if(ieta.Contains("-")) curvename = "curve5050m";

  TGraphErrors *h_curve5050 = new TGraphErrors(noffset,offset[0],mean_charge212[0], offset[1], mean_charge212[1]);
  
  //Aesthetics
  h_curve5050->SetMarkerStyle(20);
  h_curve5050->SetTitle("iEta="+ieta+", iPhi="+iphi+", Depth="+depth);
  h_curve5050->GetXaxis()->SetLabelSize(0.0425);
  h_curve5050->GetXaxis()->SetTitle("Timing Offset [ns]");
  h_curve5050->GetXaxis()->SetTitleOffset(1.15);
  h_curve5050->GetXaxis()->SetTitleSize(0.0425);
  h_curve5050->GetYaxis()->SetLabelSize(0.0425);
  h_curve5050->GetYaxis()->SetTitle("<Q2/(Q1+Q2)>");
  h_curve5050->GetYaxis()->SetTitleOffset(1.15);
  h_curve5050->GetYaxis()->SetTitleSize(0.0425);
  
  h_curve5050->Draw("AP");
  can.SaveAs(outdir+curvename+".pdf");    
}


void getMeanErr(TString histname, TString var, TString cut, double mean_err[][noffset], int idx){

  TCanvas can;
  TH1D temp(histname, "", 20, 0, 1);
  temp.Sumw2();
  ff.Project(histname,var,cut);

  mean_err[0][idx] = temp.GetMean();
  mean_err[1][idx] = temp.GetMeanError();

  //Aesthetics
  temp.SetMarkerStyle(20);
  temp.SetTitle(histname);
  temp.GetXaxis()->SetLabelSize(0.0425);
  temp.GetXaxis()->SetTitle("Q2/(Q1+Q2)");
  temp.GetXaxis()->SetTitleOffset(1.15);
  temp.GetXaxis()->SetTitleSize(0.0425);
  temp.GetYaxis()->SetLabelSize(0.0425);
  temp.GetYaxis()->SetTitle("Entries/0.05 fC");
  temp.GetYaxis()->SetTitleOffset(1.15);
  temp.GetYaxis()->SetTitleSize(0.0425);

  temp.Draw("hist");
  can.SaveAs(outdir+histname+".pdf");
}
