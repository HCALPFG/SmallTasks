#include "TChain.h"
#include <vector>

// Branches
unsigned int   run_;
unsigned int   ls_;
unsigned int   bx_;
vector<int>   *HFDigiSubdet_;
vector<int>   *HFDigiIEta_;
vector<int>   *HFDigiIPhi_;
vector<int>   *HFDigiDepth_;
vector<float>   *HFDigiRecEnergy_;
vector<vector<float> >   *HFDigiFC_;
vector<vector<float> >   *HFDigiPedFC_;
vector<vector<float> >   *HFDigiNomFC_;

TBranch   *b_run;
TBranch   *b_ls;
TBranch   *b_bx;
TBranch   *b_HFDigiSubdet;
TBranch   *b_HFDigiIEta;
TBranch   *b_HFDigiIPhi;
TBranch   *b_HFDigiDepth;
TBranch   *b_HFDigiRecEnergy;
TBranch   *b_HFDigiFC;
TBranch   *b_HFDigiPedFC;
TBranch   *b_HFDigiNomFC;

TString addCommas(double num){
  TString result(""); result += num;
  int posdot(result.First('.'));
  if(posdot==-1) posdot = result.Length();
  for(int ind(posdot-3); ind > 0; ind -= 3)
    result.Insert(ind, ",");
  return result;
}
void h1cosmetic(TH1F* &h1, int linecolor, int linewidth, int fillcolor){
  h1->SetLineColor(linecolor);
  h1->SetLineWidth(linewidth);
  h1->SetMarkerColor(linecolor);
  h1->SetFillColor(fillcolor);
  h1->SetStats(0);
  h1->SetMinimum(0.1);
}
void h2cosmetic(TH2F* &h2, char* title, TString Xvar, TString Yvar, TString Zvar){
  h2->SetTitle(title);
  h2->SetXTitle(Xvar);
  h2->SetYTitle(Yvar);
  h2->SetZTitle(Zvar);
  h2->SetStats(0);
}
void HFTimingOne(vector< vector<int> > goodbxs, int TStoCheck, int TSadjacent, int IETA, int IPHI, int DEPTH);
vector< vector<int> > findFillScheme(TChain* ch, TString bxtype, bool printbx);
TString hoursMinSec(long seconds){
  int minutes((seconds/60)%60), hours(seconds/3600);
  TString hhmmss("");
  if(hours<10) hhmmss += "0";
  hhmmss += hours; hhmmss += ":";
  if(minutes<10) hhmmss += "0";
  hhmmss += minutes; hhmmss += ":";
  if((seconds%60)<10) hhmmss += "0";
  hhmmss += seconds%60;

  return hhmmss;
}
bool isGoodBX(int run_num, int bx, vector< vector<int> > fillscheme);
void makeTextFile(TH1F* timing[][4], int ieta);
float Max(float a, float b){ 
  return a >= b ? a : b;
}
float Min(float a, float b){ 
  return a <= b ? a : b; 
}
TString roundNumber(double num, int decimals, double denom){
  if(denom==0) return " - ";
  double neg = 1; if(num*denom<0) neg = -1;
  num /= neg*denom; num += 0.5*pow(10.,-decimals);
  long num_int = static_cast<long>(num);
  long num_dec = static_cast<long>((1+num-num_int)*pow(10.,decimals));
  TString s_dec = ""; s_dec += num_dec; s_dec.Remove(0,1);
  TString result="";
  if(neg<0) result+="-";
  result+= num_int;
  if(decimals>0) {
    result+="."; result+=s_dec;
  }
  
  TString afterdot = result;
  afterdot.Remove(0,afterdot.First(".")+1);
  for(int i=0; i<decimals-afterdot.Length(); i++)
    result += "0";
  return result;
}

void init_chain(TChain *ch){

  run_ = 0;
  ls_ = 0;
  bx_ = 0;
  HFDigiSubdet_ = 0;
  HFDigiIEta_ = 0;
  HFDigiIPhi_ = 0;
  HFDigiDepth_ = 0;
  HFDigiRecEnergy_ = 0;
  HFDigiFC_ = 0;
  HFDigiPedFC_ = 0;
  HFDigiNomFC_ = 0;

  ch->SetBranchAddress("run", &run_, &b_run);
  ch->SetBranchAddress("ls", &ls_, &b_ls);
  ch->SetBranchAddress("bx", &bx_, &b_bx);
  ch->SetBranchAddress("HFDigiSubdet", &HFDigiSubdet_, &b_HFDigiSubdet);
  ch->SetBranchAddress("HFDigiIEta", &HFDigiIEta_, &b_HFDigiIEta);
  ch->SetBranchAddress("HFDigiIPhi", &HFDigiIPhi_, &b_HFDigiIPhi);
  ch->SetBranchAddress("HFDigiDepth", &HFDigiDepth_, &b_HFDigiDepth);
  ch->SetBranchAddress("HFDigiRecEnergy", &HFDigiRecEnergy_, &b_HFDigiRecEnergy);
  ch->SetBranchAddress("HFDigiFC", &HFDigiFC_, &b_HFDigiFC);  // Note that FC is after ped subtraction while NomFC is not
  ch->SetBranchAddress("HFDigiPedFC", &HFDigiPedFC_, &b_HFDigiPedFC);
  ch->SetBranchAddress("HFDigiNomFC", &HFDigiNomFC_, &b_HFDigiNomFC);
}

