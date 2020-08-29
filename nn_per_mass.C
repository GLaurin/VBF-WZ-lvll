#include <TString.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TLine.h>
#include <TColor.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TChain.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <iostream>
#include <unordered_map>
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

/*----------------------------------------------------------------------------*/

// Initializing variables
string   idir, tmass, sdir;
TString  title, proj_str, select_weight, proj_option="";

int   nbins=50; 
float xmin=0, xmax=1;

float Nsig_ocv, Nbkg_ocv, AMS_ocv, Nsig_cv, Nbkg_cv, AMS_cv;
int   tmCV;

/*----------------------------------------------------------------------------*/

string get_file_name(int mass, string phys_model="GM") {
  string insert_str="main";
  string              file_path="OutputRoot/"+sdir+"/new_"+phys_model+"_"+insert_str;
  if      (mass==200) file_path+="MVA.450765_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m200_ntuples.root";
  else if (mass==250) file_path+="MVA.450766_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m250_ntuples.root";
  else if (mass==300) file_path+="MVA.450767_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m300_ntuples.root";
  else if (mass==350) file_path+="MVA.450768_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m350_ntuples.root";
  else if (mass==400) file_path+="MVA.450769_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m400_ntuples.root";
  else if (mass==450) file_path+="MVA.450770_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m450_ntuples.root";
  else if (mass==500) file_path+="MVA.450771_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m500_ntuples.root";
  else if (mass==600) file_path+="MVA.450772_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m600_ntuples.root";
  else if (mass==700) file_path+="MVA.450773_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m700_ntuples.root";
  else if (mass==800) file_path+="MVA.450774_MGaMcAtNloPy8EG_A14NNPDF23LO_vbfGM_sH05_H5pWZ_lvll_m800_ntuples.root";
  else if (mass==900) file_path+="MVA.305035_MGPy8_A14NNPDF30NLO_VBS_H5p_lvll_900_qcd0_ntuples.root";

  if (phys_model=="HVT") {
    file_path="OutputRoot/"+sdir+"/new_"+phys_model+"_"+insert_str;
    if      (mass== 250) file_path+="MVA.307730_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0250_ntuples.root";
    else if (mass== 300) file_path+="MVA.307731_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0300_ntuples.root";
    else if (mass== 350) file_path+="MVA.309528_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0350_ntuples.root";
    else if (mass== 400) file_path+="MVA.307732_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0400_ntuples.root";
    else if (mass== 450) file_path+="MVA.309529_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0450_ntuples.root";
    else if (mass== 500) file_path+="MVA.307733_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0500_ntuples.root";
    else if (mass== 600) file_path+="MVA.307734_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0600_ntuples.root";
    else if (mass== 700) file_path+="MVA.307735_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0700_ntuples.root";
    else if (mass== 800) file_path+="MVA.307736_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0800_ntuples.root";
    else if (mass== 900) file_path+="MVA.307737_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0900_ntuples.root";
    else if (mass==1000) file_path+="MVA.307738_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1000_ntuples.root";
  }

  if (phys_model=="QQ") {
    file_path="OutputRoot/"+sdir+"/new_"+phys_model+"_"+insert_str;
    if      (mass== 500) file_path+="MVA.307733_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0500_ntuples.root";
    else if (mass== 600) file_path+="MVA.307734_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0600_ntuples.root";
    else if (mass== 700) file_path+="MVA.307735_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0700_ntuples.root";
    else if (mass== 800) file_path+="MVA.307736_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0800_ntuples.root";
    else if (mass== 900) file_path+="MVA.307737_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m0900_ntuples.root";
    else if (mass==1000) file_path+="MVA.307738_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1000_ntuples.root";
    else if (mass==1100) file_path+="MVA.307739_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1100_ntuples.root";
    else if (mass==1200) file_path+="MVA.307740_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1200_ntuples.root";
    else if (mass==1300) file_path+="MVA.307741_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1300_ntuples.root";
    else if (mass==1400) file_path+="MVA.307742_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1400_ntuples.root";
    else if (mass==1500) file_path+="MVA.307743_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1500_ntuples.root";
    else if (mass==1600) file_path+="MVA.307744_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1600_ntuples.root";
    else if (mass==1700) file_path+="MVA.307745_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1700_ntuples.root";
    else if (mass==1800) file_path+="MVA.307746_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1800_ntuples.root";
    else if (mass==1900) file_path+="MVA.307747_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m1900_ntuples.root";
    else if (mass==2000) file_path+="MVA.307748_MGPy8EG_A14NNPDF23LO_vbfHVT_Agv1_VzWZ_lvll_m2000_ntuples.root";
  }

  bool warning= (phys_model=="GM" and mass>900) or (phys_model=="HVT" and mass==200) or (phys_model=="QQ" and mass<500);
  if (warning) {cout<<"phys_model: "<<phys_model<<" while mass= "<<mass<<". Aborting the process."<<endl; exit(1);}

  return file_path;
}

string QQ_filename(int mass) {
  // HVT's qq-fusion ntuples' filenames per mass
  string insert_str="main";
  string               file_path="OutputRoot/"+sdir+"/new_QQ_"+insert_str;
  if      (mass==500)  file_path+="MVA.302266_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m0500_ntuples.root";
  else if (mass==600)  file_path+="MVA.302267_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m0600_ntuples.root";
  else if (mass==700)  file_path+="MVA.302268_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m0700_ntuples.root";
  else if (mass==800)  file_path+="MVA.302269_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m0800_ntuples.root";
  else if (mass==900)  file_path+="MVA.302270_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m0900_ntuples.root";
  else if (mass==1000) file_path+="MVA.302271_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m1000_ntuples.root";
  else if (mass==1100) file_path+="MVA.302272_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m1100_ntuples.root";
  else if (mass==1200) file_path+="MVA.302273_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m1200_ntuples.root";
  else if (mass==1300) file_path+="MVA.302274_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m1300_ntuples.root";
  else if (mass==1400) file_path+="MVA.302275_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m1400_ntuples.root";
  else if (mass==1500) file_path+="MVA.302276_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m1500_ntuples.root";
  else if (mass==1700) file_path+="MVA.302278_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m1700_ntuples.root";
  else if (mass==1800) file_path+="MVA.302279_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m1800_ntuples.root";
  else if (mass==1900) file_path+="MVA.302280_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m1900_ntuples.root";
  else if (mass==2000) file_path+="MVA.302281_MadGraphPythia8EvtGen_A14NNPDF23LO_HVT_Agv1_VcWZ_lvll_m2000_ntuples.root";
  return file_path;
}

int get_color(int mass) {
  // Color to use for each mass
  int color = kBlack;
  if      (mass==250 or mass==1000 or mass==1800)  color=kGray+2;
  else if (mass==300 or mass==1100 or mass==1900)  color=kMagenta;
  else if (mass==400 or mass==1200 or mass==2000)  color=kBlue; 
  else if (mass==500 or mass==1300)                color=kCyan+2;
  else if (mass==600 or mass==1400)                color=kGreen+2; 
  else if (mass==700 or mass==1500)                color=kYellow+2;
  else if (mass==800 or mass==1600)                color=kOrange;
  else if (mass==900 or mass==1700)                color=kRed;
  return color;
}

TH1F* get_bkg_hist(TString phys_model="GM") {
  TChain* chain = new TChain("nominal");
  TString ins_str="main";
  //chain->Add("OutputRoot/new_GM_"+ins_str+"MVA.364253_Sherpa_222_NNPDF30NNLO_lllv_ntuples.root");

  chain->Add(((TString)"OutputRoot/")+sdir.data()+"/new_"+phys_model+"_"+ins_str+"MVA.361292_MGaMcAtNloPy8EG_NNPDF30LO_A14NNPDF23LO_WZ_lvll_FxFx_ntuples.root");
  chain->Add(((TString)"OutputRoot/")+sdir.data()+"/new_"+phys_model+"_"+ins_str+"MVA.364284_Sherpa_222_NNPDF30NNLO_lllvjj_EW6_ntuples.root");

  TH1F* hist = new TH1F("bkg",title,nbins,xmin,xmax);
  chain->Project(hist->GetName(),proj_str,select_weight,proj_option);

  return hist;
}

TH1F* get_hist(int mass,TString phys_model="GM", bool qqplot=false) {

  TH1F* hist;
  if (mass>0) {
    string fname=get_file_name(mass,phys_model.Data());
    if (qqplot==true) fname = QQ_filename(mass);
    
    TFile* f = TFile::Open(fname.data(),"read");
    TTree* t = (TTree*)f->Get("nominal");
    
    TString histName = "mass"+TString::Itoa(mass,10);
    //std::cout<<histName<<std::endl;

    hist = new TH1F(histName ,title,nbins,xmin,xmax);
    if (phys_model!="QQ") select_weight += "*(abs(Weight)<10)";
    t->Project(hist->GetName(),proj_str,select_weight,proj_option);
  }
  else hist = get_bkg_hist(phys_model.Data());

  hist->SetLineColor(get_color(mass));

  return hist;
}

float AMS_old(float s, float b, bool debug=false) {
  // Older version of AMS the includes some strange terms in its log
  if (s<=0 or b<=0) return 0;

  float br = 0.00001;// #KM: systematic unc?
  float sigma=sqrt(b+br);
  float n=s+b+br;
  float radicand = 2 *( n * log (n*(b+br+sigma)/(b*b+n*sigma+br))-b*b/sigma*log(1+sigma*(n-b)/(b*(b+br+sigma))));

  float ams= 0;
  if (radicand < 0) std::cout<<"AMS: radicand is negative. Returning 0."<<std::endl;
  else       ams= sqrt(radicand);

  if (debug) std::cout<<"s, b="<<s<<"\t"<<b<<", ams="<<ams<<std::endl;

  return ams;
}

float AMS(float s, float b, bool debug=false, float br=0) {
  float rad = 2*((s+b+br)*log(1+(s/(b+br)))-s);
  if (rad>0 && rad<1000) return sqrt(rad);
  else {cout << "AMS: radicand is negative. Returning 0" << endl; return 0;}
}

TH1F* get_significance_hist(TH1F* h_sig, TH1F* h_bkg, float sf, bool is_tm=false) {
  // Significance histogram
  TH1F* significance = new TH1F("significance",title,nbins,xmin,xmax);
  significance->SetTitle("Significance for yield / 140fb-1");
  h_sig->Scale(sf);
  h_bkg->Scale(sf);

  // Filling histogram
  const int nb = nbins;
  float Nsig[nb], Nbkg[nb];
  for( int i=0; i<nbins; i++) {
    Nsig[i] = h_sig->Integral(i,h_sig->GetNbinsX());
    Nbkg[i] = h_bkg->Integral(i,h_bkg->GetNbinsX());
    significance->SetBinContent(i,AMS(Nsig[i],Nbkg[i]));
  }
  
  // Optimal cut vlue
  int ocv_bin = significance->GetMaximumBin();
  Nsig_ocv = Nsig[ocv_bin];  
  Nbkg_ocv = Nbkg[ocv_bin];
  AMS_ocv  = significance->GetBinContent(ocv_bin);

  // Training mass' optimal cut value
  if (is_tm) tmCV = ocv_bin;
  Nsig_cv = Nsig[tmCV];
  Nbkg_cv = Nbkg[tmCV];
  AMS_cv  = significance->GetBinContent(tmCV);

  return significance;
}

void nn_per_mass(string idir="", string tmass="", TString phys_model="GM", bool drawCB=true, TString varname="pSignal", bool norm2yield=true) {

  // Subdirectory in OutputRoot containing the decorated samples
  sdir = idir+'/'+tmass+"/";

  // Possible variables to plot
  if      (varname == "pSignal"     ) title="NN output : "+tmass, proj_str=varname, nbins = 50, xmin =0, xmax = 1;
  else if (varname == "M_WZ"        ) title=varname, proj_str=varname, nbins = 25, xmin =0, xmax = 1500;
  else if (varname == "M_jj"        ) title=varname, proj_str=varname, nbins = 50, xmin =0, xmax = 1500;
  else if (varname == "ZetaLep"     ) title=varname, proj_str=varname, nbins = 50, xmin =-3.5, xmax = 3.5;
  else if (varname == "DY_jj"       ) title=varname, proj_str=varname, nbins = 50, xmin =0, xmax = 10;
  else if (varname == "Deta_jj"     ) title=varname, proj_str=varname, nbins = 50, xmin =0, xmax = 10;
  else if (varname == "Dphi_jj"     ) title=varname, proj_str=varname, nbins = 50, xmin =0, xmax = TMath::Pi();
  else if (varname == "Meff"        ) title=varname, proj_str= "Pt_W+Pt_Z+Jet1Pt+Jet2Pt+Met", nbins = 50, xmin =0, xmax = 2000;
  else if (varname == "dEta_WZ"     ) title=varname, proj_str= "abs(Eta_Z-Eta_W)", nbins = 50, xmin =0, xmax = 5;
  else if (varname == "Eta_Z"       ) title=varname, proj_str= varname, nbins = 50, xmin =-3, xmax = 3;
  else if (varname == "Eta_W"       ) title=varname, proj_str= varname, nbins = 50, xmin =-3, xmax = 3;
  else if (varname == "Pt_Z"        ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 1000;
  else if (varname == "Pt_W"        ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 1000;
  else if (varname == "Met"         ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 500;
  else if (varname == "Njets"       ) title=varname, proj_str= varname, nbins = 5 , xmin =2, xmax = 7;
  else if (varname == "NBjets"      ) title=varname, proj_str= varname, nbins = 5 , xmin =0, xmax = 5;
  else if (varname == "PtBalanceZ"  ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 2;
  else if (varname == "PtBalanceW"  ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 2;
  else if (varname == "ZetaLep"     ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 10;
  else if (varname == "Jet1E"       ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 1000;
  else if (varname == "Jet2E"       ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 500;
  else if (varname == "Jet3E"       ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 300;
  else if (varname == "Jet1Pt"      ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 300;
  else if (varname == "Jet2Pt"      ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 300;
  else if (varname == "Jet3Pt"      ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 200;
  else if (varname == "Jet1Eta"     ) title=varname, proj_str= varname, nbins = 50, xmin =-5, xmax = 5;
  else if (varname == "Jet2Eta"     ) title=varname, proj_str= varname, nbins = 50, xmin =-5, xmax = 5;
  else if (varname == "Jet3Eta"     ) title=varname, proj_str= varname, nbins = 50, xmin =-5, xmax = 5;
  else if (varname == "Lep1Pt"      ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 500;
  else if (varname == "Lep2Pt"      ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 500;
  else if (varname == "Lep3Pt"      ) title=varname, proj_str= varname, nbins = 50, xmin =0, xmax = 300;
  else if (varname == "Lep1Eta"     ) title=varname, proj_str= varname, nbins = 50, xmin =-2.5, xmax = 2.5;
  else if (varname == "Lep2Eta"     ) title=varname, proj_str= varname, nbins = 50, xmin =-2.5, xmax = 2.5;
  else if (varname == "Lep3Eta"     ) title=varname, proj_str= varname, nbins = 50, xmin =-2.5, xmax = 2.5;       
  else    proj_option="norm"; //normalize to 1

  // List of the masses and list size and half-size
  vector<int> masses{0,250,300,400,500,600,700,800,900};
  if (phys_model=="QQ") masses={0,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1700,1800,1900,2000};
  const int ms = masses.size();
  int      hms = masses.size()/2+1;

  // NN output's figure
  TCanvas* c1 = new TCanvas ("tmass", "title", 800, 400);
  c1->Divide(2,1);

  //Plot 2 - second half of the masses
  c1->cd(2);
  auto legend2 = new TLegend(0.7,0.7,0.9,0.9);
  legend2->SetHeader("Mass (GeV)","C"); 
  legend2->SetFillStyle(0); 
  legend2->SetLineWidth(0); 
  legend2->SetNColumns(2);

  // Plot 1 - first half of the masses
  c1->cd(1);
  auto legend1 = new TLegend(0.7,0.7,0.9,0.9);
  legend1->SetHeader("Mass (GeV)","C"); 
  legend1->SetFillStyle(0); 
  legend1->SetLineWidth(0); 
  legend1->SetNColumns(2);

  // Initializing histograms of the Output's (NN, NN's bkg, CB and CB's bkg)
  unordered_map<int,TH1F*> hists, hists_bkg, hists_cb, hists_bkg_cb;

  // Initializing legend and histogram
  auto legend=legend1;
  char smass[3];
  TH1F* hist;

  // Plotting one mass at a time
  for (auto mass : masses) {

    // Re-initializing selection
    select_weight = "(M_jj>100)";
    if (phys_model=="QQ") select_weight="1";
    if (norm2yield) select_weight += "*WeightNormalized";

    // Separating the figure
    if (mass==masses[hms]) {
      legend->Draw();
      gStyle->SetOptStat(0);
      if (varname=="pSignal" and norm2yield) gPad->SetLogy();
      c1->cd(2); 
      legend=legend2; 
    
      //Plotting background on second figure
      hist = get_hist(0,phys_model.Data());
      TString option = "hist";
      hist->Draw(option);
      sprintf(smass, "%s", "bkg");
      legend->AddEntry(hist,smass,"f");
    }

    // Filling the histograms
    if (mass != 0) select_weight += Form("*(M_WZ>(%i*0.6)*(M_WZ<(%i*1.4)))",mass,mass);
    // Background histogram
    hists_bkg[mass] = get_hist(0,phys_model.Data());
    // Current mass histogram
    hist = get_hist(mass,phys_model.Data());
    hists[mass] = hist;

    if (mass != 0) {
      // Cut-based selection and histogram
      select_weight = Form("(M_jj>500)*(Deta_jj>3.5)*(M_WZ>(%i*0.6)*(M_WZ<(%i*1.4)))",mass,mass);
      if (norm2yield) select_weight += "*WeightNormalized";
      hists_bkg_cb[mass] = get_hist(0,phys_model.Data());
      hists_cb[mass] = get_hist(mass,phys_model.Data());
    }

    // Drawing the histogram
    TString option="same hist";
    if (mass==0) option="hist";
    hist->Draw(option);

    // Filling the legend
    if (mass != 0) sprintf(smass, "%i", mass);
    else           sprintf(smass, "%s", "bkg");
    legend->AddEntry(hist,smass,"f");
  }
  if (varname=="pSignal" and norm2yield) gPad->SetLogy();
  gStyle->SetOptStat(0);
  legend->Draw();

  // Saving the NN output figure
  string imagePath = "ControlPlots/"+idir+"/NN_output/"+varname.Data() + (idir!="" ? "_"+idir : "") + (tmass!="" ? "_"+tmass : "");
  c1->SaveAs((imagePath+".png" ).data());
  c1->SaveAs((imagePath+".root").data());

  // Further analysis for pSignal
  if (not (norm2yield and varname=="pSignal")) return;

  // Initializing AMS's figure
  auto c2 = new TCanvas("c2","title",800,400);
  c2->Divide(2,1);
  c2->cd(1);

  // Figure scale
  float sf= 1;

  // Initializing CB vs NN CSV file
  ofstream cnfile("ControlPlots/"+idir+"/NN_output/CbvsNN"+ (idir!="" ? "_"+idir : "") + (tmass!="" ? "_"+tmass : "")+".csv");
  map<int, float> sig_CB, bkg_CB, ams_CB, sig_NN_ocv, bkg_NN_ocv, ams_NN_ocv, sig_NN_cv, bkg_NN_cv, ams_NN_cv;

  // Setting the training mass' optimal cut value (tmCV)
  if (tmass!="mMulti") auto significance = get_significance_hist(hists[stoi(tmass.substr(1,4))], hists_bkg[stoi(tmass.substr(1,4))], sf, true);

  for (auto mass : masses) {
    if (mass==0) continue;

    // First line of the CSV file : masses
    cnfile << "," << mass << "," << mass;
    if (tmass!="mMulti") cnfile << "," << mass;

    // Separating the figure in two
    if (mass==masses[hms]) {
      legend1->Draw();
      if (phys_model=="QQ") gPad->SetLogy(), gStyle->SetOptStat(0);
      c2->cd(2);
    }

    // Signal, background and significance for NN for the two cuts
    auto significance = get_significance_hist(hists[mass],hists_bkg[mass],sf);
    sig_NN_ocv[mass] = Nsig_ocv;
    bkg_NN_ocv[mass] = Nbkg_ocv;
    ams_NN_ocv[mass] = AMS_ocv;
    sig_NN_cv[mass]  = Nsig_cv;
    bkg_NN_cv[mass]  = Nbkg_cv;
    ams_NN_cv[mass]  = AMS_cv;

    // Drawing significance curve
    TString option="same hist";
    if (mass==masses[1] || mass==hms) option="hist";
    significance->SetLineColor(get_color(mass));
    significance->Draw(option);

    // Comparing with Cut-based
    if (drawCB) {
      float Nsig = hists_cb[mass]->Integral();
      float Nbkg = hists_bkg_cb[mass]->Integral();
      float cb_ams = AMS(Nsig, Nbkg);

      // Drawing CB significance
      TLine *line = new TLine(0,cb_ams,1,cb_ams);
      line->SetLineColor(get_color(mass));
      line->SetLineStyle(7);
      line->Draw("same hist");

      // Filling the maps
      bkg_CB[mass] = Nbkg;
      sig_CB[mass] = Nsig;
      ams_CB[mass] = cb_ams;
    significance->SetMaximum(30);     // Significance's y-maximum - may have to be changed
    if (phys_model=="QQ") significance->SetMaximum(2), significance->SetMinimum(1e-3), gPad->SetLogy();
    }
  }
  legend2->Draw();
  if (phys_model=="QQ") gPad->SetLogy(), gStyle->SetOptStat(0);

  // Saving the significance figure
  string signPath = "ControlPlots/"+idir+"/NN_output/significance" + (idir!="" ? "_"+idir : "") + (tmass!="" ? "_"+tmass : "");
  c2->SaveAs((signPath+".png" ).data());
  c2->SaveAs((signPath+".root").data());

  // Plotting the significance per mass and method
  auto c3 = new TCanvas("c3","title",500,400);
  auto ams_legend = new TLegend(0.75,0.78,0.9,0.9);
  map<int, float> ams_TD;
  TString ams_curve_opt, label, ams_title = "Significance per mass and method - "+tmass;

  for (int l; l<3; l++) {
    auto ams_curve = new TGraph();

    // Selecting the AMS to draw and its options
    if (l==0)                         {ams_TD = ams_NN_ocv; ams_curve_opt = "APL"; ams_curve->SetMarkerStyle(20); label="NN opt";}
    else if (l==1 && tmass!="mMulti") {ams_TD = ams_NN_cv;  ams_curve_opt = "PL";  ams_curve->SetMarkerStyle(21); label="NN tcv";}
    else if (l==2)                    {ams_TD = ams_CB;     ams_curve_opt = "PL";  ams_curve->SetMarkerStyle(23); label="Cut-Based";}
    else continue;
    ams_curve->SetLineColor(2+l);
    ams_curve->SetMarkerColor(2+l);
    ams_curve->SetTitle(ams_title);

    // Setting the value for each AMS
    int k = 0;
    for (auto mass : masses) {
      if (mass == 0) continue;
      ams_curve->SetPoint(k, mass, ams_TD[mass]);
      k++;
    }

    // Figure's y-limits
    ams_curve->SetMaximum(25);
    ams_curve->SetMinimum(-0.8);
    if (tmass!="mMulti" && phys_model=="HVT" && (stoi(tmass.substr(1,4))>400)) {c3->SetLogy(); ams_curve->SetMinimum(10E-3);}
    if (phys_model=="QQ") {ams_curve->SetMaximum(2); ams_curve->SetMinimum(-0.1);}

    // Drawing the curve and updating the legend
    ams_curve->Draw(ams_curve_opt);
    ams_legend->AddEntry(ams_curve,label,"lep");
  }
  ams_legend->Draw();

  // Saving the figure
  c3->SaveAs(("ControlPlots/"+idir+"/NN_output/AMS_curves" + (idir!="" ? "_"+idir : "") + (tmass!="" ? "_"+tmass : "") + ".png").data());
  c3->SaveAs(("ControlPlots/"+idir+"/NN_output/AMS_curves" + (idir!="" ? "_"+idir : "") + (tmass!="" ? "_"+tmass : "") + ".root").data());

  // HVT's qq samples' pSignal figure
  auto c4 = new TCanvas("c4","title",800,400);
  c4->Divide(2,1);
  c4->cd(1);

  for (auto mass : masses) {

    // Selection
    select_weight = "1";
    if (norm2yield) select_weight = "WeightNormalized";
    select_weight += Form("*(M_WZ>(%i*0.6)*(M_WZ<(%i*1.4)))",mass,mass);

    //Separating the curves on 2 figures
    if (mass==masses[hms]) {
      legend1->Draw();
      gStyle->SetOptStat(0);
      gPad->SetLogy();
      c4->cd(2);
    }

    // Current mass histogram
    hist = get_hist(mass,phys_model.Data(),true);

    // Drawing the curve
    TString option="same hist";
    if (mass==0) option="hist";
    hist->Draw(option);
  }
  legend2->Draw();
  gStyle->SetOptStat(0);
  gPad->SetLogy();

  // Saving QQ figure
  c4->SaveAs(("ControlPlots/"+idir+"/NN_output/QQ_pSignal" + (idir!="" ? "_"+idir : "") + (tmass!="" ? "_"+tmass : "") + ".png").data());
  c4->SaveAs(("ControlPlots/"+idir+"/NN_output/QQ_pSignal" + (idir!="" ? "_"+idir : "") + (tmass!="" ? "_"+tmass : "") + ".root").data());

  // Filling CSV file of number of signal and bakcground, and significance by method
  cnfile << endl;
  for (int j; j<4; j++) {
    if (j==0) cnfile << ",";
    if (j==1) cnfile << "Signal,";
    if (j==2) cnfile << "Background,";
    if (j==3) cnfile << "Significance,";
    for (auto mass : masses) {
      if (mass == 0) continue;
      if (j==0) cnfile << "CB,NN_ocv" << (tmass!="mMulti" ? ",NN_tcv" : "");
      if (j==1) { cnfile << sig_CB[mass] << "," << sig_NN_ocv[mass];
        if (tmass!="mMulti") cnfile << "," << sig_NN_cv[mass];
        if (mass!=masses[ms-1]) cnfile << ",";
      }
      if (j==2) { cnfile << bkg_CB[mass] << "," << bkg_NN_ocv[mass];
        if (tmass!="mMulti") cnfile << "," << bkg_NN_cv[mass];
        if (mass!=masses[ms-1]) cnfile << ",";
      }
      if (j==3) { cnfile << ams_CB[mass] << "," << ams_NN_ocv[mass];
        if (tmass!="mMulti") cnfile << "," << ams_NN_cv[mass];
        if (mass!=masses[ms-1]) cnfile << ",";
      }
    }
    cnfile << endl;
  }
  cnfile.close();

  return;
 
}

