#include <Riostream.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom.h>
#include <TChain.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <fstream>
//Macro used to Generate all Glauber variables
//The only difference from GetEccenNpart.C is that this one chop the peripheral Glauber events, such that the overall efficiency is 93%.

//separate centrality class according to bbcq;
//one can do with npart if there is no BBC input.

//if you use a different nn cross section, then change it accordingly. It is only used to calculate TAB.
const double sigmann=42;

void GetCentrality(TH1* hbin, int npercent, double* bbcq, TH1* hbbcq, double* abbcq);

//find the centrality bin
//
int FindBin(int npercent, double* bbcq, double bbcqsum);


static float efficiency = 1.0;

void
GetEccenNpart_new(char* xname = "ile")
{
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);

  char name[100];

  TChain* ch = new TChain("t");
  int n;

  //you need a file list of root ntuple: fxxxx.lis
  //
  sprintf(name,"f%s.lis",xname);
  cout << "open the filelist : " << name << endl;
  ifstream fin(name);

  string filename;
  while(!fin.eof()){
    fin >> filename;
    if(filename.empty()) continue;
    cout << filename << endl;

    ch->AddFile(filename.c_str());
  }
  fin.close();

  double b;
  int ncoll;
  int npart;
  int nHitBbc_n, nHitBbc_s;
  double qBBC_n, qBBC_s;
  double vertex, ecc_std, ecc_rp, ecc_part,e4;
  double r_oll,r_geo,r_arith;
  ch->SetBranchAddress("b",            &b        );
  ch->SetBranchAddress("vertex",       &vertex   );
  ch->SetBranchAddress("ncoll",        &ncoll    );
  ch->SetBranchAddress("npart",        &npart    );
  ch->SetBranchAddress("ecc_std",      &ecc_std  );
  ch->SetBranchAddress("ecc_rp",       &ecc_rp  );
  ch->SetBranchAddress("ecc_part",     &ecc_part );
  ch->SetBranchAddress("r_ollitra",    &r_oll );
  ch->SetBranchAddress("r_geo",        &r_geo );
  ch->SetBranchAddress("r_arith",      &r_arith );
  ch->SetBranchAddress("e4",           &e4       );
  ch->SetBranchAddress("qBBC_n",       &qBBC_n   );
  ch->SetBranchAddress("qBBC_s",       &qBBC_s   );
  ch->SetBranchAddress("nHitBbc_n",    &nHitBbc_n);
  ch->SetBranchAddress("nHitBbc_s",    &nHitBbc_s);

  char fname[100];
  sprintf(fname,"rootfile/g%s_cent.root",xname);
  TFile* fout = new TFile(fname,"recreate");//rootfile to save distributions

  TH1* hbbcq    = new TH1D("hbbcq","hbbcq",12000,-0.5,2999.5);      hbbcq->Sumw2();
  TH1* hbbcqall = new TH1D("hbbcqall","hbbcqall",12000,-0.5,2999.5);hbbcqall->Sumw2();
  TH1* hbbcqeff = new TH1D("hbbcqeff","hbbcqeff",12000,-0.5,2999.5);hbbcqeff->Sumw2();

  n = ch->GetEntries();
  cout << "events: "<< n << endl;

  //First loop determines the BBC trigger efficiency
  //as function of bbcq;
  for(int i=0; i<n; i++){
    ch->GetEntry(i);

    if(nHitBbc_n>=2&&nHitBbc_s>=2){
      hbbcq->Fill( (qBBC_n+qBBC_s)/100. );//divide by 100 to fit in 0-3000 range
    }

    hbbcqall->Fill( (qBBC_n+qBBC_s)/100. );
  }
  hbbcqeff->Divide(hbbcq,hbbcqall);

  efficiency = hbbcq->Integral()/hbbcqall->Integral();
  cout << "efficiency : " << efficiency << endl;

  //hbin contains the integrated fraction starting from 0 bbc charge, including the BBC trigger efficiency
  TH1* hbin = new TH1D("hbin","hbin",hbbcq->GetNbinsX(),-0.5,2999.5);
  TH1* hbbcqscale = (TH1*)hbbcq->Clone("hbbcqscale");
  hbbcqscale->Scale(1.0/hbbcq->Integral());
  for(int i=1; i<=hbbcqscale->GetNbinsX(); i++){
    hbin->SetBinContent(i,hbbcqscale->Integral(1,i));
  }
  //Following two lines defines the array of cuts and average bbc charge
  //for centrality percentitle in 5%, 10% and 20% steps
  double bbcq5[21],bbcq10[11],bbcq20[6];
  double abbcq5[20],abbcq10[10],abbcq20[5];

  //calculate the various variables for 5% step
  GetCentrality(hbin,20,bbcq5,hbbcq,abbcq5);
  cout << endl << endl;

  //calculate the various variables for 10% step
  GetCentrality(hbin,10,bbcq10,hbbcq,abbcq10);
  cout << endl << endl;

  //calculate the various variables for 20% step
  GetCentrality(hbin,5,bbcq20,hbbcq,abbcq20);
  cout << endl << endl;

  cout << " Find cuts for all the centralities " << endl << endl;

  const int nval=9;//number of variables to fill
  char* centname[3] = {"5pStep","10pStep","20pStep"};
  char* varname[nval] = {"npart","ncoll","b","standard_ecc","reactionplane_ecc","participant_ecc","R_Ollitraut","R_Geo", "R_arith"};
  double vup[nval] = {499.5,2999.5,19.995,1,1,1,4,4,4};
  double vlo[nval] = {-0.5,-0.5,-0.005,-1,-1,-1,0,0,0};
  int vNb[nval] = {500,3000,2000,200,200,200,400,400,400};

  //initialize the histograms which are used to fill the distribution of variables for each centrality
  TH1* hvar[3][nval][20];
  for(int i=0; i<3; i++){
    int NC = 0;
    if(i==0) NC = 20;
    else if(i==1) NC = 10;
    else if(i==2) NC = 5;
    for(int ivar=0; ivar<nval; ivar++){
      for(int icen=0; icen<NC; icen++){
	sprintf(name,"hvar_%s_%s_cent%d",centname[i],varname[ivar],icen);
	hvar[i][ivar][icen] = new TH1D(name,name,vNb[ivar],vlo[ivar],vup[ivar]);
      }
    }
  }

  double qbbcsum;
  for(int i=0; i<n; i++){
    if(i%1000000==0) cout << i << endl;
    ch->GetEntry(i);

    if(!(nHitBbc_n>=2&&nHitBbc_s>=2)) continue;//BBC trigger condition

    qbbcsum = (qBBC_n+qBBC_s)/100.;

    int centbin5  = FindBin(20,bbcq5,qbbcsum);
    int centbin10 = FindBin(10,bbcq10,qbbcsum);
    int centbin20 = FindBin(5,bbcq20,qbbcsum);

    if(centbin5==-1) continue;
    if(centbin10==-1) continue;
    if(centbin20==-1) continue;

    //find the weight according to the corresponding average efficiency.
    double weight = hbbcqeff->GetBinContent(hbbcqeff->FindBin(qbbcsum));
    //5 percent step
    //
    hvar[0][0][centbin5]->Fill(npart,weight);
    hvar[0][1][centbin5]->Fill(ncoll,weight);
    hvar[0][2][centbin5]->Fill(b,weight);
    hvar[0][3][centbin5]->Fill(ecc_std,weight);
    hvar[0][4][centbin5]->Fill(ecc_rp,weight);
    hvar[0][5][centbin5]->Fill(ecc_part,weight);
    hvar[0][6][centbin5]->Fill(r_oll,weight);
    hvar[0][7][centbin5]->Fill(r_geo,weight);
    hvar[0][8][centbin5]->Fill(r_arith,weight);

    //10 percent step
    //
    hvar[1][0][centbin10]->Fill(npart,weight);
    hvar[1][1][centbin10]->Fill(ncoll,weight);
    hvar[1][2][centbin10]->Fill(b,weight);
    hvar[1][3][centbin10]->Fill(ecc_std,weight);
    hvar[1][4][centbin10]->Fill(ecc_rp,weight);
    hvar[1][5][centbin10]->Fill(ecc_part,weight);
    hvar[1][6][centbin10]->Fill(r_oll,weight);
    hvar[1][7][centbin10]->Fill(r_geo,weight);
    hvar[1][8][centbin10]->Fill(r_arith,weight);


    //20 percent step
    //
    hvar[2][0][centbin20]->Fill(npart,weight);
    hvar[2][1][centbin20]->Fill(ncoll,weight);
    hvar[2][2][centbin20]->Fill(b,weight);
    hvar[2][3][centbin20]->Fill(ecc_std,weight);
    hvar[2][4][centbin20]->Fill(ecc_rp,weight);
    hvar[2][5][centbin20]->Fill(ecc_part,weight);
    hvar[2][6][centbin20]->Fill(r_oll,weight);
    hvar[2][7][centbin20]->Fill(r_geo,weight);
    hvar[2][8][centbin20]->Fill(r_arith,weight);
  }

  //get mean and RMS values for the variables
  float var[3][nval+1][20];
  float rms[3][nval+1][20];
  for(int i=0; i<3; i++){
    int NC = 0;
    if(i==0) NC = 20;
    else if(i==1) NC = 10;
    else if(i==2) NC = 5;
    for(int icen=0; icen<NC; icen++){
      for(int ivar=0; ivar<nval; ivar++){
	var[i][ivar][icen] = hvar[i][ivar][icen]->GetMean();
	rms[i][ivar][icen] = hvar[i][ivar][icen]->GetRMS();
      }
      var[i][nval][icen] = var[i][1][icen]/sigmann;
      rms[i][nval][icen] = rms[i][1][icen]/sigmann;
    }
  }
  //save to file
  for(int ivar=0; ivar<4; ivar++){
    for(int icen=0; icen<16; icen++){
      cout<<var[0][ivar][icen]<<",";
    }
    cout<<var[2][ivar][4]<<",";
    cout<<endl;
  }

  cout.precision(4);

  sprintf(name,"5pstepresults/t%s.txt",xname);
  ofstream f5(name);
  cout << " Bin % &  npart &  ncoll &  b & ecc_std & ecc_rp & ecc_part & r_ollitrau & T_{AB}\\\\\\hline" << endl;
  for(int icen=0; icen<19; icen++){
    for(int ivar=0; ivar<7; ivar++){
      f5 << var[0][ivar][icen] << " ";
    }
    f5 << var[0][nval][icen] << " ";
    f5 <<endl<<  " (";
    /* for(int ivar=0; ivar<7; ivar++){
      f5 << rms[0][ivar][icen] << " ";
    }
    f5 << rms[0][nval][icen] << " ";
    f5 <<  ")"<< endl;
    */
    cout << icen*5 << "-" << icen*5+5;
    for(int ivar=0; ivar<7; ivar++){
      cout <<" & " <<var[0][ivar][icen] ;
    }
    cout <<" & " <<var[0][nval][icen] ;
    cout <<"\\\\"<<endl;
    for(int ivar=0; ivar<7; ivar++){
      if(ivar==0)cout <<" & ("<<rms[0][ivar][icen];
      else cout <<" & "<<rms[0][ivar][icen];
    }
    cout <<" & "<<rms[0][nval][icen];
    cout <<  ")\\\\\\hline" << endl;

    // printf(" & %3.1f & %4.1f & %2.2f & %1.3f & %1.3f & %1.3f & %1.3f \\\\\n",
    //	   var[0][0][icen],var[0][1][icen],var[0][2][icen],var[0][3][icen],var[0][4][icen],var[0][5][icen],var[0][6][icen]);
    //printf(" (& %3.1f & %4.1f & %2.2f & %1.3f & %1.3f & %1.3f & %1.3f) \\\\\\hline\n",
    //       rms[0][0][icen],rms[0][1][icen],rms[0][2][icen],rms[0][3][icen],rms[0][4][icen],rms[0][5][icen],rms[0][6][icen]);
  }
  f5.close();

  cout << endl << endl;
  sprintf(name,"10pstepresults/t%s.txt",xname);
  ofstream f10(name);
  cout <<endl<< " Bin % &  npart &  ncoll &  b & ecc_std & ecc_rp & ecc_part & r_ollitrau & T_{AB}\\\\\\hline" << endl;
  for(int icen=0; icen<10; icen++){
    for(int ivar=0; ivar<7; ivar++){
      f10 << var[1][ivar][icen] << " ";
    }
    f10 << var[1][nval][icen] << " ";
    f10 <<endl<<" (";
    /*for(int ivar=0; ivar<7; ivar++){
      f10 << rms[1][ivar][icen] << " ";
    }
    f10 << rms[1][nval][icen] << " ";
    f10 <<  ")" << endl;
    */
    cout << icen*10 << "-" << icen*10+10;
    for(int ivar=0; ivar<7; ivar++){
      cout <<" & " <<var[1][ivar][icen] ;
    }
    cout <<" & " <<var[1][nval][icen] ;
    cout <<"\\\\"<<endl;
    for(int ivar=0; ivar<7; ivar++){
      if(ivar==0)cout <<" & ("<<rms[1][ivar][icen];
      else cout <<" & "<<rms[1][ivar][icen];
    }
    cout <<" & "<<rms[1][nval][icen];
    cout <<  ")\\\\\\hline" << endl;

    //printf(" & %3.1f & %4.1f & %2.2f & %1.3f & %1.3f & %1.3f & %1.3f \\\\\n",
    //       var[1][0][icen],var[1][1][icen],var[1][2][icen],var[1][3][icen],var[1][4][icen],var[1][5][icen],var[1][6][icen]);
    //printf(" (& %3.1f & %4.1f & %2.2f & %1.3f & %1.3f & %1.3f & %1.3f) \\\\\\hline\n",
    //       rms[1][0][icen],rms[1][1][icen],rms[1][2][icen],rms[1][3][icen],rms[1][4][icen],rms[1][5][icen],rms[1][6][icen]);

  }
  f10.close();

  cout << endl << endl;
  sprintf(name,"20pstepresults/t%s.txt",xname);
  ofstream f20(name);
  cout <<endl<< " Bin % &  npart &  ncoll &  b & ecc_std & ecc_rp & ecc_part & r_ollitrau & T_{AB}\\\\\\hline" << endl;
  for(int icen=0; icen<5; icen++){
    for(int ivar=0; ivar<7; ivar++){
      f20 << var[2][ivar][icen] << " ";
    }
    f20 << var[2][nval][icen] << " ";
    f20 << endl<<  " (";

    /*for(int ivar=0; ivar<7; ivar++){
      f20 << rms[2][ivar][icen] << " ";
    }
    f20 << rms[2][nval][icen] << " ";
    f20 <<  ")"<< endl;
    */
    cout << icen*20 << "-" << icen*20+20;
    for(int ivar=0; ivar<7; ivar++){
      cout <<" & " <<var[2][ivar][icen] ;
    }
    cout <<" & " <<var[2][nval][icen] ;
    cout <<"\\\\"<<endl;
    for(int ivar=0; ivar<7; ivar++){
      if(ivar==0)cout <<" & ("<<rms[2][ivar][icen];
      else cout <<" & "<<rms[2][ivar][icen];
    }
    cout <<" & "<<rms[2][nval][icen];
    cout <<  ")\\\\\\hline" << endl;
    //printf(" & %3.1f & %4.1f & %2.2f & %1.3f & %1.3f & %1.3f & %1.3f \\\\\n",
    //       var[2][0][icen],var[2][1][icen],var[2][2][icen],var[2][3][icen],var[2][4][icen],var[2][5][icen],var[2][6][icen]);
    //printf(" (& %3.1f & %4.1f & %2.2f & %1.3f & %1.3f & %1.3f & %1.3f) \\\\\\hline\n",
    //       rms[2][0][icen],rms[2][1][icen],rms[2][2][icen],rms[2][3][icen],rms[2][4][icen],rms[2][5][icen],rms[2][6][icen]);

  }
  f20.close();

  fout->Write();
  fout->Close();

  return;
}

int
FindBin(int ncentbin, double* bbcq, double bbcqsum)
{

  for(int i=0; i<ncentbin; i++){

    if(bbcqsum<=bbcq[i] && bbcqsum>bbcq[i+1]){
      return i;
    }

  }

  return -1;
}


//hbbcq and hbin are the differential and integrated distribution for BBC charge after triggering condition
//ncentbin is the number of centrality bins. bbcq is the array containing the centrality bin boundary.
//abbcq is array containing the average bbc charge in each centrality bin
//centrality bin is searched starting from largest bbc charge
void
GetCentrality(TH1* hbin, int ncentbin, double* bbcq, TH1* hbbcq, double* abbcq)
{

  int ibegin = hbin->GetNbinsX();
  int iend   = 1;

  bbcq[0] = 2999.5; bbcq[ncentbin] = -0.5;

  //fraction is the fraction of events for each centrality bin
  //MinBias trigger efficiency of 93% is used for Run7
  //so we rescale the total integral of hbbcq to 0.93
  //100/ncentbin = bin width. 1/0.93 is the scale factor

  double fraction = 1.0/(100*efficiency)*(100.0/ncentbin);

  for(int ij=1; ij<=ncentbin; ij++)
    {
      double target = 1.0-ij*fraction;

      if(ncentbin==20){
	if(ij==19) target = 0.01273;//1-0.93/efficiency
	if(ij==20) continue;
      }else{
	if(ij==ncentbin)  target = 0.01273;
      }

      bool found = false;
      for(int k=ibegin; k>=iend;k--)
	{
	  double bc = hbin->GetBinContent(k);
	  double bd = hbin->GetBinCenter(k);
	  if(fabs(bc-target)<0.002) found = true;

	  if(found){
	    ibegin = k-1;
	    bbcq[ij] = bd;
	    break;
	  }
	}//k
      cout << "bbcq[" << ij << "] = " << bbcq[ij] << endl;
    }
}
