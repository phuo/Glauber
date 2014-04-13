/*************************************************************************
  > File Name: drawE.C
 ************************************************************************/
#include<iostream>
#include<fstream>
#include<vector>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1.h"
#include "TH2.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TEllipse.h"
#include "TStyle.h"
#include "ECut.C"

using namespace std;

enum{ ETOT=6, };

int findEbin(int centb, double ecc, int ihar){
    int steps=50;
    int bin=-1;
    for(int is=0; is<steps;is++){
        if(ecc<ecuts[centb][ihar][is]){ bin = is; break; }
        //if(ecc<ecuts[centb][ihar][is]){ bin = is;cout<<"ihar "<<ihar<<"  "<<is<<"  "<<ecc<<"   "<<ecuts[centb][ihar][is]<<endl; break; }
    }
    if(bin<0) {cout<<"sssswrong  "<<ecc<<"  "<<ihar<<endl;}
    return bin;
}
void drawE(){
    TChain* t = new TChain("t");
    char dfilelist[100];
    char name[100];
    sprintf(dfilelist,"outlist.txt");
    ifstream lis(dfilelist);
    int cnt=0;
    int dnt=0;
    while(!lis.eof())
    {
        string filename;
        lis >> filename;
        if(cnt>=0&&cnt<100)
        {
            cout << filename << endl;
            if(!filename.empty()) t->Add(filename.c_str());
            dnt++;
        }
        cnt++;
    }

    double fAngGl[ETOT],feGl[ETOT];
    double frAngGl[ETOT],freGl[ETOT];

    double pfAngGl[ETOT],pfeGl[ETOT];
    double pfrAngGl[ETOT],pfreGl[ETOT];

    double b;
    int ncoll;
    int npart,npartproj,nparttarg;
    int nHitBbc_n, nHitBbc_s, BBCTrig;
    int Centbin;
    double qBBC_n, qBBC_s, pTrigBBC;
    double vertex, ecc_std,ecc_rp, ecc_part,ecc_partr,r_ollitra,r_geo,r_arith,r_ollitrar,r_geor,r_arithr,e4;


    float nux[600] ;
    float nuy[600] ;
    int  st_part[600] ;

    t->SetBranchAddress("b",            &b);
    t->SetBranchAddress("Centbin",   &Centbin);
    t->SetBranchAddress("ncoll",        &ncoll);
    t->SetBranchAddress("npart",        &npart);
    /*    t->SetBranchAddress("npartproj",    &npartproj);
          t->SetBranchAddress("nparttarg",    &nparttarg);
          t->SetBranchAddress("ecc_std",      &ecc_std);
          t->SetBranchAddress("ecc_rp",       &ecc_rp);
          t->SetBranchAddress("ecc_part",     &ecc_part);

          t->SetBranchAddress("e_gl",         feGl);
          t->SetBranchAddress("ang_gl",    fAngGl);
          t->SetBranchAddress("re_gl",       freGl);
          t->SetBranchAddress("rang_gl", frAngGl);

          t->SetBranchAddress("nux", &nux);
          t->SetBranchAddress("nuy", &nuy);
          t->SetBranchAddress("st_part",&st_part);

          t->SetBranchAddress("pe_gl",        pfeGl);
          t->SetBranchAddress("pang_gl",      pfAngGl);
          */
    t->SetBranchAddress("pre_gl",       pfreGl);
    t->SetBranchAddress("prang_gl",     pfrAngGl);

    cout<<t->GetEntries()/1000000<<endl;

    TH1* hevent[7];
    for(int ic=0;ic<7; ic++){
   sprintf(name, "hevent%d", ic); 
    hevent[ic]= new TH1D(name, name,10, 0-0.5, 50-0.5);
    }


    for(int iev=0; iev<t->GetEntries(); iev++){
        t->GetEntry(iev);
        if(iev%10000000==0) cout<<iev<<endl;
        if(Centbin!=0) continue;
        int Ebin[6] = {-1,-1,-1,-1,-1,-1};
        int cutID[7] = {0,0,0,0,0,0,0};
        for(int ih=0;ih<6; ih++){ //iHar bin
            Ebin[ih] = findEbin(Centbin, pfreGl[ih], ih);
        }
        if(Ebin[5]<10 && Ebin[1]<10 && Ebin[2]<10 && Ebin[3]<10 && Ebin[4]<10) cutID[0] = 1;
        if(Ebin[0]<10 && Ebin[5]<10 && Ebin[2]<10 && Ebin[3]<10 && Ebin[4]<10) cutID[1] = 1;
        if(Ebin[0]<10 && Ebin[1]<10 && Ebin[5]<10 && Ebin[3]<10 && Ebin[4]<10) cutID[2] = 1;
        if(Ebin[0]<10 && Ebin[1]<10 && Ebin[2]<10 && Ebin[5]<10 && Ebin[4]<10) cutID[3] = 1;
        if(Ebin[0]<10 && Ebin[1]<10 && Ebin[2]<10 && Ebin[3]<10 && Ebin[5]<10) cutID[4] = 1;
        if(Ebin[0]<10 && Ebin[1]<10 && Ebin[2]<10 && Ebin[3]<10 && Ebin[4]<10) cutID[5] = 1;
        cutID[6] = 1;

        for(int i=0; i<6; i++){
            if(cutID[i] == 1) hevent[i]->Fill(Ebin[i]);
        }
    }
    //hevent->DrawClone();
    for(int i=0;i<6; i++){
    hevent[i]->Scale(1.0/hevent[i]->GetEntries());

        for(int ib=0; ib<10; ib++){
        cout<<hevent[i]->GetBinContent(ib+1)<<"  ";
        }
        cout<<endl;
    }
    /*  gStyle->SetOptStat(0);

        TLatex ptext;

        ptext.SetNDC(1);
        ptext.SetTextFont(43);
        ptext.SetTextSize(15.5);
        ptext.SetTextColor(2);


        gStyle->SetOptTitle(0);
        TH2* hist = new TH2D("hist", "hist", 104, -13, 13, 104, -13, 13);
        TCanvas* c1 = new TCanvas("c1","c1", 400*3, 400*2);
        c1->Divide(3,2);
        int ccnt=0;
        double rg[] = {0.12,0.2,0.2,0.2,0.2,0.3,0.3};
        double rg2[] = {0.05, 0.05, 0.05, 0.05,0.05, 0.05};

        for(int iev=0; iev<t->GetEntries(); iev++){
        t->GetEntry(iev);
        if(iev%100000==0) cout<<iev<<endl;

        if(ccnt==6) break;
        if(b>6) continue;
        if(fabs(pfreGl[ccnt] - rg[ccnt]) >0.02) continue;
        int badd=0;
        for(int ij=0;ij<6; ij++){
        if(ij==ccnt) continue;
        if(pfreGl[ij] > rg2[ij]) badd=1;
        }
        if(badd==1) continue;
        cout<<b<<"  "<<pfreGl[ccnt]<<"  "<<fabs(pfreGl[ccnt] - rg[ccnt])<<"  "<< rg[ccnt]<<"  "<<ccnt<<endl;
        c1->cd(ccnt+1);
        hist->Draw();
        Double_t r = 0.5*sqrt(64/TMath::Pi()/10.);
        TEllipse e;

// cout<<nux->size()<<endl;
for(int in=0; in<600; in++){
e.SetLineColor(kBlack);
e.SetFillColor(0);
e.SetLineWidth(1);
e.SetFillStyle(0);
if(abs(st_part[in])!=1) continue;
e.DrawEllipse(nux[in],nuy[in],r,r,0,360,0,"same");
    //   cout<<nux[in]<<endl;
    }

    TEllipse e1;
    for(int in=0; in<600; in++){
    e1.SetLineColor(kBlack);
    if(st_part[in]==0) continue;
    if(abs(st_part[in])==2)  {
    e1.SetLineStyle(1);
    e1.SetFillStyle(0);
    //      e1.SetFillColor(kGreen);
    if(st_part[in]==2)  e1.SetLineColor(kBlue);
    if(st_part[in]==-2)  e1.SetLineColor(kRed);
    e1.DrawEllipse(nux[in],nuy[in],r,r,0,360,0,"same");
    }
    //   cout<<nux[in]<<endl;
    }
    sprintf(name, "b = %4.2f", b);
    ptext.DrawLatex(0.8, 0.67,name);
    for(int i=0; i<6; i++){
    sprintf(name, "#epsilon_{%d}=%4.2f", i+1, pfreGl[i]);
    ptext.DrawLatex(0.8, 0.7+i*0.04,name);
    }
    e.SetFillColor(0);
    e.SetFillStyle(0);
    e.SetLineColor(1);
e.SetLineStyle(2);
e.SetLineWidth(1);
e.DrawEllipse(b/2,0,6.62,6.62,0,360,0,"same");
e.DrawEllipse(-b/2,0,6.62,6.62,0,360,0,"same");

ccnt++;
}
*/

}
