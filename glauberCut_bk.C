#include <cmath>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TProfile.h"
#include "TChain.h"
#include "TRandom.h"

using namespace std;
static const double PI = acos(-1.0);

enum{
    NHAR=6,
    NCENT=20,
    NSTEP=50,
    ETOT=6,
};


void glauberCut(int from=0, int to=100){

    char name[200];

    TChain* t = new TChain("t");
    char dfilelist[100];
    sprintf(dfilelist,"filelist.txt");
    ifstream lis(dfilelist);
    int cnt=0;
    int dnt=0;
    while(!lis.eof())
    {
        string filename;
        lis >> filename;
        if(cnt>=from&&cnt<to)
        {
            cout << filename << endl;
            if(!filename.empty()) t->Add(filename.c_str());
            dnt++;
        }
        cnt++;
    }

    double b;
    int ncoll;
    int npart,npartproj,nparttarg;
    int nHitBbc_n, nHitBbc_s, BBCTrig;
    int Centbin;
    double qBBC_n, qBBC_s, pTrigBBC;
    double vertex, ecc_std,ecc_rp, ecc_part,ecc_partr,r_ollitra,r_geo,r_arith,r_ollitrar,r_geor,r_arithr,e4;

    double e_gl[ETOT],ang_gl[ETOT];
    double re_gl[ETOT],rang_gl[ETOT];

    double pe_gl[ETOT],pang_gl[ETOT];
    double pre_gl[ETOT],prang_gl[ETOT];

    // TBranch *bnux = 0;
    // TBranch *bnuy = 0;
    // TBranch *bst_part = 0;
    std::vector<float> *nux = 0;
    std::vector<float> *nuy = 0;
    std::vector<int> *st_part = 0;

    t->SetBranchAddress("b",            &b   );
    t->SetBranchAddress("Centbin",      &Centbin);
    t->SetBranchAddress("ncoll",        &ncoll);
    t->SetBranchAddress("npart",        &npart);
    t->SetBranchAddress("npartproj",    &npartproj);
    t->SetBranchAddress("nparttarg",    &nparttarg);
    t->SetBranchAddress("ecc_std",      &ecc_std  );
    t->SetBranchAddress("ecc_rp",       &ecc_rp );
    t->SetBranchAddress("ecc_part",     &ecc_part );

    t->SetBranchAddress("e_gl",         e_gl );   //r2
    t->SetBranchAddress("ang_gl",       ang_gl );
    t->SetBranchAddress("re_gl",        re_gl );  //r3,2,3,4,5,6
    t->SetBranchAddress("rang_gl",      rang_gl);

    t->SetBranchAddress("pe_gl",         pe_gl ); //no only Conside Npart
    t->SetBranchAddress("pang_gl",       pang_gl);
    t->SetBranchAddress("pre_gl",        pre_gl);
    t->SetBranchAddress("prang_gl",      prang_gl);

    //t->SetBranchAddress("nux",  &nux);
    //t->SetBranchAddress("nuy",  &nuy);
    //t->SetBranchAddress("st_part",  &st_part);

    int nevents = t->GetEntries();
    cout<<"will run "<<nevents/1000000<<"M  events"<<endl;


    TH1D* hecc [NCENT][NHAR];
    TH1I* hcent[NCENT];
    for(int icent=0; icent<NCENT; icent++){
        sprintf(name,"hcent_%.2d", icent);
        hcent[icent] = new TH1I(name, name ,NCENT, 0-0.5, NCENT-0.5);
    }
    for(int icent=0; icent<NCENT; icent++){
        for(int ihar=0; ihar<NHAR; ihar++){
            sprintf(name, "hecc_ic%.2d_ih%d", icent, ihar);
            hecc[icent][ihar] = new TH1D(name, name, 100000, 0, 1);
        }
    }

    for(int iev=0; iev<nevents; iev++){
        int byte = t->GetEntry(iev);
        if(iev%600000==0) cout<<iev<<endl;
        if(byte==0) continue;
        if(Centbin<0 ||Centbin>19) continue;
        hcent[Centbin] ->Fill(Centbin);
        for(int ihar=0; ihar<NHAR; ihar++){
            hecc[Centbin][ihar] ->Fill(pre_gl[ihar]);
        }
        /*
           for (UInt_t j = 0; j < vpx->size(); ++j) {
           h->Fill(vpx->at(j));
           }
           */
    } //end of events

    float steps[NSTEP] = {0};
    float Ecut[NCENT][NHAR][NSTEP];

    for(int is=0; is<NSTEP-1; is++){
        steps[is] = (is+1)*1.0/NSTEP;
    }
    steps[NSTEP-1] = 1.01;

    for(int icent=0; icent<NCENT; icent++){
        for(int ihar=0; ihar<NHAR; ihar++){
            TH1* htmp;
            sprintf(name,"scale_%d_%d", icent, ihar);
            htmp = (TH1*)hecc[icent][ihar] ->Clone(name);
            htmp-> Scale(1.0/htmp->GetEntries());
            double sum = 0;
            int cnt = 0;
            for(int ib=1; ib<=htmp->GetNbinsX(); ib++){
                sum+= htmp->GetBinContent(ib);
                if(sum>steps[cnt]) { Ecut[icent][ihar][cnt] = htmp->GetBinCenter(ib); cnt++;}
            }
            Ecut[icent][ihar][NSTEP-1] = 1.0;
            cout<<icent<<"  "<<ihar<<"  "<<cnt<<"   "<<sum<<endl;
        }
    }

    ofstream tout;
    sprintf(name,"ECut.txt");
    tout.open(name);
    for(int cent=0; cent<NCENT; cent++){

    tout<<"centb "<<cent<<endl;
    for(int ihar=0; ihar<NHAR; ihar++){
        for(int iq=0; iq<NSTEP; iq++){
            if(iq==0) tout<<"{"<<Ecut[cent][ihar][iq]<<", ";
            else if(iq<NSTEP-1) tout<<Ecut[cent][ihar][iq]<<", ";
            else if(iq ==NSTEP-1) tout<<Ecut[cent][ihar][iq]<<"}, "<<endl;
        }
    }
    }
    

    tout.close();

    TFile* fout = new TFile("glauber_cut.root","recreate");
    for(int icent=0; icent<NCENT; icent++){
        for(int ihar=0; ihar<NHAR; ihar++){
            hecc[icent][ihar] ->Write();
        }
    }
    for(int icent=0; icent<NCENT; icent++){
        hcent[icent] ->Write();
    }
    fout->Close();

}
