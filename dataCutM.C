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
    NORD=12,
    MAXTRK=100000,
    NETA = 20,

    RDET=3,
    CDET=3,
    NA = 1,  //6 for ihar=2 or 3 + inclusive
    NQ=10,
};
TH1* hcent;

TH1* hqFcal[CDET][NHAR];
TH2* hQPsi[NA][NHAR];

static TFile* fileout=0;
void Init(int cent0, int cent1){
    char name[200];
    sprintf(name,"output/cut%d_%d.root", cent0,  cent1);
    cout<<name<<endl;
    fileout = new TFile(name, "recreate");

    for(int id=0; id<NA; id++){
        for(int ih=0; ih<NHAR; ih++){
            sprintf(name,"hQPsi_cut%d_har%d", id, ih+1);
            hQPsi[id][ih] = new TH2D(name, name, 1000, -PI, PI, 1000, -1, 1);
            hQPsi[id][ih] ->Sumw2();
        }
    }

    for(int idet=0; idet<CDET; idet++){
        for(int ihar=0; ihar<NHAR; ihar++){
            sprintf(name,"hq_det%d_har%d",idet, ihar+1);
            hqFcal[idet][ihar] = new TH1D(name,name, 100000, 0.0,0.5);
            hqFcal[idet][ihar] ->Sumw2();
        }
    }

    sprintf(name, "cent%.2d_%.2d", cent0, cent1);
    hcent = new TH1I(name,name, 100, 0-0.5,100-0.5);
}


void dataCutM(int cent0=20, int cent1=25, int from=0, int to=34){

    char name[200];
    int bad, centrality;
    float trk_Q   [RDET][NHAR];
    float trk_Psi [RDET][NHAR];
    float cal_Q   [CDET][NHAR];
    float cal_Psi [CDET][NHAR];

    TChain* tree = new TChain("tree");
    char dfilelist[100];
    sprintf(dfilelist,"outlist.txt");
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
            if(!filename.empty()) tree->Add(filename.c_str());
            dnt++;
        }
        cnt++;
    }

    tree->SetBranchAddress("bad",  &bad);
    tree->SetBranchAddress("centrality",   &centrality);
    tree->SetBranchAddress("trk_Q",   trk_Q   );
    tree->SetBranchAddress("trk_Psi", trk_Psi );
    tree->SetBranchAddress("cal_Q",   cal_Q   );
    tree->SetBranchAddress("cal_Psi", cal_Psi );

    float Trk_Q  [RDET][NHAR];
    float Trk_Psi[RDET][NHAR];
    float Cal_Q  [CDET][NHAR];
    float Cal_Psi[CDET][NHAR];

    Init(cent0, cent1);

    cout<<"Run Centrality "<<cent0<<"   "<<cent1<<endl;
    int nevents = tree->GetEntries();
    cout<<"nevents  "<<nevents<<endl;
       // for(int iev=0; iev<2000000; iev++){
    float maxqcut[CDET][NHAR] = {{0}};
        for(int ihar=0; ihar<NHAR; ihar++){ for(int idet=0; idet<CDET; idet++){ maxqcut[idet][ihar] = 0; }}
        for(int iev=0; iev<nevents; iev++){
        tree->GetEntry(iev);
        if(iev%1000000==0) cout<<iev<<endl;

        if(centrality>=cent1 || centrality<cent0) continue;

        hcent->Fill(centrality);
        for(int idet = 0; idet<RDET; idet++){
            for(int ihar=0; ihar<NHAR; ihar++){
                Trk_Q[idet][ihar]   = trk_Q[idet][ihar];
                Trk_Psi[idet][ihar] = trk_Psi[idet][ihar];
            }
        }

        for(int idet = 0; idet<CDET; idet++){
            for(int ihar=0; ihar<NHAR; ihar++){
                Cal_Q[idet][ihar]   = cal_Q[idet][ihar];
                Cal_Psi[idet][ihar] = cal_Psi[idet][ihar];
            }
        }
        double dQ[]   = {0,0,0,0,0,0,0};
        double dPsi[] = {0,0,0,0,0,0,0};

        for(int ih=0; ih<NHAR;ih++){
            dQ[ih]   = Cal_Q[0][ih] - Cal_Q[1][ih];
            dPsi[ih] = (Cal_Psi[0][ih] - Cal_Psi[1][ih]);
            dPsi[ih] = atan2(sin(dPsi[ih]), cos(dPsi[ih]));
        }

        for(int ihar=0; ihar<NHAR; ihar++){
            for(int idet=0; idet<CDET; idet++){
                if(Cal_Q[idet][ihar] >maxqcut[idet][ihar]) maxqcut[idet][ihar] = Cal_Q[idet][ihar];
                hqFcal[idet][ihar] ->Fill(Cal_Q[idet][ihar]);
            }
            hQPsi[NA-1][ihar] ->Fill(dPsi[ihar],dQ[ihar]);
        }

    }//end of events

    TH1* htmp;
    float steps[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.01};
    float qcut[CDET][NHAR][NQ];
    for(int idet=0; idet<CDET; idet++){
        for(int ihar=0; ihar<NHAR; ihar++){
            sprintf(name,"scale_%d_%d", idet, ihar);
            htmp = (TH1*)hqFcal[idet][ihar] ->Clone(name);
            htmp-> Scale(1.0/htmp->GetEntries());
            htmp-> SetTitle(name);
            double sum = 0;
            int cnt = 0;
            for(int ib=1; ib<=htmp->GetNbinsX(); ib++){
                sum+= htmp->GetBinContent(ib);
                if(sum>steps[cnt]) { qcut[idet][ihar][cnt] = htmp->GetBinCenter(ib); cnt++;}
            }
            qcut[idet][ihar][9] = maxqcut[idet][ihar] + 0.00001;
            cout<<idet<<"  "<<ihar<<"  "<<cnt<<"   "<<sum<<endl;
        }
    }

    ofstream tout;
    sprintf(name,"cutInfo1_%d_%d.txt", cent0, cent1);
    tout.open(name);
   tout<<"cent "<<cent0<<"--"<<cent1<<endl;
    for(int idet=0; idet<CDET; idet++){
        tout<<"IDET = "<<idet<<endl;
        for(int ihar=0; ihar<NHAR; ihar++){
            for(int iq=0; iq<NQ; iq++){
                if(iq==0) tout<<"{"<<qcut[idet][ihar][iq]<<", ";
                else if(iq<9) tout<<qcut[idet][ihar][iq]<<", ";
                else if(iq ==9) tout<<qcut[idet][ihar][iq]<<"}, "<<endl;
            }
        }
    }
    tout.close();
    fileout->Write();
    fileout->Close();

    }
