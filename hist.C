  // http://www.hepforge.org/downloads/tglaubermc
#include <vector>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TSystem.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TNtuple.h>
#include <TNamed.h>
#include <TTime.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TString.h>
#include <TEllipse.h>
#endif

#ifndef _runglauber_
#if !defined(__CINT__) || defined(__MAKECINT__)
#define _runglauber_
#endif

int hist(){


   TH1D Ncoll_cent[10];
   char name[100];
   for (int i=0;i<10;i++)
   {
      sprintf (name,"centrality_bin_%d",i);
      Ncoll_cent[i]= new TH1D (name,name,100,0,1400);
      cout<<"yes"<<endl;
   }

  TFile *f = new TFile ("result.root","recreat");
  for (i=0;i<10;i++) Ncoll_cent[i]->Write();
  f->Close();
  cout<<"It is done"<<endl;

}
