#include<iostream>
#include"TFile.h"
using namespace std;
void run(int nev, int lis){
    char name[100];
    sprintf(name,"/direct/phenix+plhf/mzhou/phuo/glauber2/glauber_Pb_Pb_%.2d.root", lis);
    //sprintf(name,"output/glauber_Pb_Pb_%.2d.root", lis);
    gSystem->Load("runglauber_saveprof4_C.so");
    runAndSaveNucleonsProf(nev,"Pb","Pb", 72, 0.4, name, 1,4.00,1.4);
    cout<<name<<endl;
}
