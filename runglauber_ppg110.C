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
//---------------------------------------------------------------------------------
void runAndSaveNtuple(Int_t n,
                      Text_t *sysA="Au",
                      Text_t *sysB="Au",
                      Double_t signn=42,
                      Double_t mind=-1,
                      const char *fname="glau_auau_ntuple.root");

//---------------------------------------------------------------------------------
//R. Wei: this function is our main routine

int runmode =0;
void runAndSaveNucleons(Int_t n,                    
                        Text_t *sysA="Au",           
                        Text_t *sysB="Au",           
                        Double_t signn=42,           
                        Double_t mind=-1,
                        const char *fname="glau_auau_nucleons.root",
			int mode=1,double mu=4.00,double k=1.4);

void runAndSaveNucleonsEbe(Int_t n,                    
			   Text_t *sysA,           
			   Text_t *sysB,           
			   Double_t signn,
			   Double_t mind,
			   const char *fname);

//---------------------------------------------------------------------------------
//R. Wei: This code is copied from PHENIX Glauber, but extended to 
//take into account the vertex dependence of BBC response.
class TNBDPar : public TNamed
{
private:
  TF1 *nmu0;
  TF1 *nmu1;
  TF1 *smu0;
  TF1 *smu1;
  TF1 *nka0;
  TF1 *nka1;
  TF1 *ska0;
  TF1 *ska1;

public:
  TNBDPar(){
    //specifically for Run7;
    //Sasha Milov provided this vertex dependent NBD (k,mu) functions
    //
    double nmu0_coef[6] = { 5.11315e+00, 6.13945e-02, 8.17016e-03, 5.63883e-04, 1.65567e-05, 1.75983e-07};
    double nmu1_coef[6] = { 5.10197e+00, 7.51575e-02,-1.84043e-03,-6.05494e-05, 3.38720e-06,-4.08911e-08};
    double smu0_coef[6] = { 4.88533e+00,-7.33741e-02,-1.64011e-03, 2.41641e-05, 8.15102e-07,-3.82866e-09};
    double smu1_coef[6] = { 4.89916e+00,-3.63039e-02, 5.67051e-03,-5.15690e-04, 1.89094e-05,-2.46071e-07};
    double nka0_coef[6] = { 1.61234e+00,-7.71948e-02,-1.47291e-02,-1.17055e-03,-4.34974e-05,-6.01865e-07};
    double nka1_coef[6] = { 1.65186e+00,-1.97867e-02, 3.28761e-03,-9.53016e-05,-9.35367e-07, 5.09596e-08};
    double ska0_coef[6] = { 1.61980e+00, 1.59291e-02, 1.67241e-03,-8.28718e-05,-9.18541e-06,-1.86663e-07};
    double ska1_coef[6] = { 1.59072e+00, 7.03501e-02,-1.42325e-02, 1.08322e-03,-3.72544e-05, 4.77481e-07};

    nmu0 = new TF1("nmu0","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5",-30,0);
    nmu1 = new TF1("nmu1","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5",0,30);
    smu0 = new TF1("smu0","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5",-30,0);
    smu1 = new TF1("smu1","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5",0,30);
    nka0 = new TF1("nka0","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5",-30,0);
    nka1 = new TF1("nka1","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5",0,30);
    ska0 = new TF1("ska0","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5",-30,0);
    ska1 = new TF1("ska1","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5",0,30);    

    nmu0->SetParameters(nmu0_coef);
    nmu1->SetParameters(nmu1_coef);
    smu0->SetParameters(smu0_coef);
    smu1->SetParameters(smu1_coef);
    
    nka0->SetParameters(nka0_coef);
    nka1->SetParameters(nka1_coef);
    ska0->SetParameters(ska0_coef);
    ska1->SetParameters(ska1_coef);
  }
  
  ~TNBDPar(){
    delete nmu0;
    delete nmu1;
    delete smu0;
    delete smu1;
    delete nka0;
    delete nka1;
    delete ska0;
    delete ska1;    
  }


  void GetParameters(double& vertex, double& mu_north, double& mu_south, double& k_north, double& k_south  ){
    if(vertex<=0){
      mu_north = nmu0->Eval(vertex);
      mu_south = smu0->Eval(vertex);
      k_north  = nka0->Eval(vertex);
      k_south  = ska0->Eval(vertex);
    }else{
      mu_north = nmu1->Eval(vertex);
      mu_south = smu1->Eval(vertex);
      k_north  = nka1->Eval(vertex);
      k_south  = ska1->Eval(vertex);
    }
    return;
  }

  ClassDef(TNBDPar,1)
};


//---------------------------------------------------------------------------------
//! random numbers from negative binomial distribution
class TRanNBD : public TNamed
{
  double mu, k;    // parameters for negative binomial distr.
  TH1D   hNBD;     // histogram for negative binomial distr.
  int    nNBDBins; // number of bins for hNBD

  void init();

public:
  TRanNBD(double nbd_mu = 4.04, double nbd_k = 2.54);
  void   SetParameters(double, double);
  double getValue(int);
  int    getRandom();

  TNBDPar nbdpar;

  ClassDef(TRanNBD,1)
};


//---------------------------------------------------------------------------------
class TGlauNucleon : public TNamed
{
   private:
      Double32_t fX;            //Position of nucleon
      Double32_t fY;            //Position of nucleon
      Double32_t fZ;            //Position of nucleon
      vector<double> fcX;       //R.Wei, added to store the position of collision.
      vector<double> fcY;       //R.Wei, added to store the position of collision.
      vector<double> fdx;       //R.Wei, added to store the distance of collision.
      vector<double> fdy;       //R.Wei, added to store the distance of collision.
      Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
      Int_t      fNColl;        //Number of binary collisions

public:
  TGlauNucleon() : fX(0), fY(0), fZ(0), fcX(0), fcY(0), fdx(0), fdy(0), fInNucleusA(0), fNColl(0) {}
  virtual   ~TGlauNucleon() {}
  
  void       Collide()            {fNColl++;}
  Int_t      GetNColl()     const {return fNColl;}
  Double_t   GetX()         const {return fX;}
  Double_t   GetY()         const {return fY;}
  Double_t   GetZ()         const {return fZ;}
  Bool_t     IsInNucleusA() const {return fInNucleusA;}
  Bool_t     IsInNucleusB() const {return !fInNucleusA;}
  Bool_t     IsSpectator()  const {return !fNColl;}
  Bool_t     IsWounded()    const {return fNColl;}
  void       Reset()              {fNColl=0; ResetcXY();}
  void       SetInNucleusA()      {fInNucleusA=1;}
  void       SetInNucleusB()      {fInNucleusA=0;}
  void       SetXYZ(Double_t x, Double_t y, Double_t z) {fX=x; fY=y; fZ=z;}
  void       SetX(Double_t x) {fX=x;}
  void       SetY(Double_t y) {fY=y;}
  void       SetZ(Double_t z) {fZ=z;}
  //R. Wei: new functions 
  vector<double>    GetcX()         const {return fcX;}
  vector<double>    GetcY()         const {return fcY;}
  vector<double>    GetdX()         const {return fdx;}
  vector<double>    GetdY()         const {return fdy;}
  void       AddcXY(Double_t x, Double_t y) {fcX.push_back(x); fcY.push_back(y);}
  void       AdddXY(Double_t x, Double_t y) {fdx.push_back(x); fdy.push_back(y);}
  void       ResetcXY() {fcX.clear(); fcY.clear();}
  void       ResetdXY() {fdx.clear(); fdy.clear();}
  void       SetcX(int i, double x) {fcX[i] = x;}
  void       SetcY(int i, double y) {fcY[i] = y;}
  void       SetdX(int i, double x) {fdx[i] = x;}
  void       SetdY(int i, double y) {fdy[i] = y;}
  
  ClassDef(TGlauNucleon,1);
};

//---------------------------------------------------------------------------------
class TGlauNucleus : public TNamed
{
   private:
      Int_t      fN;          //Number of nucleons
      Double_t   fR;          //Parameters of function
      Double_t   fA;          //Parameters of function
      Double_t   fW;          //Parameters of function
      Double_t   fMinDist;    //Minimum separation distance
      Int_t      fF;          //Type of radial distribution
      Int_t      fTrials;     //Store trials needed to complete nucleus
      TF1*       fFunction;   //Probability density function rho(r)
      TObjArray* fNucleons;   //Array of nucleons

      void       Lookup(Text_t* name);

   public:
      TGlauNucleus(Text_t* iname="Au", Int_t iN=0, Double_t iR=0, Double_t ia=0, Double_t iw=0, TF1* ifunc=0);
      virtual ~TGlauNucleus();

      using      TObject::Draw;
      void       Draw(Double_t xs, Int_t col);
      Int_t      GetN()             const {return fN;}
      Double_t   GetR()             const {return fR;}
      Double_t   GetA()             const {return fA;}
      Double_t   GetW()             const {return fW;}
      TObjArray *GetNucleons()      const {return fNucleons;}
      Int_t      GetTrials()        const {return fTrials;}
      void       SetN(Int_t in)           {fN=in;}
      void       SetR(Double_t ir);
      void       SetA(Double_t ia);
      void       SetW(Double_t iw);
      void       SetMinDist(Double_t min) {fMinDist=min;}
      void       ThrowNucleons(Double_t xshift=0.);

      ClassDef(TGlauNucleus,1)
};

//---------------------------------------------------------------------------------
class TGlauberMC : public TNamed
{
private:
  TGlauNucleus fANucleus;       //Nucleus A
  TGlauNucleus fBNucleus;       //Nucleus B
  Double_t     fXSect;          //Nucleon-nucleon cross section
  TObjArray*   fNucleonsA;      //Array of nucleons in nucleus A
  TObjArray*   fNucleonsB;      //Array of nucleons in nucleus B
  Int_t        fAN;             //Number of nucleons in nucleus A
  Int_t        fBN;             //Number of nucleons in nucleus B
  TNtuple*     fnt;             //Ntuple for results (created, but not deleted)
  Double_t     fMeanX2;         //<x^2> of wounded nucleons
  Double_t     fMeanY2;         //<y^2> of wounded nucleons
  Double_t     fMeanXY;         //<xy> of wounded nucleons
  Double_t     fMeanXParts;     //<x> of wounded nucleons
  Double_t     fMeanYParts;     //<x> of wounded nucleons
  Double_t     fMeanXSystem;    //<x> of all nucleons
  Double_t     fMeanYSystem;    //<x> of all nucleons  
  Double_t     fMeanX_A;        //<x> of nucleons in nucleus A
  Double_t     fMeanY_A;        //<x> of nucleons in nucleus A
  Double_t     fMeanX_B;        //<x> of nucleons in nucleus B
  Double_t     fMeanY_B;        //<x> of nucleons in nucleus B
  Double_t     fB_MC;           //Impact parameter (b)
  Int_t        fEvents;         //Number of events with at least one collision
  Int_t        fTotalEvents;    //All events within selected impact parameter range
  Double_t     fBMin;           //Minimum impact parameter to be generated
  Double_t     fBMax;           //Maximum impact parameter to be generated
  Int_t        fMaxNpartFound;  //Largest value of Npart obtained
  Int_t        fNpart;          //Number of wounded (participating) nucleons in current event
  Int_t        fNcoll;          //Number of binary collisions in current event
  Double_t     fSx2;            //Variance of x of wounded nucleons
  Double_t     fSy2;            //Variance of y of wounded nucleons
  Double_t     fSxy;            //Covariance of x and y of wounded nucleons

  Bool_t       CalcResults(Double_t bgen);
  Bool_t       Calc1(Double_t bgen);
  Bool_t       Calc2(Double_t bgen);
  Bool_t       CalcEvent(Double_t bgen);
  //J.Jia NBD input parameter
  int     fnbdmode;
  double  fmu, fk;
  //R. Wei: Added three different eccentricities.
  double       fSDeccen;        //standard eccentricty, such as the one used in PHENIX.
  double       fRPeccen;        //reaction plane eccentricty, consider shift but without rotation.
  double       fPeccen;         //participant eccentricity, calculated with rotation of framework
  
  //R. Wei: Variables to calculate the eccentricities \epsilon_2 and \epsilon_4 in rotated frame 
  Double_t     fSx2r;            //Variance of x of wounded nucleons in rotated
  Double_t     fSy2r;            //Variance of y of wounded nucleons in rotated
  Double_t     fMeanXPartsr;     //<x> of wounded nucleons
  Double_t     fMeanYPartsr;     //<x> of wounded nucleons
  Double_t     feccenr;          //the RP eccentricity calculated in rotated frame should = participant eccentricity
  double       feccen4r;         //the fourth order eccentricity \epsilon_{4} in the rotated frame;
  double       fSx4r;
  double       fSy4r;
  double       fSx2y2r;

  //transverse size of the overlap a la nucl-th/0508009:
  //Harmonic mean of x y variance (except no factor of 2)
  //1/R^2 =  1/<x^2> + 1/<y^2>, characterize the average gradiant or expantion rate.
  double       fR_Ollitra;     //in reaction plane frame
  double       fR_Ollitrar;    //in rotated frame
  //transverse size calcualted from the area: Geometric mean of x y variance
  //R^2=sqrt(<x^2><y^2>)
  double       fR_Geo;         //in reaction plane frame
  double       fR_Geor;        //in rotated frame
  //transeverse size as arithmatic mean of the x, y variance
  //2R^2=<x^2>+<y^2>
  double       fR_Arith;       //in reaction plane frame
  double       fR_Arithr;      //in rotated frame

  Double_t     fang;           // rotate angle;

  //R. Wei: following needed to simulate the BBC response
  TRanNBD nbd;
  double alpha;
  int fNpartProj;
  int fNpartTarg;
  void   SetPhenixBBCNbdParams(const double&, const double& k);
  double vertex;
  TF1* fvertex;

public:
  TGlauberMC(Text_t* NA = "Au", Text_t* NB = "Au", Double_t xsect = 42);
  virtual     ~TGlauberMC() {Reset();}
  
  void         Draw(Option_t* option);
  Double_t     GetB()               const {return fB_MC;}
  Double_t     GetBMin()            const {return fBMin;}
  Double_t     GetBMax()            const {return fBMax;}
  Int_t        GetNcoll()           const {return fNcoll;}
  Int_t        GetNpart()           const {return fNpart;}
  Int_t        GetNpartProj()       const {return fNpartProj;}//J.Jia
  Int_t        GetNpartTarg()       const {return fNpartTarg;}//J.Jia
  Int_t        GetNpartFound()      const {return fMaxNpartFound;}
  TNtuple*     GetNtuple()          const {return fnt;}
  TObjArray   *GetNucleons();

  Double_t     GetTotXSect()        const;
  Double_t     GetTotXSectErr()     const;
  Bool_t       NextEvent(Double_t bgen=-1);
  void         Reset()                    {delete fnt; fnt=0; delete fvertex;}
  void         Run(Int_t nevents);
  void         SetBmin(Double_t bmin)      {fBMin = bmin;}
  void         SetBmax(Double_t bmax)      {fBMax = bmax;}
  void         SetMinDistance(Double_t d)  {fANucleus.SetMinDist(d); fBNucleus.SetMinDist(d);}
  
  //R. Wei: New functions
  void         SetNBDMode(int mode)      {fnbdmode = mode;}
  void         SetMu(double mu)          {fmu = mu;}
  void         SetK(double k)            {fk = k;}

  double GetAng()         {return fang;}     //clockwise is positive;
  double GetSDEccen()     {return fSDeccen;} //standard eccentricty, such as the one used in PHENIX.
  double GetRPEccen()     {return fRPeccen;} //reaction plane eccentricty
  double GetPEccen()      {return fPeccen;}  //participant eccentricity
  double GetrEccen()      {return feccenr;}  //rotated eccentricity = participant eccen
  double GetrEccen4()     {return feccen4r;} //rotated e4
  double GetOllitraSize() {return fR_Ollitra;}
  double GetGeoSize()     {return fR_Geo;}
  double GetArithSize()   {return fR_Arith;}
  double GetOllitraSizer(){return fR_Ollitrar;}
  double GetGeoSizer()    {return fR_Geor;}
  double GetArithSizer()  {return fR_Arithr;}
  void   SetPhenixBBCMultScalingPower(const double&);
  void   GetPhenixBBCInfo(int&, int&, double&,double&,double&,int&);
  double GetVertex() {return vertex;}


  static 
  void       PrintVersion()             {cout << "TGlauberMC " << Version() << endl;}
  static 
  const char *Version()                 {return "v1.1 (adapted to PHENIX)";}
  
  ClassDef(TGlauberMC,1)
};
//---------------------------------------------------------------------------------
void runAndSaveNtuple(Int_t n,
                      Text_t *sysA,
                      Text_t *sysB,
                      Double_t signn,
                      Double_t mind,
                      const char *fname)
{
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);
  mcg->Run(n);
  TNtuple  *nt=mcg->GetNtuple();
  TFile out(fname,"recreate",fname,9);
  if(nt) nt->Write();
  out.Close();
}

float bcent[21],parcent[21];
int centbin;
int getCentBin(double b) {
  for(int i=0; i<20; i++)
    {
      if(b>=bcent[i] && b<bcent[i+1]) return i;
    }
  return -1;
}

void GetCentrality() {
  //Au+Au
  bcent[0] = 0;
  bcent[1] = 3.28;
  bcent[2] = 4.65;
  bcent[3] = 5.7;
  bcent[4] = 6.58;
  bcent[5] = 7.36;
  bcent[6] = 8.07;
  bcent[7] = 8.72;
  bcent[8] = 9.32;
  bcent[9] = 9.89;
  bcent[10] = 10.43;
  bcent[11] = 10.94;
  bcent[12] = 11.42;
  bcent[13] = 11.89;
  bcent[14] = 12.34;
  bcent[15] = 12.77;
  bcent[16] = 13.2;
  bcent[17] = 13.64;
  bcent[18] = 14.15;
  bcent[19] = 14.83;
  bcent[20] = 20;
}
// version 1 of the runglauber, used for ppg110, save event by event summary information
TH2*parta[20],*partb[20];
TH2*colla[20],*collb[20];
TH1*eccst[20],*eccp[20];
TH1* pevt;
void runAndSaveNucleons(Int_t n,                    
                        Text_t *sysA,           
                        Text_t *sysB,           
                        Double_t signn,
                        Double_t mind,
                        const char *fname,
			int mode,double mu,double k)
{
   double b;
   int ncoll;
   int npart,npartproj,nparttarg;
   int nHitBbc_n, nHitBbc_s, BBCTrig;
   double qBBC_n, qBBC_s, pTrigBBC;
   double vertex, ecc_std,ecc_rp, ecc_part,ecc_partr,r_ollitra,r_geo,r_arith,r_ollitrar,r_geor,r_arithr,e4;

   //Random generator
   runmode =0;
   GetCentrality();
   TTime* time = new TTime();
   long int seed = long(time);
   seed = abs(seed);
   cout << "seed = " << seed << endl;
   delete gRandom;
   gRandom = new TRandom3(seed);

   TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn);
   mcg->SetMinDistance(mind);
   mcg->SetNBDMode(mode);
   mcg->SetMu(mu);
   mcg->SetK(k);
   double _alpha  = 1.0;
   mcg->SetPhenixBBCMultScalingPower(_alpha);
   char namet[200];
   TFile *out =new TFile(fname,"recreate");
   for(int i=0;i<20;i++){
     sprintf(namet,"parta_%d",i);
     parta[i] = new TH2D(namet,namet,200,-10,10,200,-10,10);
     sprintf(namet,"partb_%d",i);
     partb[i] = new TH2D(namet,namet,200,-10,10,200,-10,10);
     sprintf(namet,"colla_%d",i);
     colla[i] = new TH2D(namet,namet,200,-10,10,200,-10,10);
     sprintf(namet,"collb_%d",i);
     collb[i] = new TH2D(namet,namet,200,-10,10,200,-10,10);
     sprintf(namet,"eccst_%d",i);
     eccst[i] = new TH1D(namet,namet,100,-1,1);
     sprintf(namet,"eccp_%d",i);
     eccp[i] = new TH1D(namet,namet,100,-1,1);
   }
   pevt = new TH1D("pevt","pevt",21,-0.5,20.5);
   TTree* t = new TTree("t","event tree");
   t->Branch("b",            &b,        "b/D");
   t->Branch("vertex",       &vertex,   "vertex/D");
   t->Branch("ncoll",        &ncoll,    "ncoll/I");
   t->Branch("npart",        &npart,    "npart/I");
   t->Branch("npartproj",    &npartproj,"npartproj/I");
   t->Branch("nparttarg",    &nparttarg,"nparttarg/I");
   t->Branch("ecc_std",      &ecc_std,  "ecc_std/D");
   t->Branch("ecc_rp",       &ecc_rp,   "ecc_rp/D");
   t->Branch("ecc_part",     &ecc_part, "ecc_part/D");
   t->Branch("ecc_partr",    &ecc_partr,"ecc_partr/D");
   t->Branch("r_ollitra",    &r_ollitra,"r_ollitra/D");
   t->Branch("r_geo",        &r_geo,    "r_geo/D");
   t->Branch("r_arith",      &r_arith,  "r_arith/D");
   t->Branch("r_ollitrar",   &r_ollitrar,"r_ollitrar/D");
   t->Branch("r_geor",       &r_geor,    "r_geor/D");
   t->Branch("r_arithr",     &r_arithr,  "r_arithr/D");
   t->Branch("e4",           &e4,       "e4/D");
   t->Branch("qBBC_n",       &qBBC_n,   "qBBC_n/D");
   t->Branch("qBBC_s",       &qBBC_s,   "qBBC_s/D");
   t->Branch("nHitBbc_n",    &nHitBbc_n,"nHitBbc_n/I");
   t->Branch("nHitBbc_s",    &nHitBbc_s,"nHitBbc_s/I");


   for(Int_t ievent=0;ievent<n;ievent++){
     
     if(ievent%1000 == 0 ) cout << ievent << endl;

     //get an event with at least one collision
     while(!mcg->NextEvent()) {}
 
     mcg->GetPhenixBBCInfo(nHitBbc_n,nHitBbc_s,qBBC_n,qBBC_s,pTrigBBC,BBCTrig);     
     //if(!BBCTrig) continue;  // one can cut in the analysis macro

     
     //access, save and (if wanted) print out nucleons
     TObjArray* nucleons = mcg->GetNucleons();
     if(!nucleons) continue;
     //Original Phobos code output all nucleons(next line), it takes too much disk space so we comment it out here.
     //if(out) nucleons->Write(Form("nucleonarray%d",ievent),TObject::kSingleKey); 

     b         = mcg->GetB();      
     ncoll     = mcg->GetNcoll();
     npart     = mcg->GetNpart();
     npartproj = mcg->GetNpartProj();
     nparttarg = mcg->GetNpartTarg();
     vertex    = mcg->GetVertex();
     ecc_std   = mcg->GetSDEccen();
     ecc_rp    = mcg->GetRPEccen();
     ecc_part  = mcg->GetPEccen();
     ecc_partr = mcg->GetrEccen();
     r_ollitra = mcg->GetOllitraSize();
     r_geo     = mcg->GetGeoSize();
     r_arith   = mcg->GetArithSize();
     r_ollitrar= mcg->GetOllitraSizer();
     r_geor    = mcg->GetGeoSizer();
     r_arithr  = mcg->GetArithSizer();
     e4        = mcg->GetrEccen4();
     eccst[centbin]->Fill(ecc_rp);
     eccp[centbin]->Fill(ecc_part);
     t->Fill();
   }

   out->Write();
   out->Close();
   cout << " done " << endl;
   
}
// version 2 of the runglauber, used for ppg110, save all participants and collisions 
static int cevt[20]={0};

void runAndSaveNucleonsEbe(Int_t n,                    
			   Text_t *sysA,           
			   Text_t *sysB,           
			   Double_t signn,
			   Double_t mind,
			   const char *fname){
  runmode =1;
  GetCentrality();
  TTime* time = new TTime();
  long int seed = long(time);
  seed = abs(seed);
  cout << "seed = " << seed << endl;
  delete gRandom;
  gRandom = new TRandom3(seed);
  TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn);
  mcg->SetMinDistance(mind);

  TFile* out =new TFile(fname,"recreate");
  int icen, ncoll, npart;
  float b, rotang;
  float partX[1000],partY[1000],collX[10000],collY[10000],dx[10000],dy[10000];
  char namet[200];
  for(int i=0;i<20;i++){
    sprintf(namet,"parta_%d",i);
    parta[i] = new TH2D(namet,namet,200,-10,10,200,-10,10);
    sprintf(namet,"partb_%d",i);
    partb[i] = new TH2D(namet,namet,200,-10,10,200,-10,10);
    sprintf(namet,"colla_%d",i);
    colla[i] = new TH2D(namet,namet,200,-10,10,200,-10,10);
    sprintf(namet,"collb_%d",i);
    collb[i] = new TH2D(namet,namet,200,-10,10,200,-10,10);
  }
  pevt = new TH1D("pevt","pevt",21,-0.5,20.5);
  TTree* t = new TTree("t","event tree");
  t->Branch("b",            &b,        "b/F");
  t->Branch("icen",         &icen,     "icen/I");
  t->Branch("rotang",       &rotang,   "rotang/F");
  t->Branch("ncoll",        &ncoll,    "ncoll/I");
  t->Branch("npart",        &npart,    "npart/I");
  t->Branch("partX",        &partX,    "partX[npart]/F");
  t->Branch("partY",        &partY,    "partY[npart]/F");
  t->Branch("collX",        &collX,    "collX[ncoll]/F");
  t->Branch("collY",        &collY,    "collY[ncoll]/F");
  t->Branch("dx",           &dx,       "dx[ncoll]/F");
  t->Branch("dy",           &dy,       "dy[ncoll]/F");
  for(Int_t ievent=0;ievent<n;ievent++){
    for(int i=0;i<20;i++){
      if(cevt[i]<1000000) goto cont;
    }
    cout<<"finished"<<endl;
    break;
  cont:
    if(ievent%10000 == 0 ) {
      cout << ievent << endl;
      for(int i=0;i<20;i++){
	if(cevt[i]<1000000){
	  cout<<i<<" "<<cevt[i]<<endl;
	  break;
	}
      }     
    }
    
    while(!mcg->NextEvent()) {}
    

    TObjArray* nucleons = mcg->GetNucleons();
    int cp = 0;
    int cc = 0;
    for(int i=0; i<nucleons->GetSize(); i++){
      TGlauNucleon *nucl = (TGlauNucleon*)(nucleons->At(i));
      if(nucl->IsSpectator()) continue;
      
      partX[cp] = nucl->GetX();
      partY[cp] = nucl->GetY();
      cp++;
      if(nucl->IsInNucleusA()){
	vector<double> cx = nucl->GetcX();
	vector<double> cy = nucl->GetcY();
	vector<double> ddx = nucl->GetdX();
	vector<double> ddy = nucl->GetdY();
	for(size_t j=0; j<cx.size(); j++){
	  collX[cc] = cx[j];
	  collY[cc] = cy[j];
	  dx[cc]    = ddx[j];
	  dy[cc]    = ddy[j];
	  cc++;
	}
      }
    }
    icen   = centbin;
    b      = mcg->GetB();
    rotang = mcg->GetAng();
    npart  = mcg->GetNpart();
    ncoll  = cc;
    
    //cout << cc  <<  "  " << ncoll << endl;
    t->Fill();
  }
  out->Write();
  out->Close();
  cout << " done " << endl;
}


//R. Wei: adapted from Klaus's Glauber code
//---------------------------------------------------------------------------------
ClassImp(TRanNBD)
//---------------------------------------------------------------------------------

//! random numbers from negative binomial distribution
TRanNBD::TRanNBD(double nbd_mu, double nbd_k) {

  mu = nbd_mu;
  k  = nbd_k;

  nNBDBins = 101;
  hNBD = TH1D("hNBD","NBD histogram",nNBDBins,-0.5,nNBDBins-0.5);

  init();
}


//! fill NBD histogram, a la. input from Sasha, stop at 100. Since mean is around 4-6, this is sufficient.
void TRanNBD::init() {
  hNBD.Reset();
  for (int i=0; i<nNBDBins; i++) {
    hNBD.SetBinContent(i+1,getValue(i));
  } 
}

//! Set parameters for negative binomial distr.
void TRanNBD::SetParameters(double nbd_mu, double nbd_k) {
  mu = nbd_mu;
  k  = nbd_k;
  
  init();
}

//! Return the value of the negative binomial distribution
//! with parameters (mu, k) and for given n.
//! Function taken from Sasha Milov.(Reference?)
double TRanNBD::getValue(int n) {

  double F;
  double f;
  
  F  = TMath::Gamma(n + k) / TMath::Gamma(n + 1.);
  F /= TMath::Gamma(k);
  f  = n * TMath::Log(mu/k) - (n + k) * TMath::Log(1.0 + mu/k);
  f  = TMath::Exp(f);
  F *= f;
  
  return F;  
}

//! sample negative binomial distribution
int TRanNBD::getRandom() {
  return int(hNBD.GetRandom() + 0.499999999);
}


//---------------------------------------------------------------------------------
ClassImp(TGlauNucleus)
//---------------------------------------------------------------------------------

TGlauNucleus::TGlauNucleus(Text_t* iname, Int_t iN, Double_t iR, Double_t ia, Double_t iw, TF1* ifunc) : 
   fN(iN),fR(iR),fA(ia),fW(iw),fMinDist(-1),
   fF(0),fTrials(0),fFunction(ifunc),
   fNucleons(0)
{
   if (fN==0) {
      cout << "Setting up nucleus " << iname << endl;
      Lookup(iname);
   }
}

TGlauNucleus::~TGlauNucleus()
{
   if (fNucleons) {
      delete fNucleons;
   }
   delete fFunction;
}

void TGlauNucleus::Draw(Double_t xs, Int_t col)
{
   Double_t r = 0.5*sqrt(xs/TMath::Pi()/10.);
   TEllipse e;
   e.SetLineColor(col);
   e.SetFillColor(0);
   e.SetLineWidth(1);

   for (Int_t i = 0;i<fNucleons->GetEntries();++i) {
      TGlauNucleon* gn = (TGlauNucleon*) fNucleons->At(i);
      e.SetLineStyle(1);
      if (gn->IsSpectator()) e.SetLineStyle(3);
      e.DrawEllipse(gn->GetX(),gn->GetY(),r,r,0,360,0,"");
   }
}

void TGlauNucleus::Lookup(Text_t* name)
{
   SetName(name);

   //Default values for nucleus density function
   //
   if      (TString(name) == "p")    {fN = 1;   fR = 0.6;   fA = 0;      fW =  0;      fF = 0;}
   else if (TString(name) == "d")    {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 1;}
   else if (TString(name) == "dh")   {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 3;}
   else if (TString(name) == "dhh")  {fN = 2;   fR = 0.01;  fA = 0.5882; fW =  0;      fF = 4;}
   else if (TString(name) == "O")    {fN = 16;  fR = 2.608; fA = 0.513;  fW = -0.051;  fF = 1;}
   else if (TString(name) == "Si")   {fN = 28;  fR = 3.34;  fA = 0.580;  fW = -0.233;  fF = 1;}
   else if (TString(name) == "S")    {fN = 32;  fR = 2.54;  fA = 2.191;  fW =  0.16;   fF = 2;}
   else if (TString(name) == "Ca")   {fN = 40;  fR = 3.766; fA = 0.586;  fW = -0.161;  fF = 1;}
   else if (TString(name) == "Ni")   {fN = 58;  fR = 4.309; fA = 0.517;  fW = -0.1308; fF = 1;}
   else if (TString(name) == "Cu")   {fN = 63;  fR = 4.2;   fA = 0.596;  fW =  0;      fF = 1;}
   else if (TString(name) == "W")    {fN = 186; fR = 6.58;  fA = 0.480;  fW =  0;      fF = 1;}
   else if (TString(name) == "Au")   {fN = 197; fR = 6.38;  fA = 0.535;  fW =  0;      fF = 1;}
   else if (TString(name) == "Pb")   {fN = 208; fR = 6.62;  fA = 0.546;  fW =  0;      fF = 1;}
   else if (TString(name) == "U")    {fN = 238; fR = 6.81;  fA = 0.6;    fW =  0;      fF = 1;}
   else {
      cout << "Could not find nucleus " << name << endl;
      return;
   }

   switch (fF)
   {
      case 0: // Proton
         fFunction = new TF1("prot","x*x*exp(-x/[0])",0,10);
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
	fFunction = new TF1(name,"x*x*(1+[2]*(x/[0])**2)/(1+exp((x-[0])/[1]))",0,20);//J.Jia change from 15 to 20 as in PHENIX
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 2: // 3pG
         fFunction = new TF1("3pg","x*x*(1+[2]*(x/[0])**2)/(1+exp((x**2-[0]**2)/[1]**2))",0,15);
         fFunction->SetParameters(fR,fA,fW);
         break;
      case 3: // Hulthen
         fFunction = new TF1("f3","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,10);
         fFunction->SetParameters(1/4.38,1/.85);
         break;
      case 4: // Hulthen HIJING
         fFunction = new TF1("f4","x*x*([0]*[1]*([0]+[1]))/(2*pi*(pow([0]-[1],2)))*pow((exp(-[0]*x)-exp(-[1]*x))/x,2)",0,20);
         fFunction->SetParameters(2/4.38,2/.85);
         break;
      default:
         cerr << "Could not find function type " << fF << endl;
         return;
   }
}

void TGlauNucleus::SetR(Double_t ir)
{
   fR = ir;
   switch (fF)
   {
      case 0: // Proton
         fFunction->SetParameter(0,fR);
         break;
      case 1: // 3pF
         fFunction->SetParameter(0,fR);
         break;
      case 2: // 3pG
         fFunction->SetParameter(0,fR);
         break;
   }
}

void TGlauNucleus::SetA(Double_t ia)
{
   fA = ia;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(1,fA);
         break;
      case 2: // 3pG
         fFunction->SetParameter(1,fA);
         break;
   }
}

void TGlauNucleus::SetW(Double_t iw)
{
   fW = iw;
   switch (fF)
   {
      case 0: // Proton
         break;
      case 1: // 3pF
         fFunction->SetParameter(2,fW);
         break;
      case 2: // 3pG
         fFunction->SetParameter(2,fW);
         break;
   }
}

void TGlauNucleus::ThrowNucleons(Double_t xshift)
{
   if (fNucleons==0) {
      fNucleons=new TObjArray(fN);
      fNucleons->SetOwner();
      for(Int_t i=0;i<fN;i++) {
	 TGlauNucleon *nucleon=new TGlauNucleon(); 
	 fNucleons->Add(nucleon); 
      }
   } 
   
   fTrials = 0;

   Double_t sumx=0;       
   Double_t sumy=0;       
   Double_t sumz=0;       

   Bool_t hulthen = (TString(GetName())=="dh");
   if (fN==2 && hulthen) { //special treatmeant for Hulten

      Double_t r = fFunction->GetRandom()/2;
      Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
      Double_t ctheta = 2*gRandom->Rndm() - 1 ;
      Double_t stheta = sqrt(1-ctheta*ctheta);
     
      TGlauNucleon *nucleon1=(TGlauNucleon*)(fNucleons->At(0));
      TGlauNucleon *nucleon2=(TGlauNucleon*)(fNucleons->At(1));
      nucleon1->Reset();
      nucleon1->SetXYZ(r * stheta * cos(phi) + xshift,
		       r * stheta * sin(phi),
		       r * ctheta);
      nucleon2->Reset();
      nucleon2->SetXYZ(-nucleon1->GetX() + 2*xshift,
		       -nucleon1->GetY(),
		       -nucleon1->GetZ());

      fTrials = 1;
      return;
   }

   for (Int_t i = 0; i<fN; i++) {
      TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
      nucleon->Reset();
      while(1) {
         fTrials++;
	 //J.Jia. I checked that the randomization routine for distributing nucleons are equivalent to PHENIX
         Double_t r = fFunction->GetRandom();
         Double_t phi = gRandom->Rndm() * 2 * TMath::Pi() ;
         Double_t ctheta = 2*gRandom->Rndm() - 1 ;
         Double_t stheta = TMath::Sqrt(1-ctheta*ctheta);
         Double_t x = r * stheta * cos(phi) + xshift;
         Double_t y = r * stheta * sin(phi);      
         Double_t z = r * ctheta;      
         nucleon->SetXYZ(x,y,z);
         if(fMinDist<0) break;
         Bool_t test=1;
         for (Int_t j = 0; j<i; j++) {
            TGlauNucleon *other=(TGlauNucleon*)fNucleons->At(j);
            Double_t xo=other->GetX();
            Double_t yo=other->GetY();
            Double_t zo=other->GetZ();
            Double_t dist = TMath::Sqrt((x-xo)*(x-xo)+
					(y-yo)*(y-yo)+
					(z-zo)*(z-zo));
	    
            if(dist<fMinDist) {
               test=0;
               break;
            }
         }
         if (test) break; //found nucleuon outside of mindist
      }
           
      sumx += nucleon->GetX();
      sumy += nucleon->GetY();
      sumz += nucleon->GetZ();
   }
      
   if(1) { 
     //J.Jia, The following code make sure the center-of-mass at (xshift,0,0).
     //Part of the fluctuation between b and Ncoll, b and Npart has been takenout.
     //I found it leading to a narrower width for Npart, Ncoll if we select centrality
     //by cutting on b (the mean is still very close)
     //However it will not affect the relation between BBCCharge, Ncoll and Npart 
     //since their correlation are unchanged by this recentering procedure.
     //So far I found this is the only difference from PHENIX code.
      sumx = sumx/fN;  
      sumy = sumy/fN;  
      sumz = sumz/fN;  
      for (Int_t i = 0; i<fN; i++) {
         TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleons->At(i));
         nucleon->SetXYZ(nucleon->GetX()-sumx+xshift,//J.Jia Sumx already contains the xshift. so need to cancel it.
                         nucleon->GetY()-sumy,
                         nucleon->GetZ()-sumz);
      }
   }
}

//---------------------------------------------------------------------------------
ClassImp(TGlauberMC)
//---------------------------------------------------------------------------------

TGlauberMC::TGlauberMC(Text_t* NA, Text_t* NB, Double_t xsect) :
   fANucleus(NA),fBNucleus(NB),
   fXSect(0),fNucleonsA(0),fNucleonsB(0),fnt(0),
   fMeanX2(0),fMeanY2(0),fMeanXY(0),fMeanXParts(0),
   fMeanYParts(0),fMeanXSystem(0),fMeanYSystem(0),  
   fMeanX_A(0),fMeanY_A(0),fMeanX_B(0),fMeanY_B(0),fB_MC(0),
   fEvents(0),fTotalEvents(0),fBMin(0),fBMax(0),fMaxNpartFound(0),
   fNpart(0),fNcoll(0),fSx2(0),fSy2(0),fSxy(0),
   fNpartProj(0),fNpartTarg(0)
{
   fBMin = 0;
   fBMax = 20;
   fXSect = xsect;
   
   TString name(Form("Glauber_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
   TString title(Form("Glauber %s+%s Version",fANucleus.GetName(),fBNucleus.GetName()));
   SetName(name);
   SetTitle(title);

   //R.Wei: simple vertex distribution parameterization based on data from Run7
   fvertex = new TF1("fvertex","exp(-(x+0.35)*(x+0.35)/2.0/23.9/23.9)",-30,30);
   //set default NBD parameter
   fnbdmode=1;//do not use vertex dependence.
   fmu = 3.99;//from An461, efficiency 94.5%
   fk  = 1.4;//from An461
}

Bool_t TGlauberMC::CalcEvent(Double_t bgen)
{
   // prepare event
   fANucleus.ThrowNucleons(-bgen/2.);
   fNucleonsA = fANucleus.GetNucleons();
   fAN = fANucleus.GetN();
   for (Int_t i = 0; i<fAN; i++) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      nucleonA->SetInNucleusA();
   }
   fBNucleus.ThrowNucleons(bgen/2.);
   fNucleonsB = fBNucleus.GetNucleons();
   fBN = fBNucleus.GetN();
   for (Int_t i = 0; i<fBN; i++) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      nucleonB->SetInNucleusB();
   }

   // "ball" diameter = distance at which two balls interact
   Double_t d2 = (Double_t)fXSect/(TMath::Pi()*10); // in fm^2

   // for each of the A nucleons in nucleus B
   for (Int_t i = 0; i<fBN; i++) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      for (Int_t j = 0 ; j < fAN ;j++) {
	 TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(j));
         Double_t dx = nucleonB->GetX()-nucleonA->GetX();
         Double_t dy = nucleonB->GetY()-nucleonA->GetY();
         Double_t dij = dx*dx+dy*dy;
         if (dij < d2) {
            nucleonB->Collide();
            nucleonA->Collide();

	    //R.Wei: save the average x,y location of the collision.
	    nucleonA->AddcXY( 
			     (nucleonA->GetX()+nucleonB->GetX())/2.0,
			     (nucleonA->GetY()+nucleonB->GetY())/2.0
			     );
	    nucleonB->AddcXY( 
			     (nucleonA->GetX()+nucleonB->GetX())/2.0,
			     (nucleonA->GetY()+nucleonB->GetY())/2.0
			     );
	    nucleonA->AdddXY(dx,  dy);
	    nucleonB->AdddXY(-dx,-dy);
         }
      }
   }
   
  return CalcResults(bgen);
}

//! set parameters for negative binomial distribution used 
//! in the simulation of the BBC signal
void TGlauberMC::SetPhenixBBCNbdParams(const double& mu, const double& k) {
  nbd.SetParameters(mu, k); // Au+Au at 200 GeV
}


//! set power alpha that determines scaling of BBC hit mult. with N_part
void TGlauberMC::SetPhenixBBCMultScalingPower(const double& a) {
  alpha = a;
}


Bool_t TGlauberMC::CalcResults(Double_t bgen)
{
  if(runmode==0){
    return Calc1(bgen);
  }
  if(runmode==1){
    return Calc2(bgen);
  }
  return false;
}
Bool_t TGlauberMC::Calc2(Double_t bgen){
  //! shift the profile center to center of mass of overlap
  fNpart=0;
  fNcoll=0;
  fMeanX2=0;
  fMeanY2=0;
  fMeanXY=0;
  fMeanXParts=0;
  fMeanYParts=0;
  fSx2=0;
  fSy2=0;
  fSxy=0;

  //
  // first step is to calculate the center of mass
  //   
  for (Int_t i = 0; i<fAN; i++) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    Double_t xA=nucleonA->GetX();
    Double_t yA=nucleonA->GetY();
    
    if(nucleonA->IsWounded()) {
      parta[centbin]->Fill(xA,yA);
      fNpart++;
      fMeanXParts  += xA;
      fMeanYParts  += yA;
      fMeanX2 += xA * xA;
      fMeanY2 += yA * yA;
      fMeanXY += xA * yA;
      vector<double> xx = nucleonA->GetcX();
      vector<double> yy = nucleonA->GetcY();
      vector<double> dx = nucleonA->GetdX();
      vector<double> dy = nucleonA->GetdY();
      for(size_t j=0;j<xx.size();j++ ){
	double dist = dx[j]*dx[j]+dy[j]*dy[j];
	colla[centbin]->Fill(xx[j],yy[j],0.5*exp(-dist/0.64));
      }
    }
  }
  for (Int_t i = 0; i<fBN; i++) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Double_t xB=nucleonB->GetX();
    Double_t yB=nucleonB->GetY();
    
    if(nucleonB->IsWounded()) {
      parta[centbin]->Fill(xB,yB);
      fNpart++;
      fMeanXParts  += xB;
      fMeanYParts  += yB;
      fMeanX2 += xB * xB;
      fMeanY2 += yB * yB;
      fMeanXY += xB * yB;
      vector<double> xx = nucleonB->GetcX();
      vector<double> yy = nucleonB->GetcY();
      vector<double> dx = nucleonB->GetdX();
      vector<double> dy = nucleonB->GetdY();
      for(size_t j=0;j<xx.size();j++ ){
	double dist = dx[j]*dx[j]+dy[j]*dy[j];
	colla[centbin]->Fill(xx[j],yy[j],0.5*exp(-dist/0.64));
      }
      
      //add ncoll
      fNcoll += nucleonB->GetNColl();
    }
  }
  
  if(fNpart<1) return false;
  pevt->Fill(centbin);
 
  fMeanXParts /= fNpart;
  fMeanYParts /= fNpart;
  fMeanX2 /= fNpart;
  fMeanY2 /= fNpart;
  fMeanXY /= fNpart;
   
  fSx2=fMeanX2-(fMeanXParts*fMeanXParts);
  fSy2=fMeanY2-(fMeanYParts*fMeanYParts);
  fSxy=fMeanXY-(fMeanXParts*fMeanYParts);
  
  //participant plane angle
  //
  double l = 0.5*(fSy2+fSx2+sqrt(pow(fSy2-fSx2,2)+4*fSxy*fSxy));
  fang = atan((l-fSy2)/fSxy);
  
  //
  //shift and rotate the center of mass to (0,0);
  //double x, y;
  for (Int_t i = 0; i<fAN; i++) {
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
    Double_t xA=nucleonA->GetX();
    Double_t yA=nucleonA->GetY();
    
    if(nucleonA->IsWounded()) {
      xA -= fMeanXParts;
      yA -= fMeanYParts;      
      nucleonA->SetX(xA);
      nucleonA->SetY(yA);
      double xrot = xA*cos(fang) - yA*sin(fang);
      double yrot = xA*sin(fang) + yA*cos(fang);
      partb[centbin]->Fill(xrot,yrot);
      vector<double> xx = nucleonA->GetcX();
      vector<double> yy = nucleonA->GetcY();
      vector<double> dx = nucleonA->GetdX();
      vector<double> dy = nucleonA->GetdY();
      for(size_t j=0;j<xx.size();j++ ){
	double dist = dx[j]*dx[j]+dy[j]*dy[j];
	xA=xx[j]-fMeanXParts;yA=yy[j]-fMeanYParts;
	xrot = xA*cos(fang) - yA*sin(fang);
	yrot = xA*sin(fang) + yA*cos(fang);
	collb[centbin]->Fill(xrot,yrot,0.5*exp(-dist/0.64));
	//collb[centbin]->Fill(xrot,yrot,0.5);
      }
    }
  }
  for (Int_t i = 0; i<fBN; i++) {
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
    Double_t xB=nucleonB->GetX();
    Double_t yB=nucleonB->GetY();
    
    if(nucleonB->IsWounded()) {
      xB -= fMeanXParts;
      yB -= fMeanYParts;      
      nucleonB->SetX(xB);
      nucleonB->SetY(yB);
      double xrot = xB*cos(fang) - yB*sin(fang);
      double yrot = xB*sin(fang) + yB*cos(fang);
      partb[centbin]->Fill(xrot,yrot);
      vector<double> xx = nucleonB->GetcX();
      vector<double> yy = nucleonB->GetcY();
      vector<double> dx = nucleonB->GetdX();
      vector<double> dy = nucleonB->GetdY();
      for(size_t j=0;j<xx.size();j++ ){
	double dist = dx[j]*dx[j]+dy[j]*dy[j];
	xB=xx[j]-fMeanXParts;yB=yy[j]-fMeanYParts;
	xrot = xB*cos(fang) - yB*sin(fang);
	yrot = xB*sin(fang) + yB*cos(fang);
	collb[centbin]->Fill(xrot,yrot,0.5*exp(-dist/0.64));
      }
    }
  }

  //rotate and shift ncoll
  for(int i=0; i<fAN; i++){     
    TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));    
    for(int j=0; j<(int)((nucleonA->fcX).size()); j++){
      double xc = nucleonA->fcX[j];
      double yc = nucleonA->fcY[j];           
      xc -= fMeanXParts;
      yc -= fMeanYParts;
      nucleonA->SetcX(j,xc);
      nucleonA->SetcY(j,yc);
      //distance do not need to change      
    }
  }

  for(int i=0; i<fBN; i++){     
    TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));    
    for(int j=0; j<(int)((nucleonB->fcX).size()); j++){
      double xc = nucleonB->fcX[j];
      double yc = nucleonB->fcY[j];           
      xc -= fMeanXParts;
      yc -= fMeanYParts;
      nucleonB->SetcX(j,xc);
      nucleonB->SetcY(j,yc);      
    }
  }

  return true;
}

Bool_t TGlauberMC::Calc1(Double_t bgen){

   // calc results for the given event
   fNpart=0;
   fNcoll=0;
   fMeanX2=0;
   fMeanY2=0;
   fMeanXY=0;
   fMeanXParts=0;
   fMeanYParts=0;
   fMeanXSystem=0;
   fMeanYSystem=0;
   fMeanX_A=0;
   fMeanY_A=0;
   fMeanX_B=0;
   fMeanY_B=0;
   fNpartProj = 0;
   fNpartTarg = 0;
   //R.Wei new variables
   fSx2=0;
   fSy2=0;
   fSxy=0;
   fSx2r=0;
   fSy2r=0;
   fMeanXPartsr=0;
   fMeanYPartsr=0;
   fSx4r=0;
   fSy4r=0;
   fSx2y2r=0;
   //J.Jia set default value to be some invalid number when there is no collisions
   fSDeccen = -1;
   fRPeccen = -1;
   fPeccen  = -1;
   feccenr  = -1;
   feccen4r = -1;
   fR_Ollitra  =-1;
   fR_Ollitrar =-1;
   fR_Geo      =-1;
   fR_Geor     =-1;
   fR_Arith    =-1;
   fR_Arithr   =-1;

   for (Int_t i = 0; i<fAN; i++) {
      TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
      Double_t xA=nucleonA->GetX();
      Double_t yA=nucleonA->GetY();
      fMeanXSystem  += xA;
      fMeanYSystem  += yA;
      fMeanX_A  += xA;
      fMeanY_A  += yA;

      if(nucleonA->IsWounded()) {
         parta[centbin]->Fill(xA,yA);
         fNpart++;
         fMeanXParts  += xA;
         fMeanYParts  += yA;
         fMeanX2 += xA * xA;
         fMeanY2 += yA * yA;
         fMeanXY += xA * yA;
	 vector<double> xx = nucleonA->GetcX();
	 vector<double> yy = nucleonA->GetcY();
	 vector<double> dx = nucleonA->GetdX();
	 vector<double> dy = nucleonA->GetdY();
	 for(size_t j=0;j<xx.size();j++ ){
	   double dist = dx[j]*dx[j]+dy[j]*dy[j];
	   colla[centbin]->Fill(xx[j],yy[j],0.5*exp(-dist/0.64));
	   //colla[centbin]->Fill(xx[j],yy[j],0.5);
	 }
	 fNpartProj++;
      }
   }

   for (Int_t i = 0; i<fBN; i++) {
      TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
      Double_t xB=nucleonB->GetX();
      Double_t yB=nucleonB->GetY();
      fMeanXSystem  += xB;
      fMeanYSystem  += yB;
      fMeanX_B  += xB;
      fMeanY_B  += yB;

      if(nucleonB->IsWounded()) {
         parta[centbin]->Fill(xB,yB);
         fNpart++;
         fMeanXParts  += xB;
         fMeanYParts  += yB;
         fMeanX2 += xB * xB;
         fMeanY2 += yB * yB;
         fMeanXY += xB * yB;
         vector<double> xx = nucleonB->GetcX();
         vector<double> yy = nucleonB->GetcY();
	 vector<double> dx = nucleonB->GetdX();
	 vector<double> dy = nucleonB->GetdY();
         for(size_t j=0;j<xx.size();j++ ){
	   double dist = dx[j]*dx[j]+dy[j]*dy[j];
	   colla[centbin]->Fill(xx[j],yy[j],0.5*exp(-dist/0.64));
           //colla[centbin]->Fill(xx[j],yy[j],0.5);
         }

	 fNcoll += nucleonB->GetNColl();

	 fNpartTarg++;
      }
   }


   if (fNpart>0) {
      fMeanXParts /= fNpart;
      fMeanYParts /= fNpart;
      fMeanX2 /= fNpart;
      fMeanY2 /= fNpart;
      fMeanXY /= fNpart;
   }

   if(fAN+fBN>0) {
      fMeanXSystem /= (fAN + fBN);
      fMeanYSystem /= (fAN + fBN);
   }
   if(fAN>0) {
      fMeanX_A /= fAN;
      fMeanY_A /= fAN;
   }

   if(fBN>0) {
      fMeanX_B /= fBN;
      fMeanY_B /= fBN;
   }
   if(fNpart>0){//J.Jia: calculation only make sense when there is collision
     fSx2=fMeanX2-(fMeanXParts*fMeanXParts);
     fSy2=fMeanY2-(fMeanYParts*fMeanYParts);
     fSxy=fMeanXY-fMeanXParts*fMeanYParts;
     
     //R. Wei: following code added to calculate the new variables
     fSDeccen = (fMeanY2-fMeanX2)/(fMeanY2+fMeanX2);
     fRPeccen = (fSy2-fSx2)/(fSy2+fSx2);
     fPeccen  = sqrt(pow(fSy2-fSx2,2)+4*fSxy*fSxy)/(fSx2+fSy2);
   
     double l = 0.5*(fSy2+fSx2+sqrt(pow(fSy2-fSx2,2)+4*fSxy*fSxy));
     fang = atan((l-fSy2)/fSxy);// Rotation angle
     
     //calculate the eccentricity in rotated frame, it should equal to fPeccen
     fSx2r = fSx2*cos(fang)*cos(fang)+fSy2*sin(fang)*sin(fang)-fSxy*sin(2.*fang);
     fSy2r = fSx2*sin(fang)*sin(fang)+fSy2*cos(fang)*cos(fang)+fSxy*sin(2.*fang);
     feccenr = (fSy2r-fSx2r)/(fSy2r+fSx2r);

     //fill the participants in the rotated frame
     for (Int_t i = 0; i<fAN; i++) {
       TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
       if(nucleonA->IsWounded()) {
	 Double_t xA=nucleonA->GetX();
	 Double_t yA=nucleonA->GetY();
	 double xrot = (xA-fMeanXParts)*cos(fang) - (yA-fMeanYParts)*sin(fang);
	 double yrot = (xA-fMeanXParts)*sin(fang) + (yA-fMeanYParts)*cos(fang);
         partb[centbin]->Fill(xrot,yrot);
         vector<double> xx = nucleonA->GetcX();
         vector<double> yy = nucleonA->GetcY();
	 vector<double> dx = nucleonA->GetdX();
	 vector<double> dy = nucleonA->GetdY();
         for(size_t j=0;j<xx.size();j++ ){
	   double dist = dx[j]*dx[j]+dy[j]*dy[j];
	   xA=xx[j];yA=yy[j];
	   xrot = (xA-fMeanXParts)*cos(fang) - (yA-fMeanYParts)*sin(fang);
	   yrot = (xA-fMeanXParts)*sin(fang) + (yA-fMeanYParts)*cos(fang);
           collb[centbin]->Fill(xrot,yrot,0.5*exp(-dist/0.64));
           //collb[centbin]->Fill(xrot,yrot,0.5);
         }
       }
     }
     for (Int_t i = 0; i<fBN; i++) {
       TGlauNucleon *nucleonB=(TGlauNucleon*)(fNucleonsB->At(i));
       if(nucleonB->IsWounded()) {
	 Double_t xB=nucleonB->GetX();
	 Double_t yB=nucleonB->GetY();
	 double xrot = (xB-fMeanXParts)*cos(fang) - (yB-fMeanYParts)*sin(fang);
	 double yrot = (xB-fMeanXParts)*sin(fang) + (yB-fMeanYParts)*cos(fang);
         partb[centbin]->Fill(xrot,yrot);
         vector<double> xx = nucleonB->GetcX();
         vector<double> yy = nucleonB->GetcY();
	 vector<double> dx = nucleonB->GetdX();
	 vector<double> dy = nucleonB->GetdY();
         for(size_t j=0;j<xx.size();j++ ){
	   double dist = dx[j]*dx[j]+dy[j]*dy[j];
           xB=xx[j];yB=yy[j];
	   xrot = (xB-fMeanXParts)*cos(fang) - (yB-fMeanYParts)*sin(fang);
	   yrot = (xB-fMeanXParts)*sin(fang) + (yB-fMeanYParts)*cos(fang);
           collb[centbin]->Fill(xrot,yrot,0.5*exp(-dist/0.64));
           //collb[centbin]->Fill(xrot,yrot,0.5);
         }
       }
     }
     pevt->Fill(centbin);
     //calculate the e4 in rotated frame;
     fMeanXPartsr = fMeanXParts*cos(fang) -fMeanYParts*sin(fang);
     fMeanYPartsr = fMeanXParts*sin(fang) +fMeanYParts*cos(fang);
     for (Int_t i = 0; i<fAN; i++) {
       TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleonsA->At(i));
       Double_t x=nucleon->GetX()*cos(fang)-nucleon->GetY()*sin(fang);
       Double_t y=nucleon->GetX()*sin(fang)+nucleon->GetY()*cos(fang);;
       if(nucleon->IsWounded()) {
         fSx4r   += pow(x- fMeanXPartsr,4);
         fSy4r   += pow(y- fMeanYPartsr,4);
         fSx2y2r += pow((x- fMeanXPartsr)*(y- fMeanYPartsr),2);
       }
     }
     for (Int_t i = 0; i<fBN; i++) {
       TGlauNucleon *nucleon=(TGlauNucleon*)(fNucleonsB->At(i));
       Double_t x=nucleon->GetX()*cos(fang)-nucleon->GetY()*sin(fang);
       Double_t y=nucleon->GetX()*sin(fang)+nucleon->GetY()*cos(fang);;
       if(nucleon->IsWounded()) {
         fSx4r   += pow( x- fMeanXPartsr,4);
         fSy4r   += pow( y- fMeanYPartsr,4);
         fSx2y2r += pow((x- fMeanXPartsr)*(y- fMeanYPartsr),2);
       }
     } 
     fSx4r   /= fNpart;
     fSy4r   /= fNpart;
     fSx2y2r /= fNpart;
   
     feccen4r = 1-8*fSx2y2r/(fSx4r+fSy4r+2*fSx2y2r);
     
     //transver size in reaction plane frame. 
     fR_Ollitra = 1.0/(sqrt(1.0/fSx2+1.0/fSy2));
     fR_Geo     = pow(fSx2*fSy2,0.25);
     fR_Arith   = sqrt((fSx2+fSy2)*0.5);
     //transver size in rotated frame. 
     fR_Ollitrar = 1.0/(sqrt(1.0/fSx2r+1.0/fSy2r));
     fR_Geor     = pow(fSx2r*fSy2r,0.25);
     fR_Arithr   = sqrt((fSx2r+fSy2r)*0.5);

     //R. Wei: end of change
     /////////////////////////////////////////////////////////////////////////////////   
   }

   fB_MC = bgen;

   fTotalEvents++;
   if (fNpart>0) fEvents++;

   if (fNpart==0) return kFALSE;
   if (fNpart > fMaxNpartFound) fMaxNpartFound = fNpart;

   return kTRUE;
}

void TGlauberMC::Draw(Option_t* /*option*/)
{
   fANucleus.Draw(fXSect, 2);
   fBNucleus.Draw(fXSect, 4);

   TEllipse e;
   e.SetFillColor(0);
   e.SetLineColor(1);
   e.SetLineStyle(2);
   e.SetLineWidth(1);
   e.DrawEllipse(GetB()/2,0,fBNucleus.GetR(),fBNucleus.GetR(),0,360,0);
   e.DrawEllipse(-GetB()/2,0,fANucleus.GetR(),fANucleus.GetR(),0,360,0);
}

Double_t TGlauberMC::GetTotXSect() const
{
   return (1.*fEvents/fTotalEvents)*TMath::Pi()*fBMax*fBMax/100;
}

Double_t TGlauberMC::GetTotXSectErr() const
{
   return GetTotXSect()/TMath::Sqrt((Double_t)fEvents) * 
      TMath::Sqrt(Double_t(1.-fEvents/fTotalEvents));
}

TObjArray *TGlauberMC::GetNucleons() 
{
   if(!fNucleonsA || !fNucleonsB) return 0;
   fNucleonsA->SetOwner(0);
   fNucleonsB->SetOwner(0);
   TObjArray *allnucleons=new TObjArray(fAN+fBN);
   allnucleons->SetOwner();
   for (Int_t i = 0; i<fAN; i++) {
      allnucleons->Add(fNucleonsA->At(i));
   }
   for (Int_t i = 0; i<fBN; i++) {
      allnucleons->Add(fNucleonsB->At(i));
   }
   return allnucleons;
}

Bool_t TGlauberMC::NextEvent(Double_t bgen)
{
  //   if(bgen<0) 
  bgen = TMath::Sqrt((fBMax*fBMax-fBMin*fBMin)*gRandom->Rndm()+fBMin*fBMin);
  fB_MC = bgen;
  
  centbin = getCentBin(bgen);
  if(centbin<0) centbin=0;
  if(runmode==1){
    if(cevt[centbin]>1000000) return false;
    cevt[centbin]++;
  }
  vertex = fvertex->GetRandom();

  return CalcEvent(bgen);
}

void TGlauberMC::Run(Int_t nevents)
{
   TString name(Form("nt_%s_%s",fANucleus.GetName(),fBNucleus.GetName()));
   TString title(Form("%s + %s (x-sect = %d mb)",fANucleus.GetName(),fBNucleus.GetName(),(Int_t) fXSect));
   if (fnt == 0) {
      fnt = new TNtuple(name,title,
                        "Npart:Ncoll:B:MeanX:MeanY:MeanX2:MeanY2:MeanXY:VarX:VarY:VarXY:MeanXSystem:MeanYSystem:MeanXA:MeanYA:MeanXB:MeanYB");
      fnt->SetDirectory(0);
   }

   for (int i = 0;i<nevents;i++) {

      while(!NextEvent()) {}

      Float_t v[18];
      v[0]  = GetNpart();
      v[1]  = GetNcoll();
      v[2]  = fB_MC;
      v[3]  = fMeanXParts;
      v[4]  = fMeanYParts;
      v[5]  = fMeanX2;
      v[6]  = fMeanY2;
      v[7]  = fMeanXY;
      v[8]  = fSx2;
      v[9]  = fSy2;
      v[10] = fSxy;
      v[11] = fMeanXSystem;
      v[12] = fMeanYSystem;
      v[13] = fMeanX_A;
      v[14] = fMeanY_A;
      v[16] = fMeanX_B;
      v[17] = fMeanY_B;

      fnt->Fill(v);

      if (!(i%50)) cout << "Event # " << i << " x-sect = " << GetTotXSect() << " +- " << GetTotXSectErr() << " b        \r" << flush;  
   }
   cout << endl << "Done!" << endl;
}

//R. Wei: add code for BBC.
//! Returns the simulated BBC charge signal and the number of BBC hits, this is taken from Klaus PHENIX glauber code (reference?)
void TGlauberMC::GetPhenixBBCInfo(int& nHitBbc_n, int& nHitBbc_s, 
                                   double& qBBC_n, double& qBBC_s,
                                   double& ptrig_aa, int& BBCTrig) {

  double mu_north, mu_south, k_north, k_south;
  if(fnbdmode==0){  
    nbd.nbdpar.GetParameters(vertex,mu_north,mu_south,k_north,k_south);
  }else{
    mu_north = fmu; mu_south =fmu;
    k_north  = fk; k_south   =fk;
  }

  SetPhenixBBCNbdParams(mu_north,k_north);
  nHitBbc_n = 0;
  for (int i=0; i <fNpartProj; i++) {
    nHitBbc_n += nbd.getRandom();
  }

  SetPhenixBBCNbdParams(mu_south,k_south);
  nHitBbc_s = 0;  
  for (int i=0; i <fNpartTarg; i++) {
    nHitBbc_s += nbd.getRandom();
  }
  
  // get BBC trigger decisions
  BBCTrig = 0;
  if (nHitBbc_n>1 && nHitBbc_s>1) BBCTrig = 1;

  //cout << "NpartProj = " << fNpartProj <<", nHitBbc_n = " << nHitBbc_n << "   ";
  //cout << "NpartTarg = " << fNpartTarg <<", nHitBbc_s = " << nHitBbc_s << endl;

  const double ldmean  = 75; // pC
  const double ldsigma = 10; // pC

  const double qmax    = 2*ldmean;

  // determine BBC charged signal without consideration of
  // of multiple photomultiplier hits (default)
  qBBC_n = 0.;
  for(int i=0; i<nHitBbc_n; i++) {
    double q = gRandom->Landau(ldmean,ldsigma);
    if (q<0) q = 0; // take care of ROOT bug that produces
                    // -inf as result of Landau sampling 
    if (q<qmax) qBBC_n += q;
    else qBBC_n += qmax;
  }

  qBBC_s = 0.;
  for(int i=0; i<nHitBbc_s; i++) {
    double q = gRandom->Landau(ldmean,ldsigma);
    if (q<0) q = 0; // take care of ROOT bug that produces
                    // -inf as result of Landau sampling 
    if (q<qmax) qBBC_s += q;
    else qBBC_s += qmax;
  }

  ptrig_aa = 0.;
}


#endif
