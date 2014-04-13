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

#include "ECut_full.C"
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
void runAndSaveNucleons(Int_t n,
        Text_t *sysA="Au",
        Text_t *sysB="Au",
        Double_t signn=42,
        Double_t mind=-1,
        const char *fname="glau_auau_nucleons.root",
        int mode=1,double mu=4.00,double k=1.4);
//This function also save the transverse profile of the nuclear overlap
void runAndSaveNucleonsProf(Int_t n,
        Text_t *sysA="Pb",
        Text_t *sysB="Pb",
        Double_t signn=64,
        Double_t mind=0.4,
        const char *fname="glau_auau_nucleons.root",
        int mode=1,double mu=4.00,double k=1.4);


float bcent[21];
float bave[21];
int   centmode;
int   collmode;
int getCentBin(double b)
{

    for(int i=0; i<20; i++)
    {
        if(b>=bcent[i] && b<bcent[i+1]) return i;
    }
    return -1;

}

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

    void
GetCentrality()
{
    if(collmode==0){
        //Pb+Pb
        bcent[0] = 0;
        bcent[1] = 3.317;
        bcent[2] = 4.687;
        bcent[3] = 5.739;
        bcent[4] = 6.627;
        bcent[5] = 7.409;
        bcent[6] = 8.117;
        bcent[7] = 8.767;
        bcent[8] = 9.373;
        bcent[9] = 9.943;
        bcent[10] = 10.479;
        bcent[11] = 10.991;
        bcent[12] = 11.479;
        bcent[13] = 11.947;
        bcent[14] = 12.399;
        bcent[15] = 12.835;
        bcent[16] = 13.267;
        bcent[17] = 13.713;
        bcent[18] = 14.217;
        bcent[19] = 14.907;
        bcent[20] = 20;
    }
    if(collmode==1){
        //Cu+Au
        cout<<" Cu+Au "<<endl;
        bcent[0] = 0;
        bcent[1] = 2.79;
        bcent[2] = 3.96;
        bcent[3] = 4.86;
        bcent[4] = 5.61;
        bcent[5] = 6.28;
        bcent[6] = 6.88;
        bcent[7] = 7.43;
        bcent[8] = 7.94;
        bcent[9] = 8.43;
        bcent[10] = 8.88;
        bcent[11] = 9.32;
        bcent[12] = 9.73;
        bcent[13] = 10.13;
        bcent[14] = 10.51;
        bcent[15] = 10.89;
        bcent[16] = 11.28;
        bcent[17] = 11.7;
        bcent[18] = 12.18;
        bcent[19] = 12.87;
        bcent[20] = 20;
        bave[0] = 1.86407;
        bave[1] = 3.41421;
        bave[2] = 4.43099;
        bave[3] = 5.24769;
        bave[4] = 5.95781;
        bave[5] = 6.5903;
        bave[6] = 7.1646;
        bave[7] = 7.69294;
        bave[8] = 8.19276;
        bave[9] = 8.66078;
        bave[10] = 9.10694;
        bave[11] = 9.53071;
        bave[12] = 9.93612;
        bave[13] = 10.3264;
        bave[14] = 10.7051;
        bave[15] = 11.0896;
        bave[16] = 11.4909;
        bave[17] = 11.9354;
        bave[18] = 12.4999;
        bave[19] = 13.6975;
    }
    if(collmode==2){
        //Zr+Au
        cout<<" Zr+Au "<<endl;
        bcent[0] = 0;
        bcent[1] = 2.88;
        bcent[2] = 4.12;
        bcent[3] = 5.05;
        bcent[4] = 5.83;
        bcent[5] = 6.53;
        bcent[6] = 7.15;
        bcent[7] = 7.71;
        bcent[8] = 8.26;
        bcent[9] = 8.75;
        bcent[10] = 9.23;
        bcent[11] = 9.69;
        bcent[12] = 10.12;
        bcent[13] = 10.54;
        bcent[14] = 10.95;
        bcent[15] = 11.33;
        bcent[16] = 11.72;
        bcent[17] = 12.14;
        bcent[18] = 12.62;
        bcent[19] = 13.29;
        bcent[20] = 20;
        bave[0] = 1.92443;
        bave[1] = 3.54382;
        bave[2] = 4.60936;
        bave[3] = 5.44982;
        bave[4] = 6.19073;
        bave[5] = 6.85519;
        bave[6] = 7.43473;
        bave[7] = 7.98787;
        bave[8] = 8.50987;
        bave[9] = 8.99597;
        bave[10] = 9.4694;
        bave[11] = 9.91246;
        bave[12] = 10.3366;
        bave[13] = 10.7491;
        bave[14] = 11.1476;
        bave[15] = 11.5379;
        bave[16] = 8.03764;
        bave[17] = 12.3758;
        bave[18] = 12.933;
        bave[19] = 14.0942;
    }
    if(collmode==3){
        //I+Au
        cout<<" I+Au "<<endl;
        bcent[0] = 0;
        bcent[1] = 3.02;
        bcent[2] = 4.32;
        bcent[3] = 5.32;
        bcent[4] = 6.14;
        bcent[5] = 6.85;
        bcent[6] = 7.5;
        bcent[7] = 8.1;
        bcent[8] = 8.67;
        bcent[9] = 9.21;
        bcent[10] = 9.7;
        bcent[11] = 10.19;
        bcent[12] = 10.65;
        bcent[13] = 11.07;
        bcent[14] = 11.48;
        bcent[15] = 11.88;
        bcent[16] = 12.28;
        bcent[17] = 12.72;
        bcent[18] = 13.21;
        bcent[19] = 13.88;
        bcent[20] = 20;
        bave[0] = 2.01163;
        bave[1] = 3.70697;
        bave[2] = 4.84519;
        bave[3] = 5.74302;
        bave[4] = 6.5081;
        bave[5] = 7.18319;
        bave[6] = 7.80587;
        bave[7] = 8.39716;
        bave[8] = 8.94877;
        bave[9] = 9.46321;
        bave[10] = 9.95129;
        bave[11] = 10.4288;
        bave[12] = 10.8651;
        bave[13] = 11.2812;
        bave[14] = 11.6865;
        bave[15] = 12.0845;
        bave[16] = 12.4976;
        bave[17] = 12.9621;
        bave[18] = 13.5241;
        bave[19] = 14.6727;
    }
    if(collmode==4){
        //Dy+Au
        cout<<" Zy+Au "<<endl;
        bcent[0] = 0;
        bcent[1] = 3.15;
        bcent[2] = 4.49;
        bcent[3] = 5.5;
        bcent[4] = 6.36;
        bcent[5] = 7.1;
        bcent[6] = 7.79;
        bcent[7] = 8.41;
        bcent[8] = 8.98;
        bcent[9] = 9.53;
        bcent[10] = 10.06;
        bcent[11] = 10.55;
        bcent[12] = 11.02;
        bcent[13] = 11.46;
        bcent[14] = 11.89;
        bcent[15] = 12.3;
        bcent[16] = 12.72;
        bcent[17] = 13.15;
        bcent[18] = 13.64;
        bcent[19] = 14.32;
        bcent[20] = 20;
        bave[0] = 2.08476;
        bave[1] = 3.86239;
        bave[2] = 5.00562;
        bave[3] = 5.94357;
        bave[4] = 6.74536;
        bave[5] = 7.45618;
        bave[6] = 8.11235;
        bave[7] = 8.70312;
        bave[8] = 9.2646;
        bave[9] = 9.80496;
        bave[10] = 10.3112;
        bave[11] = 10.7892;
        bave[12] = 11.2454;
        bave[13] = 11.6818;
        bave[14] = 12.1013;
        bave[15] = 12.5126;
        bave[16] = 12.9373;
        bave[17] = 13.3938;
        bave[18] = 13.9606;
        bave[19] = 15.1262;
    }
}

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
        vector<double> fdij;        //J.Jia, tranverse distance between two nucleons.
        Bool_t     fInNucleusA;   //=1 from nucleus A, =0 from nucleus B
        Int_t      fNColl;        //Number of binary collisions

    public:
        TGlauNucleon() : fX(0), fY(0), fZ(0), fcX(0), fcY(0), fInNucleusA(0), fNColl(0) {}
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
        //R. Wei: new functions
        vector<double>    GetcX()         const {return fcX;}
        vector<double>    GetcY()         const {return fcY;}
        vector<double>    Getdij()        const {return fdij;}
        void       AddcXY(Double_t x, Double_t y, Double_t z) {fcX.push_back(x); fcY.push_back(y); fdij.push_back(z);}
        void       ResetcXY() {fcX.clear(); fcY.clear(); fdij.clear();}

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
        enum{
            ETOT=6,
        };
        double fAngGl[ETOT],feGl[ETOT];
        double frAngGl[ETOT],freGl[ETOT];

        double pfAngGl[ETOT],pfeGl[ETOT];
        double pfrAngGl[ETOT],pfreGl[ETOT];

       float nux[500];
       float nuy[500];
       int st_part[500];

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
        //P.huo

        //      vector<float>    GetNux()         const {return nux;}
        //      vector<float>    GetNuy()         const {return nuy;}
        //      vector<float>    GetSt_part()     const {return st_part;}

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

//---------------------------------------------------------------------------------
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
    //
    TTime* time = new TTime();
    long int seed = long(time);
    seed = abs(seed);
    cout << "seed = " << seed << "cross section "<<signn<<endl;
    delete gRandom;
    gRandom = new TRandom3(seed);

    TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn);
    mcg->SetMinDistance(mind);
    mcg->SetNBDMode(mode);
    mcg->SetMu(mu);
    mcg->SetK(k);
    double _alpha  = 1.0;
    mcg->SetPhenixBBCMultScalingPower(_alpha);

    TFile *out =new TFile(fname,"recreate");

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

        if(ievent%10000 == 0 ) cout << ievent << endl;

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

        t->Fill();
    }

    out->Write();
    out->Close();
    cout << " done " << endl;

}

//---------------------------------------------------------------------------------
void runAndSaveNucleonsProf(Int_t n,
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
    int Centbin;
    double qBBC_n, qBBC_s, pTrigBBC;
    double vertex, ecc_std,ecc_rp, ecc_part,ecc_partr,r_ollitra,r_geo,r_arith,r_ollitrar,r_geor,r_arithr,e4;

        int cutID[6] = {-1,-1,-1,-1,-1,-1};
    centmode=0;
    collmode =0;
    if((TString(sysA) == "Cu")) collmode=1;
    if((TString(sysA) == "Zr")) collmode=2;
    if((TString(sysA) == "I"))  collmode=3;
    if((TString(sysA) == "Dy")) collmode=4;

    GetCentrality();

    //Random generator
    //
    TTime* time = new TTime();
    long int seed = long(time);
    seed = abs(seed);
    cout << "seed = " << seed << "cross section "<<signn<<endl;
    delete gRandom;
    gRandom = new TRandom3(seed);

    TGlauberMC *mcg=new TGlauberMC(sysA,sysB,signn);
    mcg->SetMinDistance(mind);
    mcg->SetNBDMode(mode);
    mcg->SetMu(mu);
    mcg->SetK(k);
    double _alpha  = 1.0;
    mcg->SetPhenixBBCMultScalingPower(_alpha);

    TFile *out =new TFile(fname,"recreate");

    TTree* t = new TTree("t","event tree");

    //  t->Branch("seed",        &seed,     "seed/I");
    t->Branch("b",            &b,       "b/D");
    t->Branch("Centbin",   &Centbin, "Centbin/I");
    t->Branch("ncoll",        &ncoll,    "ncoll/I");
    t->Branch("npart",        &npart,    "npart/I");
       t->Branch("pre_gl",        &mcg->pfreGl,  "pre_gl[6]/D");
          t->Branch("prang_gl",      &mcg->pfrAngGl,"prang_gl[6]/D");
    t->Branch("nux",  &mcg->nux, "nux[500]/F");
    t->Branch("nuy",  &mcg->nuy, "nuy[500]/F");
    t->Branch("st_part",  &mcg->st_part,"st_part[500]/I");
    t->Branch("cutID",  cutID,"cutID[6]/I");
//    t->Branch("pre_gl",        &mcg->pfreGl,  "pre_gl[6]/D");
    //  t->Branch("prang_gl",      &mcg->pfrAngGl,"prang_gl[6]/D");
    /*    t->Branch("npartproj",    &npartproj,"npartproj/I");
          t->Branch("nparttarg",    &nparttarg,"nparttarg/I");
          t->Branch("ecc_std",      &ecc_std,  "ecc_std/D");
          t->Branch("ecc_rp",       &ecc_rp,   "ecc_rp/D");
          t->Branch("ecc_part",     &ecc_part, "ecc_part/D");

          t->Branch("e_gl",         &mcg->feGl,  "e_gl[6]/D");
          t->Branch("ang_gl",       &mcg->fAngGl,"ang_gl[6]/D");
          t->Branch("re_gl",         &mcg->freGl,  "re_gl[6]/D");
          t->Branch("rang_gl",       &mcg->frAngGl,"rang_gl[6]/D");

          t->Branch("pe_gl",         &mcg->pfeGl,  "pe_gl[6]/D");
          t->Branch("pang_gl",       &mcg->pfAngGl,"pang_gl[6]/D");
          t->Branch("pre_gl",        &mcg->pfreGl,  "pre_gl[6]/D");
          t->Branch("prang_gl",      &mcg->pfrAngGl,"prang_gl[6]/D");
          t->Branch("nux",  &mcg->nux);
          t->Branch("nuy",  &mcg->nuy);
          t->Branch("st_part",  &mcg->st_part);
          */
    cout << " booking the histograms " << endl;
    TH1* hImpactAll = new TH1D("hImpactAll","Impact parameter (all collisions)",2000,-0.005,19.995);
    TH1* hImpactIn  = new TH1D("hImpactIn","Impact parameter (Ncoll >= 1)",2000,-0.005,19.995);
    TH1* hNpart     = new TH1D("hNpart","dN/dN?part! (Ncoll >= 1)",500,-0.5,499.5);
    TH1* hNcoll     = new TH1D("hNcoll","dN/dN?coll! (Ncoll >= 1)",1000,-0.5,999.5);
    char name[200];
 /*   TH2* hpartXY[6][10][6];//Cent Ebin Har
    for(int ic=0;ic<6; ic++){
        for(int ib=0; ib<10;  ib++){
            for(int ih=0; ih<6; ih++){
                sprintf(name,"hpartXY_c%d_b%d_h%d",ic, ib, ih);
                hpartXY[ic][ib][ih] = new TH2D(name,name,400,-10,10,400,-10,10);
                hpartXY[ic][ib][ih]->Sumw2();
            }
        }
    }
*/
    double neve[20] = {0};

    TH1* hevent[6];
    for(int ic=0;ic<6; ic++){
   sprintf(name, "hevent%d", ic); 
    hevent[ic]= new TH1D(name, name,10, 0-0.5, 50-0.5);
    hevent[ic]->Sumw2();
    }
    for(int ic=0;ic<6;ic++){
    cout<<"Cent "<<ic<<"  cut"<<endl;
        for(int ih=0; ih<6; ih++){
        cout<<ecuts[ic][ih][9]<<"  ";
        }
        cout<<endl;
    }

    for(Int_t ievent=0;ievent<n;ievent++){

        if(ievent%100000 == 0 ) cout << ievent << endl;
        double frecc[6], freccAg[6];
        for(int ih=0; ih<6; ih++){
            frecc[ih] = -0.1; freccAg[ih] = 0;
            frecc[ih] = mcg->pfreGl[ih]; freccAg[ih] = mcg->pfrAngGl[ih];
        }
        int badd=0;
        for(int ih=0; ih<6; ih++){ if(frecc[ih]<0) badd=1; }
        if(badd==1) continue;

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

        int centbin = getCentBin(b);
        Centbin = centbin;
        //if(ievent%10000 == 0)cout<<"HI "<<b<<" "<<centbin<<endl;
        if(centbin<0) continue;
        if(centbin>=6) continue;
        neve[centbin]++;
        hImpactAll->Fill(b);
        if(ncoll>=1) hImpactIn->Fill(b);
        hNpart->Fill(npart);
        hNcoll->Fill(ncoll);
        //himp[centbin]->Fill(b);

        //cout<<"centbin  "<<centbin<<"*******************************"<<endl;
        int Ebin[6] = {-1,-1,-1,-1,-1,-1};
        for(int ih=0;ih<6; ih++){ //iHar bin
            Ebin[ih] = findEbin(centbin, frecc[ih], ih);
        }
        for(int id=0; id<6; id++) cutID[id] = -1;
        if(Ebin[5]<10 && Ebin[1]<10 && Ebin[2]<10 && Ebin[3]<10 && Ebin[4]<10) cutID[0] = 1;
        if(Ebin[0]<10 && Ebin[5]<10 && Ebin[2]<10 && Ebin[3]<10 && Ebin[4]<10) cutID[1] = 1;
        if(Ebin[0]<10 && Ebin[1]<10 && Ebin[5]<10 && Ebin[3]<10 && Ebin[4]<10) cutID[2] = 1;
        if(Ebin[0]<10 && Ebin[1]<10 && Ebin[2]<10 && Ebin[5]<10 && Ebin[4]<10) cutID[3] = 1;
        if(Ebin[0]<10 && Ebin[1]<10 && Ebin[2]<10 && Ebin[3]<10 && Ebin[5]<10) cutID[4] = 1;
        if(Ebin[0]<10 && Ebin[1]<10 && Ebin[2]<10 && Ebin[3]<10 && Ebin[4]<10) cutID[5] = 1;
       int ppd = 0;
        for(int id=0; id<6; id++){if(cutID[id]==1) ppd=1; }
        if(ppd==0) continue;
        t->Fill();

        for(int i=0; i<6; i++){
            if(cutID[i] == 1) hevent[i]->Fill(Ebin[i]);
        }
       /*

        double xxm=0.0, yym=0.0;
        int count=0;
        for(int i=0; i<nucleons->GetSize(); i++){
            TGlauNucleon* tgn = (TGlauNucleon*)nucleons->At(i);
            if(tgn->IsSpectator()) continue;
            double x = tgn->GetX();
            double y = tgn->GetY();
            xxm+=x;
            yym+=y;
            count++;
        }
        if(count!=npart) cout<<"wired"<<endl;
        xxm/=count; yym/=count;


        for(int ih=0; ih<6; ih++){
            if(cutID[ih]!=1) continue;
            int cnt =(int)(Ebin[ih]/5);
            if(cnt<0) continue;
            double psi =mcg->pfrAngGl[ih]*(ih+1);
            for(int i=0; i<nucleons->GetSize(); i++){
                TGlauNucleon* tgn = (TGlauNucleon*)nucleons->At(i);
                if(tgn->IsSpectator()) continue;
                double x = tgn->GetX() - xxm;
                double y = tgn->GetY() - yym;
                double x1 = x*cos(psi) + y*sin(psi);
                double y1 =-x*sin(psi) + y*cos(psi);
                hpartXY[centbin][cnt][ih]->Fill(x1, y1); // in the recenter-rotated frame
            }
        }
        */
    }
    out->cd();
    t->Write();
    hNpart->Write();
    hNcoll->Write();
    hImpactAll->Write();
//    for(int ic=0;ic<6; ic++){ for(int ib=0; ib<10;  ib++){ for(int ih=0; ih<6; ih++){ hpartXY[ic][ib][ih]->Write(); } } }
//    for(int ic=0;ic<6; ic++){ for(int ib=0; ib<10;  ib++){ for(int ih=0; ih<6; ih++){cout<<ic<<" "<<ib<<" "<<ih<<" "<<hpartXY[ic][ib][ih]->GetEntries()<<endl;; } } }
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

    else if (TString(name) == "Zr")   {fN = 91;  fR = 4.81;   fA = 0.546;  fW =  0;      fF = 1;}
    else if (TString(name) == "I")    {fN = 127; fR = 5.38;   fA = 0.546;  fW =  0;      fF = 1;}
    else if (TString(name) == "Dy")   {fN = 162; fR = 5.83;   fA = 0.546;  fW =  0;      fF = 1;}

    else if (TString(name) == "W")    {fN = 186; fR = 6.58;  fA = 0.480;  fW =  0;      fF = 1;}
    else if (TString(name) == "Au")   {fN = 197; fR = 6.38;  fA = 0.535;  fW =  0;      fF = 1;}
    else if (TString(name) == "Pb")   {fN = 208; fR = 6.62;  fA = 0.546;  fW =  0;      fF = 1;}
    else if (TString(name) == "U")    {fN = 238; fR = 6.81;  fA = 0.6;    fW =  0;      fF = 1;
        cout<<"Pb+Pb"<<endl;
    }
    else {
        cout << "Could not find nucleus " << name << endl;
        return;
    }
    // new parameterization
    // (from textbook: Povh, Rith, Scholz, Zetsche:
    //  Teilchen und Kerne)
    //ws_radius = 1.07*pow(A,1./3.);
    //ws_diff   = 0.546;

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
    float RA = fANucleus.GetR();
    float RB = fBNucleus.GetR();


    // prepare event
    fANucleus.ThrowNucleons(-bgen*RA/(RA+RB));
    fNucleonsA = fANucleus.GetNucleons();
    fAN = fANucleus.GetN();
    for (Int_t i = 0; i<fAN; i++) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
        nucleonA->SetInNucleusA();
    }
    fBNucleus.ThrowNucleons(bgen*RB/(RA+RB));
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
                        (nucleonA->GetY()+nucleonB->GetY())/2.0, dij
                        );
                nucleonB->AddcXY(
                        (nucleonA->GetX()+nucleonB->GetX())/2.0,
                        (nucleonA->GetY()+nucleonB->GetY())/2.0, dij
                        );
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
    double frcos[ETOT] = {0};
    double frsin[ETOT] = {0};
    double frr[ETOT]={0};
    double fcos[ETOT] = {0};
    double fsin[ETOT] = {0};
    double fr[ETOT]={0};
    double xm=0.0, ym=0.0,xmc=0.0, ymc=0.0,w=0;

    double a = (1.0-0.14)/2.0;
    double b = 0.14;
    //double a =1,b=0;

    for (Int_t i = 0; i<fAN; i++) {
        TGlauNucleon *nucleonA=(TGlauNucleon*)(fNucleonsA->At(i));
        Double_t xA=nucleonA->GetX();
        Double_t yA=nucleonA->GetY();
        fMeanXSystem  += xA;
        fMeanYSystem  += yA;
        fMeanX_A  += xA;
        fMeanY_A  += yA;

        if(nucleonA->IsWounded()) {
            fNpart++;
            fMeanXParts  += xA;
            fMeanYParts  += yA;
            fMeanX2 += xA * xA;
            fMeanY2 += yA * yA;
            fMeanXY += xA * yA;
            xm += xA;
            ym += yA;

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
            fNpart++;
            fMeanXParts  += xB;
            fMeanYParts  += yB;
            fMeanX2 += xB * xB;
            fMeanY2 += yB * yB;
            fMeanXY += xB * yB;
            fNcoll += nucleonB->GetNColl();
            xm += xB;
            ym += yB;

            fNpartTarg++;
        }
    }
    int ncc=0;
    for (Int_t i = 0; i<fBN; i++) {
        TGlauNucleon *nu=(TGlauNucleon*)(fNucleonsB->At(i));
        if(nu->IsWounded()) {
            ncc+=nu->GetNColl();
            for(int nc=0;nc<nu->GetNColl();nc++){
                xmc += nu->GetcX()[nc];
                ymc += nu->GetcY()[nc];
            }
        }
    }
    if(ncc!=fNcoll){
        cout<<"wrong"<<endl;
        exit(1);
    }
    w=a*fNpart+b*fNcoll;
    xm  = a*xm + b*xmc;
    ym  = a*ym + b*ymc;
    xm /= w;
    ym /= w;

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

    /*
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
        for (Int_t j = 0; j<fAN; j++) {
            TGlauNucleon *nu=(TGlauNucleon*)(fNucleonsA->At(j));
            if(!nu->IsWounded()) continue;
            double x=nu->GetX();
            double y=nu->GetY();
            double r=0;
            double rr = sqrt(pow(x-xm,2)+pow(y-ym,2));
            double p = atan2(y-ym,x-xm);
            for(int i=0; i<ETOT; i++){
                if(i==0) r = pow(rr,3)*a;
                else r = pow(rr,i+1)*a;
                frcos[i] += r*cos((i+1)*p);
                frsin[i] += r*sin((i+1)*p);
                frr[i]   += r;
            }
            rr = (pow(x-xm,2)+pow(y-ym,2))*a;
            for(int i=0; i<ETOT; i++){
                fcos[i] += rr*cos((i+1)*p);
                fsin[i] += rr*sin((i+1)*p);
                fr[i]   += rr;
            }
        }
        for (Int_t j = 0; j<fBN; j++) {
            TGlauNucleon *nu=(TGlauNucleon*)(fNucleonsB->At(j));
            if(!nu->IsWounded()) continue;

            double x=nu->GetX();
            double y=nu->GetY();
            double r=0;
            double rr = sqrt(pow(x-xm,2)+pow(y-ym,2));
            double p = atan2(y-ym,x-xm);
            for(int i=0; i<ETOT; i++){
                if(i==0) r = pow(rr,3)*a;
                else r = pow(rr,i+1)*a;
                frcos[i] += r*cos((i+1)*p);
                frsin[i] += r*sin((i+1)*p);
                frr[i]   += r;
            }
            rr = (pow(x-xm,2)+pow(y-ym,2))*a;
            for(int i=0; i<ETOT; i++){
                fcos[i] += rr*cos((i+1)*p);
                fsin[i] += rr*sin((i+1)*p);
                fr[i]   += rr;
            }
            //ncoll
            for(int nc=0;nc<nu->GetNColl();nc++){
                double x = nu->GetcX()[nc];
                double y = nu->GetcY()[nc];
                double r=0;
                double rr = sqrt(pow(x-xm,2)+pow(y-ym,2));
                double p = atan2(y-ym,x-xm);
                for(int i=0; i<ETOT; i++){
                    if(i==0) r = pow(rr,3)*b;
                    else r = pow(rr,i+1)*b;
                    frcos[i] += r*cos((i+1)*p);
                    frsin[i] += r*sin((i+1)*p);
                    frr[i]   += r;
                }
                rr = (pow(x-xm,2)+pow(y-ym,2))*b;
                for(int i=0; i<ETOT; i++){
                    fcos[i] += rr*cos((i+1)*p);
                    fsin[i] += rr*sin((i+1)*p);
                    fr[i]   += rr;
                }
            }
        }
        for(int i=0; i<ETOT; i++){
            frAngGl[i] = (atan2(frsin[i],frcos[i])+TMath::Pi())/double(i+1);//a.la. 1003.0194
            if(frAngGl[i]<-TMath::Pi()/(i+1)) frAngGl[i]+=2*TMath::Pi()/(i+1);//round it to be center around zero
            if(frAngGl[i]> TMath::Pi()/(i+1)) frAngGl[i]-=2*TMath::Pi()/(i+1);//round it to be center around zero
            freGl[i] = sqrt(frcos[i]*frcos[i]+frsin[i]*frsin[i])/frr[i];

            fAngGl[i] = (atan2(fsin[i],fcos[i])+TMath::Pi())/double(i+1);//a.la. 1003.0194
            if(fAngGl[i]<-TMath::Pi()/(i+1)) fAngGl[i]+=2*TMath::Pi()/(i+1);//round it to be center around zero
            if(fAngGl[i]> TMath::Pi()/(i+1)) fAngGl[i]-=2*TMath::Pi()/(i+1);//round it to be center around zero
            feGl[i] = sqrt(fcos[i]*fcos[i]+fsin[i]*fsin[i])/fr[i];
        }

    }else{
        for(int i=0; i<ETOT; i++){
            frAngGl[i] = 0;       freGl[i] =0;
            fAngGl[i] = 0;       feGl[i] =0;
        }
    }
    */
    //*** P. Huo add
    double pfrcos[ETOT] = {0};
    double pfrsin[ETOT] = {0};
    double pfrr[ETOT]={0};
    double pfcos[ETOT] = {0};
    double pfsin[ETOT] = {0};
    double pfr[ETOT]={0};
    double px=0.0, py=0.0, pxm=0.0, pym=0.0;

    int ncount=0;
    for (Int_t j = 0; j<fAN; j++) {
        TGlauNucleon *nu=(TGlauNucleon*)(fNucleonsA->At(j));
        if(!nu->IsWounded()) continue;
        double x=nu->GetX();
        double y=nu->GetY();
        px+=x;
        py+=y;
        ncount++;
    }

    for (Int_t i = 0; i<fBN; i++) {
        TGlauNucleon *nu=(TGlauNucleon*)(fNucleonsB->At(i));
        if(!nu->IsWounded()) continue;
        double x=nu->GetX();
        double y=nu->GetY();
        px+=x;
        py+=y;
        ncount++;
    }

    if(ncount!=fNpart) {cout<<"wired"<<endl;}

    for(int i=0; i<ETOT; i++){
        pfAngGl[i] = -1; pfeGl[i] = -1;
        pfrAngGl[i] = -1; pfreGl[i] =-1;
    }
    if(ncount>1){
        double pxm = px/ncount;
        double pym = py/ncount;

        for (Int_t j = 0; j<fAN; j++) {
            TGlauNucleon *nu=(TGlauNucleon*)(fNucleonsA->At(j));
            if(!nu->IsWounded()) continue;
            double x=nu->GetX();
            double y=nu->GetY();
            double r=0;
            double rr = sqrt(pow(x-pxm,2)+pow(y-pym,2));
            double p = atan2(y-pym,x-pxm);
            for(int i=0; i<ETOT; i++){
                if(i==0) r = pow(rr,3);
                else r = pow(rr,i+1);
                pfrcos[i] += r*cos((i+1)*p);
                pfrsin[i] += r*sin((i+1)*p);
                pfrr[i]   += r;
            }
            rr = (pow(x-xm,2)+pow(y-ym,2));
            for(int i=0; i<ETOT; i++){
                pfcos[i] += rr*cos((i+1)*p);
                pfsin[i] += rr*sin((i+1)*p);
                pfr[i]   += rr;
            }
        }

        for (Int_t j = 0; j<fBN; j++) {
            TGlauNucleon *nu=(TGlauNucleon*)(fNucleonsB->At(j));
            if(!nu->IsWounded()) continue;
            double x=nu->GetX();
            double y=nu->GetY();
            double r=0;
            double rr = sqrt(pow(x-pxm,2)+pow(y-pym,2));
            double p = atan2(y-pym,x-pxm);
            for(int i=0; i<ETOT; i++){
                if(i==0) r = pow(rr,3);
                else r = pow(rr,i+1);
                pfrcos[i] += r*cos((i+1)*p);
                pfrsin[i] += r*sin((i+1)*p);
                pfrr[i]   += r;
            }
            rr = (pow(x-xm,2)+pow(y-ym,2));
            for(int i=0; i<ETOT; i++){
                pfcos[i] += rr*cos((i+1)*p);
                pfsin[i] += rr*sin((i+1)*p);
                pfr[i]   += rr;
            }
        }
        for(int i=0; i<ETOT; i++){
            if(pfr[i]<=0||pfrr[i]<=0) continue;
            pfrAngGl[i] = (atan2(pfrsin[i],pfrcos[i])+TMath::Pi())/double(i+1);//a.la. 1003.0194
            if(pfrAngGl[i]<-TMath::Pi()/(i+1)) pfrAngGl[i]+=2*TMath::Pi()/(i+1);//round it to be center around zero
            if(pfrAngGl[i]> TMath::Pi()/(i+1)) pfrAngGl[i]-=2*TMath::Pi()/(i+1);//round it to be center around zero
            pfreGl[i] = sqrt(pfrcos[i]*pfrcos[i]+pfrsin[i]*pfrsin[i])/pfrr[i];

            pfAngGl[i] = (atan2(pfsin[i],pfcos[i])+TMath::Pi())/double(i+1);//a.la. 1003.0194
            if(pfAngGl[i]<-TMath::Pi()/(i+1)) pfAngGl[i]+=2*TMath::Pi()/(i+1);//round it to be center around zero
            if(pfAngGl[i]> TMath::Pi()/(i+1)) pfAngGl[i]-=2*TMath::Pi()/(i+1);//round it to be center around zero
            pfeGl[i] = sqrt(pfcos[i]*pfcos[i]+pfsin[i]*pfsin[i])/pfr[i];
        }
    }else{
        for(int i=0; i<ETOT; i++){
            pfAngGl[i] = 0; pfeGl[i] = 0;
            pfrAngGl[i] = 0; pfreGl[i] =0;
        }
    }
      for(int i=0;i<500; i++){
        nux[i] = 0; nuy[i] = 0; st_part[i] = 0;
    }

    for (Int_t j = 0; j<fAN; j++) {
        TGlauNucleon *nu=(TGlauNucleon*)(fNucleonsA->At(j));
        int id;
        if(!nu->IsWounded()) {id=-1;}
        else{id=-2;}
        nux[j] = nu->GetX();
        nuy[j] = nu->GetY();
        st_part[j] = id;
    }

    for (Int_t i = 0; i<fBN; i++) {
        TGlauNucleon *nu=(TGlauNucleon*)(fNucleonsB->At(i));
        int id;
        if(!nu->IsWounded()) {id=1;}
        else{id=2;}
        nux[i+fAN] = nu->GetX();
        nuy[i+fAN] = nu->GetY();
        st_part[i+fAN] = id;
    }

    //****P. Huo end of change

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
    if(centmode){
        int centbin = getCentBin(bgen);
        if(centbin<0) centbin=0;
        bgen = bave[centbin];
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
