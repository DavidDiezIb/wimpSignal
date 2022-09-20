#include "TF1.h"
#include "TF2.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "TMath.h"

#include <cstdio>
#include <cstdlib>

using namespace RooFit;

Double_t fFv(Double_t *x,Double_t *par); //Velocity distribution function f(v)*v^2
Double_t WIMPvdistr(Double_t E,Double_t *par); //Integral of the velocity distribution
Double_t fHelm_f(Double_t T,Int_t A); //Helm nuclear form factor
Double_t fMr(Double_t Mx,Double_t Mn); //WIMP-nucleus reduced mass
Double_t fWIMPrate(Double_t *x, Double_t *par); //WIMP event rate in Counts/keV/tonne/year
Double_t fWIMPrateKKD(Double_t *x, Double_t *par); //WIMP event rate in Counts/keV/kg/day



void WIMPsignalOnly() {
   /*// --- Observable --- //
   RooRealVar rawAna_CalibratedThresholdIntegral("rawAna_CalibratedThresholdIntegral","Recoil energy (keV)", 4.9, 40.9);

   // --- Build WIMP PDF --- //
   RooRealVar sigmaX("X sigma(cm^-2)(e-45)", "sigma_x (cm^{-2})", 0.1);//, 0.00001, 100000); //pow(10, -45), pow(10, -50), pow(10, -40)); // WIMP-nucleon cross section.
   RooRealVar massX("X mass (GeV)", "WIMP mass(GeV/c^2)", 200.); //WIMP mass . Will be used as a variable parameter.
   
   TF3 **tfWIMPrate; tfWIMPrate = new TF3 *[1];
   *tfWIMPrate = new TF3("Event rate [keV-1 ton-1 year-1]", fWIMPrate, 4.9, 40.9);  //Function generated from WIMP signal
      
   //RooAbsReal* pdfProj = tfWIMPrate[0].createProjection(RooArgSet(rawAna_CalibratedThresholdIntegral)) 

   RooAbsPdf **pdfWIMPrate; pdfWIMPrate = new RooAbsPdf *[1];
   *pdfWIMPrate = bindPdf(tfWIMPrate[0], rawAna_CalibratedThresholdIntegral, sigmaX, massX);  //PDF generated from WIMP signal function 
   
   // --- Number of expected WIMP events --- //   
   RooAbsReal* signal = bindFunction(tfWIMPrate[0], rawAna_CalibratedThresholdIntegral, sigmaX, massX); //WIMP signal function
   TF1* sigTF1 = signal->asTF(rawAna_CalibratedThresholdIntegral, RooArgList(sigmaX, massX) );
   signal->Print();
   cout << sigTF1->Integral(4.9, 40.9) << endl;
      
   // --- Time factor --- //
   Double_t Tfactor = 1.; // Xenon1T - Trex-Dm simulations
   //Double_t Tfactor = 1016.08/(24*365); //Trex-DM data neon 4bar
   
   // --- Mass factor --- //
   Double_t Mfactor = 1.; // Xenon1T
   //Double_t Mfactor = 66.3/pow(10,6);  //Trex-DM data neon 4bar
   //Double_t Mfactor = 10000./pow(10,6); //Trex-DM simulation 10kg
   //Double_t Mfactor = 300./pow(10,6); //Trex-DM simulation 0.3kg
   
   RooAbsReal* fracInt = pdfWIMPrate[0]->createIntegral(rawAna_CalibratedThresholdIntegral);  //Number of expected events per tonne and year
   RooAbsReal* funcInt = signal->createIntegral(rawAna_CalibratedThresholdIntegral );
   //Double_t numberExpectedEvents = Mfactor*Tfactor*fracInt->getVal();
   cout << "Expected events from pdfWIMPrate[0]: " << fracInt->getVal() << endl;
   cout << "Signal integral from bindFunction(tfWIMPrate[0]: " << funcInt->getVal() << endl;
   *//*
   //// Read data from file ////
   TFile DataFile("BackgroundCalibratedWithCutsFullRange(1309to1343)ThresholdInt.root");
   TTree *tree = (TTree*) DataFile.Get("Tree");
   RooDataSet data("data","dataset with calibrated threshold integral", RooArgSet(rawAna_CalibratedThresholdIntegral), Import(*tree)); 
   RooDataHist* Hdata = data.binnedClone();   

   // --- Generate a toyMC sample from composite PDF --- //
   RooDataSet *dataW = pdfWIMPrate[0]->generate(rawAna_CalibratedThresholdIntegral, 350);
   pdfWIMPrate[0]->fitTo(*dataW); // Perform extended ML fit of composite PDF to toy data
    */
    /*
   // --- Plot --- //
   TCanvas *c_WIMPs=new TCanvas();
   //c_WIMPs->SetLogx();
   //c_WIMPs->SetLogy();
   //rawAna_CalibratedThresholdIntegral.setBins(100);
   RooPlot* Myframe = rawAna_CalibratedThresholdIntegral.frame();
   //data.plotOn(Myframe);
   //Hdata->plotOn(Myframe);
   //dataW->plotOn(Myframe);
   pdfWIMPrate[0]->plotOn(Myframe);
   signal->plotOn(Myframe, RooFit::LineColor(kRed));
   
   //Myframe->SetMinimum(pow(10,-13));
   //Myframe->SetMaximum(pow(10,4));
   //Myframe->SetTitle("Background and 10 GeV 2.5 10^(39)cm^2 WIMP in Neon (1016h, 66.3g)"); // Y sacale in events/keV/year/ton
   Myframe->Draw();
   //c_WIMPs->WaitPrimitive();
   */
   //gStyle->SetPalette(kOcean);
   TCanvas *c=new TCanvas();
   //c->SetLogx();
   //c->SetLogy();
   //c->SetXTitle("Er(keV)");
   //c->SetYTitle("Counts/keV/tonne/year");
   
   
   int const nMasses = 2;
   array<double,nMasses> masses = {10,100}; //5, 10, 30, 50, 100
   double crossSection = 1; // In units of 10^(-45) cm^2
   
   // Event rate in Counts/keV/tonne/year
   TF1 *f1 = new TF1("Rate in counts/keV/tonne/year",fWIMPrate,0.01, 40, 10000);
   f1->SetParameters(crossSection,1.0); //// --- Cross section (E-45), WIMP mass (GeV) --- ////
   //Draw signal
   f1->SetMinimum(0);
   f1->SetMaximum(10);
   f1->DrawCopy("");
   
   //array<double,nMasses> colors = {"kRed", "kRed +10", "kOrange + 10", "kOrange"};
   for(int i=0; i<masses.size(); i++){
     //f1->GetYaxis()->SetRangeUser(pow(10,-7),pow(10,4));
     f1->SetParameters(crossSection,masses[i]);
     f1->SetLineColor(kOrange+i);
     f1->DrawCopy("same");
   }
   
   
   // Event rate in Counts/keV/kg/day
   TCanvas *c2=new TCanvas();
   //c2->SetLogx();
   c2->SetLogy();
   
   TF1 *f2 = new TF1("Rate in counts/keV/kg/day",fWIMPrateKKD,0.01, 50., 10000);
   f2->SetParameters(crossSection,1.0); //// --- Cross section (E-45), WIMP mass (GeV) --- ////
   f2->SetLineColor(kTeal);
   f2->SetMinimum(pow(10,-8));
   f2->SetMaximum(pow(10,-4));
   f2->DrawCopy("");
   //f2->GetHistogram()->GetXaxis()->SetTitle("Er(keV)");
   //f2->GetHistogram()->GetYaxis()->SetTitle("Counts/keV/tonne/year");
   
   for(int i=0; i<masses.size(); i++){
     //f2->GetYaxis()->SetRangeUser(pow(10,-7),pow(10,4));
     //f2->SetAxisRange(pow(10,-7),pow(10,4),"Y");
     f2->SetParameters(crossSection,masses[i]);
     f2->SetLineColor(kTeal+1+i);
     f2->DrawCopy("same");
   }

}



///////////// Functions //////////////////

//Helm form factor
/*Double_t fHelm_f(Double_t T,Int_t A){
    // T = recoil energy in keV
    //Modulus of the Helm Factor
    //Defined as in J.D.LewinP.F.Smith/AstroparticlePhysics6(1996)87-112
    //Femtometers-Effective radius of the target nucleus
    Double_t s=0.9; //Femtometers-Skin thickness of the nucleus
    Double_t R0=sqrt(pow(1.23*std::cbrt(A)-0.60,2)+(7./3.)*pow(TMath::Pi(),2)*0.52*0.52-5*s*s);
    Double_t q=6.92*pow(10,-3)*sqrt(A*T); // Este factor numerico tiene que venir de dividir entre h barra las distancias en fm y de sqrt(2Mn) -> sqrt(2Mn)sqrt(A*Er)/h (hbar = 197)
    //hc = 197 Mev fm = 197*1000 keV fm
    // sqrt(2*0.931)/0.197 = 6.927 Juraría que para que todo esté en keV sobra ese 10^-3. Si lo quito queda demasiado pequeña la señal
    Double_t qR_0=q*R0;
    
    return pow(3*(sin(qR_0)-qR_0*cos(qR_0))/pow(qR_0,3), 2)*exp(-q*q*s*s); //  /2.
}*/
/*
//Helm form factor translating python version
Double_t fHelm_f(Double_t E,Int_t A){
    // E = recoil energy in keV
    Double_t h = 197.3; //MeV fm (hc)
    Double_t nucleon_mass = 0.938; //GeV/c2
    Double_t Mn = nucleon_mass * (double)A;
    Double_t s = 0.9/h; //Femtometers-Skin thickness of the nucleus
    Double_t R = 1.23*std::cbrt(A)-0.6;
    Double_t R0 = sqrt(pow(R,2)+7.*pow(TMath::Pi(),2)*0.52*0.52/3.-5*h*h*s*s) /h;
    Double_t q = sqrt(2*Mn*E); 
    
    return exp(-q*q*s*s)*pow(3*(sin(q*R0)-q*R0*cos(q*R0))/pow(q*R0,3), 2);
}*/

// Nuclear radius parametrization
Double_t fR0(Int_t A){
    Double_t s=0.9; //Femtometers-Skin thickness of the nucleus
    Double_t R0=sqrt(pow(1.23*std::cbrt(A)-0.60,2)+(7./3.)*pow(TMath::Pi(),2)*0.52*0.52-5*s*s);
    return R0;
}

// First Bessel fuction. X in degrees, converted to radians inside this function.
Double_t fBessel1(Double_t x){
    Double_t rad = x*TMath::Pi()/180;
    //return sin(rad)/(rad*rad)-cos(rad)/rad;
    return sin(x)/(x*x)-cos(x)/x;
}

//Helm form factor
Double_t fHelm_f(Double_t E,Int_t A){
    // E = recoil energy in keV
    //Modulus of the Helm Factor
    //Defined as in J.D.LewinP.F.Smith/AstroparticlePhysics6(1996)87-112
    //Femtometers-Effective radius of the target nucleus
    Double_t s=0.9; //Femtometers-Skin thickness of the nucleus
    //Double_t R0=sqrt(pow(1.23*pow(A,1/3)-0.60,2)+(7./3.)*pow(TMath::Pi(),2)*0.52*0.52-5*s*s);
    Double_t q=6.927*pow(10,-3)*sqrt(A*E); //*pow(10,-3) Este factor numerico tiene que venir de dividir entre h barra las distancias en fm y de sqrt(2Mn) -> sqrt(2Mn)sqrt(A*Er)/h (hbar = 197)
    //hc = 197 Mev fm = 197*1000 keV fm
    // sqrt(2*0.931)/0.197 = 6.927 Juraría que para que todo esté en keV sobra ese 10^-3. Si lo quito queda demasiado pequeña la señal. -> No hay que quitarlo, tiene que ver con poner hbar en Mev en lugar de keV qeu es para lo qeu se introduce ese 10^-3, para pasar de keV a Mev en hbar dividiendo: sqrt(A*T)*sqrt(2*Mn)/hbar(en keV)
    Double_t qR_0=q*fR0(A);
    
    //return pow(3*(sin(qR_0)-qR_0*cos(qR_0))/pow(qR_0,3), 2)*exp(-q*q*s*s); //  /2. TMath::Pi()/180 para pasar a radianes
    return pow(3*fBessel1(qR_0)/qR_0, 2)*exp(-q*q*s*s); //  /2. TMath::Pi()/180 para pasar a radianes
}

//Reduced mass
Double_t fMr(Double_t M1,Double_t M2){
    return M1*M2/(M1+M2);
}

//Velocity distribution function f(v)/v
Double_t fFv(Double_t *x,Double_t *par){
    Double_t sigmav=220.*sqrt(1./2.);//v_0/sqrt(2)(km/s) Dispersión en v: sqrt(3/2)*v0 para Javier Garcia Garza, sqrt(1/2)*v0 para Victor
    Double_t v0=220.;//v_0(km/s)
    Double_t vlab=232.;//Vlab(km/s)
    Double_t vesc=544.;//vesc(km/s)
    Double_t Nesc=ROOT::Math::erf(vesc/v0)-(2/sqrt(TMath::Pi()))*(vesc/v0)*exp(-vesc*vesc/(v0*v0));
    Double_t xx=x[0];
    Double_t xmax = std::min( 1., (vesc*vesc - vlab*vlab - xx*xx)/(2*vlab*xx) );
    //Nesci.e.correctioninthenormalizationduetofinitevesc.
    //DefinedatLewisetal.(2.2)
    
    //return 2./(Nesc*Vlab*pow(2*TMath::Pi()*sigmav*sigmav,1./2.))*exp(-Vlab*Vlab/(2*sigmav*sigmav))*(exp(-(xx-Vlab)*xx/(2*sigmav*sigmav))-exp(-(xx+Vlab)*xx/(2*sigmav*sigmav))); //Victor
    //return exp(-xx*xx/(2*sigmav*sigmav))/(xx*sigmav*sqrt(2*TMath::Pi())) ; //// Mía, algo más grande la seccion eficaz estimada.
    return 1./xx/Nesc/(vlab*v0*pow(TMath::Pi(),1./2.)) * ( exp(-(xx-vlab)*(xx-vlab)/(v0*v0))-exp(-(xx*xx+vlab*vlab+2*xx*vlab*xmax)/(v0*v0)) ); // Variant including xmax. vf(v) in the local frame (geen in https://github.com/JelleAalbers/wimprates/blob/master/notebooks/Checks%2C%20plots.ipynb) 
    
    /*Parameter list:
    par[0]=sigmav //WIMP velocity dispersion
    par[1]=v0 //WIMP local circular velocity
    par[2]=Vlab //Laboratory velocity
    par[3]=vesc //Scape velocity
    par[4]=Nesc //Correction to the normalization due to velocity cut off
    */
}


//Velocity distribution function vf(v) in the local frame (geen in https://github.com/JelleAalbers/wimprates/blob/master/notebooks/Checks%2C%20plots.ipynb)
/*Double_t fFv(Double_t xx){
    Double_t sigmav=220.*sqrt(1./2.);//v_0/sqrt(2)(km/s) Dispersión en v: sqrt(3/2)*v0 para Javier Garcia Garza, sqrt(1/2)*v0 para Victor
    Double_t v0=220.;//v_0(km/s)
    Double_t vlab=232.;//Vlab(km/s) (v_earth) in Jelle
    Double_t vesc=544.;//vesc(km/s)
    Double_t Nesc=ROOT::Math::erf(vesc/v0)-(2/sqrt(TMath::Pi()))*(vesc/v0)*exp(-vesc*vesc/(v0*v0)); // k = erf(v_esc/v_0) - 2/np.pi**0.5 * v_esc/v_0 * np.exp(-(v_esc/v_0)**2) (Jelle) // Usa v0 y vesc
    Double_t xmax = std::min( 1., (vesc*vesc - vlab*vlab - xx*xx)/(2*vlab*xx) );
    //Nesci.e.correctioninthenormalizationduetofinitevesc.
    //DefinedatLewisetal.(2.2)

    //return exp(-xx*xx/(2*sigmav*sigmav))/(xx*sigmav*sqrt(2*TMath::Pi())) ; //// Mía, algo más grande la seccion eficaz estimada.
    //return 2./(Nesc*vlab*pow(2*TMath::Pi()*sigmav*sigmav,1./2.))*exp(-vlab*vlab/(2*sigmav*sigmav))*(exp(-(xx-vlab)*xx/(2*sigmav*sigmav))-exp(-(xx+vlab)*xx/(2*sigmav*sigmav))); //Victor
    //return xx/(vlab*v0*pow(TMath::Pi(),1./2.)) * ( exp(-(xx-vlab)*(xx-vlab)/(v0*v0))-exp(-(xx+vlab)*(xx+vlab)/(v0*v0)) ); // Victor/Jelle without normalization. Without normalization constant it adds up 1.
    return xx/Nesc/(vlab*v0*pow(TMath::Pi(),1./2.)) * ( exp(-(xx-vlab)*(xx-vlab)/(v0*v0))-exp(-(xx*xx+vlab*vlab+2*xx*vlab*xmax)/(v0*v0)) ); // Variant including xmax
}*/

//Integral of the velocity distribution
Double_t WIMPvdistr(Double_t E,Double_t *par){
    /*Double_t parVel[8];
    par[0]=A;
    par[1]=par[1]; //WIMP mass
    par[2]=Mn; //Gev/c2
    par[3]=Sigma0;
    par[4]=v0; km/s
    par[5]=Vlab; km/s
    par[6]=Vesc; km/s
    par[7]=Nesc;*/
    Double_t c=3*pow(10,5); // km/s
    
    //Based on the Maxwell-Boltzmann model. Define the lower integration limit.
    //Double_t vmin=c*sqrt(E*pow(10,-6)*pow(par[1]+par[0]*par[2],2)/(2*par[1]*par[1]*par[0]*par[2])); //E is converted into GeV. vmin sale en km/s porqeu hay un factor c
    Double_t vmin = c*sqrt(E*pow(10,-6)*par[0]*par[2]/(2*pow(fMr(par[1], par[0]*par[2]), 2))); //
    //E is converted into GeV. vmin is in km/s
    if(par[6]+par[5]<=vmin){return pow(10,-25);}
    
    //Define integrand. Angular part of the integration is already calculated
    TF1 *f=new TF1("fv",fFv,vmin,par[6]+par[5],0);
    //Wrap the function
    ROOT::Math::WrappedTF1 wIntegrand(*f);
    //Create the Integrator
    ROOT::Math::GaussLegendreIntegrator ig; //Gauss-Legendre quadrature
    
    //Set parameters of the integration
    ig.SetFunction(wIntegrand);
    ig.SetRelTolerance(0.000001);
    
    return ig.Integral(vmin,par[6]+par[5]); //10^-5 factor adjust km/s to cm/s
    
    /*Parameter list:
    par[0]=A
    par[1]=Mx
    par[2]=Mn
    par[3]=sigmav //WIMP velocity dispersion
    par[4]=v0 //WIMP local circular velocity
    par[5]=Vlab //Laboratory velocity
    par[6]=vesc //Escape velocity
    par[7]=Nesc //Correction to the normalization due to velocity cut off
    */
}

//WIMP rate
Double_t fWIMPrate(Double_t *x, Double_t *par){
    //Calculates the WIMP differential rate as a function of recoil energy
    // Parece que aquí x[0]=recoilenergy , par[0]=WIMPcrosssection, par[1]=WIMPmass
    
    //INPUTPARAMETERS/////////////////////////////////////////////////////////////
    //Atomic and mass numbers for the element under study
    //Xenon Z=54,N=77;
    //Argon Z=20,N=20;
    //Neon Z=10,N=10;
    //Hidrogen Z=1,N=1
    //Carbon Z=6,N=6
    //Germanium Z=32,N=32
    
    Int_t Z=54,N=77; 
    Int_t A=N+Z;
    Double_t Mn=0.938; //Nucleon mass(GeV)
    
    //Quantity of material
    Double_t Weight=pow(10,6); // One tonne in grames
    Double_t NAtoms=6.02214*pow(10,23)*Weight/A; //# of atoms in the detector
    
    //Quantity of time
    Double_t Time=3600*24*365; //Seconds, rate will be given in years. One year in seconds
    
    //Dark matter distribution characterization
    Double_t DMDensity=0.3; //(GeV/c^2/cm^3)
    Double_t Sigma0=220./sqrt(2); //WIMP velocity dispersion(km/s)
    Double_t v0=220.; //WIMP mean velocity(km/s)
    Double_t Vlab=232.; //Laboratory velocity(km/s)
    Double_t Vesc=544.; //Scape velocity from the galaxy(km/s)
    Double_t Nesc=ROOT::Math::erf(544./220)-(2/sqrt(TMath::Pi()))*(544./220.)*exp(-544.*544/(220*220));
    //Nesc i.e. correction in the normalization due to finite vesc. Defined at Lewisetal.(2.2)
    
    Double_t parVel[8];
    parVel[0]=A;
    parVel[1]=par[1]; //WIMP mass
    parVel[2]=Mn;
    parVel[3]=Sigma0;
    parVel[4]=v0;
    parVel[5]=Vlab;
    parVel[6]=Vesc;
    parVel[7]=Nesc;
    Double_t c=3*pow(10,5); // km/s
    ///////////////////////////////////////////////////////////////////////////////
    Double_t xx=x[0];
    //return  pow(10,-1)*c*c*(Time*NAtoms*A*Mn*DMDensity*A*A*(par[0]*pow(10, -45))/(2.*par[1]*pow(fMr(par[1],A*Mn),2))*pow(fHelm_f(xx,A),2)*WIMPvdistr(xx,parVel)); // Rate in counts/tonne/year/keV
    return c*c*Time*NAtoms*A*Mn
    *DMDensity/(2.*par[1]*pow(fMr(par[1],A*Mn),2))
    *A*A*(par[0]*pow(10, -45))*pow(1+par[1]/Mn,2)*pow(1+par[1]/(A*Mn),-2)
    *fHelm_f(xx,A)
    *WIMPvdistr(xx,parVel); // Factor de forma ya elevado al cuadrado en la definicion // *pow(1+par[1]/Mn,2)*pow(1+par[1]/(A*Mn),-2) = pow(fMr(mx, A*Mn),2)/pow(fMr(mx, Mn),2) //pow(10,-1)*
}

Double_t fWIMPrateKKD(Double_t *x, Double_t *par){
    return fWIMPrate(x,par)/365000; // rate in days and kg instead of years and tonnes
}
