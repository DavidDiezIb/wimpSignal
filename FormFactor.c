#include "TF1.h"
#include "TF2.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "TMath.h"
//#include "/media/storage/home/davidp/git/units/units/units.hpp"
using namespace std;
/*using namespace units;
using namespace units::domains;
using namespace units::precise;
using namespace units::constants;
using namespace units::precise::energy;
*/

#include <cstdio>
#include <cstdlib>
#include <string>




Double_t fHelm_f(Double_t T,Int_t A); //Helm nuclear form factor
Double_t fR0(Int_t A); // Nuclear radius parametrization
Double_t fHelm_fq(Double_t q,Int_t A); // In q units (fm^-1)(?)
Double_t fHelm_fp(Double_t E,Int_t A);
Double_t fBessel1(Double_t x);



void FormFactor() {
    //Atomic and mass numbers for the element under study
    //Xenon Z=54,N=77;
    //Argon Z=20,N=20;
    //Neon Z=10,N=10;
    //Hidrogen Z=1,N=1
    //Carbon Z=6,N=6
    //Germanium Z=32,N=32
    
    Int_t Z=54,N=77; 
    Int_t A=N+Z;
    
    cout << "R0: " << fR0(A) << endl;
    
    double xmin = 10;
    double xmax = 10000;
    int n = 100000;
    double step = (xmax-xmin)/(double)n;
    double x[n], y[n], yp[n];
    
    for (int i=0;i<n;i++) {
      x[i] = xmin +i*step;
      y[i] = fHelm_f(x[i], A);
      yp[i] = fHelm_fp(x[i], A);
    }
    
    
    Double_t h = 197.3; //MeV fm (hc)
    Double_t nucleon_mass = 0.938; //GeV/c2
    Double_t Mn = nucleon_mass * (double)A;
    Double_t s = 0.9/h; //Femtometers-Skin thickness of the nucleus
    Double_t R = 1.23*pow((double)A,1/3)-0.6;
    Double_t R0 = sqrt(pow(R,2)+7.*pow(TMath::Pi(),2)*0.52*0.52/3.-5*h*h*s*s) /h;
    
    /*
    // Using units package:
    precise_measurement hc = 197.3 *units::precise::mega*eV * units::precise::femto*units::precise::m; //MeV fm (hc)
    precise_measurement nucleon_mass = 0.938 *units::precise::giga*eV / (constants::c*constants::c); //GeV/c2
    precise_measurement Mn = nucleon_mass * (double)A;
    Double_t s = 0.9; //Femtometers-Skin thickness of the nucleus
    Double_t R = 1.23*std::cbrt(A)-0.6;
    Double_t R0 = sqrt(pow(R,2)+7.*pow(TMath::Pi(),2)*0.52*0.52/3.-5*s*s); // /hc
    //R0=sqrt(pow(1.23*std::cbrt(A)-0.60,2)+(7./3.)*pow(TMath::Pi(),2)*0.52*0.52-5*s*s);
 */
    cout << s << endl;
    cout << std::cbrt(A) << endl;
    cout << R << endl;
    cout << R0 << endl;
    //cout << R0.value();
    //cout << to_string(R0.units()) << endl;
    //cout << (sqrt(2*Mn*200)).value() << (sqrt(2*Mn*200)).units() << endl;
    cout << fHelm_fp(200,A) << endl;
    
    //// FF in energy units. Two different definitions. ////
    TCanvas *c=new TCanvas();
    c->SetLogx();
    c->SetLogy();
    auto g = new TGraph(n,x,y);
    g->SetLineColor(kMagenta+3);
    g->SetLineWidth(4);
    g->SetMinimum(pow(10,-10));
    g->SetMaximum(10.);
    g->Draw("");
    
    auto gp = new TGraph(n,x,yp);
    gp->SetLineColor(kGreen+3);
    gp->SetLineWidth(4);
    //gp->SetMinimum(pow(10,-12));
    //gp->SetMaximum(10.);
    gp->Draw("same");
    
    double qmin = 0;
    double qmax = 5;
    int nq = 100;
    double qstepsize = (qmax-qmin)/(double)nq;
    double xq[nq], yq[nq];
    
    for (int i=0;i<nq;i++) {
      xq[i] = qmin +i*qstepsize;
      yq[i] = fHelm_fq(xq[i], A);
    }
    
    //// FF in q units ////
    TCanvas *cq=new TCanvas();
    //cq->SetLogx();
    cq->SetLogy();
    auto gq = new TGraph(nq,xq,yq);
    gq->SetLineColor(kBlue+2);
    gq->SetLineWidth(4);
    gq->SetMinimum(pow(10,-12));
    gq->SetMaximum(10.);
    gq->Draw("");
    
    cout << "F^2(3000keV): " << fHelm_f(3000, A) << endl;
    cout << "sin(30): " << sin(30*TMath::Pi()/180) << endl;
    

}

///////////// Functions //////////////////
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

//Helm form factor squared translating python version. // Jelle Aalbers
Double_t fHelm_fp(Double_t E,Int_t A){
    // E = recoil energy in keV
    //precise_measurement E_keV = E *units::precise::kilo*eV;
    //precise_measurement hc = 197327 *units::precise::kilo*eV * units::precise::femto*units::precise::m; //keV fm (hc)
    Double_t hc = 197327; //keV fm (hc)
    //Double_t nucleon_mass = 0.938; // GeV/c2
    Double_t nucleon_mass_keV = 0.938*pow(10,6); // keV/c2
    //precise_measurement nucleon_mass = 0.938 *units::precise::giga*eV / (constants::c*constants::c);
    Double_t Mn = nucleon_mass_keV * (double)A; //keV/c2
    //precise_measurement Mn = nucleon_mass * A; //keV/c2
    Double_t s = 0.9; //Femtometers-Skin thickness of the nucleus
    //Double_t R = 1.23*std::cbrt(A)-0.6; // Problems with cubic root, use this std::cbrt() function instead
    //Double_t R0 = sqrt(pow(R,2)+7.*pow(TMath::Pi(),2)*0.52*0.52/3.-5*h*h*s*s);
    Double_t q = sqrt(2*Mn*E); 
    
    //return exp(-q*q*s*s/(hc*hc))*pow(3*(sin(q*R0/ (hc*hc))-q*R0*cos(q*R0))/pow(q*R0,3), 2);
    return exp(-q*q*s*s/(hc*hc))*9*fBessel1(q*fR0(A)/hc)*fBessel1(q*fR0(A)/hc)/pow(q*fR0(A)/hc, 2);
    //return exp(-(q*q*s*s/(hc*hc)).value())*9*fBessel1((q*fR0(A)/hc).value())*fBessel1((q*fR0(A)/hc).value())/pow((q*fR0(A)/hc).value(), 2);
    
    
}

/*
/// Version python

def F2(E,A): #Squared form factor. Radius formula Ciaran's thesis
    # Esta expresión del factor de forma está sacada de J.Gracia Garza, pero los parámetros R1, R, s provienen del trabajo de C. O`Hare.
    
    h=197.3 #MeV fm (hc)
    Mn = nucleon_mass * A

    s = 0.9/h #*10**(-15)
    R = ((1.23 * A**(1/3))-0.6)  #*10**(-15)
    R1 = math.sqrt(R**2 +(7*pi**2*0.52**2)/3- 5*h*h*s**2)/h ### Esta es la representación que coincide con la gráfica de Ciaran

    q=math.sqrt(2*Mn*E)

    return exp(-(q*s)**2)*(3*(sin(q*R1)-q*R1*cos(q*R1))/(q*R1)**3)**2  #(exp(-(q**2)*(s**2)))*(3*(sin(q*R1)/((q**2)*(R1**2))-(cos(q*R1)/(q*R1)))/(q*R1))**2
*/

//Helm form factor in q units
Double_t fHelm_fq(Double_t q,Int_t A){
    // T = recoil energy in keV
    //Modulus of the Helm Factor
    //Defined as in J.D.LewinP.F.Smith/AstroparticlePhysics6(1996)87-112
    //Femtometers-Effective radius of the target nucleus
    Double_t s=0.9; //Femtometers-Skin thickness of the nucleus
    //Double_t R0=sqrt(pow(1.23*pow(A,1/3)-0.60,2)+(7./3.)*pow(TMath::Pi(),2)*0.52*0.52-5*s*s);
    //Double_t q=6.927*pow(10,-3)*sqrt(A*T); //*pow(10,-3) Este factor numerico tiene que venir de dividir entre h barra las distancias en fm y de sqrt(2Mn) -> sqrt(2Mn)sqrt(A*Er)/h (hbar = 197)
    //hc = 197 Mev fm = 197*1000 keV fm
    // sqrt(2*0.931)/0.197 = 6.927 Juraría que para que todo esté en keV sobra ese 10^-3. Si lo quito queda demasiado pequeña la señal
    Double_t qR_0=q*fR0(A);
    
    return pow(3*(sin(qR_0)-qR_0*cos(qR_0))/pow(qR_0,3), 2)*exp(-q*q*s*s); //  /2.
}