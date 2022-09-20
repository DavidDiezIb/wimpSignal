#include "TF1.h"
#include "TF2.h"
#include "Math/WrappedTF1.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/GaussLegendreIntegrator.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "TMath.h"

#include <cstdio>
#include <cstdlib>
#include <algorithm> 
#include <random>
#include <chrono>
#include <iostream>

Double_t fMaxwell(Double_t x, Double_t a); // Maxwell Boltzmann distribution
Double_t fFv(Double_t xx); // Distribución de velocidades
Double_t fFvInt(); // Integral


void VelocityHaloDistribution() {
    //Atomic and mass numbers for the element under study
    //Xenon Z=54,N=77;
    //Argon Z=20,N=20;
    //Neon Z=10,N=10;
    //Hidrogen Z=1,N=1
    //Carbon Z=6,N=6
    //Germanium Z=32,N=32
    
    Int_t Z=10,N=10; 
    Int_t A=N+Z;
    
    Double_t v0=220.;//v_0(km/s)
    Double_t vesc=544.;//vesc(km/s)
    Double_t vlab=232.;//Vlab(km/s) (v_earth) in Jelle
    
    double xmin = 0.1;
    double xmax = 800;
    int n = 100000;
    double step = (xmax-xmin)/(double)n;
    double x[n], y[n], z[n];
    
    for (int i=0;i<n;i++) {
      x[i] = xmin +i*step;
      y[i] = fFv(x[i]);
      z[i] = fMaxwell(x[i], v0/sqrt(2)); 
    }
    
    Double_t Nesc=ROOT::Math::erf(vesc/v0)-(2/sqrt(TMath::Pi()))*(vesc/v0)*exp(-vesc*vesc/(v0*v0)); // k = erf(v_esc/v_0) - 2/np.pi**0.5 * v_esc/v_0 * np.exp(-(v_esc/v_0)**2) (Jelle) // Usa v0 y vesc
    
    cout << "Normalization constant: " << Nesc << endl;
    cout << "Integral de 0 a vesc + vlab: " << fFvInt() << endl;
    
    TCanvas *c=new TCanvas();
    //c->SetLogx();
    //c->SetLogy();
    TGraph *g = new TGraph(n,x,y);
    g->SetLineColor(kGreen+3);
    g->SetLineWidth(4);
    g->SetMinimum(0);
    g->SetMaximum(0.005);
    g->Draw("");
    
    TGraph *gp = new TGraph(n,x,z);
    gp->SetLineColor(kBlue+3);
    gp->SetLineWidth(4);
    gp->Draw("same");
    
    auto legend = new TLegend(0.55,0.7,0.9,0.9);
    legend->AddEntry(g,"Local frame","l");
    legend->AddEntry(gp,"Galactic frame","l");
    legend->Draw();
    
    TLine *lg = new TLine(vesc,0,vesc,0.0037);
    lg->SetLineStyle(2);
    lg->SetLineWidth(4);
    lg->SetLineColor(kBlue+3);
    lg->Draw();
    
    TLine *ll = new TLine(vesc+vlab,0,vesc+vlab,0.0037);
    ll->SetLineStyle(2);
    ll->SetLineWidth(4);
    ll->SetLineColor(kGreen+3);
    ll->Draw();

}

///////////// Functions //////////////////   
//Maxwell Boltzmann distribution
Double_t fMaxwell(Double_t x, Double_t a){
    return sqrt(2./TMath::Pi())*x*x*exp(-x*x/(2*a*a))/pow(a,3);
}


//Velocity distribution function f(v) in the local frame (geen in https://github.com/JelleAalbers/wimprates/blob/master/notebooks/Checks%2C%20plots.ipynb)
Double_t fFv(Double_t xx){
    //Double_t sigmav=220.*sqrt(1./2.);//v_0/sqrt(2)(km/s) Dispersión en v: sqrt(3/2)*v0 para Javier Garcia Garza, sqrt(1/2)*v0 para Victor
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
}

    
//Integral of the velocity distribution. Vesc in galactic frame transforms to Vesc+Vlab in earth frame.
Double_t fFvInt(){
    Double_t vlab=232.;//Vlab(km/s) (v_earth) in Jelle
    Double_t vesc=544.;//vesc(km/s)
    
    double xmin = 0.001;
    double xmax = vesc+vlab;
    int n = 100000;
    double step = (xmax-xmin)/(double)n;
    double add;
    
    for (int i=0;i<n;i++) {
      add += step*fFv(xmin + i*step);
    }
    return add;
    
}