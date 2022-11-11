'''Wimp rate calculation: N = t * flux * effective area
flux= v*ro/mw
eff area = Mt/mN * cross sectionN
differential cross section = cross section / Emax * F^2
Integrate over velocity distribution'''


import numpy as np
from numpy import sin, cos, pi, exp
from scipy.special import erf
from scipy.integrate import quad
from scipy.stats import poisson

import numericalunits as nu

HC = 197327. ## keV fm (hc)
Mn = 1 #0.9314941024171443 #0.938 ## GeV/c2 J.Aalbers takes Mn = 1 implicitly. Except in the Helm Form Factor where it takes Mn = 0.9314941024171443
# For more info look for Unified Atomic Mass Unit (Dalton unit): 1u =  1.66053906660(50)×10−27 kg =  931.49410242(28) MeV/c^2
MnFF = 0.9314941024171443
Mn_KEV = MnFF*10**6 ## keV/c2
S = 0.9 ## Femtometers-Skin thickness of the nucleus
C = 300000. # speed of light in km/s

RO = 0.3 ## DM density(GeV/c^2/cm^3)
VESC = 544. #* nu.km/nu.s## km/s
VEARTH = 232. #* nu.km/nu.s## km/s
SIG0 = 220. #* nu.km/nu.s## km/s

ATOMIC_WEIGHT = dict(
    Xe=131.293,
    Ar=39.948,
    Ge=72.64,
    Ne=20.1797,
    C=12.0107, 
    H=1.00797)

#A = ATOMIC_WEIGHT['Xe'] ## Nucleus mass


def bessel1(x):
    '''First Bessel function'''
    return sin(x)/(x*x)-cos(x)/x

def RN(A):
    '''Parametrization of atomic nucleous. A in atomic units'''
    return ((1.23*A**(1./3.)-0.60)**2 + (7./3.) * pi**2 * 0.52**2 - 5 * S**2)**(1./2.)

def FormFactor(E, A):
    '''Helm form factor parametrization. E in keV, A atomic weight in atomic units. 
    Factors HC to convert units. Form Factor is unitless.'''
    q = (2*Mn_KEV*A*E)**(1./2.) # In keV/c
    qs = q*S/HC
    qR = q*RN(A)/HC
    return exp(-0.5*qs**2)*3*bessel1(qR)/qR

def FormFactor2(E, A):
    '''Squared form factor'''
    return FormFactor(E,A)**2

def mr(m1, m2):
    '''Reduced mass'''
    return (m1*m2)/(m1+m2)

def Emax(v, mw, A):
    '''Maximum energy in [keV] deposited in an elastic recoil by a WIMP with velocity v. 
    v in km/s, mw in Gev/c^2, A in atomic units. Second version: A atomic weight, not conversion to energy units.
    Factor 10**6/C**2 to convert output units to keV. '''
    #return 10**6/C**2 * 2*mr(mw, Mn*A)**2 * v**2 / (Mn*A)
    return 10**6/C**2 * 2*mr(mw, A)**2 * v**2 / (A)

def sigmaNucleus(mw, sigma0, A):
    '''Cross section WIMP nucleus in [cm^2]. mw in Gev/c^2, sigma0 cross section per nucleon in cm^2, 
    A in atomic units. Second version: A atomic weight, not conversion to energy units.'''
    #return sigma0 * A**2 * mr(mw, Mn*A)**2 / mr(mw, Mn)**2
    return sigma0 * A**2 * mr(mw, A)**2 / mr(mw, 1)**2

def diffCrossSec(E, v, mw, sigma0, A):
    '''Differential cross section in energy, in [cm/keV]. E in keV, v in km/s, mw in Gev/c^2, 
    sigma0 cross section per nucleon in cm^2, A in atomic units (or atomic weight)'''
    return sigmaNucleus(mw, sigma0, A) * FormFactor2(E, A) / Emax(v, mw, A)

def flux(v, mw):
    '''Velocity in km/s times DM density in GeV/c^2/cm^3, by mw in Gev/c^2. Factor 10**5 to express velocity in cm/s. Flux in units of [counts/s/cm^2]'''
    return 10**5 * v * RO / mw

def time(T):
    '''Exposure time in years. Output: equivalent time in [seconds]'''
    return 365. * 24. * 3600. * T

def mass(MD, A):
    '''Detector mass in tonnes, A in atomic units (g/mol). Output: number of atoms in the detector. Factor 10**6 to convert tonnes in grams. N Avogadro = 6.02214*10**23'''
    return 10**6 * MD * 6.02214*10**23 / (Mn*A)

def vDist(v):
    '''Velocity distribution in Earth frame, f(v) in units of 1/(km/s), v in km/s'''
    # Normalization constant, see Lewin&Smith appendix 1a
    _w = VESC/SIG0
    k = erf(_w) - 2/np.pi**0.5 * _w * np.exp(-_w**2)  # unitless
    xmax = np.fmin( 1., (VESC**2 - VEARTH**2 - v**2)/(2*VEARTH*v) )
    return v * k / (VEARTH*SIG0*pi**0.5) * ( exp(-(v-VEARTH)**2/SIG0**2) - exp(-(v**2+VEARTH**2+2*v*VEARTH*xmax)/SIG0**2) )  ##  f(v) units / (velocity)
    # J.Aalbers: v * k ... In my opinion, as a normalization constant, should be v / k. As long as k is close to 1, this is not very relevant.

def vmin_elastic(E, mw, A):
    """Minimum WIMP velocity that can produce a recoil of energy E(keV), mw in GeV/c^2, 
    A in atomic units. Factor 10**(-6) to convert energy from keV to GeV. 
    Factor C from masess. """
    ##print(np.size(np.sqrt(Mn*A * E / (2 * mr(mw, Mn*A)**2))))
    return np.sqrt(A * E * 10**(-6) / (2 * mr(mw, A)**2)) * C # Mn = 1 as used in J.Aalbers
    #return np.sqrt(Mn*A * E * 10**(-6) / (2 * mr(mw, Mn*A)**2)) * C ##Mn*

def rate(E, mw, sigma0, A, T, MD):
    '''Integration of all components of the WIMP signal. E in keV, mw in GeV/c^2, sigma0 in cm^2, 
    A in atomic units, T in years, MD in tonnes'''
    integrand = lambda v : time(T) * flux(v, mw) * diffCrossSec(E, v, mw, sigma0, A) * mass(MD, A) * vDist(v)
    v_min = vmin_elastic(E, mw, A)
    v_max = VESC + VEARTH
    
    if v_min >= v_max:
        return 0
    
    result = quad(lambda v: integrand(v),v_min, v_max)[0]
    return result

def resolution(E, spectrum, Ereff, FWHMreff):
    ''' Convolution of energy spectrum with gaussian of variable width depending on the energy.
    E array of energyes in keV, spectrum array with the distribution function, 
    Ereff reference energy with known resolution, FWHMreff width of the gaussian resolution for Ereff (in %, eg. 21 for 21% FWHM).'''

    fwSig = (2*np.sqrt(2*np.log(2))) # FWHM = 2*sqrt(2*ln(2)) * sigma
    eRes =  np.zeros(np.size(E))

    for i in range(np.size(E)):
        fwhm = FWHMreff/100  * np.sqrt(Ereff/E[i])
        sigma = fwhm / fwSig
        gaussNorm = sigma * np.sqrt(2*np.pi)
        for j in range(np.size(E)):
            eRes[j] += spectrum[i] * exp(-((E[j]-E[i])**2)/(2*sigma**2)) / gaussNorm

    return eRes
          








#### POISSON METHOD ####

def survivalPoiss(alpha, b):
    '''Right hand side of the Poisson distribution for signal+background. Alpha in (0, 1), b number of background signals expected.'''
    if alpha > 1 or alpha < 0: print("Value for alpha outside bounds (0, 1)")
    a = 1
    nsig = 0
    while a > alpha :
        a = poisson.sf(nsig+b, b)
        nsig+=1
    return nsig