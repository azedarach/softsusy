#ifndef TWOLOOPGLUINO_H
#define TWOLOOPGLUINO_H

#include <iostream>
using namespace std;

/// Returns the two-loop part of delta in physical mgluino=m3(1 + delta)
/// All running couplings in DRbar scheme and Feynman gauge
/// g3 is running QCD gauge coupling
/// mdl, mdr etc are squark masses. "mtl" is a mostly left-handed eigenstate,
/// but the "l" is just a label really. 
/// thetab and thetat are the bottom and top mixing angles
/// ht and hb are the top and bottom Yukawa couplings
/// q is the renormalisation scale at which the corrections are called
double twoLoopGluino(double g3, double mdl, double mdr, double msl, 
		     double msr, double mbl, double mbr, double mul, 
		     double mur, double mcl, double mcr, double thetat, 
		     double thetab, double m3, double ht, double hb, double q);

#endif
