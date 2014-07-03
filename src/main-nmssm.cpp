#include <iostream>
#include "def.h"
#include "linalg.h"
#include "lowe.h"
#include "nmssmsoftsusy.h"

using namespace softsusy;



int main() {

  /// Sets format of output: 6 decimal places
  outputCharacteristics(6);
  softsusy::PRINTOUT = 0;

  /// DGS NMSSM parameters
  int sgnMu = 1;      ///< sign of mu parameter
  int numPoints = 30; ///< number of scan points
  double lambda = 0.05, kappa = 0.002, s = 0.0, xiF = 0.0, mupr = 0.0;
  double n5 = 2, mMess = 1.74e13, LAMBDA = 1.e5, cgrav = 1., xiPhi = 0.1, 
    xiT = 0.1;

  QedQcd oneset;      ///< See "lowe.h" for default definitions parameters

  /// most important Standard Model inputs: you may change these and recompile
  double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.18;
  oneset.setAlpha(ALPHAS, alphasMZ);
  oneset.setPoleMt(mtop);
  oneset.setMass(mBottom, mbmb);
  oneset.toMz();      ///< Runs SM fermion masses to MZ

  /// Set limits of tan beta scan
  double startTanb = 5.0, endTanb = 40.0, tanb = 0.;

  DoubleVector pars(6);
  pars(1) = n5; 
  pars(2) = LAMBDA;
  pars(3) = mMess; 
  pars(4) = cgrav; 
  pars(5) = xiPhi ; 
  pars(6) = xiT;
  DoubleVector nmpars(5);
  nmpars(1) = lambda; nmpars(2) = kappa; nmpars(3) = s;
  nmpars(4) = xiF; nmpars(5) = mupr;
  bool uni = false; // MGUT defined by g1(MGUT)=g2(MGUT)
  double mGutGuess=mMess;
  softsusy::Z3 = true;

  for (int i = 0; i < numPoints; i++) {
     tanb = (endTanb - startTanb) / double(numPoints) * double(i) +
        startTanb; // set tan beta ready for the scan.

     NmssmSoftsusy n;

     try {
       n.lowOrg(gmsbNMSSMBcs, mGutGuess, pars, nmpars, sgnMu, tanb,
                oneset, uni);
       //       cout << n; exit(0);
     } catch (const std::string& error) {
       n.flagProblemThrown(true);
     } catch (const char* error) {
       n.flagProblemThrown(true);
     }

     /// check the point in question is problem free: if so print the output
     if (!n.displayProblem().test()) {
        cout << tanb << ' '
             << n.displayPhys().mh0(1) << ' '
             << n.displayPhys().mh0(2) << ' '
             << n.displayPhys().mA0(1) << ' '
             << n.displayPhys().mHpm <<  ' '
	     << n.displayPhys().mu(1, 3) << ' '
	     << n.displayPhys().mu(2, 3) << '\n';
     } else {
        cout << tanb << ' ' << n.displayProblem() << '\n';
     }
  }

  return 0;
}
