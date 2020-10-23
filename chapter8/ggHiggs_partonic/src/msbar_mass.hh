#pragma once

// routine to compute ms-bar bottom mass
extern double runmass(double mass0, double api0, double apif, unsigned int order, int nf);

extern double pole_to_MSbar(double mbPole, double alphas_at_mbPole, int order); // converts pole  to MSbar mass, using value of mbPole
extern double MSbar_to_pole(double mbmb,   double alphas_as_mbmb,   int order); // converts MSbar to pole  mass, using value of mb(mb)

extern double beta0(int nf);
extern double beta1(int nf);
extern double beta2(int nf);
extern double beta3(int nf);

