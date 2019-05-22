#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

void readInput(int &N, double &L, double &xd, double &t, int &itMax, double &pL, double &uL, double &roL, double &pR, double &uR, double &roR, double &gama, double &pinf);

double soundSpeedEOS(double gama, double pinf, double p, double ro);

double roHugoniot(double gama, double p, double p0, double ro0);

double roIsentropic(double gama, double p, double p0, double ro0);

double roFanL(double gama, double roL, double uL, double cL, double S);
double uFanL(double gama, double uL, double cL, double S);
double pFanL(double gama, double pL, double uL, double cL, double S);

double roFanR(double gama, double roR, double uR, double cR, double S);
double uFanR(double gama, double uR, double cR, double S);
double pFanR(double gama, double pR, double uR, double cR, double S);


#endif