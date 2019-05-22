#include "functions.h"
#include <math.h>
#include <iostream>
#include <fstream> 
#include <ostream> 

using namespace std;


void readInput(int &N, double &L, double &xd, double &t, int &itMax, double &pL, double &uL, double &roL, double &pR, double &uR, double &roR, double &gama, double &pinf)
{
	// Purpose : to read input user data and initialize corresponding variables
	string const fileNameIn("input.txt");
	ifstream streamIn(fileNameIn.c_str());
	if (streamIn)
	{
		string line("");
		getline(streamIn,line); 
		getline(streamIn,line);
		getline(streamIn,line);
		getline(streamIn,line);
		getline(streamIn,line);   
		N = stoi(line);           // Nb mailles internes

		getline(streamIn,line);
		getline(streamIn,line);   
		L = stod(line);           // Longueur tube

		getline(streamIn,line);
		getline(streamIn,line);   
		xd = stod(line);          // Position interface

		getline(streamIn,line);
		getline(streamIn,line);

		getline(streamIn,line);
		getline(streamIn,line);   
		t = stod(line);           // Tps final

		getline(streamIn,line);
		getline(streamIn,line);   
		itMax = stoi(line);       // Nbre max. d'iterations

		getline(streamIn,line);
		getline(streamIn,line);
		getline(streamIn,line);

		getline(streamIn,line);
		getline(streamIn,line);   
		pL = stod(line);          // Pression chambre Gauche

		getline(streamIn,line);
		getline(streamIn,line);   
		uL = stod(line);          // Vit. chambre Gauche

		getline(streamIn,line);
		getline(streamIn,line);   
		roL = stod(line);         // Masse vol. chambre Gauche

		getline(streamIn,line);
		getline(streamIn,line);

		getline(streamIn,line);
		getline(streamIn,line);   
		pR = stod(line);          // Pression chambre Droite

		getline(streamIn,line);
		getline(streamIn,line);   
		uR = stod(line);          // Vit. chambre Droite

		getline(streamIn,line);
		getline(streamIn,line);   
		roR = stod(line);         // Masse vol. chambre Droite

		getline(streamIn,line);
		getline(streamIn,line);

		getline(streamIn,line);
		getline(streamIn,line);   
		gama = stod(line);        // Coeff. de Laplace

		getline(streamIn,line);
		getline(streamIn,line);   
		pinf = stod(line);        // Parametre SG



		// Rappel du fichier input
		cout << "Input file" << endl;
		cout << "---------------------------------------------" << endl;
		cout << "Space mesh + Geometry \n";
		cout << "Number of mesh                  (N)        [-] : " << N << endl;
		cout << "Shock tube length               (L)        [m] : " << L << endl;
		cout << "Interface position              (xd)       [m] : " << xd << endl;

		cout << "---------------------------------------------" << endl;
		cout << "Time + condition \n";
		cout << "Duration of simulation             (tend)  [s] : " << t << endl;
		cout << "Max of iterations Riemann solver   (-)     [-] : " << itMax << endl;

		cout << "---------------------------------------------" << endl;
		cout << "State conditions \n";
		cout << " - Left Side \n";
		cout << "Pressure                   (pL)    [Pa]    :   " << pL << endl;
		cout << "Speed                      (uL)    [m/s]   :   " << uL << endl;
		cout << "Density                    (roL)   [kg/m3] :   " << roL << endl;
		cout << " - Right Side \n";
		cout << "Pressure                   (pR)    [Pa]    :   " << pR << endl;
		cout << "Speed                      (uR)    [m/s]   :   " << uR << endl;
		cout << "Density                    (roR)   [kg/m3] :   " << roR << endl;

		cout << "---------------------------------------------" << endl;
		cout << "Fluid EOS (SG) \n";
		cout << "Gamma        (gama)   [-]  : " << gama << endl;
		cout << "Pinfini      (pinf)   [Pa] : " << pinf << endl;
		
		cout << "---------------------------------------------" << endl;
		cout << "Start of the simulation" << endl;

        streamIn.close();
	}
	else
	{
		cout << "Error : reading input.txt file" << endl;
		exit(0);
	}
}


double soundSpeedEOS(double gama, double pinf, double p, double ro)
{
	// Purpose : calculate the sound speed with Equation Of State (here for Stiffened Gas)
	double c(0.);
	c = sqrt((gama*(p+pinf))/ro);
	return c;
}



// ----- Shock wave function -----
double roHugoniot(double gama, double p, double p0, double ro0)
{
	// Purpose : calculate the density for a shocked state
	double ro(0.), num(0.), den(1.);
	num = p*(gama+1.) + (gama-1.)*p0;
	den = p*(gama-1.) + (gama+1.)*p0;
	ro = ro0*(num/den);
	return ro;
}



// ----- Rarefaction wave functions -----

double roIsentropic(double gama, double p, double p0, double ro0)
{
	// Purpose : calculate the density for an isentropic flow in case of a rarefaction wave
	double ro(0.);
	ro = ro0 * pow(p/p0, 1./gama);
	return ro;
}

double roFanL(double gama, double roL, double uL, double cL, double S)
{
	// Purpose : calculate the density inside a left fan of a rarefaction wave
	double roFanL(0.), fn(0.), G1(0.), G2(0.), G3(0.);
    
    G1 = 2./(gama+1.);
    G2 = (gama - 1.)/(gama + 1.);
    G3 = 2./(gama-1.);
	
	fn = G1 + (G2/cL)*(uL-S);
	roFanL = roL*pow(fn,G3);
	return roFanL;
}

double uFanL(double gama, double uL, double cL, double S)
{
	// Purpose : calculate the speed inside a left fan of a rarefaction wave
	double uFanL(0.), G1(0.), G2(0.);

	G1 = 2./(gama+1.);
	G2 = (gama-1.)/2.;

	uFanL = G1*(cL+G2*uL+S);
	return uFanL;
}

double pFanL(double gama, double pL, double uL, double cL, double S)
{
	// Purpose : calculate the pressure inside a left fan of a rarefaction wave
	double pFanL(0.), fn(0.), G1(0.), G2(0.), G3(0.);

    G1 = 2./(gama+1.);
    G2 = (gama - 1.)/(gama + 1.);
    G3 = (2.*gama)/(gama-1.);

    fn = G1 + (G2/cL)*(uL-S);
    pFanL = pL*pow(fn,G3);
    return pFanL;
}


double roFanR(double gama, double roR, double uR, double cR, double S)
{
	// Purpose : calculate the density inside a right fan of a rarefaction wave
	double roFanR(0.), fn(0.), G1(0.), G2(0.), G3(0.);
    
    G1 = 2./(gama+1.);
    G2 = (gama - 1.)/(gama + 1.);
    G3 = 2./(gama-1.);
	
	fn = G1 - (G2/cR)*(uR-S);
	roFanR = roR*pow(fn,G3);
	return roFanR;
}

double uFanR(double gama, double uR, double cR, double S)
{
	// Purpose : calculate the speed inside a right fan of a rarefaction wave
	double uFanR(0.), G1(0.), G2(0.);

	G1 = 2./(gama+1.);
	G2 = (gama-1.)/2.;

	uFanR = G1*(-cR+G2*uR+S);
	return uFanR;
}

double pFanR(double gama, double pR, double uR, double cR, double S)
{
	// Purpose : calculate the pressure inside a right fan of a rarefaction wave
	double pFanR(0.), fn(0.), G1(0.), G2(0.), G3(0.);

    G1 = 2./(gama+1.);
    G2 = (gama - 1.)/(gama + 1.);
    G3 = (2.*gama)/(gama-1.);

    fn = G1 - (G2/cR)*(uR-S);
    pFanR = pR*pow(fn,G3);
    return pFanR;
}

