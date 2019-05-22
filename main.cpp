#include "functions.h" // Fonctions
#include <iostream>    // Bib communication avec console
#include <string>      // Bib chaine de char
#include <cmath>       // Bib fn math
#include <math.h>	   // Bib fn math
#include <fstream>     // Bib fichier

using namespace std;

int main()
{
    int N(1000), itMax(50), iteration(0);
    double L(1.), xd(0.5), t(1.e-4); 
    double pL(1.e5), uL(0.), roL(1.2);
    double pR(1.e5), uR(0.), roR(1.2);
    double gama(1.4), pinf(0.);

    double dx(0.);

    double eps(1.e-6), error(1.);
    double pStar(1.), uStar(0.), pOld(0.);
    double G1(0.), G2(0.), G3(0.), G4(0.);
    double Al(0.), Bl(0.), Ar(0.), Br(0.);
    double cL(0.), cR(0.);
    double phiL(0.), phiR(0.);
    double dphiL(0.), dphiR(0.);
    double fp(0.), dfp(0.);

    double S(0.), sL(0.), sR(0.), cStarL(0.), cStarR(0.), sTL(0.), sHL(0.), sTR(0.), sHR(0.);

    // Reading data
    readInput(N,L,xd,t,itMax,pL,uL,roL,pR,uR,roR,gama,pinf);
    double x[N];
    double ro[N], u[N], p[N], c[N], M[N];


    // Mesh
    dx = L/N;
    for(int i = 0; i < N; i++)
    {
        x[i] = i*dx;
    }

    // Computing sound speed in initial states
    cL = soundSpeedEOS(gama, pinf, pL, roL);
    cR = soundSpeedEOS(gama, pinf, pR, roR);

    // Calculation of pStar with searching of zero by Newton-Raphson 
    pStar = (pL+pR)*0.5;
    
    G1 = gama + 1.;
    G2 = gama - 1.;
    G3 = G1/(2*gama);
    G4 = G2/(2*gama);
    
    Al = 2./(G1*roL);
    Bl = pL*G2/G1;
    Ar = 2./(G1*roR);
    Br = pR*G2/G1;

    while(error>eps)
    {	
        if (pStar>pL) // Shock 
        {
            phiL = (pStar-pL)*sqrt(Al/(pStar+Bl)); 
            dphiL = sqrt(Al/(Bl+pStar))*(1.-(pStar-pL)/(2.*(Bl+pStar))); 
        }
        else          // Rarefaction
        {
            phiL = ((2.*cL)/G2)*( pow(pStar/pL,G4)-1. );
            dphiL = (1./(roL*cL))*pow(pStar/pL,-G3);
        }

        if (pStar>pR) // Shock 
        {
            phiR = (pStar-pR)*sqrt(Ar/(pStar+Br));  
            dphiR = sqrt(Ar/(Br+pStar))*(1.-(pStar-pR)/(2.*(Br+pStar)));
        }
        else          // Rarefaction
        {
            phiR = ((2.*cR)/G2)*( pow(pStar/pR,G4)-1. );
            dphiR = (1./(roR*cR))*pow(pStar/pR,-G3);
        }

        iteration += 1;
        pOld = pStar;
        fp = phiL + phiR + uR - uL;
        dfp = dphiL + dphiR;
        pStar -= fp/dfp;
        error = 2.*fabs((pStar-pOld)/(pStar+pOld));
        if (iteration >= itMax) { cout << "Error : non convergence of Newton-Raphson for finding pStar \n"; exit(0); }
    }
    cout << "pStar : " << pStar << endl;
    cout << "Convergence after " << iteration << " iteration\n";
    

    // Calculation of uStar
    uStar = 0.5*(uL+uR) + 0.5*(phiR-phiL);
    cout << "uStar : " << uStar << endl;


    // Computing wave's speed
    if (pStar>pL) { // Left shock
        sL = uL - cL*sqrt(G3*(pStar/pL) + G4); 
    }
    else { // Left rarefaction
        cStarL = cL*pow(pStar/pL, G4); // Sound speed behind the rarefaction wave
        sTL = uStar - cStarL;
        sHL = uL - cL;
    }

    if (pStar>pR) { // Right shock
        sR = uR + cR*sqrt(G3*(pStar/pR) + G4); 
    }
    else { // Right rarefaction
        cStarR = cR*pow(pStar/pR, G4); // Sound speed behind the rarefaction wave
        sTR = uStar + cStarR;
        sHR = uR + cR;
    }


    // Sampling of the solution
    for (int i = 0; i < N; ++i)
    {
    	S = (x[i]-xd)/t;
    	if (S<=uStar)
    	{
    		if (pStar>pL) // Left Shock
    		{
    			if (S<=sL)    // Initial state
    			{
	    			ro[i] = roL;
	    			u[i] = uL;
	    			p[i] = pL;
    			}
    			else         // Shocked state
    			{
	    			ro[i] = roHugoniot(gama, pStar, pL, roL);
	    			u[i] = uStar;
	    			p[i] = pStar;
    			}
    		}

    		else // Left Rarefaction
    		{
    			if (S>sTL)  // Relaxed state after the fan
    			{
	    			ro[i] = roIsentropic(gama, pStar, pL, roL);
	    			u[i] = uStar;
	    			p[i] = pStar;
    			}
    			else if (S>sHL)    // Relaxed state inside the fan 
    			{
    				ro[i] = roFanL(gama, roL, uL, cL, S);
    				u[i] = uFanL(gama, uL, cL, S);
    				p[i] = pFanL(gama, pL, uL, cL, S);
    			}
    			else                   // Initial state
    			{
    				ro[i] = roL;
    				u[i] = uL;
    				p[i] = pL;
    			}
    		}
    	}

    	else // S>uStar	
    	{
			if (pStar>pR) // Right Shock
    		{
    			if (S>=sR) // Initial state
    			{
	    			ro[i] = roR;
	    			p[i] = pR;
	    			u[i] = uR;
    			}
    			else // Shocked state
    			{
	    			ro[i] = roHugoniot(gama, pStar, pR, roR);
	    			p[i] = pStar;
	    			u[i] = uStar;
    			}
    		}
    		else // Right Rarefaction
    		{
    			if (S<sTR)  // Relaxed state before the fan
    			{
	    			ro[i] = roIsentropic(gama, pStar, pR, roR);
	    			u[i] = uStar;
	    			p[i] = pStar;
    			}
    			else if (S<sHR)    // Relaxed state inside the fan 
    			{
    				ro[i] = roFanR(gama, roR, uR, cR, S);
    				u[i] = uFanR(gama, uR, cR, S);
    				p[i] = pFanR(gama, pR, uR, cR, S);
    			}
    			else                   // Initial state
    			{
    				ro[i] = roR;
    				u[i] = uR;
    				p[i] = pR;
    			}
    	   }
        }

        // Computing the sound speed and the mach number
        c[i] = soundSpeedEOS(gama, pinf, p[i], ro[i]);
        M[i] = u[i]/c[i];
    }

    // Writing results
    ofstream streamOut("output.txt");   
    if (streamOut) {
        for (int i = 0; i < N; ++i) {
            streamOut << x[i] << " " << ro[i] << " " << u[i] << " " << p[i] << " " << M[i] << endl;
        }
    }
    else {cout << "Error : impossible to open the output file \n"; exit(0);}
    streamOut.close();
    
    return 0;
}






