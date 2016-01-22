#include <cmath>
using namespace std;

#include "GridFunctionNd/UCLAQ_GridFunction1d.h"


#ifndef _ExactPeriodicSolution_
#define _ExactPeriodicSolution_
//
// Member functions of the ExactPeriodicSolution class returns
// the exact solution to u_t + a u_x = 0
// for x in the interval [xMin, xMax].
//
// The initial data implemented is a smoothed  "hump"
// centered in the interval 
//
// Initial u(x,0) = [16*(x-1/4)*(3/4-x) ]^3 for x in [1/4,3/4]
//                = 0                       otherwise
//
// 
class ExactPeriodicSolution
{
public :

    //
    // Constructor : Specify the interval [xMin, xMax] and the
    //               convection speed a. 
    //

    ExactPeriodicSolution(double xMin, double xMax, double a)
    {
        this->xMin        = xMin;
        this->xMax        = xMax;
        this->xInterval   = xMax-xMin;
        this->a           = a;
    }

    //
    // Given x in the interval [xMin,xMax] this function returns
    // the exact solution at time t.
    //
    double evaluateExactSolution(double x, double t)
    {
        double pStar;

        pStar = x - a*t;

        // determine relative location in the
        // first periodic interval 

        pStar = pStar - xMin;
        pStar = pStar - xInterval*floor(pStar/xInterval);
 
        // evaluate initial data at this location 

        return initialFunction(pStar);
    }

    //
    // This routine samples the exact solution at time t at
    // the nodes of an equispaced grid whose number of panels
    // is one less than the size of the input array u. 
    //
    void evaluateExactSolution(UCLAQ::GridFunction1d& u, double t)
    {
        long mPanel = u.getXpanelCount();
        double dx   = u.getHx();
        double xp;
    
        long i;

        for(i = 0; i < mPanel+1; i++)
        {
        xp = xMin + double(i)*dx;
        u(i) = evaluateExactSolution(xp,t);
        }
        return;
    }

    //
    // The initial function is a smoothed  hump centered in  
    // the middle of the domain. The specification is with respect
    // to the periodic interval [0, xInterval].
    //
    // Initial u(x,0) = [16*(x-1/4)*(3/4-x) ]^3 for x in [1/4,3/4]
    //                = 0                       otherwise
    //
    double initialFunction(double pStar)
    {
        double s = pStar/(xInterval);   // normalize to [0,1]
        if(s < 0.25) return 0.0;
        if(s > 0.75) return 0.0;
        return  pow(16.0*(s-.25)*(.75-s),3);
    }

    double xMin;        // periodic domain bounds
    double xMax;
    double xInterval; // periodic size
    double a;         // convective wave speed
};

#endif
