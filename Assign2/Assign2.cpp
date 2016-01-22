//
//#####################################################################
//                            Assign2.cpp 
//#####################################################################
//
// This program computes the exact solution to 
// 
// u_t + a u_x = 0 
//
// for xMin <= x <= xMax, 0 <= t <= tFinal with
// u(0,x) = u0(x) and u(t,x) being periodic in x 
// with period xMax-xMin (e.g. "periodic" boundary
// conditions). 
//
// The spatial mesh is uniform with M panels between
// x = xMin and x = xMax. The array indices run from 0 to M
// as indicated below: 
// 
//
//               M panels 
// xMin                             xMax
//   
// |---x---x---x---x---x---x---x---x---x
// 0   1   2   3   4   5               M  
//
// Output is printed to files of the form
//
//         uExactXXX.dat 
//
// where XXX is the output file index.
//
// Currently the format of the data is an ASCII file
// and can be read by either GNUplot or Matlab.
//

//#####################################################################
// Created for Math 269B. 
// Version: Dec 23, 2015
// Author : Chris Anderson 
//#####################################################################

#include <cmath>      // For math functions
#include <cstdio>     // For "C" style input/output; printf, fprintf, etc,
#include <iostream>   // For "C++" stream input/output
using namespace std;


#include "GridFunctionNd/UCLAQ_GridFunction1d.h"          // Provides a "uniform grid" function
#include "GridFunctionNd/UCLAQ_GridFunction1dUtility.h"   // Provides output routines

#include "ExactPeriodicSolution.h"                        // Class for the exact solution

// Function prototype for utility to construct filenames for output
// The code for the routine follows the main(..) routine

string composeFileName(string fileNamePrefix, long outputIndex);

UCLAQ::GridFunction1d forwardDiffOpP(const UCLAQ::GridFunction1d& V, double dt, double a)
{
	// Extract grid structure from input argument

	double xMin = V.getXmin();
	double xMax = V.getXmax();
	double hx = V.getHx();
	long  xPanels = V.getXpanelCount();


	UCLAQ::GridFunction1d  R(xPanels, xMin, xMax);   // For the return argument

	R.initialize(V);             // Make a copy of the input vector for the
								 // return argument. This will automatically
								 // initialize the values of R with those of V.

								 // Right end. Using periodicity to determine V(-1)

	//R(0) = V(0) - a*dt*(V(0) - V(xPanels - 1)) / (hx);

	// Interior points

	for (long i = 0; i <= xPanels - 1; i++)
	{
		R(i) = V(i) - a*dt*(V(i+1) - V(i)) / (hx);
	}

	R(xPanels) = R(0);
	return R;
}

UCLAQ::GridFunction1d backDiffOpP(const UCLAQ::GridFunction1d& V, double dt, double a)
{
	// Extract grid structure from input argument

	double xMin = V.getXmin();
	double xMax = V.getXmax();
	double hx = V.getHx();
	long  xPanels = V.getXpanelCount();


	UCLAQ::GridFunction1d  R(xPanels, xMin, xMax);   // For the return argument

	R.initialize(V);             // Make a copy of the input vector for the
								 // return argument. This will automatically
								 // initialize the values of R with those of V.

								 // Right end. Using periodicity to determine V(-1)

	R(0) = V(0) - a*dt*(V(0) - V(xPanels - 1)) / (hx);

	// Interior points

	for (long i = 1; i <= xPanels - 1; i++)
	{
		R(i) = V(i) - a*dt*(V(i) - V(i - 1)) / (hx);
	}

	R(xPanels) = R(0);
	return R;
}

UCLAQ::GridFunction1d bDiffOpP(const UCLAQ::GridFunction1d& V, double dt, double a)
{
	// Extract grid structure from input argument

	double xMin = V.getXmin();
	double xMax = V.getXmax();
	double hx = V.getHx();
	long  xPanels = V.getXpanelCount();


	UCLAQ::GridFunction1d  R(xPanels, xMin, xMax);   // For the return argument

	R.initialize(V);             // Make a copy of the input vector for the
								 // return argument. This will automatically
								 // initialize the values of R with those of V.

								 // Right end. Using periodicity to determine V(-1)

	R(0) = V(0) - a*dt*(3.0*V(0)/2.0 - 2.0*V(xPanels - 1)+0.5*V(xPanels - 2)) / (hx);
	R(1) = V(1) - a*dt*(3.0 * V(1) / 2.0 - 2.0 * V(0) + 0.5*V(xPanels - 1)) / (hx);


	// Interior points

	for (long i = 2; i <= xPanels - 1; i++)
	{
		R(i) = V(i) - a*dt*(3.0 * V(i) / 2.0 - 2.0 * V(i-1) + 0.5*V(i - 2)) / (hx);
	}

	R(xPanels) = R(0);
	return R;
}

UCLAQ::GridFunction1d centralDiffOpP(const UCLAQ::GridFunction1d& V, double dt, double a)
{
	// Extract grid structure from input argument

	double xMin = V.getXmin();
	double xMax = V.getXmax();
	double hx = V.getHx();
	long  xPanels = V.getXpanelCount();


	UCLAQ::GridFunction1d  R(xPanels, xMin, xMax);   // For the return argument

	R.initialize(V);             // Make a copy of the input vector for the
								 // return argument. This will automatically
								 // initialize the values of R with those of V.

								 // Right end. Using periodicity to determine V(-1)

	R(0) = V(0) - a*dt*(V(1) - V(xPanels - 1)) / (2.0*hx);

	// Interior points

	for (long i = 1; i <= xPanels - 1; i++)
	{
		R(i) = V(i) - a*dt*(V(i+1) - V(i - 1)) / (2.0*hx);
	}

	R(xPanels) = R(0);
	return R;
}


int main(int argc, char* argv[]) 
{ 

    double tFinal;             // final time

    double xMin    = 0.0;      // left endpoint
    double xMax    = 1.0;      // right endpoint
    double a       = -1.0;      // convection speed

    long   M;                   // number of panels
    double dx;                  // mesh spacing 

    double    time;             // total simulation time
    double      dt;             // timestep
    double  dtTemp;             // timestep for sub-stepping 

    long   outputCount;         // number of timsteps output
    double outputTimeIncrement; // time interval between output times
    double outputTime;          // next output time
    long   outputIndex;         // index of output steps

	double error;

    string       exactFileNamePrefix = "uExact";   // output filename prefix.

	string      numFileNamePrefix = "v";   // output filename prefix.

    ostringstream  fileStringStream;               // String stream creating output file titles


    cout << " Enter number of grid panels M : " << endl;
    cin  >> M;
 
    cout << " Enter timestep size dt : " << endl;
    cin >> dt;

    cout <<  "Enter Final Time : " << endl;
    cin  >>  tFinal;

    cout <<  "Number of Output Times : " << endl;
    cin  >>  outputCount;


    UCLAQ::GridFunction1d uExact(M,xMin,xMax);     // GridFunction1d to hold exact solution values

	UCLAQ::GridFunction1d v(M, xMin, xMax);     // GridFunction1d to hold numerical solution values

	UCLAQ::GridFunction1d vpls(M, xMin, xMax);

    ExactPeriodicSolution exactSoln(xMin,xMax,a);  // Exact solution class instance

    // Initialize spatial variables 

    dx = (xMax-xMin)/double(M);

    // Initialize time-stepping variables

    time                  = 0.0;
    outputTimeIncrement   = tFinal/double(outputCount);
    outputTime            = time + outputTimeIncrement;
    outputIndex           = 0;
    string  outputFileName;
  
    // Initialize solution using the exact solution evaluated at
    // t = tInitial 

    uExact.initialize(M,xMin,xMax);

    exactSoln.evaluateExactSolution(uExact,0.0);

	v.initialize(M, xMin, xMax);

	exactSoln.evaluateExactSolution(v, 0.0);

    // Instantiate utility class to use to create output

    UCLAQ::GridFunction1dUtility gUtility1d;

    // Output the initial solution profile

    outputFileName = composeFileName(exactFileNamePrefix, outputIndex);
    gUtility1d.outputToGNUplot(uExact,outputFileName);

	outputFileName = composeFileName(numFileNamePrefix, outputIndex);
	gUtility1d.outputToGNUplot(v, outputFileName);


    printf("Output %3ld Time = %10.5f  Timesteps taken : %d \n", outputIndex, time,0);
    
    // Main time-stepping loop 

    int outputFlag       = 0;
    long kStep           = 0;
    int exitFlag         = 0;                    
    while((time < tFinal)&&(exitFlag == 0))     
    {

    // Simulation evolution 

    if(time + dt < outputTime)              // next step requires no output   
    {
        //
        // Advance numerical solution : Nothing done here until Assignment 2
        //
		//vpls = centralDiffOpP(v, dt, a);
		//vpls = bDiffOpP(v, dt, a);
		//vpls = backDiffOpP(v, dt, a);
		vpls = forwardDiffOpP(v, dt, a);

        // update time, stepcount, and exact solution 
        time += dt;
        kStep++;
        exactSoln.evaluateExactSolution(uExact,time); 
		v = vpls;
    }
    else                                    // next step requires output
    {
        outputFlag = 1;                      
        dtTemp     = outputTime - time;     // substep so we hit output time exactly

        if(dtTemp > 1.0e-10)                // advance the solution to output time
        {                                   // (if necessary) 
        //
        // Advance numerical solution : Nothing done here until Assignment 2
        //

			//vpls = centralDiffOpP(v, dt, a);
			//vpls = bDiffOpP(v, dt, a);
			//vpls = backDiffOpP(v, dt, a);
			vpls = forwardDiffOpP(v, dt, a);

        // update time and exact solution 

        time += dtTemp;
        kStep++;
        exactSoln.evaluateExactSolution(uExact,time);
		v = vpls;

        }

        if(abs(time-tFinal) < 1.0e-10){exitFlag = 1;}
    }

    // Simulation output 

    if(outputFlag  == 1)
    {
        outputFlag  = 0;                    // reset flags and time markers
        outputTime += outputTimeIncrement;
        outputIndex++;      
        

        printf("Output %3ld Time = %10.5f  Timesteps taken : %ld \n", outputIndex, time, kStep);


        outputFileName = composeFileName(exactFileNamePrefix, outputIndex);
        gUtility1d.outputToGNUplot(uExact,outputFileName);

		outputFileName = composeFileName(numFileNamePrefix, outputIndex);
		gUtility1d.outputToGNUplot(v, outputFileName);


 
    }

    }

	error = 0;
	for (long i = 0; i <= M-1 ; i++)
	{
		error = error+(v(i)-uExact(i))*(v(i) - uExact(i))*dt;
	}

	error = sqrt(error);

	printf("error = %10.5f  \n ", error);

    return 0;
 }


string composeFileName(string fileNamePrefix, long outputIndex)
{
	ostringstream  fileStringStream;  // String stream creating output file titles
	fileStringStream.str("");
	fileStringStream << fileNamePrefix;
	if      (outputIndex <= 9)   {fileStringStream << "00" << outputIndex << ".dat";}
	else if(outputIndex <= 99)   {fileStringStream << "0" << outputIndex  << ".dat";}
	else                         {fileStringStream <<        outputIndex  << ".dat";}

    return fileStringStream.str();
}

