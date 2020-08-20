#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "../nr.h"
#include "../def.h"

using namespace std;

//define multiple shooting intervals for experiment expNr
void setMesh(GlobExp *ex,Glob *globs,long expNr)
{
  long i,nPoints=ex[expNr].nPoints; //number of mesh or multiple shooting points; initialized through
                                    //nms parameter   
  ex[expNr].mesh=dvector(1,nPoints);
  ex[expNr].mesh[1]= ex[expNr].fitstart; //first mesh point at start time of data
  ex[expNr].mesh[nPoints]= ex[expNr].fitend; //last mesh point at end time of data

  if (globs->noMeasurements==FALSE) 
    {         
      // there are measurements; mesh at data points
      if (ex[expNr].nMeasure < nPoints) //the # of time points (nMeasure) is smaller than the # of mesh points
	{
	  cerr << "Too few measurements to construct mesh\n";
	  exit(1);
	}
      
      //set intermediate multiple shooting points, evenly spaced in indices
      //Note: firstMeasure set to 1 in readData.cc
      
      if(ex[expNr].tms!=NULL)
	{
          //for (i=1; i <= nPoints; ++i)
          for (i=2; i < nPoints; ++i) // Use fitstart and fitend assigned above, don't use -tms options provided 
	    ex[expNr].mesh[i]=ex[expNr].tms[i];
      	}
      else
	{
	  for (i=2; i < nPoints; ++i)
	    ex[expNr].mesh[i]=ex[expNr].xMeasure[ex[expNr].firstMeasure+int(double((i-1)*(ex[expNr].nMeasure-1))/double(nPoints-1))];
	}
    }
  else 
    {        
      // evenly spaced mesh in time, this may not coincide with data points
      // fitstart and fitend set to the 1st and last time point in readData.cc
      for (i=2; i < nPoints; ++i)
	ex[expNr].mesh[i]=ex[expNr].fitstart+(ex[expNr].fitend-ex[expNr].fitstart)*double(i-1)/double(nPoints-1);
    }
}
