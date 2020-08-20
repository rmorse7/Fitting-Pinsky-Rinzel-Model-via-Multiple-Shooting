//! \b Module for reading the data file of the measurements
/*! \file  */

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ctype.h>
#include <math.h>
#include <string.h>

#include "../nr.h"
#include "../def.h"

using namespace std;

//read in data for experiment number expNr
void readData(GlobExp *ex,Glob *globs,long expNr)
{
  int lineNr=0;
  long k,l,index;
  long nLines=0,maxLeng=0;
  long nobsmax=ex[expNr].nobs; //default value read from model.cc or set as parameter
  char c,*line,*arg;
  double *x;
  double **y;
  double **sig;

  ifstream in;

  FILE *fp;

  //open the file name
  fp=fopen(ex[expNr].fileName,"r");
  if (fp==NULL) {
    cerr << "Error opening data file " << ex[expNr].fileName << "\n";
    exit(1);
  }
  
  //determine total number of lines and maximum length of the lines
  while((c=fgetc(fp))!=EOF)
        {
	  for(k=1;c!='\n';k++)
	    {
	      if(k>maxLeng)
		maxLeng=k;
	      c=fgetc(fp);
	    }
	  nLines++;
	}
#ifdef DEBUGREADDATA
  *dbg << "Total number of lines in " <<  ex[expNr].fileName << " : " << nLines << endl;
  *dbg << "Maximal length of line in " << ex[expNr].fileName << " : " << maxLeng << endl;
#endif
  fclose(fp);

  in.open(ex[expNr].fileName);

  //allocate memory
  line=new char[maxLeng+2];
  arg=new char[maxLeng+2];
  x=dvector(1,nLines); // time vector
  y=dmatrix(1,nLines,1,ex[expNr].nobs); //observations 
  sig=dmatrix(1,nLines,1,ex[expNr].nobs); // standard deviation for each observation

  //extracting the data
  while(!in.eof())
        {
	  //process one line at the time
	  in.getline(line,maxLeng+2,'\n');
	  if(isspace(line[0])==FALSE && line[0]!= '#' && line[0] != '\0' )
	    {
	      lineNr++; //process line
#ifdef DEBUGREADDATA
	      *dbg << "Extracting data line no. : " << lineNr << endl;
#endif
	      l=0;
	      index=0;
	      for(k=0;k<=strlen(line);k++)
		{
		  if(isspace(line[k])!=FALSE) //found a blank space
		    {
		      arg[l]='\0'; //terminate string
		      l=0;
		      if(index==0)  //for each line we have the time value, then the y value of each observation
			            //followed by its sigma
			x[lineNr]=atof(arg);
		      else if(index <=2*ex[expNr].nobs && (index % 2)==1)
			y[lineNr][(index+1)/2]=atof(arg);
		      else if(index <=2*ex[expNr].nobs && (index % 2)==0)
			sig[lineNr][index/2]=atof(arg);
		      index++; //increase index after each data save
		    }
		  else
		    {
		      arg[l]=line[k]; //add current character to arg string 
		      l++;
		    } //found blank space or not
		} //loop over characters

	      //save the last element of the line
	      arg[l]='\0';	
	      if(index==0)
		x[lineNr]=atof(arg);
	      else if(index <=2*ex[expNr].nobs && (index % 2)==1)
		y[lineNr][(index+1)/2]=atof(arg);
	      else if(index <=2*ex[expNr].nobs && (index % 2)==0)
		sig[lineNr][index/2]=atof(arg);
	      
	      //determ. max. complete obs.
	      if((index % 2)==1)
		index--;
	      if(index==0)
		{
		  cerr << "Data line " << lineNr << " contains no data." << endl;
		  exit(1);
		}
	      if(index/2 <=  nobsmax)
		nobsmax=index/2;
	    }
		    
	}
  in.close();

  //update the number of observations if there is a mismatch with initial default value
  ex[expNr].nobs=nobsmax;
  ex[expNr].nMeasure=lineNr; //number of time points
  if(ex[expNr].nMeasure==0)
    {
      cerr << "File " << ex[expNr].fileName << " contains no data." << endl;
      exit(1);
    }

  //check if the time is in ascending order
    for(k=2;k<=ex[expNr].nMeasure;k++)
      {
	if(x[k-1] >= x[k])
	  {
	    cerr << "Data not in temporal ascending order at dataline " << k << "\n";
	    exit(1);
	  }
      }
    //update the start and end time of the fit based on file values
    if(ex[expNr].fitstart == -1e99)
      ex[expNr].fitstart=x[1];
    if(ex[expNr].fitend ==  1e99)
      ex[expNr].fitend=x[ex[expNr].nMeasure];
    
    //for memory allocation
    index=1;
    for(k=1;k<=ex[expNr].nMeasure;k++)
      {
	if((x[k] >= ex[expNr].fitstart) && (x[k] <= ex[expNr].fitend))
	  { 
	    index++;
	  }
      }
    index--; 
  
    //the actual amount of data
    ex[expNr].nMeasure=index;

    //memory allocation, time vector, observations, sigma
    ex[expNr].xMeasure=dvector(1,ex[expNr].nMeasure);
    ex[expNr].yMeasure=dmatrix(1,ex[expNr].nMeasure,1,ex[expNr].nobs);
    ex[expNr].sigma=dmatrix(1,ex[expNr].nMeasure,1,ex[expNr].nobs);

    // copy data to data structure
    index=1;
    for(k=1;k<=lineNr;k++)
      {
	if((x[k] >= ex[expNr].fitstart) && (x[k] <= ex[expNr].fitend))
	  {
	    ex[expNr].xMeasure[index]=x[k];
	    for(l=1;l<=ex[expNr].nobs;l++)
	      {
		ex[expNr].yMeasure[index][l]=y[k][l];
		ex[expNr].sigma[index][l]=sig[k][l];
	      }
	    index++;
	  }
      }
    index--;

    //setting some variables
    ex[expNr].firstMeasure=1; //index of first measure time point
    ex[expNr].lastMeasure=ex[expNr].nMeasure; //index of last measure time point
    
    //free memory
    free_dvector(x,1,nLines);
    free_dmatrix(y,1,nLines,1,ex[expNr].nobs);
    free_dmatrix(sig,1,nLines,1,ex[expNr].nobs);
}
