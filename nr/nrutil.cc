#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void cmat2fvec(double **cmat,long n,long m,double *fvec)
{
  long i,j;
  for(i=1;i<=n;i++)
    {
      for(j=1;j<=m;j++)
	{
	  fvec[(j-1)*n+i-1]=cmat[i][j];
	}
    }
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v=0;
	unsigned long sz=(nh-nl+1);
	if(sz>0){
	    v=(double *)malloc((size_t) ((sz+NR_END)*sizeof(*v)));
	    if (!v) nrerror("allocation failure in dvector()");
	    return v-nl+NR_END;
	}
	return v;
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    return dvector(nl,nh);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v=0;
	unsigned long sz=(nh-nl+1);
	if(sz>0){
	    v=(int *)malloc((size_t) ((sz+NR_END)*sizeof(*v)));
	    if (!v) nrerror("allocation failure in ivector()");
	    return v-nl+NR_END;
	}
	return v;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v=0;
	unsigned long sz=(nh-nl+1);
	if(sz>0){
	    v=(unsigned char *)malloc((size_t) ((sz+NR_END)*sizeof(*v)));
	    if (!v) nrerror("allocation failure in dvector()");
	    return v-nl+NR_END;
	}
	return v;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v=0;
	unsigned long sz=(nh-nl+1);
	if(sz>0){
	    v=(unsigned long *)malloc((size_t) ((sz+NR_END)*sizeof(*v)));
	    if (!v) nrerror("allocation failure in dvector()");
	    return v-nl+NR_END;
	}
	return v;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m=0;
	
	if(nrow!=0){
	    m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(*m)));
	    if (!m) nrerror("allocation failure 1 in dmatrix()");
	    m += NR_END;
	    m -= nrl;
	}

	
	/* allocate rows and set pointers to them */
	if((nrow*ncol)!=0){
	    m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(**m)));
	    if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
	    m[nrl] += NR_END;
	    m[nrl] -= ncl;
	    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	}else{
	    for(i=nrl+1;i<=nrh;i++) m[i]=0;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    return dmatrix(nrl,nrh,ncl,nch);
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m=0;

	if(nrow!=0){
	    m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(*m)));
	    if (!m) nrerror("allocation failure 1 in imatrix()");
	    m += NR_END;
	    m -= nrl;
	}
	
	/* allocate rows and set pointers to them */
	if((nrow*ncol)!=0){
	    m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(**m)));
	    if (!m[nrl]) nrerror("allocation failure 2 in imatrix()");
	    m[nrl] += NR_END;
	    m[nrl] -= ncl;
	    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	}else{
	    for(i=nrl+1;i<=nrh;i++) m[i]=0;
	}
	/* return pointer to array of pointers to rows */
	return m;
}
char **cmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	char **m=0;

	if(nrow!=0){
	    m=(char **) malloc((size_t)((nrow+NR_END)*sizeof(*m)));
	    if (!m) nrerror("allocation failure 1 in cmatrix()");
	    m += NR_END;
	    m -= nrl;
	}
	
	/* allocate rows and set pointers to them */
	if((nrow*ncol)!=0){
	    m[nrl]=(char *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(**m)));
	    if (!m[nrl]) nrerror("allocation failure 2 in cmatrix()");
	    m[nrl] += NR_END;
	    m[nrl] -= ncl;
	    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	}else{
	    for(i=nrl+1;i<=nrh;i++) m[i]=0;
	}
	/* return pointer to array of pointers to rows */
	return m;
}


double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	double **m=0;

	if(nrow>0){
	    /* allocate array of pointers to rows */
	    m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	    if (!m) nrerror("allocation failure in submatrix()");
	    m += NR_END;
	    m -= newrl;
	}
	if((nrow*ncol)>0){
	    /* set pointers to rows */
	    for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m=0;
	if(nrow>0){
	/* allocate pointers to rows */
	    m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	    if (!m) nrerror("allocation failure in convert_matrix()");
	    m += NR_END;
	    m -= nrl;
	}
	/* set pointers to rows */
	if((nrow*ncol)>0){
	    m[nrl]=a-ncl;
	    for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	    /* return pointer to array of pointers to rows */
	}
	return m;
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t=0;

	if(nrow>0){
	    /* allocate pointers to pointers to rows */
	    t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	    if (!t) nrerror("allocation failure 1 in d3tensor()");
	    t += NR_END;
	    t -= nrl;
	}

	if((nrow*ncol)>0){
	    /* allocate pointers to rows and set pointers to them */
	    t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	    if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
	    t[nrl] += NR_END;
	    t[nrl] -= ncl;
	}else{
	    for(j=nrl;j<=nrh;j++)
		t[j]=0;
	}

	if((nrow*ncol*ndep)>0){
	    /* allocate rows and set pointers to them */
	    t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	    if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
	    t[nrl][ncl] += NR_END;
	    t[nrl][ncl] -= ndl;
	    for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	    for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	    }
	}else{
	    for(i=nrl;j<=nrh;j++)
		for(j=ncl;i<=nch;i++)
		    t[i][j]=0;
	}
	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
    if(v>0)
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
    if(v>0)
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
    if(v>0)
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
    if(v>0)
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
    if(v>0)
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
    if(m!=0){
	if(m[nrl]!=0)
	    free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
    }
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
    free_dmatrix(m,nrl,nrh,ncl,nch);
}


void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
    if(m>0){
	if(m[nrl]>0)
	    free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
    }
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
    if(m>0){
	if(m[nrl]>0)
	    free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
    }
}




void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
    if(b>0){
	free((FREE_ARG) (b+nrl-NR_END));
    }
}

void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
    if(b>0){
	free((FREE_ARG) (b+nrl-NR_END));
    }
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a double d3tensor allocated by d3tensor() */
{
    if(t>0){
	if(t[nrl]>0){
	    if(t[nrl][ncl]>0){
		free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	    }
	    free((FREE_ARG) (t[nrl]+ncl-NR_END));
	}
	free((FREE_ARG) (t+nrl-NR_END));
    }
}

