/* Copyright (C) 2005 The Scalable Software Infrastructure Project. All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. Neither the name of the project nor the names of its contributors 
      may be used to endorse or promote products derived from this software 
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE SCALABLE SOFTWARE INFRASTRUCTURE PROJECT
   ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
   TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
   PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE SCALABLE SOFTWARE INFRASTRUCTURE
   PROJECT BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_CONFIG_H
        #include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
        #include "lis_config_win.h"
#endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lis.h"

#undef __FUNC__
#define __FUNC__ "main"
LIS_INT main(LIS_INT argc, char* argv[])
{
        LIS_SCALAR *a,*q,*r;
	LIS_INT m,n,nn;
	LIS_INT	i,j,ii,jj,nnz,qriter;
	double time,time0;
	LIS_REAL qrerr;	

 	LIS_DEBUG_FUNC_IN;

	lis_initialize(&argc, &argv);

	if( argc < 3 )
	{
	  printf("Usage: %s m n\n", argv[0]);
	  CHKERR(1);
	}

	m  = atoi(argv[1]);
	n  = atoi(argv[2]);
	if( m<=0 || n<=0 )
	{
	  printf("m=%d <=0 or n=%d <=0\n", m,n);
	  CHKERR(1);
	}
	
	printf("\n");

	/* create arrays */

	nn = m*n;

	a = (LIS_SCALAR *)malloc(nn*nn*sizeof(LIS_SCALAR));
	q = (LIS_SCALAR *)malloc(nn*nn*sizeof(LIS_SCALAR));
	r = (LIS_SCALAR *)malloc(nn*nn*sizeof(LIS_SCALAR));

	/* define two-dimensional Laplacian */

	lis_array_set_all(nn*nn,(LIS_SCALAR)0.0,a);

	nnz = 0;
	for(ii=0;ii<nn;ii++)
	  {
	    i = ii/m;
	    j = ii - i*m;
	    if( i>0 )   { jj = ii - m; a[ii + jj * nn] = -1.0; nnz++;} 
	    if( i<n-1 ) { jj = ii + m; a[ii + jj * nn] = -1.0; nnz++;}
	    if( j>0 )   { jj = ii - 1; a[ii + jj * nn] = -1.0; nnz++;}
	    if( j<m-1 ) { jj = ii + 1; a[ii + jj * nn] = -1.0; nnz++;}
	    jj = ii; a[ii + jj * nn] = 4.0; nnz++;
	  }

	printf("matrix size = %d x %d (%d nonzero entries)\n\n", nn,nn,nnz);

	/* solve eigenproblem */

	time0 = lis_wtime();
	lis_array_qr(nn,a,q,r,&qriter,&qrerr);
	time = lis_wtime() - time0;

	printf("QR    : number of iterations = %d\n", qriter);
	printf("QR    : elapsed time         = %e sec.\n", time);
	printf("QR    :   eigensolver        = %e sec.\n", time);
#ifdef _LONG__DOUBLE
	printf("QR    : 2-norm of A(2,1)     = %Le\n\n", qrerr);
#else
	printf("QR    : 2-norm of A(2,1)     = %e\n\n", qrerr);
#endif

	lis_finalize();

	LIS_DEBUG_FUNC_OUT;
	return 0;
}


