#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <cholmod.h>

double cholmod_logdet(double *G, int *ia, int *ja, int N)
{
    static cholmod_common c;
    static cholmod_sparse *A, *Ls;
    static cholmod_factor *L ;
    
    int     i, nnz, s;
    int     *Ai, *Ap, *Lp;
    double    *Ax, *Lx;

	int	ldabsr, mblk, info;
	int job[6] = { 1, 0, 1, 0, 0, 1 };
	char matdescra[7]="GxxFxx";

    double logDetA;

	FILE	*f;

    cholmod_start (&c);
    c.final_asis = 0 ;
    c.final_super = 0 ;
    c.final_ll = 1 ;
    c.final_pack = 1 ;
    c.final_monotonic = 1 ;
	c.nmethods = 1 ;
	c.method [0].ordering = CHOLMOD_NATURAL ;
	c.postorder = 0 ;

    nnz = ia[N] - 1;
    A = cholmod_allocate_sparse(N, N, nnz, 1, 1, -1, CHOLMOD_REAL, &c);
    Ai = (int *)A->i;
    Ap = (int *)A->p;

    memcpy(A->x, G, nnz*sizeof(double));
    for (i=0; i<=N; i++) {
	    Ap[i] = ia[i] - 1;
    }

    for (i=0; i<nnz; i++) {
	    Ai[i] = ja[i] - 1;
    }

    L = cholmod_analyze (A, &c) ;
    cholmod_factorize (A, L, &c) ;

	if (c.status == CHOLMOD_NOT_POSDEF) {
		return 1234.1234;
		f = fopen("G_cholmod.mtx", "w+");
		cholmod_write_sparse(f, A, NULL, NULL, &c);
		fclose(f);
		exit(EXIT_FAILURE);
		return 0.0;
	}

    Ls = cholmod_factor_to_sparse(L, &c);

    cholmod_band_inplace(0,0,1, Ls, &c);

    
    Lx = (double *)Ls->x;
    Lp = (int *)Ls->p;
    logDetA = 0.0;

    s = 1;
    for (i=0; i<N; i++) {
        if (Lx[Lp[i]] < 0.0) {
            s = -s;
        }
        logDetA = logDetA + 2.0*log(fabs(Lx[Lp[i]]));
    }


    if (s<0) {
        fprintf(stderr, "det(G) < 0!\n");
        exit(EXIT_FAILURE);
    }

    cholmod_free_factor (&L, &c) ;
    cholmod_free_sparse (&A, &c) ;
    cholmod_free_sparse (&Ls, &c) ;

    cholmod_finish (&c) ;

    return logDetA;
}
