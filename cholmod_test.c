#include <stdio.h>
#include <time.h>
#include <cholmod.h>

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
static double tic (void) { return (clock () / (double) CLOCKS_PER_SEC) ; }
static double toc (double t) { double s = tic () ; return (CS_MAX (0, s-t)) ; }

cholmod_sparse *read_csr(const char *fname, cholmod_common *c)
{
	FILE *f;
	cholmod_sparse *A;
	int	*Ai, *Ap, i, M, N, nnz, iv;
	double	dv, *Ax;

	f = fopen(fname, "r");
	fscanf(f, "%d %d %d", &M, &N, &nnz);
		
	A = cholmod_allocate_sparse(M, N, nnz, 1, 1, -1, CHOLMOD_REAL, c);
	Ai = (int *) A->i;
	Ap = (int *) A->p;
	Ax = (double *) A->x;
	for (i=0; i<=M; i++) {
		fscanf(f, "%d", &iv);
		Ap[i] = iv - 0;
	}
	for (i=0; i<nnz; i++) {
		fscanf(f, "%d", &iv);
		Ai[i] = iv - 0;
	}
	for (i=0; i<nnz; i++) fscanf(f, "%lg", Ax + i);

	fclose(f);

	return A;
}

void factor_super_diagonal(cholmod_factor *L, double Ld[])
{
    double  *Lx;
    int *Lpi, *Lpx, *Super;
    int s, nsuper, k1, k2, psi, psend, psx, nscol, nsrow, jj;
    
    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Lx = L->x;
    Super = L->super ;
        
    for (s = 0 ; s < nsuper ; s++)
    {
        k1 = Super [s] ;
        k2 = Super [s+1] ;
        psi = Lpi [s] ;
        psend = Lpi [s+1] ;
        psx = Lpx [s] ;
        nsrow = psend - psi ;
        nscol = k2 - k1 ;
        
        
        for (jj = 0 ; jj < nscol ; jj++)
        {
            Ld[jj + k1]=Lx[psx + jj + jj*nsrow];
        }
    }
    
}


int main (int argc, char **argv)
{
    cholmod_sparse *A, *Ls ;
    cholmod_factor *L ;
    cholmod_common c ;
    double t;
	double *Ax, *Lx, *Ld, logDetA, logDetA2;
	int	*Ai, *Ap, *Lp;
	
	size_t	i;
	
    /* stype < 0 fastest option for the natural ordering */
    /* stype > 0 fastest option for factoring a permuted matrix */
    
    cholmod_start (&c) ;
	/* convert to packed LL' when done */
    c.supernodal = CHOLMOD_AUTO;
    c.nmethods = 0;
	
	A = read_csr(argv[1], &c);
	t = tic();
	
	Ai = (int *)A->i;
	Ap = (int *)A->p;
	Ax = (double *)A->x;


	for (i=0; i<=A->ncol; i++) Ap[i]--;
	for (i=0; i<Ap[A->ncol]; i++) Ai[i]--;
	
    L = cholmod_analyze (A, &c) ;		    /* analyze */
    printf("The selected method is :%d\n", c.selected);
    printf("Supernodal: %d\n", L->is_super);
    printf("L->is_ll: %d\n", L->is_ll);
    printf("L->nsuper: %d\n", L->nsuper);
    cholmod_factorize (A, L, &c) ;		    /* factorize */
    
    Ld = (double *)malloc(L->n * sizeof(double));
    factor_super_diagonal(L, Ld);
    Ls = cholmod_factor_to_sparse(L, &c);
    
    cholmod_band_inplace(0,0,1, Ls, &c);
    Lx = (double *)Ls->x;
    Lp = (int *)Ls->p;
    logDetA = 0.0;
    logDetA2 = 0.0;

    for (i=0; i<L->n; i++) {
        logDetA = logDetA + 2.0*log(fabs(Ld[i]));
        logDetA2 = logDetA2 + 2.0*log(fabs(Lx[Lp[i]]));
    }
    
    printf("logdetA = %lg, logdetA2 = %lg\n", logDetA, logDetA2);
    	
	
    cholmod_free_factor (&L, &c) ;
    cholmod_free_sparse (&A, &c) ;
    cholmod_finish (&c) ;
    return (0) ;
}
