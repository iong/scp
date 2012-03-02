#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <cholmod.h>

typedef struct cholmod_state_struct {
     cholmod_common *c;
     cholmod_sparse *A;
     cholmod_factor *L;
} cholmod_state;


void cholmod_init(cholmod_state *s, int N, int *ia, int *ja, double *G)
{
     s->c = (cholmod_common *)malloc(sizeof(cholmod_common));
     cholmod_start (s->c);
     s->c->supernodal = CHOLMOD_SUPERNODAL ;
     s->c->final_super = 1 ;

     s->A = (cholmod_sparse *) malloc(sizeof(cholmod_sparse));

     s->A->nrow = N ;
     s->A->ncol = N ;
     s->A->nzmax = ia[N] - ia[0];
     s->A->packed = 1 ;    /* default is packed (A->nz not present) */
     s->A->sorted = 1;
     s->A->stype = 1 ;
     s->A->itype = s->c->itype ;
     s->A->xtype = CHOLMOD_REAL ;
     s->A->dtype = s->c->dtype ;

     s->A->nz = NULL ;
     s->A->p = ia ;
     s->A->i = ja ;
     s->A->x = G ;
     s->A->z = NULL ;

     s->L = cholmod_analyze (s->A, s->c) ;

      printf("The selected method is :%d\n", s->c->selected);
}

static double logdet_supernodal_factor(cholmod_factor *L)
{
     double logdet;
     double  *Lx;
     int *Lpi, *Lpx, *Super;
     int s, nsuper, k1, k2, psi, psend, psx, nscol, nsrow, jj;
    
     nsuper = L->nsuper ;
     Lpi = L->pi ;
     Lpx = L->px ;
     Lx = L->x;
     Super = L->super ;

     logdet = 0.0;
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
	       logdet += log(fabs(Lx[psx + jj + jj*nsrow]));
	  }
     }
     logdet *= 2.0;

     printf ("logdet = %lg\n", logdet);
     return logdet;
}


double cholmod_logdet(cholmod_state *s, double *G)
{    
     s->A->x = G;
     s->A->xtype = CHOLMOD_REAL ;
     cholmod_factorize(s->A, s->L, s->c);

     return logdet_supernodal_factor(s->L);
}

void cholmod_cleanup(cholmod_state *s)
{
     if (s->c == NULL) return;

     cholmod_free_factor (&(s->L), s->c) ;
     free(s->A);
     cholmod_finish (s->c) ;
     free(s->c);
}
