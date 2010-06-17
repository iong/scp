#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>
#include <vecLib/vecLib.h>
#include <gsl/gsl_rng.h>
#include <xmmintrin.h>

#include "vectmath.h"

int	N, NMC, nfine, nstreams;
double	(*r)[3], (*rnew)[3], (*rmin)[3], (*a)[3], (*v)[3];
double	*U, *Ucache, *Ucacheline, *Umin, *U0;
int	*naccepted, nswap=100;
double	epsilon, sigma, sigma6, mass;
double	*kT, *beta, Tmin, Tmax;
double	dt, ddt, tanneal, tstop, R0, R0sq, Omega, OmegaSq, dUrel;
double	t;
int	minidx, Uminiter;
double	Uminall;

static double *xstep, *xsteps;
FILE	*eout;
char *inputf=NULL, *outputf, *reference_file;

gsl_rng *rng;

static inline void pol2cart(double *src, double *dest)
{
	dest[0] = (src[0] * sin(src[1])) * cos(src[2]);
	dest[1] = (src[0] * sin(src[1])) * sin(src[2]);
	dest[2] = src[0] * cos(src[1]);
}

/**
 \short \f[ f(x) = exp(-(x-x_0)^2/(2*\sigma^2))/(\sqrt{2*\pi}*\sigma) \f]
 */
static inline double gaussran(double sigma, double x0)
{
	double x1, x2, y;
	
	x1 = gsl_rng_uniform(rng);
	x2 = gsl_rng_uniform(rng);
	y = sqrt(-2 * log(x1)) * cos(2 * M_PI * x2);
	
	return y * sigma + x0;
}


static void dumpxyz(char *fname, double (*rr)[3], char *label)
{
	FILE	*fout;
	int	i;
	
	
	fout = fopen(fname, "w");
	fprintf(fout, "%d\n%s\n", N, label);
	for (i=0; i<N; i++) {
		fprintf(fout, "N %lg %lg %g\n", rr[i][0], rr[i][1], rr[i][2]);
	}
	fclose(fout);
}


static int converged(int mciter)
{
	int	i;
	
	minidx = 0;
	for (i=1; i<nstreams; i++) {
		if (Umin[i] < Umin[minidx]) {
			minidx = i;
		}
	}
	if (Umin[minidx] < Uminall) {
		Uminall = Umin[minidx];
		Uminiter = mciter;
	} else if (mciter - Uminiter > 1000000) {
		return 1;
	}
	return	0;
		
}


static void dump_Umin(int MCiter)
{
	char	label[128];
	double	lUmin;
	int	i, imin;
	
	lUmin = 1e100;
	printf("n = %d, Umin =", MCiter);
	for (i=0; i<nstreams; i++) {
		printf(" %lg,", Umin[i]);
		if (Umin[i] <= lUmin) {
			lUmin = Umin[i];
			imin = i;
		}
	}
	fputc('\n', stdout);
	
	sprintf(label, "(N2)%d Umin = %lg", Nmol, lUmin);
	dumpxyz(outputf, rmin + imin*N, label);
}




static inline double ULJ(int i, int j, int nrep)
{
	double	dr[3], dr2, y6, Uij;
	int	ioff, joff;
	
	ioff = nrep*Nmol + i;
	joff = nrep*Nmol + j;

	DOTPSUBV(dr2, dr, rnew[joff], rnew[ioff]);
	y6 = sigma6/(dr2*dr2*dr2);
	Uij = 4.0*epsilon*(y6*y6-y6);
		
	return Uij;
}


static double ljUtot_single_jump(int nrep, int j)
{
	double	Uij, dE;
	int	i;
	
	dE = 0;
	for (i=0; i<N; i++) {
		if (i==j) {
			continue;
		}
		Uij = ULJ(i, j, nrep);
		dE += Uij - Ucache[j*N+i];
		Ucacheline[i] = Uij;
	}
	return U0[nrep] + dE;
}



static double ljUtot_cache(int nrep)
{
	double	Uij, Utot;
	int	i, j;
	
	Utot=0;
	for (i=0; i<N; i++) {
		for (j=i+1; j<N; j++) {
			Uij = ULJ(i, j, nrep);
			Ucache[i*N+j] = Uij;
			Ucache[j*N+i] = Uij;
			Utot += Uij;
		}
	}
	
	return Utot;
}



static void mc_move_atom(int i, double step_radius)
{
	double	polar[3], move[3];
	
	polar[0] = gaussran(step_radius, 0);
	polar[1] = acos(2.0*gsl_rng_uniform(rng) - 1.0);
	polar[2] = 2.0*M_PI*gsl_rng_uniform(rng);
	pol2cart(polar, move);
	ADDV(rnew[i], r[i], move);
}



static void swap_ptrs(double *p1, double *p2, int len)
{
	memcpy(Ucache, p2, len*sizeof(double));
	memcpy(p2, p1, len*sizeof(double));
	memcpy(p1, Ucache, len*sizeof(double));
}



static void swap_streams()
{
	double	p;
	int	i, j;
	i = gsl_rng_uniform_int(rng, nstreams-1);
	j = i+1;
	
	p = exp((beta[i] - beta[j]) * (U0[i] - U0[j]));
	if (p > gsl_rng_uniform(rng)) {
		swap_ptrs(r[i*N], r[j*N], 3*N);
		swap_ptrs(U0+i, U0+j, 1);
	}
}


static void accept_move(int nrep, int m, double	Unew)
{
	int	moff, k;
	
	moff = m + N*nrep;
	memcpy(r[moff], rnew[moff], 3*sizeof(double));
	memcpy(Ucache + m*N, Ucacheline, N*sizeof(double));
	for (k=0; k<N; k++) {
		Ucache[k*N + m] = Ucacheline[k];
	}
	U0[nrep] = Unew;
}

void mc_one_by_one(int burnin_len, int nrep)
{
	double Unew, p,rsq, Ulmin;
	int acceptance_trials = 1000, i;
	int	j, joff;
	
	memcpy(rnew[nrep*N], r[nrep*N], 3*N*sizeof(double));
	U0[nrep] = ljUtot_cache(nrep);
	Ulmin = 1e6;
	for (i = 1; i <= burnin_len; i++) {
		//printf("j = ");
		j = gsl_rng_uniform_int(rng, N);
		joff = j+N*nrep;
		//printf("%d\n", j);
		mc_move_atom(j+N*nrep, xsteps[nrep]);
		
		DOTVP(rsq, rnew[joff], rnew[joff]);
		if (rsq > R0sq) {
			SETV(rnew[joff], r[joff]);
			continue;
		}
		DOTVP(rsq, r[joff], r[joff]);
		if (rsq > R0sq) {
			printf("Aaaaaaaaaaaarrrrrrrrrrrrrr\n");
			kill(getpid(), SIGSEGV);
		}
		
		Unew = ljUtot_single_jump(nrep, j);
		p = exp(-(Unew - U0[nrep]) * beta[nrep]);
		
		if (p > gsl_rng_uniform(rng)) {
			accept_move(nrep, j, Unew);
			naccepted[nrep]++;
		}
		if (U0[nrep] < Umin[nrep]) {
			Umin[nrep] = U0[nrep];
			memcpy(rmin[nrep*N], r[nrep*N], 3*N*sizeof(double));
		}
		if ((i % acceptance_trials) == 0) {
			printf("%d: %lg s%lg %d\n", nrep, Umin[nrep], xsteps[nrep], naccepted[nrep]);
			if (naccepted[nrep] < 0.5 * (double) acceptance_trials) {
				xsteps[nrep] = xsteps[nrep] / 1.10779652734191;
			} else {
				xsteps[nrep] = xsteps[nrep] * 1.12799165273419;
			}
			
			naccepted[nrep] = 0;
		}
	}
}


static void initial_config()
{
	double	rsq,sigma2, dr2, dr2min=1e100, dr[3];
	int	i, j;
	
	sigma2 = sigma*sigma;
	CLRV(r[0]);
	for (i=1; i<N; i++) {
		do {
			r[i][0] = R0*(gsl_rng_uniform(rng) - 0.5);
			r[i][1] = R0*(gsl_rng_uniform(rng) - 0.5);
			r[i][2] = R0*(gsl_rng_uniform(rng) - 0.5);
			DOTVP(rsq, r[i], r[i]);
			dr2min = 4*R0sq;
			for (j=0; j<i; j++) {
				DOTPSUBV(dr2, dr, r[i], r[j]);
				if (dr2<dr2min) {
					dr2min = dr2;
				}
			}
		} while ((rsq > R0sq) || (dr2min < 0.25*sigma2) );
	}
}


void load_data(char *coords)
{
        FILE    *fin;
        int  i;
        
        fin = fopen(coords, "r");
        if(!fin) {
                fprintf(stderr, "load_data()\n");
                exit(EXIT_FAILURE);
        }
        
        for (i=0; i<N; i++) {
                fscanf(fin, "%lg %lg %lg", &(r[i][0]), &(r[i][1]), &(r[i][2]));
        }
        fclose(fin);
}


static void read_options(const char fname[])
{
	char valname[100], valstring[100];
	int required;
	FILE *fin;
	
	fin = fopen(fname, "r");
	required = 0;
	while (fscanf(fin, "%s %s", valname, valstring) == 2) {
		if (strcmp(valname, "include") == 0) {
			read_options(valstring);
		}
		else if (strcmp(valname, "N") == 0) {
			N = atoi(valstring);
			required++;
		} else if (strcmp(valname, "LJ") == 0) {
			epsilon = atof(valstring);
			fscanf(fin, "%lg", &sigma);
			sigma6 = pow(sigma, 6.0);
			required++;
		} else if (strcmp(valname, "NMC") == 0) {
			NMC = atoi(valstring);
			required++;
		} else if (strcmp(valname, "Tmax") == 0) {
			Tmax = atof(valstring);
			required++;
		} else if (strcmp(valname, "Tmin") == 0) {
			Tmin = atof(valstring);
			required++;
		} else if (strcmp(valname, "nstreams") == 0) {
			nstreams = atoi(valstring);
			required++;
		} else if (strcmp(valname, "R0") == 0) {
			R0 = atof(valstring);
			R0sq = R0*R0;
			required++;
		} else if (strcmp(valname, "Omega") == 0) {
			Omega = atof(valstring);
			OmegaSq = Omega*Omega;
			required++;
		} else if (strcmp(valname, "input") == 0) {
			inputf = strdup(valstring);
			required++;
		} else if (strcmp(valname, "output") == 0) {
			outputf = strdup(valstring);
			required++;
		} 
	}
	
	fclose(fin);
}




int main (int argc, const char * argv[])
{
	char	fname[128];
	struct timeval tv;
	int	n, i;
	
	dUrel = 1e-6;
	read_options(argv[1]);
	eout = fopen("eout.dat", "w");



	//dUrel = 0.5*dUrel/(double)(N*N);
	
	r = (double (*)[3])calloc(3*N*nstreams, sizeof(double));
	rnew = (double (*)[3])calloc(3*N*nstreams, sizeof(double));
	rmin = (double (*)[3])calloc(3*N*nstreams, sizeof(double));
	kT = (double *)calloc(nstreams, sizeof(double));
	beta = (double *)calloc(nstreams, sizeof(double));
	Umin = (double *)calloc(nstreams, sizeof(double));
	U0 = (double *)calloc(nstreams, sizeof(double));
	xstep = (double *)calloc(nstreams, sizeof(double));
	xsteps = (double *)calloc(nstreams, sizeof(double));
	naccepted = (int *)calloc(nstreams, sizeof(int));
	Ucache = (double *)malloc(N*N*sizeof(double));
	Ucacheline = (double *)malloc(N*sizeof(double));

	
	gettimeofday(&tv, NULL);
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, tv.tv_usec);

	if (inputf==NULL) {
		initial_config();
	}
	else {
		load_data(inputf);
	}

	for (i=0; i<nstreams; i++) {
		xsteps[i] = xstep[i] = 0.2*sigma;
		kT[i] = Tmin*pow(Tmax/Tmin, (double)i/(double)(nstreams-1));
		beta[i] = 1.0/kT[i];
		Umin[i] = 1e100;
	}
	for (i=1; i<nstreams; i++) {
		memcpy(r[N*i], r[0], 3*N*sizeof(double));
	}

	dumpxyz("startcfg.xyz", r, "X trial");
	for (n=0; n<NMC; n += 100) {
		if(n%50000==0) dump_Umin(n);
		for (i=0; i<nstreams; i++) {
			mc_one_by_one(100, i);
		}
		swap_streams();
	}
	dump_Umin(n);
	
	fclose(eout);

    return 0;
}
