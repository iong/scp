#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>
#include <gsl/gsl_rng.h>
#include <xmmintrin.h>

#ifdef __ICC
#include <mkl_cblas.h>
#elif __APPLE__
#include <vecLib/vecLib.h>
#else
#include <cblas.h>
#endif


#include "vectmath.h"

int	N, NMC, nfine, nstreams, nT;
double	(*r)[3], (*rnew)[3], (*rmin)[3], (*f)[3], *y, *invmeff, *sqrtmeff;
double	*U, *Ucache, *Ucacheline, *Umin, *U0;
int	*naccepted, nswap=100;
double	epsilon, sigma, sigma6, mass;
double	*refkT, *refbeta, refTmin, refTmax, *kT, *beta, *Z, Tmin, Tmax;
double	dt, ddt, tanneal, tstop, R0, R0sq, Omega, OmegaSq, dUrel;
double	t;
int	minidx, Uminiter;
double	Uminall;

double	bl;
static double *xstep, *xsteps;
FILE	*eout;
char *inputf=NULL, *outputf;

gsl_rng *rng;





double sg_c[4] = {96609.488289873, 14584.62075507514, -365.460614956589, -19.5534697800036};
double sg_a[4] = {1.038252215127, 0.5974039109464, 0.196476572277834, 0.06668611771781};
int sg_n = 4;


double LJC[3] = {10998.6151589526, -0.165282225586247, -6.46198166728172};
double LJA[3] = {8.81337773201576,  0.36684190090077,   1.43757007579231};
int LJ_ngauss = 3;


/*
extern void vgwquenchspb_(double *q, double *f, double *W, double *enrg,
			double *taumax, double *taui, double *atol,
			double *rc, double *y);
extern void vgwquenchspb_(int *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
*/

//extern void vgw0_(int *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

extern void vgwinit_(double *, int *, double *, double *, double *, double *,
		     double *, double *);
extern void vgw0_(int *N, double *q0, double *Ueff, double *taumax, double *taumin, double *y);

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
	
	sprintf(label, "(N2)%d Umin = %lg", N, lUmin);
	dumpxyz(outputf, rmin + imin*N, label);
}


static inline double ULJ(int i, int j, int nrep)
{
	double	dr[3], dr2, y6, Uij;
	int	ioff, joff;
	
	ioff = nrep*N + i;
	joff = nrep*N + j;

	DOTPSUBV(dr2, dr, rnew[joff], rnew[ioff]);
	y6 = sigma6/(dr2*dr2*dr2);
	Uij = 4.0*epsilon*(y6*y6-y6);
		
	return Uij;
}

static double ljUtot(int nrep)
{
	double	Utot;
	int	i, j;
	
	Utot=0;
	for (i=0; i<N; i++) {
		for (j=i+1; j<N; j++) {
			Utot += ULJ(i, j, nrep);
		}
	}
	
	return Utot;
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


static void accept_trial(int nrep, double Unew)
{
	U0[nrep] = Unew;
	memcpy(r[nrep*N], rnew[nrep*N], 3*N*sizeof(double));
}


static void newtrial(int nrep)
{
	double	polar[3], move[3];
	int	i;
	
	for (i=1; i<N; i++) {
		polar[0] = gaussran(xstep[nrep], 0);
		polar[1] = acos(2.0*gsl_rng_uniform(rng) - 1.0);
		polar[2] = 2.0*M_PI*gsl_rng_uniform(rng);
		pol2cart(polar, move);
		ADDV(rnew[i+nrep*N], r[i+nrep*N], move);
	}
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
	
	p = exp((refbeta[i] - refbeta[j]) * (U0[i] - U0[j]));
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
		j = 1 + gsl_rng_uniform_int(rng, N-1);
		joff = j+N*nrep;
		mc_move_atom(joff, xsteps[nrep]);
		
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
		p = exp(-(Unew - U0[nrep]) * refbeta[nrep]);
		
		if (p > gsl_rng_uniform(rng)) {
			accept_move(nrep, j, Unew);
			naccepted[nrep]++;
		}
		if (U0[nrep] < Umin[nrep]) {
			Umin[nrep] = U0[nrep];
			memcpy(rmin[nrep*N], r[nrep*N], 3*N*sizeof(double));
		}
		if ((i % acceptance_trials) == 0) {
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
		} while ((rsq > R0sq) || (dr2min < 0.5*sigma2) );
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
		} else if (strcmp(valname, "mass") == 0) {
			mass = atof(valstring)/48.5086;
			required++;
		} else if (strcmp(valname, "NMC") == 0) {
			NMC = atoi(valstring);
			required++;
		} else if (strcmp(valname, "refTmax") == 0) {
			refTmax = atof(valstring);
			required++;
		} else if (strcmp(valname, "refTmin") == 0) {
			refTmin = atof(valstring);
			required++;
		} else if (strcmp(valname, "nstreams") == 0) {
			nstreams = atoi(valstring);
			required++;
		} else if (strcmp(valname, "Tmax") == 0) {
			Tmax = atof(valstring);
			required++;
		} else if (strcmp(valname, "Tmin") == 0) {
			Tmin = atof(valstring);
			required++;
		} else if (strcmp(valname, "nT") == 0) {
			nT = atoi(valstring);
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


static double heat_capacity(const int j, const double Z[], const double T[])
{
	double	left_der, right_der, Cv;
	
	left_der = 0.25*(T[j]+T[j-1])*(T[j]+T[j-1])*log(Z[j]/Z[j-1])/(T[j]-T[j-1]);
	right_der = 0.25*(T[j]+T[j+1])*(T[j]+T[j+1])*log(Z[j+1]/Z[j])/(T[j+1]-T[j]);
	Cv = (right_der - left_der) *2.0/ (T[j+1]-T[j-1]);	
	return Cv;
}

static double heat_capacity2(const int j, const double Z[], const double beta[])
{
	double	left_der, right_der, Cv;
	
	left_der = log(Z[j]/Z[j-1])/(beta[j]-beta[j-1]);
	right_der = log(Z[j+1]/Z[j])/(beta[j+1]-beta[j]);
	Cv = 2.0*beta[j]*beta[j]*(right_der - left_der)/ (beta[j+1]-beta[j-1]);	
	return Cv;
}

int main (int argc, const char * argv[])
{
	struct timeval tv;
	int	n, i;
	double	atoler, rc, taumin, imass;
	FILE	*Zout;
	
	
	dUrel = 1e-6;
	read_options(argv[1]);
	eout = fopen("eout.dat", "w");



	//dUrel = 0.5*dUrel/(double)(N*N);
	
	r = (double (*)[3])calloc(3*N*nstreams, sizeof(double));
	rnew = (double (*)[3])calloc(3*N*nstreams, sizeof(double));
	rmin = (double (*)[3])calloc(3*N*nstreams, sizeof(double));
	f = (double (*)[3])calloc(3*N, sizeof(double));
	y = (double *)calloc(1+21*N, sizeof(double));
	invmeff = (double *)calloc(9*N, sizeof(double));
	sqrtmeff = (double *)calloc(9*N, sizeof(double));

	refkT = (double *)calloc(nstreams, sizeof(double));
	refbeta = (double *)calloc(nstreams, sizeof(double));
	kT = (double *)calloc(nT, sizeof(double));
	beta = (double *)calloc((nT+1), sizeof(double));
	Z = (double *)calloc((nT+1), sizeof(double));
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
		refkT[i] = refTmin*pow(refTmax/refTmin, (double)i/(double)(nstreams-1));
		refbeta[i] = 1.0/refkT[i];
		Umin[i] = 1e100;
	}
	for (i=1; i<nstreams; i++) {
		memcpy(r[N*i], r[0], 3*N*sizeof(double));
	}
	/*
	for (i=0; i<nT; i++) {
		beta[nT-i-1] = 1.0/Tmax + (double)i/(double)(nT-1)*(1.0/Tmin-1.0/Tmax);
		kT[nT-i-1] = 1.0/beta[nT-i-1];
	}
	 */
	for (i=0; i<nT; i++) {
		kT[i] = Tmin + (double)i/(double)(nT-1)*(Tmax-Tmin);
		beta[i] = 1.0/kT[i];
	}
	beta[nT] = 0.0;
	
	bl = 10*R0;
	cblas_dscal(LJ_ngauss, epsilon, LJC, 1);
	cblas_dscal(LJ_ngauss, 1.0/(sigma*sigma), LJA, 1);
	
	atoler = 1e-4;
	taumin=1e-4;
	imass = 2.0;
	rc = 4.0*sigma;
	vgwinit_(&imass, &sg_n, sg_c, sg_a, &bl, &rc, &taumin, &atoler);

	dumpxyz("startcfg.xyz", r, "X trial");
	memset(Z, 0, nT*sizeof(double));
	Zout = fopen("Z.dat", "w");
	for (n=0; n<NMC; n += 100) {

		for (i=0; i<nstreams; i++) {
			mc_one_by_one(1000, i);
		}
		
		for (i=nT-1; i>=0; i--) {
			double	Ueff;
			vgw0_(&N, rnew[0], &Ueff, beta+i, beta+(i+1), y);
			Z[i] += exp(-beta[i+1]*Ueff + refbeta[0]*U0[0]);
		}
		swap_streams();
		if(n%10000==0) {
			dump_Umin(n);
			rewind(Zout);
			ftruncate(fileno(Zout), 0);
			for (i=1; i<(nT-2); i++) {
				fprintf(Zout, "%lg %lg %lg %lg %lg\n", kT[i],
					heat_capacity(i, Z, kT),
					heat_capacity2(i, Z, beta),
					beta[i], Z[i]);
			}
			fflush(Zout);
		}
	}
	dump_Umin(n);
	
	
	fclose(eout);
	fclose(Zout);
    return 0;
}
