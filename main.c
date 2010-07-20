#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>
#include <gsl/gsl_rng.h>
#include <xmmintrin.h>

#include <mpi.h>

#ifdef __ICC
#include <mkl_cblas.h>
#elif __APPLE__
#include <vecLib/vecLib.h>
#else
#include <cblas.h>
#endif


#include "vectmath.h"


struct silvera_goldman {
	double alpha, beta, gamma, C6, C8, C9, C10, rc, rc2;
};

int	myrank, ncpu;
int	N, nline, NMC, nstreams, nT;
double	(*r)[3], (*rnew)[3], (*rmin)[3], (*f)[3], *y, *invmeff, *sqrtmeff;
double	*U, *Ucache, *Ucacheline, *Umin, *globalUmin, *U0;
int	*naccepted, *ntrials, nswap=100;
double	epsilon, sigma, sigma6, mass;
double	*kT, *beta, Tmin, Tmax, *Z, *Ztotal;
double	t;
int	minidx, Uminiter;
double	Uminall;

double	rho, ulen, bl, bl2, Rc;
static double *xstep, *xsteps;
FILE	*eout;
char *inputf=NULL, *outputf;

gsl_rng *rng;





double gexp_c[10], gexp_a[10];
int	gexp_n;

struct silvera_goldman sg = {14.375784, 2.961819, 0.035471, 84105.640272, 417373.713331, 146845.504557, 2613703.276242, 4.402116, 19.37862};
double sg_c[3] = {31319, -282.51, -9.2781};
double sg_a[3] = {0.68069, 0.17253, 0.050958};
int sg_n = 3;


double lj_c[3] = {10998.6151589526, -0.165282225586247, -6.46198166728172};
double lj_a[3] = {8.81337773201576,  0.36684190090077,   1.43757007579231};
int lj_n = 3;


extern void vgwinit_(int *N, double *, int *, double *, double *, double *, double *,
		     double *, double *);
extern void vgw0_(double *q0, double *Ueff, double *taumax, double *taumin, double *y);

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


static void dump_acceptance_rates()
{
	double	p;
	int	i;
	FILE	*aout;

	aout = fopen("acceptance_rates.dat", "w");
	for (i=0; i<(nstreams-1); i++) {
		p = exp((beta[i] - beta[i+1]) * (U0[i] - U0[i+1]));
		fprintf(aout, "%d %lg\n", i, p);
	}
	fclose(aout);
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


static void dump_Umin(int myrank, int MCiter)
{
	char	label[128];
	double	lUmin;
	int	i, imin;
	
	MPI_Allreduce(Umin, globalUmin, nstreams, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	if (myrank==0) {
		printf("n = %d, Umin =", MCiter);
		for (i=0; i<nstreams; i++) {
			printf(" %lg,", globalUmin[i]);
		}
		fputc('\n', stdout);
	}

	lUmin = globalUmin[0];
	for (i=0; i<nstreams; i++) {
		if (globalUmin[i] <= lUmin) {
			lUmin = globalUmin[i];
			imin = i;
		}
	}
	
	if (fabs(Umin[imin] - lUmin) < 0.1) {
		sprintf(label, "(N2)%d Umin = %lg", N, lUmin);
		dumpxyz(outputf, rmin + imin*N, label);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

static void accept_trial(int nrep, double Unew)
{
	U0[nrep] = Unew;
	memcpy(r[nrep*N], rnew[nrep*N], 3*N*sizeof(double));
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
	const double	beta0=0.0;
	double Unew, p, Ulmin;
	const int acceptance_trials = 1000;
	int	i, j, joff, k;
	
	memcpy(rnew[nrep*N], r[nrep*N], 3*N*sizeof(double));
	vgw0_(rnew[0], U0+nrep, beta+nrep, &beta0, y);
	Ulmin = 1e6;
	for (i = 1; i <= burnin_len; i++) {
		ntrials[nrep]++;
		j = 1 + gsl_rng_uniform_int(rng, N-1);
		joff = j+N*nrep;
		mc_move_atom(joff, xsteps[nrep]);
		
		/* wrap around the box */
		for (k=0; k<3; k++) {
			if (rnew[joff][k] > bl2) {
				rnew[joff][k] = rnew[joff][k] - bl;
			} else if (rnew[joff][k] < -bl2) {
				rnew[joff][k] = rnew[joff][k] + bl;
			} 
		}
		
		vgw0_(rnew[0], &Unew, beta+nrep, &beta0, y);
		p = exp(-(Unew - U0[nrep]) * beta[nrep]);
		
		if (p > gsl_rng_uniform(rng)) {
			accept_move(nrep, j, Unew);
			naccepted[nrep]++;
		}
		if (U0[nrep] < Umin[nrep]) {
			Umin[nrep] = U0[nrep];
			memcpy(rmin[nrep*N], r[nrep*N], 3*N*sizeof(double));
		}
		if ((ntrials[nrep] % acceptance_trials) == 0) {
			if (naccepted[nrep] < 0.3 * (double) acceptance_trials) {
				xsteps[nrep] = xsteps[nrep] / 1.10779652734191;
			} else 
			if (naccepted[nrep] > 0.4 * (double) acceptance_trials) {
				xsteps[nrep] = xsteps[nrep] * 1.12799165273419;
			}
			
			naccepted[nrep] = 0;
		}
	}
}


static void init_unit_cell()
{
	double	x, y, z;
	int	i, j, k, idx;

	for (i=0; i<nline; i++) {
		x = (0.5 + (double)i)*ulen - 0.5*bl;
		for (j=0; j<nline; j++) {
			y = (0.5 + (double)j)*ulen - 0.5*bl;
			for (k=0; k<nline; k++) {
				z = (0.5 + (double)k)*ulen - 0.5*bl;
				idx = nline*(j + i*nline) + k;
				r[idx][0]=x;
				r[idx][1]=y;
				r[idx][2]=z;
			}
		}
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
	int	i, required;
	FILE *fin;
	
	fin = fopen(fname, "r");
	required = 0;
	while (fscanf(fin, "%s %s", valname, valstring) == 2) {
		if (strcmp(valname, "include") == 0) {
			read_options(valstring);
		}
		else if (strcmp(valname, "nline") == 0) {
			nline = atoi(valstring);
			N = nline*nline*nline;
			required++;
		} else if (strcmp(valname, "LJ") == 0) {
			epsilon = atof(valstring);
			fscanf(fin, "%lg", &sigma);
			sigma6 = pow(sigma, 6.0);

			for (i=0; i<lj_n; i++) {
				gexp_c[i] = epsilon*lj_c[i];
				gexp_a[i] = lj_a[i]/(sigma*sigma);
			}
			gexp_n = lj_n;
			Rc = 4.0*sigma;
	
			required++;
		} else if (strcmp(valname, "SG") == 0) {
			memcpy(gexp_c, sg_c, sg_n*sizeof(double));
			memcpy(gexp_a, sg_a, sg_n*sizeof(double));
			gexp_n = sg_n;
			Rc = 4.0*sg.rc;
			required++;
		} else if (strcmp(valname, "mass") == 0) {
			mass = atof(valstring)/48.5086;;
			required++;
		} else if (strcmp(valname, "NMC") == 0) {
			NMC = atoi(valstring);
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
		} else if (strcmp(valname, "rho") == 0) {
			rho = atof(valstring);
			ulen = cbrt(1.0/rho);
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

static void init_mem()
{
	r = (double (*)[3])calloc(3*N*nstreams, sizeof(double));
	rnew = (double (*)[3])calloc(3*N*nstreams, sizeof(double));
	rmin = (double (*)[3])calloc(3*N*nstreams, sizeof(double));
	f = (double (*)[3])calloc(3*N, sizeof(double));
	y = (double *)calloc(1+21*N, sizeof(double));
	invmeff = (double *)calloc(9*N, sizeof(double));
	sqrtmeff = (double *)calloc(9*N, sizeof(double));

	kT = (double *)calloc(nstreams, sizeof(double));
	beta = (double *)calloc(nstreams, sizeof(double));
	Z = (double *)calloc(nstreams, sizeof(double));
	Ztotal = (double *)calloc(nstreams, sizeof(double));
	Umin = (double *)calloc(nstreams, sizeof(double));
	globalUmin = (double *)malloc(nstreams*sizeof(double));
	U0 = (double *)calloc(nstreams, sizeof(double));
	xstep = (double *)calloc(nstreams, sizeof(double));
	xsteps = (double *)calloc(nstreams, sizeof(double));
	naccepted = (int *)calloc(nstreams, sizeof(int));
	ntrials = (int *)calloc(nstreams, sizeof(int));
	Ucache = (double *)malloc(N*N*sizeof(double));
	Ucacheline = (double *)malloc(N*sizeof(double));
}

int main (int argc,  char * argv[])
{
	struct timeval tv;
	int	n, i;
	double	atoler, taumin;
	FILE	*Zout;
	char	fname[200];
	
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
	
	read_options(argv[1]);
	init_mem();

	bl = nline * ulen;
	bl2 = 0.5*bl;
	

	gettimeofday(&tv, NULL);
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, tv.tv_usec+myrank);

	if (inputf==NULL) {
		init_unit_cell();
	}
	else {
		load_data(inputf);
	}

	for (i=0; i<nstreams; i++) {
		xsteps[i] = xstep[i] = 0.125*ulen;
		kT[i] = Tmin*pow(Tmax/Tmin, (double)i/(double)(nstreams-1));
		beta[i] = 1.0/kT[i];
		Umin[i] = 1e100;
	}
	for (i=1; i<nstreams; i++) {
		memcpy(r[N*i], r[0], 3*N*sizeof(double));
	}

	atoler = 1e-4;
	taumin=1e-2;
	vgwinit_(&N, &mass, &gexp_n, gexp_c, gexp_a, &bl, &Rc, &taumin, &atoler);
	
	if (myrank==0) {
		eout = fopen("eout.dat", "w");
		dumpxyz("startcfg.xyz", r, "X trial");
	}
		
	sprintf(fname, "Z%02d.dat", myrank);
			
	for (n=0; n<NMC; n += 1) {

		for (i=0; i<nstreams; i++) {
			mc_one_by_one(1000, i);
			Z[i] += exp(-beta[i]*U0[i]);
		}
		dump_acceptance_rates();
		swap_streams();
				
		if(n%1000==0) {
			MPI_Allreduce(Z, Ztotal, nT, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			dump_Umin(myrank, ncpu*n);
		}

		if(n%1000==0 && myrank==0) {
			Zout = fopen("Z.dat", "w");
			/*
			rewind(Zout);
			ftruncate(fileno(Zout), 0);
			 */
			for (i=1; i<(nT-2); i++) {
				fprintf(Zout, "%lg %lg %lg %lg %lg %lg\n", kT[i],
					heat_capacity(i, Ztotal, kT),
					heat_capacity2(i, Ztotal, beta),
					beta[i], Z[i], Ztotal[i]);
			}
			fflush(Zout);
			fclose(Zout);
		}
	}

	dump_Umin(myrank, ncpu*n);

	MPI_Finalize();
	return 0;
}
