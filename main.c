#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vecLib/vecLib.h>
#include <gsl/gsl_rng.h>
#include <xmmintrin.h>

#include "vectmath.h"

int	N, NMC, nfine;
double	(*r)[3], (*rnew)[3], (*rmin)[3], (*a)[3], (*v)[3], *U, Utot, Umin, Ekin;
double	epsilon, sigma, sigma6, mass;
double	kT, beta;
double	dt, ddt, tanneal, tstop, R0, R0sq, Omega, OmegaSq, dUrel;
double	t;

static double xstep;
FILE	*eout;

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


static void dump_particles(double tnow)
{
	char	fname[128];
	FILE	*fout;
	int	i;
	
	snprintf(fname, 128, "%06.0lf.dat", tnow);
	fout = fopen(fname, "w");
	for (i=0; i<N; i++) {
		fprintf(fout, "%lg %lg %g %lg %lg %g %lg %lg %g %lg\n",
			r[i][0], r[i][1], r[i][2],
			v[i][0], v[i][1], v[i][2],
			a[i][0], a[i][1], a[i][2], U[i]);
	}
	fclose(fout);
}

static void dump_Umin()
{
	FILE	*fout;
	int	i;
	
	fout = fopen("Umin_cfg.dat", "w");
	for (i=0; i<N; i++) {
		fprintf(fout, "%lg %lg %g\n",
			rmin[i][0], rmin[i][1], rmin[i][2]);
	}
	fclose(fout);
}


static double ljUtot()
{
	double	dr[3], dr2, y6;
	int	i, j;
	
	Utot=0;
	for (i=0; i<N; i++) {
		for (j=i+1; j<N; j++) {
			DOTPSUBV(dr2, dr, rnew[j], rnew[i]);
			y6 = sigma6/(dr2*dr2*dr2);
			Utot += 4.0*epsilon*(y6*y6 - y6);
		}
	}
	
	return Utot;
}


static void newtrial(double step_radius)
{
	double	polar[3], move[3];
	int	i;
	
	for (i=0; i<N; i++) {
		polar[0] = gaussran(step_radius, 0);
		polar[1] = acos(2.0*gsl_rng_uniform(rng) - 1.0);
		polar[2] = 2.0*M_PI*gsl_rng_uniform(rng);
		pol2cart(polar, move);
		ADDV(rnew[i], r[i], move);
	}
}


static void accept_trial()
{
	memcpy(r, rnew, 3*N*sizeof(double));
}


static double energy()
{
	double	E;
	int	i;
	
	//Ekin = cblas_dnrm2(3*N, v[0], 1);
	Ekin = 0;
	for (i=0; i<N; i++) {
		Ekin += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
	}
	Ekin = 0.5*mass*Ekin;
	
	E = Ekin + Utot;
	
	return E;
}


void mc_metropolis(int burnin_len)
{
	double U0, Unew, p,rsq, Ulmin;
	int acceptance_trials = 1000, i, j, naccepted = 0;
	
	memcpy(rnew, r, 3*N*sizeof(double));
	U0 = ljUtot();
	Ulmin = 1e6;
	for (i = 1; i <= burnin_len; i++) {
		if (U0 < Umin) {
			Umin = U0;
			memcpy(rmin, r, 3*N*sizeof(double));
		}
		if (U0 < Ulmin) {
			Ulmin = U0;
		}
		newtrial(xstep);
		for (j=0; j<N; j++) {
			DOTVP(rsq, rnew[j], rnew[j]);
			if (rsq > R0sq) {
				continue;
			}
		}
		Unew = ljUtot();
		p = exp(-(Unew - U0) * beta);
		
		if (p > gsl_rng_uniform(rng)) {
			accept_trial();
			U0 = Unew;
			naccepted++;
		}
		if ((i % acceptance_trials) == 0) {
			printf("%lg %lg %d\n", Ulmin, xstep, naccepted);
			if (naccepted < 0.5 * (double) acceptance_trials) {
				xstep = xstep / 1.10779652734191;
			} else {
				xstep = xstep * 1.12799165273419;
			}
			
			naccepted = 0;
		}
	}
}


static void confinement()
{
	double	rsq, f_, Uij;
	int i;
	
	for (i=0; i<N; i++) {
		DOTVP(rsq, r[i], r[i]);
		if (rsq <= R0sq) {
			continue;
		}
		f_ = -mass*OmegaSq*(1 - sqrt(R0sq/rsq));
		ADDMULVS(a[i], r[i], f_);
		Uij =  0.5*mass*OmegaSq*(R0sq + rsq - 2.0*sqrt(rsq) * R0);
		Utot += Uij;
		U[i] += Uij;
	}
}


static void directsum()
{
	double	dr[3], dr2, y6, y12, ai[3], Ui, Uij, f_, fij[3];
	int	i, j;
	
	memset(a[0], 0, 3*N*sizeof(double));
	memset(U, 0, N*sizeof(double));
	Utot=0;
	for (i=0; i<N; i++) {
		SETV(ai, a[i]);
		Ui=U[i];

		for (j=i+1; j<N; j++) {
			DOTPSUBV(dr2, dr, r[j], r[i]);
			y6 = sigma6/(dr2*dr2*dr2);
			y12 = y6*y6;
			Uij = 4.0 * epsilon * (y12 - y6);
			f_ = 4.0 * epsilon * (12 * y12 - 6 *y6) / dr2;
			MULVS(fij, dr, f_);
			ADDV(a[j], a[j], fij);
			SUBV(ai, ai, fij);
			Ui += Uij;
			U[j] += Uij;
			Utot += Uij;
		}
		SETV(a[i], ai);
		U[i] = Ui;
	}	
}


static double d3nrm2(const double r[])
{
	return sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
}

static double dr2min = 1e6;
static double single_dt(double dr2, double vsq, double ldt)
{
	double	y6, y12, f_, df2, Uij, dU, dUrel_, dx;
	
	dr2min = 1e6;
	if (dr2<dr2min) {
		dr2min = dr2;
	}
	for (dUrel_ = 1; dUrel_ > dUrel; ldt = ldt/2) {
		y6 = sigma6/(dr2*dr2*dr2);
		y12 = y6*y6;
		Uij = 4.0 * epsilon * (y12 - y6);
		f_ = 4.0 * epsilon * (6 * y6 - 12 * y12) / sqrt(dr2);
		dx = 0.5*sqrt(vsq)*ldt + 0.125*fabs(f_)/mass*ldt*ldt;
		df2 = 8.0 * f_ * epsilon * (42*y6-156*y12) / dr2;
		
		dU = 0.125*df2 * dx * ldt * ldt/mass;
		dUrel_ = fabs(dU/Uij);
	}
	return 2.0*ldt;
}

static double single_dt_c(double ldt, int i)
{
	double	f_, df2, Uij, dU, dUrel_, dx, rsq, vsq;
	
	DOTVP(rsq, r[i], r[i]);
	DOTVP(vsq, v[i], v[i]);
	
	for (dUrel_ = 1; dUrel_ > dUrel; ldt = ldt/2) {
		f_ = -mass*OmegaSq*(1 - sqrt(R0sq/rsq));
		Uij =  0.5*mass*OmegaSq*(R0sq + rsq - 2.0*sqrt(rsq) * R0);
		dx = 0.5*sqrt(vsq)*ldt + 0.125*fabs(f_)/mass*ldt*ldt;
		df2 = 2.0 * f_ * mass * OmegaSq;
		
		dU = 0.125*df2 * dx * ldt * ldt/mass;
		dUrel_ = fabs(dU/Uij);
	}
	return 2.0*ldt;
}

static double estimate_dt(double vmaxsq)
{
	double	dr[3], dr2, ldt;
	int	i, j;
	
	ldt = dt;
	for (i=0; i<N; i++) {
		for (j=i+1; j<N; j++) {
			DOTPSUBV(dr2, dr, r[j], r[i]);
			ldt = single_dt(dr2, vmaxsq, ldt);
		}
	}
	
	return ldt;
}


static void leapfrog(double trun)
{
	int	i, j, nrun;
	
	nrun = round(trun/dt);
	directsum();
	confinement();
	for (i=0; i<nrun; i++) {
		cblas_daxpy(3*N, 0.5*dt/mass, a[0], 1, v[0], 1);
		cblas_daxpy(3*N, dt, v[0], 1, r[0], 1);
		/*
		memset(a[0], 0, 3*N*sizeof(double));
		confinement();
		for (j=0; j<nfine; j++) {
			cblas_daxpy(3*N, 0.5*ddt/mass, a[0], 1, v[0], 1);
			cblas_daxpy(3*N, ddt, v[0], 1, r[0], 1);
			memset(a[0], 0, 3*N*sizeof(double));
			confinement();
			cblas_daxpy(3*N, 0.5*ddt/mass, a[0], 1, v[0], 1);
		}*/
		directsum();
		confinement();
		cblas_daxpy(3*N, 0.5*dt/mass, a[0], 1, v[0], 1);/*
		memcpy(rnew[0], a[0], 3*N*sizeof(double));
		confinement();*/
		//dump_particles(i);
		fprintf(eout, "%lg %lg %lg %lg\n", t+i*dt, energy(), Utot, Ekin);
		//memcpy(a[0], rnew[0], 3*N*sizeof(double));
	}
}


static void anneal(double kBT, double trun)
{
	int	i, imax, rsq, rsqmax=0.0;
	
	for (i=0; i<N; i++) {
		v[i][0] = gaussran(sqrt(kBT / mass), 0);
		v[i][1] = gaussran(sqrt(kBT / mass), 0);
		v[i][2] = gaussran(sqrt(kBT / mass), 0);
		DOTVP(rsq, r[i], r[i]);
		if (rsq > rsqmax) {
			rsqmax = rsq;
			imax = i;
		}
	}

	dt = estimate_dt(kBT/mass);
	ddt = M_PI/Omega/20;
	nfine = ceil(dt/ddt);
	if (nfine%2 == 1) {
		nfine++;
	}
	ddt = dt/(double) nfine;
	
	leapfrog(trun);
	//memcpy(r, rnew, 3*N*sizeof(double));
}


static void initial_config()
{
	double	rsq,sigma2, dr2, dr[3];
	int	i, j;
	
	sigma2 = sigma*sigma;
	for (i=0; i<N; i++) {
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
			mass = atof(valstring)*1836;
			required++;
			
		} else if (strcmp(valname, "tanneal") == 0) {
			tanneal = atof(valstring);
			required++;
		} else if (strcmp(valname, "dt") == 0) {
			dt = atof(valstring);
			required++;
		} else if (strcmp(valname, "tstop") == 0) {
			tstop = atof(valstring);
			required++;
		} else if (strcmp(valname, "NMC") == 0) {
			NMC = atoi(valstring);
			required++;
		} else if (strcmp(valname, "T") == 0) {
			kT = atof(valstring);
			beta = 1.0/kT;
			required++;
		} else if (strcmp(valname, "R0") == 0) {
			R0 = atof(valstring);
			R0sq = R0*R0;
			required++;
		} else if (strcmp(valname, "Omega") == 0) {
			Omega = atof(valstring);
			OmegaSq = Omega*Omega;
			required++;
		} else if (strcmp(valname, "dUrel") == 0) {
			dUrel = atof(valstring);
			required++;
		}
	}
	
	fclose(fin);
}




int main (int argc, const char * argv[])
{
	double	kTanneal;
	dUrel = 1e-3;
	read_options(argv[1]);
	eout = fopen("eout.dat", "w");


	xstep = 0.2*sigma;
	dUrel = 0.5*dUrel/(double)(N*N);
	
	r = (double (*)[3])calloc(3*N, sizeof(double));
	rnew = (double (*)[3])calloc(3*N, sizeof(double));
	rmin = (double (*)[3])calloc(3*N, sizeof(double));
	v = (double (*)[3])calloc(3*N, sizeof(double));
	a = (double (*)[3])calloc(3*N, sizeof(double));
	U = (double *)calloc(N, sizeof(double));
	
	rng = gsl_rng_alloc(gsl_rng_mt19937);

	initial_config();
	Umin = 1e100;
	
	
	//_mm_setcsr( _MM_MASK_MASK &~ (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO) );
	kTanneal = 100*epsilon;
	anneal(kTanneal, sigma/sqrt(kTanneal/mass));
	dump_particles(0);
	for (t=0; t<tstop; t+= tanneal) {
		mc_metropolis(NMC);
		dump_Umin();
		anneal(kT, tanneal);
		printf("t = %lg, Umin = %lg\n", t, Umin);
	}
	
	fclose(eout);

    return 0;
}
