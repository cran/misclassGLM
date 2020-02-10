/*
 * misclassGAMhelper.c
 *
 *  Created on: 29.11.2012
 *      Author: sdl
 */


#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>
#include <R_ext/Applic.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

struct dataset {
	double *y;
	double *x;
	double *setm;
	double *w;
	int *px;
	int *ms;
	int *pm;
	int *n;
	int *ell;
};

double _cfgaussValidation(int p, double *par, void *ex) {
	double *beta = par;

	struct dataset *ds = (struct dataset*)ex;

	int n = *(ds->n);
	double *y = ds->y;
	double *x = ds->x;
	double *setm = ds->setm;
	double *w = ds->w;
	int px = *(ds->px);
	int ms = *(ds->ms);
	int pm = *(ds->pm);

	double tmp,tmp2,tmp3,sigma2,lnsqrt2pi;
	unsigned i,j,d;

	sigma2 = -2*par[px+pm+1]*par[px+pm+1];
	lnsqrt2pi = log(sqrt(2*M_PI)*par[px+pm+1]);

	double ret = 0.0;

	for(i=0; i<n; i++)	{
		tmp = 0;

		tmp2 = y[i]-beta[0];
		for(d=1; d<=px; d++) {
			tmp2 -= beta[d] * x[(d-1)*n+i];
		}

		for(j=0; j<ms; j++) {
			tmp3 = tmp2;
			for(d=px+1; d<=px+pm; d++) {
				tmp3 -= beta[d] * setm[(d-px-1)*ms+j];
			}
			tmp += exp(tmp3*tmp3/sigma2) * w[j*n+i];
		}
		ret += log(tmp);
	}
	ret -= n*lnsqrt2pi;

	ret = -1 * ret;

	return(ret);
}

extern SEXP cfgaussValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM);

SEXP cfgaussValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM)
{
	double *y=NUMERIC_POINTER(Y);
	double *x=NUMERIC_POINTER(X);
	double *setm=NUMERIC_POINTER(SetM);
	double *w=NUMERIC_POINTER(W);
	int n=INTEGER(GET_DIM(X))[0];
	int px=INTEGER(GET_DIM(X))[1];
	int ms=INTEGER(GET_DIM(SetM))[0];
	int pm=INTEGER(GET_DIM(SetM))[1];
	double *beta=NUMERIC_POINTER(Beta);

	SEXP Ret;

	struct dataset ds;
	ds.n = &n;
	ds.y = y;
	ds.x = x;
	ds.setm = setm;
	ds.w = w;
	ds.px = &px;
	ds.ms = &ms;
	ds.pm = &pm;

	void *ex = &ds;

	PROTECT(Ret=NEW_NUMERIC(1));
	double *ret = NUMERIC_POINTER(Ret);

	*ret = _cfgaussValidation(px+pm+2, beta, ex);

	UNPROTECT(1);
	return(Ret);
}

double lpredict(double *x, double *setm, int ms, double *beta, int px, int pm, int j) {
  unsigned d;
  double ret;

  ret = beta[0];

  for(d=px+1; d<=px+pm; d++) {
    ret += beta[d] * setm[(d-px-1)*ms+j];
  }

  return(ret);
}

double lpredicti(double *y, double *x, double *beta, int px, int pm, int i, int n) {
  unsigned d;
  double ret;

  ret = y[i];
  for(d=1; d<=px; d++) {
    ret -= beta[d] * x[(d-1)*n+i];
  }

  return(ret);
}

double summ(double *exppredicts, int n, int ms, int i) {
  unsigned j;
  double ret;

  ret = 0;
  for(j=0; j<ms; j++) {
    ret += exppredicts[j*n+i];
  }

  return(ret);
}

void _cggaussValidation(int p, double *par, double *ret, void *ex) {

	double *beta = par;
	double tmp,sigma2,sigma2a,tempsumm;
	double *predicts,*exppredicts;
	unsigned i,j,k;

	struct dataset *ds = (struct dataset*)ex;

	int n = *(ds->n);
	double *y = ds->y;
	double *x = ds->x;
	double *setm = ds->setm;
	double *w = ds->w;
	int px = *(ds->px);
	int ms = *(ds->ms);
	int pm = *(ds->pm);

	predicts = (double *) calloc(n*ms, sizeof(double));
	if (NULL == predicts) {
		error("not enough memory");
	}
	exppredicts = (double *) calloc(n*ms, sizeof(double));
	if (NULL == exppredicts) {
			error("not enough memory");
		}

	for(k=0; k<px+pm+2; k++) {
		ret[k] = 0;
	}

	sigma2a = par[px+pm+1]*par[px+pm+1];
	sigma2 = -2*sigma2a;



	for(j=0; j<ms; j++) {
		tmp = lpredict(x,setm,ms,beta,px,pm,j);
		for(i=0; i<n; i++) {
			predicts[j*n+i] = lpredicti(y,x,beta,px,pm,i,n)-tmp;
			exppredicts[j*n+i] = exp(predicts[j*n+i]*predicts[j*n+i]/sigma2) * w[j*n+i];
		}
	}

	for(i=0; i<n; i++) {
		tempsumm = summ(exppredicts,n,ms,i);
		tmp = 0;
		for(j=0; j<ms; j++) {
			tmp += predicts[j*n+i] * exppredicts[j*n+i];
		}
		tmp /= tempsumm;
		ret[0] += tmp;
		for(k=1; k<=px; k++) {
			ret[k] += tmp * x[(k-1)*n+i];
		}
		for(k=px+1; k<=px+pm; k++) {
			tmp = 0;
			for(j=0; j<ms; j++) {
				tmp += predicts[j*n+i] * exppredicts[j*n+i] * setm[(k-px-1)*ms+j];
			}
			ret[k] += tmp/tempsumm;
		}
		tmp = 0;
		for(j=0; j<ms; j++) {
			tmp += (predicts[j*n+i]*predicts[j*n+i]/sigma2a-1) * exppredicts[j*n+i];
		}
		ret[px+pm+1] += tmp/tempsumm;
	}

	for(k=0; k<px+pm+1; k++) {
		ret[k] /= sigma2a;
	}
	ret[px+pm+1] /= par[px+pm+1];

	for(k=0; k<px+pm+2; k++) {
		ret[k] = -1 * ret[k];
	}

	free(exppredicts);
	free(predicts);
}


extern SEXP cggaussValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM);

SEXP cggaussValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM)
{
	double *y=NUMERIC_POINTER(Y);
	double *x=NUMERIC_POINTER(X);
	double *setm=NUMERIC_POINTER(SetM);
	double *w=NUMERIC_POINTER(W);
	int n=INTEGER(GET_DIM(X))[0];
	int px=INTEGER(GET_DIM(X))[1];
	int ms=INTEGER(GET_DIM(SetM))[0];
	int pm=INTEGER(GET_DIM(SetM))[1];
	double *beta=NUMERIC_POINTER(Beta);

	SEXP Ret;
	PROTECT(Ret=NEW_NUMERIC(px+pm+2));

	double *ret = NUMERIC_POINTER(Ret);

	struct dataset ds;
	ds.n = &n;
	ds.y = y;
	ds.x = x;
	ds.setm = setm;
	ds.w = w;
	ds.px = &px;
	ds.ms = &ms;
	ds.pm = &pm;

	void *ex = &ds;

	_cggaussValidation(px+pm+2, beta, ret, ex);

	UNPROTECT(1);

	return(Ret);
}

double _cflogitValidation(int p, double *par, void *ex) {
	double tmp,tmp2,tmp3,ret;
	unsigned i,j,d;

	struct dataset *ds = (struct dataset*)ex;

	int n = *(ds->n);
	double *y = ds->y;
	double *x = ds->x;
	double *setm = ds->setm;
	double *w = ds->w;
	int px = *(ds->px);
	int ms = *(ds->ms);
	int pm = *(ds->pm);
	double *beta = par;

	ret = 0;

	for(i=0; i<n; i++)	{
		tmp = 0;

		tmp2 = beta[0];
		for(d=1; d<=px; d++) {
			tmp2 += beta[d] * x[(d-1)*n+i];
		}

		for(j=0; j<ms; j++) {
			tmp3 = tmp2;
			for(d=px+1; d<=px+pm; d++) {
				tmp3 += beta[d] * setm[(d-px-1)*ms+j];
			}

			tmp3 = exp(tmp3);
			tmp3 = tmp3/(1+tmp3);
			if (y[i]) {
				tmp += tmp3 * w[j*n+i];
			} else {
				tmp += (1-tmp3) * w[j*n+i];
			}
		}
		ret += log(tmp);
	}

	ret = -1 * ret;

	return(ret);
}

extern SEXP cflogitValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM);

SEXP cflogitValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM)
{
	double *y=NUMERIC_POINTER(Y);
	double *x=NUMERIC_POINTER(X);
	double *setm=NUMERIC_POINTER(SetM);
	double *w=NUMERIC_POINTER(W);
	int n=INTEGER(GET_DIM(X))[0];
	int px=INTEGER(GET_DIM(X))[1];
	int ms=INTEGER(GET_DIM(SetM))[0];
	int pm=INTEGER(GET_DIM(SetM))[1];
	double *beta=NUMERIC_POINTER(Beta);

	struct dataset ds;
	ds.n = &n;
	ds.y = y;
	ds.x = x;
	ds.setm = setm;
	ds.w = w;
	ds.px = &px;
	ds.ms = &ms;
	ds.pm = &pm;

	void *ex = &ds;

	SEXP Ret;
	PROTECT(Ret=NEW_NUMERIC(1));

	double *ret = NUMERIC_POINTER(Ret);

	*ret = _cflogitValidation(px+pm+1, beta, ex);

	UNPROTECT(1);

	return(Ret);
}

void _cglogitValidation(int p, double *par, double *ret, void *ex) {
	double tmp,tmp2,tmp3,lik;
	double *eta,*tmp4;
	unsigned i,j,d;

	struct dataset *ds = (struct dataset*)ex;

	int n = *(ds->n);
	double *y = ds->y;
	double *x = ds->x;
	double *setm = ds->setm;
	double *w = ds->w;
	int px = *(ds->px);
	int ms = *(ds->ms);
	int pm = *(ds->pm);

	tmp4 = (double *) calloc(p, sizeof(double));
	if (NULL == tmp4) {
		error("not enough memory");
	}
	eta = (double *) calloc(ms, sizeof(double));
	if (NULL == eta) {
		error("not enough memory");
	}

	for(j=0; j<p; j++) {
		ret[j] = 0;
	}

	for(i=0; i<n; i++)	{
		lik = 0;
		tmp2 = par[0];
		for(d=1; d<=px; d++) {
			tmp2 += par[d] * x[(d-1)*n+i];
		}

		for(j=0; j<ms; j++) {
			tmp3 = tmp2;
			for(d=px+1; d<=px+pm; d++) {
				tmp3 += par[d] * setm[(d-px-1)*ms+j];
			}

			eta[j] = exp(tmp3);
			tmp3 = eta[j]/(1+eta[j]);

			if (y[i]) {
				lik += tmp3 * w[j*n+i];
			} else {
				lik += (1-tmp3) * w[j*n+i];
			}
		}

		for(d=0; d<p; d++) {
			tmp4[d]=0;
		}

		for(j=0; j<ms; j++) {
			if (y[i]) {
				tmp = eta[j] / ((1+eta[j])*(1+eta[j])) * w[j*n+i];
			} else {
				tmp = -eta[j] / ((1+eta[j])*(1+eta[j])) * w[j*n+i];
			}

			tmp4[0] += tmp;
			for(d=0; d<px; d++) {
				tmp4[d+1] += x[d*n+i] * tmp;
			}
			for(d=0; d<pm; d++) {
				if(setm[d*ms+j]) tmp4[d+px+1] += setm[d*ms+j] * tmp;
			}
		}

		for(d=0; d<p; d++) {
			ret[d] += tmp4[d]/lik;
		}
	}

	for(i=0; i<p; i++) {
		ret[i] = -1 * ret[i];
	}

	free(tmp4);
	free(eta);
}

extern SEXP cglogitValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM);

SEXP cglogitValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM)
{
	double *y=NUMERIC_POINTER(Y);
	double *x=NUMERIC_POINTER(X);
	double *setm=NUMERIC_POINTER(SetM);
	double *w=NUMERIC_POINTER(W);
	int n=INTEGER(GET_DIM(X))[0];
	int px=INTEGER(GET_DIM(X))[1];
	int ms=INTEGER(GET_DIM(SetM))[0];
	int pm=INTEGER(GET_DIM(SetM))[1];
	double *beta=NUMERIC_POINTER(Beta);

	struct dataset ds;
	ds.n = &n;
	ds.y = y;
	ds.x = x;
	ds.setm = setm;
	ds.w = w;
	ds.px = &px;
	ds.ms = &ms;
	ds.pm = &pm;

	void *ex = &ds;

	SEXP Ret;
	PROTECT(Ret=NEW_NUMERIC(px+pm+1));

	double *ret = NUMERIC_POINTER(Ret);

	_cglogitValidation(px+pm+1, beta, ret, ex);

	UNPROTECT(1);

	return(Ret);
}



double _cfmlogitValidation(int p, double *par, void *ex) {
  double tmp,tmph;
  double *eta,*etah;
  unsigned i,j,d,l;

  struct dataset *ds = (struct dataset*)ex;

  int n = *(ds->n);
  double *y = ds->y;
  double *x = ds->x;
  double *setm = ds->setm;
  double *w = ds->w;
  int px = *(ds->px);
  int ms = *(ds->ms);
  int pm = *(ds->pm);
  int ell = *(ds->ell);

  double *beta = par;
  p = px+pm+1;

  eta = (double *) calloc(ell, sizeof(double));
  if (NULL == eta) {
    error("not enough memory");
  }
  etah = (double *) calloc(ell, sizeof(double));
  if (NULL == etah) {
    error("not enough memory");
  }
  double ret;

  ret = 0;

  for(i=0; i<n; i++) {
    tmp = 0;

    for (l=0; l<ell; l++) {
      etah[l] = beta[l*p];
      for(d=1; d<=px; d++) {
        etah[l] += beta[l*p+d] * x[(d-1)*n+i];
      }
    }

    for(j=0; j<ms; j++) {
      tmph = 1;
      for (l=0; l<ell; l++) {
        eta[l] = etah[l];
        for(d=px+1; d<=px+pm; d++) {
          eta[l] += beta[l*p+d] * setm[(d-px-1)*ms+j];
        }
        eta[l] = exp(eta[l]);
        tmph += eta[l];
      }
      for (l=0; l<ell; l++) {
        if(y[l*n+i]) tmp += eta[l]/tmph * w[j*n+i];
      }
      if(y[ell*n+i]) tmp += w[j*n+i]/tmph;
    }
    ret += log(tmp);
  }

  ret = -1 * ret;

  free(eta);
  free(etah);

  return(ret);
}

extern SEXP cfmlogitValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM);

SEXP cfmlogitValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM)
{
  double *y=NUMERIC_POINTER(Y);
  double *x=NUMERIC_POINTER(X);
  double *setm=NUMERIC_POINTER(SetM);
  double *w=NUMERIC_POINTER(W);
  int ell=INTEGER(GET_DIM(Y))[1] - 1;
  int n=INTEGER(GET_DIM(X))[0];
  int px=INTEGER(GET_DIM(X))[1];
  int ms=INTEGER(GET_DIM(SetM))[0];
  int pm=INTEGER(GET_DIM(SetM))[1];
  double *beta=NUMERIC_POINTER(Beta);

  struct dataset ds;
  ds.n = &n;
  ds.y = y;
  ds.x = x;
  ds.setm = setm;
  ds.w = w;
  ds.px = &px;
  ds.ms = &ms;
  ds.pm = &pm;
  ds.ell = &ell;

  void *ex = &ds;

  SEXP Ret;
  PROTECT(Ret=NEW_NUMERIC(1));

  double *ret = NUMERIC_POINTER(Ret);
  *ret = _cfmlogitValidation((px+pm+1)*ell, beta, ex);

  UNPROTECT(1);
  return(Ret);
}

void _cgmlogitValidation(int p, double *par, double *ret, void *ex) {
  double lik,tmp3;
  double *eta,*etah,*tmp2,*tmp,*etah2;
  unsigned i,j,k,d,l,o;

  double *beta = par;

  struct dataset *ds = (struct dataset*)ex;

  int n = *(ds->n);
  double *y = ds->y;
  double *x = ds->x;
  double *setm = ds->setm;
  double *w = ds->w;
  int px = *(ds->px);
  int ms = *(ds->ms);
  int pm = *(ds->pm);
  int ell = *(ds->ell);

  p = px + pm + 1;
  eta = (double *) calloc(ell*ms, sizeof(double));
  if (NULL == eta) error("not enough memory");
  etah = (double *) calloc(ell, sizeof(double));
  if (NULL == etah) error("not enough memory");
  etah2 = (double *) calloc(ell*ms, sizeof(double));
  if (NULL == etah2) error("not enough memory");
  tmp = (double *) calloc(ms, sizeof(double));
  if (NULL == tmp) error("not enough memory");
  tmp2 = (double *) calloc(ell*p, sizeof(double));
  if (NULL == tmp2) error("not enough memory");

  for(j=0; j<ell*p; j++) {
    ret[j] = 0;
  }

  for(j=0; j<ms; j++) {
    for(d=px+1; d<=px+pm; d++) {
      if (setm[(d-px-1)*ms+j]) {
        for (l=0; l<ell; l++) {
          etah2[l*ms+j] += beta[l*p+d] * setm[(d-px-1)*ms+j];
        }
      }
    }
  }

  for(i=0; i<n; i++) {
    lik = 0;

    for (l=0; l<ell; l++) {
      etah[l] = beta[l*p];
      for(d=1; d<=px; d++) {
        etah[l] += beta[l*p+d] * x[(d-1)*n+i];
      }
    }

    for(j=0; j<ms; j++) {
      tmp[j] = 1;
      for (l=0; l<ell; l++) {
        eta[l*ms+j] = exp(etah[l]+etah2[l*ms+j]);
        tmp[j] += eta[l*ms+j];
      }
      for (l=0; l<ell; l++) {
        if(y[l*n+i]) lik += eta[l*ms+j]/tmp[j] * w[j*n+i];
      }
      if(y[ell*n+i]) lik += w[j*n+i]/tmp[j];
    }

    for(j=0; j<ell*p; j++) {
      tmp2[j]=0;
    }

    for(o=0; o<=ell; o++) {
      if (y[o*n+i]) {
        for(l=0; l<ell; l++) {
          for (k=0; k<ms; k++) {
            if (o==l) {
              tmp3 = ((eta[l*ms+k] * (tmp[k] - eta[l*ms+k])) / (tmp[k] * tmp[k])) * w[k*n+i];
            } else {
              if (o<ell) {
                tmp3 = ((-eta[l*ms+k] * eta[o*ms+k]) / (tmp[k] * tmp[k])) * w[k*n+i];
              } else {
                tmp3 = (-eta[l*ms+k] / (tmp[k] * tmp[k])) * w[k*n+i];
              }
            }
            tmp2[l*p] += tmp3;
            for (j=0; j<px; j++) {
              tmp2[l*p+1+j] += x[j*n+i] * tmp3;
            }
            for (j=0; j<pm; j++) {
              if (setm[j*ms+k]) tmp2[l*p+px+1+j] += setm[j*ms+k] * tmp3;
            }
          }
        }
      }
    }
    for(j=0; j<ell*p; j++) {
      ret[j] += tmp2[j]/lik;
    }
  }

  for(i=0; i<ell*p; i++) {
    ret[i] = -1 * ret[i];
  }

  free(eta);
  free(etah);
  free(etah2);
  free(tmp);
  free(tmp2);
}

extern SEXP cgmlogitValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM);

SEXP cgmlogitValidation(SEXP Beta, SEXP Y, SEXP X, SEXP W, SEXP SetM)
{
  double *y=NUMERIC_POINTER(Y);
  double *x=NUMERIC_POINTER(X);
  double *setm=NUMERIC_POINTER(SetM);
  double *w=NUMERIC_POINTER(W);
  int ell=INTEGER(GET_DIM(Y))[1] - 1;
  int n=INTEGER(GET_DIM(X))[0];
  int px=INTEGER(GET_DIM(X))[1];
  int ms=INTEGER(GET_DIM(SetM))[0];
  int pm=INTEGER(GET_DIM(SetM))[1];
  double *beta=NUMERIC_POINTER(Beta);

  int p = px + pm + 1;

  SEXP Ret;
  PROTECT(Ret=NEW_NUMERIC(ell*p));

  double *ret = NUMERIC_POINTER(Ret);

  struct dataset ds;
  ds.n = &n;
  ds.y = y;
  ds.x = x;
  ds.setm = setm;
  ds.w = w;
  ds.px = &px;
  ds.ms = &ms;
  ds.pm = &pm;
  ds.ell = &ell;

  void *ex = &ds;

  _cgmlogitValidation(p*ell, beta, ret, ex);

  UNPROTECT(1);

  return(Ret);
}

extern SEXP get_link_mlogit(SEXP Link, SEXP Link_d_eta, SEXP Link_d2_eta, SEXP Y, SEXP Eta, SEXP Theta);

SEXP get_link_mlogit(SEXP Link, SEXP Link_d_eta, SEXP Link_d2_eta, SEXP Y, SEXP Eta, SEXP Theta)
{
  unsigned j, l, idx, k1;
  double tmp, tmp2, linktmp, *expeta;

  double *link = NUMERIC_POINTER(Link);
  double *link_d_eta = NUMERIC_POINTER(Link_d_eta);
  double *link_d2_eta = NUMERIC_POINTER(Link_d2_eta);
  double *y = NUMERIC_POINTER(Y);
  unsigned k = GET_LENGTH(Y);
  k1 = k-1;

  double *eta = NUMERIC_POINTER(Eta);

  expeta = (double *) calloc(k, sizeof(double));
  if (NULL == expeta) error("not enough memory");

  tmp = 1;

  idx = k1;
  for (j=0; j<k1; j++) {
    expeta[j] = exp(eta[j]);
    tmp += expeta[j];
    if (y[j]>0.5) idx = j;
  }

  if (idx==k1) {
    *link = 1;
  } else {
    *link = expeta[idx];
  }
  *link /= tmp;

  linktmp = *link / tmp;

  for (j=0; j<k1; j++) {
    link_d_eta[j] = -expeta[j] * linktmp;
  }
  if (idx<k1) link_d_eta[idx] += expeta[idx] / tmp;

  for (j=0; j<k1; j++) {
    link_d2_eta[j * k1 + j] = 2/tmp * expeta[j] * expeta[j] - expeta[j];
    link_d2_eta[j * k1 + j] *= linktmp;
  }
  for (j=0; j<k1; j++) {
    for (l=0; l<j; l++) {
      tmp2 = 2/tmp * expeta[j] * expeta[l] * linktmp;
      link_d2_eta[j * k1 + l] = tmp2;
      link_d2_eta[l * k1 + j] = tmp2;
    }
  }

  if (idx<k1) {
    tmp2 = tmp * tmp;
    link_d2_eta[idx * k1 + idx] += expeta[idx] / tmp;

    for (j=0; j<k1; j++) {
      tmp = expeta[idx] * expeta[j] / tmp2;
      link_d2_eta[idx * k1 + j] -= tmp;
      link_d2_eta[j * k1 + idx] -= tmp;
    }
  }

  free(expeta);

  return(R_NilValue);
}

extern SEXP get_Gdeta(SEXP G_d_eta, SEXP G_d_eta_2, SEXP Eta_u, SEXP Y, SEXP W, SEXP link_derivs, SEXP Theta, SEXP K);

SEXP get_Gdeta(SEXP G_d_eta, SEXP G_d_eta_2, SEXP Eta_u, SEXP Y, SEXP W, SEXP link_derivs, SEXP Theta, SEXP K)
{
  unsigned lk, i, u, h, g, r, n, uni, nlk;
  double *tmp, *tmp2, *tmp_h, *tmp2_h, *link_d2_eta_h, nenner, temp, nenner2;

  SEXP Link, Link_d_eta, Link_d2_eta, R_fcall, Y_tmp, Eta_tmp;

  double *gdeta = NUMERIC_POINTER(G_d_eta);
  double *gdeta2 = NUMERIC_POINTER(G_d_eta_2);
  double *etau = NUMERIC_POINTER(Eta_u);
  double *y = NUMERIC_POINTER(Y);
  double *w = NUMERIC_POINTER(W);
  int l = INTEGER(GET_DIM(Eta_u))[0];
  int k = INTEGER_VALUE(K);

  if (k>1) {
    n = INTEGER(GET_DIM(Y))[0];
  } else {
    n = GET_LENGTH(Y);
  }

  if(!isFunction(link_derivs)) error("link_derivs must be a function");

  lk = l * k;
  nlk = n * lk;

  tmp = (double *) calloc(lk, sizeof(double));
  if (NULL == tmp) error("not enough memory");
  tmp2 = (double *) calloc(lk * lk, sizeof(double));
  if (NULL == tmp2) error("not enough memory");

  PROTECT(Link = NEW_NUMERIC(1));
  PROTECT(Link_d_eta = NEW_NUMERIC(k));
  PROTECT(Link_d2_eta = allocMatrix(REALSXP, k, k));
  if (k>1) {
    PROTECT(Y_tmp = NEW_NUMERIC(k + 1));
  } else {
    PROTECT(Y_tmp = NEW_NUMERIC(k));
  }
  PROTECT(Eta_tmp = NEW_NUMERIC(k));

  double *y_tmp = NUMERIC_POINTER(Y_tmp);
  double *eta_tmp = NUMERIC_POINTER(Eta_tmp);
  double *link = NUMERIC_POINTER(Link);
  double *link_d_eta = NUMERIC_POINTER(Link_d_eta);
  double *link_d2_eta = NUMERIC_POINTER(Link_d2_eta);

  PROTECT(R_fcall = allocVector(LANGSXP,7));
  SETCAR(R_fcall, link_derivs);

  for(i=0; i<n; i++) {
    nenner = 0.0;

    for (h=0; h<k; h++) {
      y_tmp[h] = y[h * n + i];
    }
    if (k>1) y_tmp[k] = y[k * n + i];

    for(u=0; u<l; u++) {
      uni = u * n + i;
      SETCADR(R_fcall, Link);
      SETCADDR(R_fcall, Link_d_eta);
      SETCADDDR(R_fcall, Link_d2_eta);
      SETCAD4R(R_fcall, Y_tmp);

      for (h=0; h<k; h++) {
        eta_tmp[h] = etau[h * n * l + i * l + u];
      }

      SETCAD4R(CDR(R_fcall), Eta_tmp);
      SETCAD4R(CDDR(R_fcall), Theta);

      EVAL(R_fcall);

      nenner += *link * w[uni];

      r = u * k;
      tmp_h = &tmp[r];
      for(g=0; g<k; g++) {
        tmp_h[g] = link_d_eta[g] * w[uni];
        tmp2_h = &tmp2[(r + g) * lk + r];
        link_d2_eta_h = &link_d2_eta[g * k];
        for(h=0; h<k; h++) {
          tmp2_h[h] = link_d2_eta_h[h] * w[uni];
        }
      }
    }

    nenner2 = nenner * nenner;

    for (h=0; h<lk; h++) {
      gdeta[h * n + i] = tmp[h] / nenner;
      gdeta2[h * nlk + h * n + i] = tmp2[h * lk + h] / nenner - (tmp[h] * tmp[h] / nenner2);
    }

    for (h=0; h<lk; h++) {
      tmp2_h = &tmp2[h * lk];
      for (g=0; g<h; g++) {
        temp = tmp2_h[g] / nenner - (tmp[g] * tmp[h] / nenner2); // TODO: stimmt das so?
        gdeta2[h * nlk + g * n + i] = temp;
        gdeta2[g * nlk + h * n + i] = temp;
      }
    }
  }

  UNPROTECT(6);

  free(tmp);
  free(tmp2);

  return(R_NilValue);
}


static const R_CallMethodDef callMethods[] = {
  {"cfgaussValidation", (DL_FUNC) &cfgaussValidation, 5},
  {"cggaussValidation", (DL_FUNC) &cggaussValidation, 5},
  {"cflogitValidation", (DL_FUNC) &cflogitValidation, 5},
  {"cglogitValidation", (DL_FUNC) &cglogitValidation, 5},
  {"cfmlogitValidation", (DL_FUNC) &cfmlogitValidation, 5},
  {"cgmlogitValidation", (DL_FUNC) &cgmlogitValidation, 5},
  {"get_link_mlogit", (DL_FUNC)  &get_link_mlogit, 6},
  {"get_Gdeta", (DL_FUNC) &get_Gdeta, 8},
  {NULL, NULL, 0}
};

void R_init_misclassGLM(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

