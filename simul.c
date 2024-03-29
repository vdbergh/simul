#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>
#include <inttypes.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
#define MAX_THREADS 128
pthread_t threads[MAX_THREADS];

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define NPROC_COMMAND "nproc"
int nproc(void) {
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  int nproc_;
  FILE *f = popen(NPROC_COMMAND, "r");
  if (f == NULL) {
    return 1;
  }
  read = getline(&line, &len, f);
  pclose(f);
  nproc_ = read > 0 ? atoi(line) : 1;
  free(line);
  return nproc_;
}

#define phi(x) gsl_sf_erf_Z(x)
#define Phi(x) (1 - gsl_sf_erf_Q(x))

typedef struct {
  int funcalls;
  int iterations;
  int error_num;
} stats_t;

#define SIGNERR -1
#define CONVERR -2

typedef double (*callback_type)(double, void *);
// extern double brentq(callback_type f, double xa, double xb, double xtol,
// double rtol, int iter, stats_t *stats, void *args);

typedef struct {
  callback_type f;
  int funcalls;
  void *args;
} callback_wrapper_t;

double callback_wrapper(double x, void *args) {
  callback_wrapper_t *c = (callback_wrapper_t *)(args);
  c->funcalls++;
  return (*(c->f))(x, c->args);
}

double brentq(callback_type f, double xa, double xb, double xtol, double rtol,
              int iter, stats_t *stats, void *args) {
  int status;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  gsl_function F;
  callback_wrapper_t c;
  c.f = f;
  c.funcalls = 0;
  c.args = args;
  stats->error_num = CONVERR;
  stats->iterations = 0;
  F.function = callback_wrapper;
  F.params = (void *)(&c);
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_set_error_handler_off();
  status = gsl_root_fsolver_set(s, &F, xa, xb);
  if (status == GSL_EINVAL) {
    stats->error_num = SIGNERR;
  } else {
    assert(status == GSL_SUCCESS);
    status = GSL_CONTINUE;
  }
  for (int i = 0; i < iter && status == GSL_CONTINUE; i++) {
    stats->iterations++;
    status = gsl_root_fsolver_iterate(s);
    assert(status == GSL_SUCCESS);
    r = gsl_root_fsolver_root(s);
    xa = gsl_root_fsolver_x_lower(s);
    xb = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(xa, xb, xtol, rtol);
    if (status == GSL_SUCCESS) {
      stats->error_num = 0;
      break;
    }
    assert(status == GSL_CONTINUE);
  }
  stats->funcalls = c.funcalls;
  gsl_root_fsolver_free(s);
  return r;
}

#define N 5
#define H0 0
#define H1 1
#define CONTINUE 2

/*
    Probality distributions (generally having a name starting with
    "pdf") are represented by an array[2*N] consisting of pairs
    (ai,pi), i=1,...,N.  It is usually assumed that the ai are strictly
    ascending and p1>0, pN>0.
*/

double L_(double x) { return 1 / (1 + pow(10, -x / 400.0)); }

double Linv(double s) { return -400 * log10(1 / s - 1); }

void disp(double pdf[]) {
  int i;
  printf("[");
  for (i = 0; i < N - 1; i++) {
    printf("(%f,%f), ", pdf[2 * i], pdf[2 * i + 1]);
  }
  i = N - 1;
  printf("(%f,%f)]\n", pdf[2 * i], pdf[2 * i + 1]);
}

void muvar(double pdf_in[], double *mu, double *var) {
  int i;
  double a, p;
  double epsilon = 1e-6;
  double sum = 0.0;
  double sum2 = 0.0;
  *mu = 0.0;
  for (i = 0; i < N; i++) {
    a = pdf_in[2 * i];
    p = pdf_in[2 * i + 1];
    assert(-epsilon <= p);
    assert(p <= 1 + epsilon);
    sum += p;
    *mu += a * p;
    sum2 += a * a * p;
  }
  assert(fabs(sum - 1) < epsilon);
  *var = sum2 - (*mu) * (*mu);
}

double ucp(double mu, double var, int batch, double delta) {
  double mu_ = batch * mu;
  double var_ = batch * var;
  double sigma_ = sqrt(var_);
  double delta_z = (delta - mu_) / sigma_;
  double o, p;
  if (delta <= 0) {
    p = 0.0;
  } else if (delta_z >= 5) {
    p = 1.0;
  } else if (delta_z <= -5) {
    p = delta / mu_;
  } else {
    o = sigma_ * (phi(delta_z) - delta_z * (1 - Phi(delta_z)));
    p = delta / (o + delta);
  }
  // assert(p>=0 && p<=1);
  return p;
}

double f(double x, void *args) {
  int i;
  double a, p;
  double sum = 0.0;
  double *pdf_in = (double *)(args);
  for (i = 0; i < N; i++) {
    a = pdf_in[2 * i];
    p = pdf_in[2 * i + 1];
    sum += p * a / (1 + x * a);
  }
  return sum;
}

double secular(double pdf_in[]) {
  stats_t stats = {0, 0, 0};
  int i;
  double t, v, w, l, u, epsilon, xtol, rtol, x;
  v = pdf_in[0];
  w = pdf_in[0];
  for (i = 2; i < 2 * N; i += 2) {
    t = pdf_in[i];
    if (t < v) {
      v = t;
    }
    if (t > w) {
      w = t;
    }
  }
  assert(v * w < 0);
  l = -1 / w;
  u = -1 / v;
  epsilon = 1e-9;
  xtol = 2e-12; // scipy defaults
  rtol = 8.881784197001252e-16;
  x = brentq(f, l + epsilon, u - epsilon, xtol, rtol, 1000, &stats,
             (void *)pdf_in);
  assert(stats.error_num == 0);
  return x;
}

void MLE_expected(double pdf_in[], double s, double pdf_out[]) {
  /*
    This function computes the maximum likelood estimate for a
    discrete distribution with expectation value s, given an observed
    (i.e. empirical) distribution pdf_in.

    The theory behind this function can be found in the online
    document

    http://hardy.uhasselt.be/Fishtest/support_MLE_multinomial.pdf

    (see Proposition 1.1).
  */
  double x, p, a, mu, var;
  int i;
  double pdf1[2 * N];
  for (i = 0; i < N; i++) {
    a = pdf_in[2 * i];
    p = pdf_in[2 * i + 1];
    pdf1[2 * i] = a - s;
    pdf1[2 * i + 1] = p;
  }
  x = secular(pdf1);
  for (i = 0; i < N; i++) {
    a = pdf_in[2 * i];
    p = pdf_in[2 * i + 1];
    pdf_out[2 * i] = a;
    pdf_out[2 * i + 1] = p / (1 + x * (a - s));
  }
  muvar(pdf_out, &mu, &var); /* for validation */
  assert(fabs(s - mu) < 1e-6);
}

void uniform(double pdf_in[], double pdf_out[]) {
  double Ninv = 1 / ((double)N);
  int i;
  double a;
  for (i = 0; i < N; i++) {
    a = pdf_in[2 * i];
    pdf_out[2 * i] = a;
    pdf_out[2 * i + 1] = Ninv;
  }
}

void MLE_t_value(double pdf_in[], double ref, double s, double pdf_out[]) {
  /*
    See https://hardy.uhasselt.be/Fishtest/normalized_elo_practical.pdf
    Section 4.1
  */
  double pdf_[2 * N];
  double mu, var, a, p, x, m, d;
  int i, j;
  uniform(pdf_in, pdf_out);
  for (i = 0; i < 10; i++) {
    memcpy(pdf_, pdf_out, 2 * N * sizeof(double));
    muvar(pdf_out, &mu, &var);
    double sigma = sqrt(var);
    double pdf1[2 * N];
    for (j = 0; j < N; j++) {
      a = pdf_in[2 * j];
      p = pdf_in[2 * j + 1];
      pdf1[2 * j] = a - ref - s * sigma * (1 + pow((mu - a) / sigma, 2)) / 2;
      pdf1[2 * j + 1] = p;
    }
    x = secular(pdf1);
    for (j = 0; j < N; j++) {
      pdf_out[2 * j + 1] = pdf_in[2 * j + 1] / (1 + x * pdf1[2 * j]);
    }
    m = 0;
    for (j = 0; j < N; j++) {
      d = fabs(pdf_[2 * j + 1] - pdf_out[2 * j + 1]);
      if (d > m) {
        m = d;
      }
    }
    if (m < 1e-9) {
      break;
    }
  }
  muvar(pdf_out, &mu, &var);
  assert(fabs(s - (mu - ref) / sqrt(var)) < 1e-5);
}

double myrand(uint64_t *prng) {
  /*
    https://nuclear.llnl.gov/CNP/rng/rngman/node4.html
  */
  uint64_t a = UINT64_C(2862933555777941757);
  uint64_t b = UINT64_C(3037000493);
  uint64_t current = *prng;
  *prng = a * (*prng) + b;
  return current / pow(2, 64);
}

void jump(uint64_t *prng) {
  /* do 2^48 steps */
  uint64_t a = UINT64_C(3311271626024157185);
  uint64_t b = UINT64_C(8774982398954700800);
  *prng = a * (*prng) + b;
}

double pick(uint64_t *prng, double pdf[]) {
  double x = myrand(prng);
  int i;
  double s = 0.0;
  for (i = 0; i < N; i++) {
    s += pdf[2 * i + 1];
    if (x < s) {
      return pdf[2 * i];
    }
  }
  /* we should not get here */
  return pdf[2 * (N - 1)];
}

double LLR_expected(double pdf_in[], double s0, double s1) {
  double pdf0[2 * N], pdf1[2 * N];
  double p, p0, p1;
  double sum = 0.0;
  int i;
  MLE_expected(pdf_in, s0, pdf0);
  MLE_expected(pdf_in, s1, pdf1);
  for (i = 0; i < N; i++) {
    p = pdf_in[2 * i + 1];
    p0 = pdf0[2 * i + 1];
    p1 = pdf1[2 * i + 1];
    sum += p * log(p1 / p0);
  }
  return sum;
}

void LLRjumps_expected(double pdf_in[], double s0, double s1,
                       double pdf_out[]) {
  double pdf0[2 * N], pdf1[2 * N];
  double p, p0, p1;
  int i;
  MLE_expected(pdf_in, s0, pdf0);
  MLE_expected(pdf_in, s1, pdf1);
  for (i = 0; i < N; i++) {
    p = pdf_in[2 * i + 1];
    p0 = pdf0[2 * i + 1];
    p1 = pdf1[2 * i + 1];
    pdf_out[2 * i] = log(p1 / p0);
    pdf_out[2 * i + 1] = p;
  }
}

double LLR_t_value(double pdf_in[], double ref, double s0, double s1) {
  double pdf0[2 * N], pdf1[2 * N];
  double p, p0, p1;
  double sum = 0.0;
  int i;
  MLE_t_value(pdf_in, ref, s0, pdf0);
  MLE_t_value(pdf_in, ref, s1, pdf1);
  for (i = 0; i < N; i++) {
    p = pdf_in[2 * i + 1];
    p0 = pdf0[2 * i + 1];
    p1 = pdf1[2 * i + 1];
    sum += p * log(p1 / p0);
  }
  return sum;
}

void LLRjumps_t_value(double pdf_in[], double ref, double s0, double s1,
                      double pdf_out[]) {
  double pdf0[2 * N], pdf1[2 * N];
  double p, p0, p1;
  int i;
  MLE_t_value(pdf_in, ref, s0, pdf0);
  MLE_t_value(pdf_in, ref, s1, pdf1);
  for (i = 0; i < N; i++) {
    p = pdf_in[2 * i + 1];
    p0 = pdf0[2 * i + 1];
    p1 = pdf1[2 * i + 1];
    pdf_out[2 * i] = log(p1 / p0);
    pdf_out[2 * i + 1] = p;
  }
}

double LLR_alt(double pdf_in[], double s0, double s1) {
  /*
    This function computes the approximate generalized log likelihood ratio
    (divided by N) for s=s1 versus s=s0 where pdf is an empirical distribution
    and s is the expectation value of the true distribution.

    http://hardy.uhasselt.be/Fishtest/support_MLE_multinomial.pdf
  */
  int i;
  double p, v, r0 = 0.0, r1 = 0.0;
  for (i = 0; i < N; i++) {
    p = pdf_in[2 * i + 1];
    v = pdf_in[2 * i];
    r0 += p * (v - s0) * (v - s0);
    r1 += p * (v - s1) * (v - s1);
  }
  return 1 / 2.0 * log(r0 / r1);
}

void regularize(int results_in[], double results_out[]) {
  /*
     Replace zeros with a small value to avoid division by
     zero issues.
  */
  double epsilon = 1e-4;
  int i;
  for (i = 0; i < N; i++) {
    if (results_in[i] == 0) {
      results_out[i] = epsilon;
    } else {
      results_out[i] = (double)(results_in[i]);
    }
  }
}

void results_to_pdf(int results_in[], double *count, double pdf_out[]) {
  double results_out[N];
  int i;
  *count = 0.0;
  regularize(results_in, results_out);
  for (i = 0; i < N; i++) {
    *count += results_out[i];
  }
  for (i = 0; i < N; i++) {
    pdf_out[2 * i] = i / (N - 1.0);
    pdf_out[2 * i + 1] = results_out[i] / (*count);
  }
}

double LLR_logistic(double s0, double s1, int results_in[]) {
  /*
    This function computes the generalized log-likelihood ratio for
    "results_in" which should be an array of length 5 containing the
    frequencies of the game pairs LL,LD+DL,LW+DD+WL,DW+WD,WW.
  */
  double pdf_out[2 * N];
  double count;
  results_to_pdf(results_in, &count, pdf_out);
  return count * LLR_expected(pdf_out, s0, s1);
}

double LLR_normalized(double nt0, double nt1, int results_in[]) {
  /*
    Like LLR_logistic but using normalized t-values.
    See Section 4.1 in
    http://hardy.uhasselt.be/Fishtest/normalized_elo_practical.pdf
  */
  double pdf_out[2 * N];
  double count;
  double t0, t1;
  double sqrt2 = sqrt(2);
  if (N == 3) {
    t0 = nt0;
    t1 = nt1;
  } else if (N == 5) {
    t0 = nt0 * sqrt2;
    t1 = nt1 * sqrt2;
  } else {
    assert(0);
  }
  results_to_pdf(results_in, &count, pdf_out);
  return count * LLR_t_value(pdf_out, 0.5, t0, t1);
}

double LLR_normalized_alt(double nt0, double nt1, int results_in[]) {
  /*
    See Section 4.2 in
    http://hardy.uhasselt.be/Fishtest/normalized_elo_practical.pdf
  */
  double pdf_out[2 * N];
  double count;
  double mu, var, sigma_pg, games, nt;
  results_to_pdf(results_in, &count, pdf_out);
  muvar(pdf_out, &mu, &var);
  if (N == 5) {
    sigma_pg = sqrt(2 * var);
    games = 2 * count;
  } else if (N == 3) {
    sigma_pg = sqrt(var);
    games = count;
  } else {
    assert(0);
  }
  nt = (mu - 0.5) / sigma_pg;
  return (games / 2.0) *
         log((1 + (nt - nt0) * (nt - nt0)) / (1 + (nt - nt1) * (nt - nt1)));
}

/*
  We use the BayesElo model to generate realistic pentanomial
  frequencies. Therefore, our logistic input parameters have to be
  converted to the BayesElo model. Strategy:

  - Convert (draw_ratio, bias) to (draw_elo, advantage).

  - Determine the Elo of the BayesElo model in such a way that the score
  as calculated using pentanomial probabilities derived from this Elo,
  corresponds to the given logistic/normalized Elo. This requires numerically
  solving a suitable equation.
*/

void proba_to_bayeselo(double P[], double *belo, double *drawelo) {
  /*
    Takes a probability: P[2], P[0]
    Returns elo, drawelo.
  */
  assert(0 < P[2] && P[2] < 1 && 0 < P[0] && P[0] < 1);
  *belo = 200 * log10(P[2] / P[0] * (1 - P[0]) / (1 - P[2]));
  *drawelo = 200 * log10((1 - P[0]) / P[0] * (1 - P[2]) / P[2]);
}

void ldw_calc(double belo, double de, double ldw[]) {
  ldw[2] = 1 / (1 + pow(10, (-belo + de) / 400));
  ldw[0] = 1 / (1 + pow(10, (belo + de) / 400));
  ldw[1] = 1 - ldw[2] - ldw[0];
}

void pent_calc(double belo, double draw_elo, double advantage, double pdf[]) {
  double ldw1[3], ldw2[3];
  int i, j, k;
  ldw_calc(belo + advantage, draw_elo, ldw1);
  ldw_calc(belo - advantage, draw_elo, ldw2);
  for (i = 0; i < 2 * N; i++) {
    pdf[i] = 0.0;
  }
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      k = i + j;
      pdf[2 * k] = k / 4.0;
      pdf[2 * k + 1] += ldw1[i] * ldw2[j];
    }
  }
}

typedef struct q {
  double s;
  double draw_elo;
  double advantage;
} q_t;

double g(double belo, void *args) { /* logistic Elo */
  double pdf[2 * N];
  q_t *qs = (q_t *)(args);
  double mu, var;
  pent_calc(belo, qs->draw_elo, qs->advantage, pdf);
  muvar(pdf, &mu, &var);
  return (qs->s) - mu;
}

double elo_to_belo(double elo, double draw_elo, double advantage) {
  double epsilon = 1e-9;
  double s = L_(elo);
  q_t qs = {s, draw_elo, advantage};
  stats_t stats = {0, 0, 0};
  double belo = brentq(g, -1000, 1000, epsilon, epsilon, 1000, &stats, &qs);
  assert(stats.error_num == 0);
  return belo;
}

double h(double belo, void *args) { /* normalized Elo, assumes pentanomial */
  double pdf[2 * N];
  q_t *qs = (q_t *)(args);
  double mu, var;
  pent_calc(belo, qs->draw_elo, qs->advantage, pdf);
  muvar(pdf, &mu, &var);
  return (qs->s) - (mu - 1 / 2.0) / sqrt(2 * var);
}

const double nelo_divided_by_nt = 347.43558552260146; // 800/log(10)

double nelo_to_belo(double nelo, double draw_elo, double advantage) {
  double epsilon = 1e-9;
  double s = nelo / nelo_divided_by_nt;
  q_t qs = {s, draw_elo, advantage};
  stats_t stats = {0, 0, 0};
  double belo = brentq(h, -1000, 1000, epsilon, epsilon, 1000, &stats, &qs);
  assert(stats.error_num == 0);
  return belo;
}

void be_data(double draw_ratio, double bias, double *draw_elo,
             double *advantage) {
  double P[3];
  double bias_s = L_(bias);
  P[2] = bias_s - draw_ratio / 2;
  P[1] = draw_ratio;
  P[0] = 1.0 - P[1] - P[2];
  proba_to_bayeselo(P, advantage, draw_elo);
}

/*
  End of BayesElo conversion.
*/

#define ELO_LOGISTIC 0
#define ELO_NORMALIZED 1

void simulate(uint64_t *prng, double alpha, double beta, double elo0,
              double elo1, int elo_model, double pdf[], int batch,
              int overshoot, int *status, int *duration, int *invalid) {
  int results[N] = {0, 0, 0, 0, 0};
  double LA = log(beta / (1 - alpha));
  double LB = log((1 - beta) / alpha);
  double l;
  double LLR_;
  double min_LLR = 0.0;
  double max_LLR = 0.0;
  double sq0 = 0.0;
  double sq1 = 0.0;
  double o0 = 0.0;
  double o1 = 0.0;
  double score0;
  double score1;
  int i;

  if (elo_model == ELO_LOGISTIC) {
    score0 = L_(elo0);
    score1 = L_(elo1);
  } else if (elo_model == ELO_NORMALIZED) {
    score0 = elo0 / nelo_divided_by_nt;
    score1 = elo1 / nelo_divided_by_nt;
  } else {
    assert(0);
  }

  *duration = 0;
  *status = CONTINUE;
  *invalid = 0;

  while (1) {
    (*duration)++;
    l = pick(prng, pdf);
    results[(int)(4.0 * l + 0.001)]++; /* excess of caution */
    if ((*duration) % batch != 0) {
      continue;
    }
    double count;
    double pdf2[2 * N];
    double pdf3[2 * N];
    double mu;
    double var;
    results_to_pdf(results, &count, pdf2);
    if (elo_model == ELO_LOGISTIC) {
      LLRjumps_expected(pdf2, score0, score1, pdf3);
    } else if (elo_model == ELO_NORMALIZED) {
      LLRjumps_t_value(pdf2, 0.5, score0, score1, pdf3);
    } else {
      assert(0);
    }
    muvar(pdf3, &mu, &var);
    LLR_ = count * mu;
    if (overshoot == 0) {
      if (LLR_ > LB) {
        *status = H1;
      } else if (LLR_ < LA) {
        *status = H0;
      }
    } else if (overshoot == 1) {
      /*
        https://hardy.uhasselt.be/Fishtest/stochastic_stopping.pdf
       */
      double p_upper = ucp(mu, var, batch, LB - LLR_);
      double p_lower = ucp(-mu, var, batch, LLR_ - LA);
      double u = myrand(prng);
      if (u >= p_upper) {
        *status = H1;
      } else if (u >= p_lower) {
        *status = H0;
      }
    } else if (overshoot == 2) {
      /*
         Dynamic overshoot correction using
         Siegmund - Sequential Analysis - Corollary 8.33.
      */
      if (LLR_ > max_LLR) {
        sq1 += (LLR_ - max_LLR) * (LLR_ - max_LLR);
        max_LLR = LLR_;
        o1 = sq1 / LLR_ / 2;
      }
      if (LLR_ < min_LLR) {
        sq0 += (LLR_ - min_LLR) * (LLR_ - min_LLR);
        min_LLR = LLR_;
        o0 = -sq0 / LLR_ / 2;
      }
      if (LLR_ > LB - o1) {
        *status = H1;
      } else if (LLR_ < LA + o0) {
        *status = H0;
      }
    } else {
      assert(0);
    }
    if (*status != CONTINUE) {
      /* The GSPRT does not work well with very low outcome values */
      for (i = 1; i < N - 1; i++) {
        if (results[i] == 0) {
          *invalid = 1;
          break;
        }
      }
      break;
    }
  }
}

typedef struct sim {
  /* identical for every thread */
  double alpha;
  double beta;
  double elo0;
  double elo1;
  double *pdf;
  int batch;
  int overshoot;
  int elo_model;
  /* in/out data */
  uint64_t prng;
  volatile int stop;
  volatile int count;
  volatile int pass;
  volatile double total_duration;
  volatile int invalid;
} sim_t;

void *sim_function(void *args) {
  sim_t *sim_;
  sim_ = (sim_t *)(args);
  uint64_t prng;
  pthread_mutex_lock(&mutex);
  jump(&(sim_->prng));
  prng = sim_->prng;
  pthread_mutex_unlock(&mutex);
  assert(sim_->stop == 0 || sim_->stop == 1);
  while (!sim_->stop) {
    int status;
    int duration;
    int invalid;
    simulate(&prng, sim_->alpha, sim_->beta, sim_->elo0, sim_->elo1,
             sim_->elo_model, sim_->pdf, sim_->batch, sim_->overshoot, &status,
             &duration, &invalid);
    duration *= 2;
    pthread_mutex_lock(&mutex);
    sim_->total_duration += duration;
    if (status == H1) {
      (sim_->pass)++;
    }
    (sim_->count)++;
    sim_->invalid += invalid;
    pthread_mutex_unlock(&mutex);
  }
  return NULL;
}

void usage() {
  printf(
      "simul [-h] [--alpha ALPHA] [--beta BETA] [--elo0 ELO0] [--elo1 ELO1] "
      "[--elo ELO] [--draw_ratio DRAW_RATIO] [--bias BIAS] [--ovcor ALGORITHM] "
      "[--threads THREADS] [--truncate TRUNCATE] [--batch BATCH] "
      "[--elo_model ELO_MODEL] "
      "[--seed SEED]\n");
}

void disp_elo_models(double pdf[], double pdf0[], double pdf1[], double belo,
                     double belo0, double belo1) { /* Assumes pentanomial */
  double mu, var, mu0, var0, mu1, var1, elo, elo0, elo1, nelo, nelo0, nelo1;
  double eps = 1e-6; /* to suppress spurious -0.0000 */
  muvar(pdf, &mu, &var);
  muvar(pdf0, &mu0, &var0);
  muvar(pdf1, &mu1, &var1);
  elo = Linv(mu) + eps;
  elo0 = Linv(mu0) + eps;
  elo1 = Linv(mu1) + eps;
  nelo = (mu - 0.5) / sqrt(2 * var) * nelo_divided_by_nt + eps;
  nelo0 = (mu0 - 0.5) / sqrt(2 * var0) * nelo_divided_by_nt + eps;
  nelo1 = (mu1 - 0.5) / sqrt(2 * var1) * nelo_divided_by_nt + eps;
  belo = belo + eps;
  belo0 = belo0 + eps;
  belo1 = belo1 + eps;
  printf("Elo         Logistic          Normalized            Bayes\n");
  printf("===         ========          ==========            =====\n");
  printf("Elo0      %10.5f          %10.5f       %10.5f                   \n",
         elo0, nelo0, belo0);
  printf("Elo1      %10.5f          %10.5f       %10.5f                   \n",
         elo1, nelo1, belo1);
  printf("Elo       %10.5f          %10.5f       %10.5f                   \n",
         elo, nelo, belo);
}

int main(int argc, char **argv) {
  double alpha = 0.05, beta = 0.05, elo0 = 0.0, elo1 = 5.0, elo = 0.0,
         draw_ratio = 0.61, bias = 0.0;
  int batch = 1;
  double p, ci;
  double belo, belo0, belo1, draw_elo, advantage;
  double pdf[2 * N], pdf0[2 * N], pdf1[2 * N];
  int num_threads = nproc();
  int overshoot = 2;
  int elo_model = ELO_LOGISTIC;
  int i;
  sim_t sim_;
  double av_duration;
  int truncate = 0;
  uint64_t seed = (uint64_t)time(0);
  for (i = 1; i <= argc - 1; i++) {
    if (strcmp(argv[i], "-h") == 0) {
      usage();
      return 0;
    } else if (strcmp(argv[i], "--alpha") == 0) {
      if (i < argc - 1) {
        alpha = atof(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--beta") == 0) {
      if (i < argc - 1) {
        beta = atof(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--elo0") == 0) {
      if (i < argc - 1) {
        elo0 = atof(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--elo1") == 0) {
      if (i < argc - 1) {
        elo1 = atof(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--elo") == 0) {
      if (i < argc - 1) {
        elo = atof(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--draw_ratio") == 0) {
      if (i < argc - 1) {
        draw_ratio = atof(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--bias") == 0) {
      if (i < argc - 1) {
        bias = atof(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--threads") == 0) {
      if (i < argc - 1) {
        num_threads = atoi(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--truncate") == 0) {
      if (i < argc - 1) {
        truncate = atoi(argv[i + 1]);
        if (truncate == 0) { /* hack */
          truncate = 1;
        }
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--batch") == 0) {
      if (i < argc - 1) {
        batch = atoi(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--elo_model") == 0) {
      if (i < argc - 1) {
        if (strcmp(argv[i + 1], "logistic") == 0) {
          elo_model = ELO_LOGISTIC;
        } else if (strcmp(argv[i + 1], "normalized") == 0) {
          elo_model = ELO_NORMALIZED;
        } else {
          usage();
          return 0;
        }
        i += 1;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--seed") == 0) {
      if (i < argc - 1) {
        seed = strtoull(argv[i + 1], NULL, 0);
        i++;
      } else {
        usage();
        return 0;
      }
    } else if (strcmp(argv[i], "--ovcor") == 0) {
      if (i < argc - 1) {
        overshoot = atoi(argv[i + 1]);
        i++;
      } else {
        usage();
        return 0;
      }
    } else {
      usage();
      return 0;
    }
  }
  if (alpha <= 0 || alpha >= 1 || beta <= 0 || beta >= 1 || elo1 <= elo0 ||
      !(overshoot == 0 || overshoot == 1 || overshoot == 2) || batch <= 0) {
    usage();
    return 0;
  }
  num_threads = MIN(num_threads, MAX_THREADS);
  num_threads = MAX(num_threads, 1);
  if (draw_ratio / 2 >= MIN(L_(bias), 1 - L_(bias))) {
    printf("The bias and the draw_ratio are not compatible.\n");
    return 0;
  }
  printf("Design parameters\n");
  printf("=================\n");
  printf("alpha      = %8.4f\nbeta       = %8.4f\nelo0       = %8.4f\n"
         "elo1       = %8.4f\nelo        = %8.4f\ndraw_ratio = %8.4f\n"
         "bias       = %8.4f\n"
         "ovcor      = %3d\nthreads    = %3d\ntruncate   =   %d\n"
         "batch      = %3d\n"
         "elo_model  =   %s\n"
         "seed       =   %" PRIu64 "\n\n",
         alpha, beta, elo0, elo1, elo, draw_ratio, bias, overshoot, num_threads,
         truncate, batch, elo_model == ELO_LOGISTIC ? "logistic" : "normalized",
         seed);
  be_data(draw_ratio, bias, &draw_elo, &advantage);
  if (elo_model == ELO_LOGISTIC) {
    belo = elo_to_belo(elo, draw_elo, advantage);
    belo0 = elo_to_belo(elo0, draw_elo, advantage);
    belo1 = elo_to_belo(elo1, draw_elo, advantage);
  } else if (elo_model == ELO_NORMALIZED) {
    belo = nelo_to_belo(elo, draw_elo, advantage);
    belo0 = nelo_to_belo(elo0, draw_elo, advantage);
    belo1 = nelo_to_belo(elo1, draw_elo, advantage);
  } else {
    assert(0);
  }
  pent_calc(belo, draw_elo, advantage, pdf);
  pent_calc(belo0, draw_elo, advantage, pdf0);
  pent_calc(belo1, draw_elo, advantage, pdf1);
  printf("BayesElo\n");
  printf("========\n");
  printf("draw_elo   = %8.4f\nadvantage  = %8.4f\n"
         "probs      =  [%f, %f, %f, %f, %f]\n\n",
         draw_elo, advantage, pdf[1], pdf[3], pdf[5], pdf[7], pdf[9]);
  disp_elo_models(pdf, pdf0, pdf1, belo, belo0, belo1);
  printf("\n");
  sim_.alpha = alpha;
  sim_.beta = beta;
  sim_.elo0 = elo0;
  sim_.elo1 = elo1;
  sim_.pdf = pdf;
  sim_.batch = batch;
  sim_.stop = 0;
  sim_.count = 0;
  sim_.pass = 0;
  sim_.total_duration = 0.0;
  sim_.invalid = 0;
  sim_.overshoot = overshoot;
  sim_.prng = seed;
  sim_.elo_model = elo_model;

  for (i = 0; i < num_threads; i++) {
    pthread_create(&(threads[i]), NULL, sim_function, (void *)(&sim_));
  }
  while (1) {
    sleep(2);
    if (sim_.count == 0) {
      continue;
    }
    p = sim_.pass / (sim_.count + 0.0);
    av_duration = sim_.total_duration / sim_.count;
    /* 3 sigma */
    ci = 3 * sqrt(p * (1 - p)) / sqrt(sim_.count);
    printf("sims=%d pass=%f[%f,%f] length=%.1f\n", sim_.count, p, p - ci,
           p + ci, av_duration);
    fflush(stdout);
    if (truncate != 0 && sim_.count >= truncate) {
      sim_.stop = 1;
      for (i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
      }
      break;
    }
  }
  return 0;
}
