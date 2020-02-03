#include <pthread.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
#define MAX_THREADS 16
pthread_t threads[MAX_THREADS];

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

typedef struct {
    int funcalls;
    int iterations;
    int error_num;
} stats_t;

#define SIGNERR -1
#define CONVERR -2

typedef double (*callback_type)(double,void*);

/* This is a slightly modified version of brentq provided by SciPy */

/* Written by Charles Harris charles.harris@sdl.usu.edu */

/*
  At the top of the loop the situation is the following:
    1. the root is bracketed between xa and xb
    2. xa is the most recent estimate
    3. xp is the previous estimate
    4. |fp| < |fb|
  The order of xa and xp doesn't matter, but assume xp < xb. Then xa lies to
  the right of xp and the assumption is that xa is increasing towards the root.
  In this situation we will attempt quadratic extrapolation as long as the
  condition
  *  |fa| < |fp| < |fb|
  is satisfied. That is, the function value is decreasing as we go along.
  Note the 4 above implies that the right inequlity already holds.
  The first check is that xa is still to the left of the root. If not, xb is
  replaced by xp and the interval reverses, with xb < xa. In this situation
  we will try linear interpolation. That this has happened is signaled by the
  equality xb == xp;
  The second check is that |fa| < |fb|. If this is not the case, we swap
  xa and xb and resort to bisection.
*/

double
brentq(callback_type f, double xa, double xb, double xtol, double rtol,
       int iter, stats_t *stats, void* args)
{
    double xpre = xa, xcur = xb;
    double xblk = 0., fpre, fcur, fblk = 0., spre = 0., scur = 0., sbis;
    /* the tolerance is 2*delta */
    double delta;
    double stry, dpre, dblk;
    int i;

    fpre = (*f)(xpre, args);
    fcur = (*f)(xcur, args);
    stats->funcalls = 2;
    if (fpre*fcur > 0) {
        stats->error_num = SIGNERR;
        return 0.;
    }
    if (fpre == 0) {
        return xpre;
    }
    if (fcur == 0) {
        return xcur;
    }

    stats->iterations = 0;
    for (i = 0; i < iter; i++) {
        stats->iterations++;
        if (fpre*fcur < 0) {
            xblk = xpre;
            fblk = fpre;
            spre = scur = xcur - xpre;
        }
        if (fabs(fblk) < fabs(fcur)) {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;

            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }

        delta = (xtol + rtol*fabs(xcur))/2;
        sbis = (xblk - xcur)/2;
        if (fcur == 0 || fabs(sbis) < delta) {
            return xcur;
        }

        if (fabs(spre) > delta && fabs(fcur) < fabs(fpre)) {
            if (xpre == xblk) {
                /* interpolate */
                stry = -fcur*(xcur - xpre)/(fcur - fpre);
            }
            else {
                /* extrapolate */
                dpre = (fpre - fcur)/(xpre - xcur);
                dblk = (fblk - fcur)/(xblk - xcur);
                stry = -fcur*(fblk*dblk - fpre*dpre)
                    /(dblk*dpre*(fblk - fpre));
            }
            if (2*fabs(stry) < MIN(fabs(spre), 3*fabs(sbis) - delta)) {
                /* good short step */
                spre = scur;
                scur = stry;
            } else {
                /* bisect */
                spre = sbis;
                scur = sbis;
            }
        }
        else {
            /* bisect */
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur; fpre = fcur;
        if (fabs(scur) > delta) {
            xcur += scur;
        }
        else {
            xcur += (sbis > 0 ? delta : -delta);
        }

        fcur = (*f)(xcur, args);
        stats->funcalls++;
    }
    stats->error_num = CONVERR;
    return xcur;
}

#define N 5
#define H0 0
#define H1 1
#define CONTINUE 2

double L_(double x){
  return 1/(1+pow(10,-x/400.0));
}

void muvar(double pdf_in[], double *mu, double *var){
  int i;
  double a,p;
  double epsilon=1e-6;
  double sum=0.0;
  double sum2=0.0;
  *mu=0.0;
  for(i=0;i<N;i++){
    a=pdf_in[2*i];
    p=pdf_in[2*i+1];
    assert(-epsilon<=p);
    assert(p<=1+epsilon);
    sum+=p;
    *mu+=a*p;
    sum2+=a*a*p;
  } 
  assert(fabs(sum-1)<epsilon);
  *var=sum2-(*mu)*(*mu);
}

typedef struct p {
  double *pdf_in;
  double s;
} p_t;

double f(double x, void *args){
  int i;
  double a,p;
  double sum=0.0;
  p_t *ps=(p_t*)(args);
  double s=ps->s;
  double *pdf_in=ps->pdf_in;
  for(i=0;i<N;i++){
    a=pdf_in[2*i];
    p=pdf_in[2*i+1];
    sum+=p*(a-s)/(1+x*(a-s));
  }
  return sum;
}

void MLE(double pdf_in[], double s, double pdf_out[]){
  /*
    This function computes the maximum likelood estimate for a
    discrete distribution with expectation value s, given an observed
    (i.e. empirical) distribution pdf.
    
    pdf is an array[2*N] consisting of pairs (ai,pi), i=1,...,N.  It
    is assumed that that the ai are strictly ascending, a1<s<aN and
    p1>0, pN>0.
    
    The theory behind this function can be found in the online 
    document

    http://hardy.uhasselt.be/Fishtest/support_MLE_multinomial.pdf
    
    (see Proposition 1.1).
  */
  double epsilon=1e-9;
  double v,w,l,u,x,p,a,mu,var;
  stats_t stats={0,0,0};
  int i;
  p_t ps;
  v=pdf_in[0];
  w=pdf_in[2*(N-1)];
  l=-1/(w-s);
  u=1/(s-v);
  ps.s=s;
  ps.pdf_in=pdf_in;
  x=brentq(f,l+epsilon,u-epsilon,epsilon,epsilon,1000,&stats,&ps);
  assert(stats.error_num==0);
  for(i=0;i<N;i++){
    a=pdf_in[2*i];
    p=pdf_in[2*i+1];
    pdf_out[2*i]=a;
    pdf_out[2*i+1]=p/(1+x*(a-s));
  }
  muvar(pdf_out,&mu,&var); /* for validation */
  assert(fabs(s-mu)<1e-6);
}

double pick(double pdf[]){
  double x = (float)rand()/(float)(RAND_MAX);
  int i;
  double s=0.0;
  for(i=0;i<N;i++){
    s+=pdf[2*i+1];
    if(x<s){
      return pdf[2*i];
    }
  }
  /* we should not get here */
  return pdf[2*N];
}

double LLR(double pdf_in[], double s0, double s1){
  double pdf0[2*N], pdf1[2*N];
  double p,p0,p1;
  double sum=0.0;
  int i;
  MLE(pdf_in,s0,pdf0);
  MLE(pdf_in,s1,pdf1);
  for(i=0;i<N;i++){
    p=pdf_in[2*i+1];
    p0=pdf0[2*i+1];
    p1=pdf1[2*i+1];
    sum+=p*log(p1/p0);
  }
  return sum;
}

void regularize(int results_in[], double results_out[]){
  /* 
     Replace zeros with a small value to avoid division by
     zero later on.
  */
  double epsilon=1e-3;
  int i;
  for(i=0; i<N; i++){
    if(results_in[i]==0){
      results_out[i]=epsilon;
    }else{
      results_out[i]=(double)(results_in[i]);
    }
  }
}

void results_to_pdf(int results_in[], double *count, double pdf_out[]){
  double results_out[N];
  int i;
  *count=0.0;
  regularize(results_in,results_out);
  for(i=0;i<N;i++){
    *count+=results_out[i];
  }
  for(i=0;i<N;i++){
    pdf_out[2*i]=i/(N-1.0);
    pdf_out[2*i+1] =results_out[i]/(*count);
  }
}
    
void LLR_logistic(double s0, double s1, int results_in[], double *LLR_){
  /*
    This function computes the generalized log-likelihood ratio for
    "results" which should be an array if length 5 containing the
    frequencies of the game pairs LL,LD+DL,LW+DD+WL,DW+WD,WW.
  */
  double pdf_out[2*N];
  double count;
  results_to_pdf(results_in, &count, pdf_out);
  *LLR_=count*LLR(pdf_out,s0,s1);
}

/*
  We use the BayesElo model to generate realistic pentanomial
  frequencies. Therefore, our logistic input parameters have to be
  converted to the BayesElo model. Strategy:

  - Convert the draw_ratio to draw_elo, assuming equal strength
  engines.

  - Convert the bias expressed in logistic Elo to the advantage
  parameter in BayesElo using the standard scale factor. 

  - Determine the Elo of the BayesElo model in such a way that the score
  as calculated using pentanomial probabilities derived from this Elo,
  corresponds to the given logistic Elo. This requires numerically solving
  a suitable equation.
*/

double draw_elo_calc(double draw_ratio){
  return 400*(log(1/((1-draw_ratio)/2.0)-1)/log(10));
}

double scale(double de){
  return 4*pow(10,(-de/400))/pow(1+pow(10,(-de/400)),2);
}

void ldw_calc(double belo, double de, double ldw[]){
  ldw[2]=1/(1+pow(10,(-belo+de)/400));
  ldw[0]=1/(1+pow(10,(belo+de)/400));
  ldw[1]=1-ldw[2]-ldw[0];
}

void pent_calc(double belo, double draw_elo, double advantage, double pdf[]){
  double ldw1[3],ldw2[3];
  int i,j,k;
  ldw_calc(belo+advantage,draw_elo,ldw1);
  ldw_calc(belo-advantage,draw_elo,ldw2);
  for(i=0;i<2*N;i++){
    pdf[i]=0.0;
  }
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      k=i+j;
      pdf[2*k]=k;
      pdf[2*k+1]+=ldw1[i]*ldw2[j];
    }
  }
}

typedef struct q {
  double s;
  double draw_elo;
  double advantage;
} q_t;

double g(double belo, void *args){
  double pdf[2*N];
  q_t *qs=(q_t*)(args);
  double mu,var;
  pent_calc(belo,qs->draw_elo,qs->advantage,pdf);
  muvar(pdf,&mu,&var);
  return 4*(qs->s)-mu;
}

double elo_to_belo(double elo, double draw_elo, double advantage){
  double epsilon=1e-9;
  double s=L_(elo);
  q_t qs={s,draw_elo,advantage};
  stats_t stats={0,0,0};
  double belo=brentq(g,-1000,1000,epsilon,epsilon,1000,&stats,&qs);
  assert(stats.error_num==0);
  return belo;
}

void be_data(double draw_ratio, double bias, double *draw_elo, double *advantage){
  *draw_elo=draw_elo_calc(draw_ratio);
  *advantage=bias/scale(*draw_elo);
}

/*
  End of BayesElo conversion.
*/

void simulate(double alpha,double beta,double elo0,double elo1,double pdf[],
	      int overshoot,int *status,int *duration,int *invalid){
  int results[N]={0,0,0,0,0};
  double LA=log(beta/(1-alpha));
  double LB=log((1-beta)/alpha);
  int l;
  double LLR_;
  double min_LLR=0.0;
  double max_LLR=0.0;
  double sq0=0.0;
  double sq1=0.0;
  double o0=0.0;
  double o1=0.0;
  double score0=L_(elo0);
  double score1=L_(elo1);
  int i;
  
  *duration=0;
  *status=CONTINUE;
  *invalid=0;

  while(1){
    (*duration)++;
    l=pick(pdf);   /* no locking as rand() is marked as thread safe in the man page */
    results[l]++;
    LLR_logistic(score0,score1,results, &LLR_);
    /* 
       Dynamic overshoot correction using
       Siegmund - Sequential Analysis - Corollary 8.33.
    */
    if(LLR_>max_LLR){
      sq1+=(LLR_-max_LLR)*(LLR_-max_LLR);
      max_LLR=LLR_;
      o1=sq1/LLR_/2;
    }
    if(LLR_<min_LLR) {
      sq0+=(LLR_-min_LLR)*(LLR_-min_LLR);
      min_LLR=LLR_;
      o0=-sq0/LLR_/2;
    }
    assert(overshoot==0 || overshoot==1);
    if(!overshoot){
       o0=0;o1=0;
    }
    if(LLR_>LB-o1){
      *status=H1;
    }else if(LLR_ < LA+o0){
      *status=H0;
    }
    if(*status!=CONTINUE){
      /* The GSPRT does not work well with very low outcome values */
      for(i=1;i<N-1;i++){
	if(results[i]==0){
	  *invalid=1;
	  break;
	}
      }
      break;
    }
  }
}

typedef struct sim {
  double alpha;
  double beta;
  double elo0;
  double elo1;
  double *pdf;
  int    overshoot;
  /* out data */
  int    count;
  int    pass;
  double total_duration;
  int    invalid;
} sim_t;

void *sim_function(void *args){
  while(1){
    int status;
    int duration;
    int invalid;
    sim_t *sim_;
    sim_=(sim_t *)(args);
    simulate(sim_->alpha,sim_->beta,sim_->elo0,sim_->elo1,sim_->pdf,
	     sim_->overshoot,&status,&duration,&invalid);
    duration*=2;
    pthread_mutex_lock(&mutex);
    sim_->total_duration+=duration;
    if(status==H1){
      (sim_->pass)++;
    }
    (sim_->count)++;
    sim_->invalid+=invalid;
    pthread_mutex_unlock(&mutex);
  }
}

void usage(){
  printf("simul [-h] [--alpha ALPHA] [--beta BETA] [--elo0 ELO0] [--elo1 ELO1] "
	 "[--draw_ratio DRAW_RATIO] [--bias BIAS] [--overshoot OVERSHOOT] "
	 "[--threads THREADS]\n");
}

int main(int argc, char **argv){
  double alpha=0.05,beta=0.05,elo0=0.0,elo1=5.0,elo=0.0,draw_ratio=0.61,bias=0.0;
  double p,ci;
  double belo,belo0,belo1,draw_elo,advantage;
  double pdf[2*N];
  int num_threads=1;
  int overshoot=1;
  int i;
  sim_t sim_;
  double av_duration;
  double inv;
  for(i=1;i<=argc-1;i++){
    if(strcmp(argv[i],"-h")==0){
      usage();
      return 0;
    }else if(strcmp(argv[i],"--alpha")==0){
      if(i<argc-1){
	alpha=strtof(argv[i+1],NULL);
	i++;
      }else{
	usage();
	return 0;
      }
    }else if(strcmp(argv[i],"--beta")==0){
      if(i<argc-1){
	beta=strtof(argv[i+1],NULL);
	i++;
      }else{
	usage();
	return 0;
      }
    }else if(strcmp(argv[i],"--elo0")==0){
      if(i<argc-1){
	elo0=strtof(argv[i+1],NULL);
	i++;
      }else{
	usage();
	return 0;
      }
    }else if(strcmp(argv[i],"--elo1")==0){
      if(i<argc-1){
	elo1=strtof(argv[i+1],NULL);
	i++;
      }else{
	usage();
	return 0;
      }
    }else if(strcmp(argv[i],"--elo")==0){
      if(i<argc-1){
	elo=strtof(argv[i+1],NULL);
	i++;
      }else{
	usage();
	return 0;
      }
    }else if(strcmp(argv[i],"--draw_ratio")==0){
      if(i<argc-1){
	draw_ratio=strtof(argv[i+1],NULL);
	i++;
      }else{
	usage();
	return 0;
      }
    }else if(strcmp(argv[i],"--bias")==0){
      if(i<argc-1){
	bias=strtof(argv[i+1],NULL);
	i++;
      }else{
	usage();
	return 0;
      }
    }else if(strcmp(argv[i],"--threads")==0){
      if(i<argc-1){
	num_threads=(int)(strtof(argv[i+1],NULL));
	i++;
      }else{
	usage();
	return 0;
      }
    }else if(strcmp(argv[i],"--overshoot")==0){
      if(i<argc-1){
	overshoot=(int)(strtof(argv[i+1],NULL));
	i++;
      }else{
	usage();
	return 0;
      }
    }else{
      usage();
      return 0;
    }
  }
  if(alpha<=0||alpha>=1||beta<=0||beta>=1||elo1<=elo0||!(overshoot==0 || overshoot==1)){
    usage();
    return 0;
  }
  srand(time(0));
  num_threads=MIN(num_threads,MAX_THREADS);
  num_threads=MAX(num_threads,1);
  printf("Design parameters\n");
  printf("=================\n");
  printf("alpha      = %8.4f\nbeta       = %8.4f\nelo0       = %8.4f\n"
	 "elo1       = %8.4f\nelo        = %8.4f\ndraw_ratio = %8.4f\n"
	 "bias       = %8.4f\novershoot  = %3d\nthreads    = %3d\n\n",
	 alpha,beta,elo0,elo1,elo,draw_ratio,bias,overshoot,num_threads);
  be_data(draw_ratio,bias,&draw_elo,&advantage);
  belo=elo_to_belo(elo,draw_elo,advantage);
  belo0=elo_to_belo(elo0,draw_elo,advantage);
  belo1=elo_to_belo(elo1,draw_elo,advantage);
  pent_calc(belo,draw_elo,advantage,pdf);
  printf("BayesElo\n");
  printf("========\n");
  printf("belo0      = %8.4f\nbelo1      = %8.4f\nbelo       = %8.4f\n"
	 "draw_elo   = %8.4f\nadvantage  = %8.4f\n"
	 "probs      =  [%f, %f, %f, %f, %f]\n\n",
	 belo0,belo1,belo,draw_elo,advantage,
	 pdf[1],pdf[3],pdf[5],pdf[7],pdf[9]);
  sim_.alpha=alpha;
  sim_.beta=beta;
  sim_.elo0=elo0;
  sim_.elo1=elo1;
  sim_.pdf=pdf;
  sim_.pass=0;
  sim_.count=0;
  sim_.total_duration=0.0;
  sim_.invalid=0;
  sim_.overshoot=overshoot;
  
  for(i=0;i<num_threads;i++){
    pthread_create(&(threads[i]), NULL, sim_function, (void*) (&sim_));
  }
  while(1){
    sleep(2);
    if(sim_.count==0){
      continue;
    }
    p=sim_.pass/(sim_.count+0.0);
    inv=sim_.invalid/(sim_.count+0.0);
    av_duration=sim_.total_duration/sim_.count;
    /* 3 sigma */
    ci=3*sqrt(p*(1-p))/sqrt(sim_.count);
    printf("sims=%d pass=%f[%f,%f] length=%.1f zeros=%f\n",sim_.count,p,p-ci,p+ci,av_duration,inv);
    fflush(stdout);
  }
}


