## A multi-threaded pentanomial simulator in C

1. The simulator implements a GSPRT for the `pentanomial model`.
The theoretical basis for the GSPRT is:

http://stat.columbia.edu/~jcliu/paper/GSPRT_SQA3.pdf

2. The theoretical basis for the pentanomial model is:

http://hardy.uhasselt.be/Fishtest/support_MLE_multinomial.pdf

3. Compile command:

```gcc -Wall -O3 simul.c -lm -lpthread -o simul```

4. Sample usage:

```
$ ./simul -h
simul [-h] [--alpha ALPHA] [--beta BETA] [--elo0 ELO0] [--elo1 ELO1] [--draw_ratio DRAW_RATIO] [--bias BIAS] [--overshoot OVERSHOOT] [--threads THREADS]

$ ./simul --elo0 0 --elo1 2 --elo 0 --draw_ratio 0.74 --bias 300 --threads 4

Design parameters
=================
alpha      =   0.0500
beta       =   0.0500
elo0       =   0.0000
elo1       =   2.0000
elo        =   0.0000
draw_ratio =   0.7400
bias       = 300.0000
overshoot  =   1
threads    =   4

BayesElo
========
belo0      =   0.0000
belo1      =   8.6881
belo       =   0.0000
draw_elo   = 330.2304
advantage  = 663.1300
probs      =  [0.002855, 0.109373, 0.775545, 0.109373, 0.002855]

sims=570 pass=0.036842[0.013172,0.060512] length=19175.6 zeros=0.000000
sims=1203 pass=0.045719[0.027652,0.063786] length=19308.8 zeros=0.000000
sims=1714 pass=0.046674[0.031389,0.061960] length=19297.2 zeros=0.000000
sims=2329 pass=0.045513[0.032557,0.058470] length=19224.8 zeros=0.000000
sims=2911 pass=0.046032[0.034380,0.057684] length=19343.3 zeros=0.000000
sims=3461 pass=0.046518[0.035779,0.057258] length=19371.2 zeros=0.000000
sims=4107 pass=0.047236[0.037305,0.057167] length=19310.9 zeros=0.000000
sims=4612 pass=0.047918[0.038483,0.057354] length=19345.4 zeros=0.000000
sims=5217 pass=0.049645[0.040624,0.058667] length=19418.5 zeros=0.000000
sims=5787 pass=0.049248[0.040715,0.057782] length=19496.3 zeros=0.000000
sims=6321 pass=0.049834[0.041623,0.058045] length=19452.8 zeros=0.000000
sims=6932 pass=0.049625[0.041800,0.057450] length=19502.3 zeros=0.000000
sims=7473 pass=0.050047[0.042480,0.057614] length=19453.0 zeros=0.000000
sims=8091 pass=0.050550[0.043243,0.057857] length=19436.3 zeros=0.000000
sims=8694 pass=0.050495[0.043450,0.057540] length=19436.7 zeros=0.000000
sims=9246 pass=0.049751[0.042968,0.056535] length=19408.0 zeros=0.000000
...

$ ./simul --elo0 0 --elo1 2 --elo 2 --draw_ratio 0.74 --bias 300 --threads 4

Design parameters
=================
alpha      =   0.0500
beta       =   0.0500
elo0       =   0.0000
elo1       =   2.0000
elo        =   2.0000
draw_ratio =   0.7400
bias       = 300.0000
overshoot  =   1
threads    =   4

BayesElo
========
belo0      =   0.0000
belo1      =   8.6881
belo       =   8.6881
draw_elo   = 330.2304
advantage  = 663.1300
probs      =  [0.002698, 0.104042, 0.775328, 0.114912, 0.003019]

sims=605 pass=0.952066[0.926011,0.978122] length=19535.3 zeros=0.000000
sims=1193 pass=0.958089[0.940684,0.975494] length=19209.9 zeros=0.000000
sims=1745 pass=0.957593[0.943121,0.972065] length=19496.7 zeros=0.000000
sims=2370 pass=0.955696[0.943016,0.968376] length=19358.2 zeros=0.000000
sims=2921 pass=0.956522[0.945202,0.967842] length=19205.8 zeros=0.000000
sims=3485 pass=0.954950[0.944409,0.965490] length=19267.0 zeros=0.000000
sims=4039 pass=0.953949[0.944055,0.963843] length=19383.4 zeros=0.000000
sims=4571 pass=0.953621[0.944289,0.962952] length=19415.8 zeros=0.000000
sims=5186 pass=0.951600[0.942660,0.960541] length=19397.1 zeros=0.000000
sims=5752 pass=0.952017[0.943562,0.960471] length=19395.0 zeros=0.000000
sims=6358 pass=0.952973[0.945008,0.960937] length=19300.8 zeros=0.000000
sims=6979 pass=0.952285[0.944631,0.959940] length=19340.8 zeros=0.000000
sims=7488 pass=0.952190[0.944793,0.959587] length=19334.5 zeros=0.000000
sims=8048 pass=0.950795[0.943562,0.958028] length=19361.2 zeros=0.000000
sims=8635 pass=0.950434[0.943427,0.957441] length=19348.8 zeros=0.000000
sims=9198 pass=0.950968[0.944213,0.957722] length=19335.6 zeros=0.000000
sims=9839 pass=0.949995[0.943403,0.956587] length=19324.7 zeros=0.000000
...
```
5. All inputs for the simulator are in logistic Elo.

6. The `bias` parameter is a proxy for the `RMS bias` of the opening book.
The RMS bias is the Root Mean Square of the biases of the openings in
the book where the bias of an opening is defined as the conversion to
Elo (using the standard logistic formula) of the expected score for
white between engines of "equal strength". Explicitly the RMS bias is
the square root of the average of the squares of the biases expressed
in Elo.

7. To perform simulations we need a method to obtain `realistic
pentanomial frequencies`. To this end we use the BayesElo
model. Therefore, our logistic input parameters `draw_ratio, bias,
elo` have to be converted to the BayesElo model. We follow the
following strategy:

  * Convert the `draw_ratio` to `draw_elo`, assuming equal strength
engines.

  * Convert the `bias` expressed in logistic Elo to the `advantage`
parameter in BayesElo using the standard scale factor. 

  * Determine the Elo for the BayesElo model (`belo`) in such a way
that the expected score as calculated using the pentanomial
probabilities (`probs`) derived from `belo, draw_elo, advantage`,
corresponds to the given logistic Elo (`elo`). This requires
numerically solving a suitable equation.

8. The theoretical guarantees of the GSPRT are asymptotic.
It does not work so well for very low outcome values.
The proportion of zero outcomes for LD, LW+DD, DW is contained in the
`zeros` field.

9. We perform dynamic overshoot correction using `Siegmund -
Sequential Analysis - Corollary 8.33.` This is currently not yet
implemented in Fishtest as it requires keeping a bit of state.