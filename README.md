## A multi-threaded pentanomial simulator in C

1. The simulator implements a GSPRT for the `pentanomial model`.
It follows the [Fishtest implementation](https://github.com/glinscott/fishtest). The simulator
is not particularly optimized (in order not to obscure the
code too much) but it is fast enough to simulate tens of millions
of tests with realistic bounds (e.g. {0,2} in logistic Elo)
in a reasonable time frame. This is enough to measure pass
probabilities with 4 decimal digits.

2. The theoretical basis for the GSPRT is:

http://stat.columbia.edu/~jcliu/paper/GSPRT_SQA3.pdf

3. The theoretical basis for the pentanomial model is:

http://hardy.uhasselt.be/Fishtest/support_MLE_multinomial.pdf

4. The results of this simulator can be compared with the
[SPRT calculator](https://tests.stockfishchess.org/html/SPRTcalculator.html) used by
[Fishtest](https://tests.stockfishchess.org/tests).

5. Compile command:

```gcc -Wall -O3 simul.c -lm -lpthread -o simul```

6. Sample usage:

```
$ ./simul -h
simul [-h] [--alpha ALPHA] [--beta BETA] [--elo0 ELO0] [--elo1 ELO1] [--draw_ratio DRAW_RATIO] [--bias BIAS] [--noovcor] [--threads THREADS] [--truncate TRUNCATE]

$  ./simul --elo0 0 --elo1 2 --draw_ratio 0.56 --elo 2 --bias 90 --threads 4
Design parameters
=================
alpha      =   0.0500
beta       =   0.0500
elo0       =   0.0000
elo1       =   2.0000
elo        =   2.0000
draw_ratio =   0.5600
bias       =  90.0000
ovcor      =   1
threads    =   4
truncate   =   0

BayesElo
========
belo0      =   0.0000
belo1      =   3.2144
belo       =   3.2144
draw_elo   = 252.5476
advantage  = 142.4833
probs      =  [0.031425, 0.242521, 0.442462, 0.250299, 0.033293]

sims=192 pass=0.921875[0.863772,0.979978] length=63663.6 zeros=0.000000
sims=382 pass=0.929319[0.889980,0.968658] length=62359.4 zeros=0.000000
sims=590 pass=0.940678[0.911502,0.969854] length=61227.3 zeros=0.000000
sims=780 pass=0.950000[0.926589,0.973411] length=62116.7 zeros=0.000000
sims=973 pass=0.954779[0.934795,0.974763] length=61057.4 zeros=0.000000
sims=1166 pass=0.955403[0.937268,0.973538] length=61900.8 zeros=0.000000
sims=1352 pass=0.953402[0.936205,0.970599] length=62035.8 zeros=0.000000
sims=1549 pass=0.952227[0.935970,0.968485] length=61941.6 zeros=0.000000
sims=1782 pass=0.951740[0.936509,0.966970] length=61145.8 zeros=0.000000
sims=1962 pass=0.951070[0.936460,0.965681] length=60964.2 zeros=0.000000
sims=2167 pass=0.949239[0.935092,0.963385] length=61111.7 zeros=0.000000
sims=2367 pass=0.950993[0.937681,0.964305] length=60992.7 zeros=0.000000
sims=2572 pass=0.952177[0.939554,0.964800] length=60554.6 zeros=0.000000
sims=2771 pass=0.950559[0.938205,0.962914] length=60848.6 zeros=0.000000
...

$ ./simul --elo0 0 --elo1 2 --draw_ratio 0.56 --elo 0 --bias 90 --threads 4
Design parameters
=================
alpha      =   0.0500
beta       =   0.0500
elo0       =   0.0000
elo1       =   2.0000
elo        =   0.0000
draw_ratio =   0.5600
bias       =  90.0000
ovcor      =   1
threads    =   4
truncate   =   0

BayesElo
========
belo0      =   0.0000
belo1      =   3.2144
belo       =   0.0000
draw_elo   = 252.5476
advantage  = 142.4833
probs      =  [0.032347, 0.246400, 0.442505, 0.246400, 0.032347]

sims=223 pass=0.062780[0.014050,0.111511] length=56532.5 zeros=0.000000
sims=403 pass=0.059553[0.024187,0.094920] length=60041.7 zeros=0.000000
sims=590 pass=0.059322[0.030146,0.088498] length=61325.5 zeros=0.000000
sims=812 pass=0.057882[0.033297,0.082467] length=60446.0 zeros=0.000000
sims=983 pass=0.055951[0.033960,0.077942] length=60883.8 zeros=0.000000
sims=1190 pass=0.053782[0.034163,0.073400] length=61182.4 zeros=0.000000
sims=1388 pass=0.054035[0.035829,0.072240] length=61321.8 zeros=0.000000
sims=1579 pass=0.051298[0.034643,0.067953] length=61134.4 zeros=0.000000
sims=1789 pass=0.053661[0.037678,0.069645] length=60998.3 zeros=0.000000
sims=1989 pass=0.054299[0.039055,0.069542] length=60610.1 zeros=0.000000
sims=2201 pass=0.054521[0.040002,0.069039] length=60354.0 zeros=0.000000
sims=2414 pass=0.052610[0.038978,0.066241] length=60246.9 zeros=0.000000
sims=2607 pass=0.052551[0.039440,0.065661] length=60097.2 zeros=0.000000
sims=2818 pass=0.052874[0.040228,0.065521] length=60006.6 zeros=0.000000
...
```
7. The `bias` parameter is a proxy for the `RMS bias` of the opening
book.  The RMS bias is the Root Mean Square of the biases of the
openings in the book where the bias of an opening is defined as the
conversion to Elo (using the standard logistic formula) of the
expected score for white between engines of "equal
strength". Explicitly the RMS bias is the square root of the average
of the squares of the biases expressed in Elo. In the simulation we
assume that every opening has the same bias. One may show that in
first approximation this is correct.

8. The command line arguments for the simulator are in logistic Elo but to perform
simulations we need a method to obtain `realistic pentanomial
frequencies`. To this end invoke the [BayesElo
model](https://www.remi-coulom.fr/Bayesian-Elo/#theory) internally. Therefore,
our logistic input parameters `draw_ratio, bias, elo` have to be
converted to the BayesElo model. We follow the following strategy:

  * Convert `(draw_ratio, bias)` to `(draw_elo, advantage)`.

  * Determine the Elo for the BayesElo model (`belo`) in such a way
that the expected score as calculated using the pentanomial
probabilities derived from `belo, draw_elo, advantage`,
corresponds to the given logistic Elo (`elo`). This requires
numerically solving a suitable equation.

9. The theoretical guarantees of the GSPRT are asymptotic.
It does not work so well for very low outcome values.
The proportion of zero outcomes for LD, LW+DD, DW is contained in the
`zeros` output field.

10. We perform `dynamic overshoot correction` using `Siegmund -
Sequential Analysis - Corollary 8.33`. For a rough introduction
see http://hardy.uhasselt.be/Fishtest/dynamic_overshoot_correction.pdf.

## Summary of parameters

| Parameter | Explanation |
| --------- | ----------- |
| elo       | Actual Elo (difference) |
| elo0      | H0 corresponds to Elo=Elo0 |
| elo1      | H1 corresponds to Elo=Elo1 |
| alpha     | Pass probability if H0 is true |
| beta     | Fail probability if H1 is true |
| draw_ratio | Draw ratio between equal strength engines (taking into account the opening book) |
| bias       | A proxy for `RMS bias` |
| noovcor       | A flag to disable dynamic overshoot correction (for testing) |
| threads       | Simultaneous runs |
| truncate      | Stop the simulation after this many runs |


