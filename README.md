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

5. Sample usage:

```
$ ./simul --elo_model normalized --draw_ratio 0.95 --elo1 5 --elo 2.5

Design parameters
=================
alpha      =   0.0500
beta       =   0.0500
elo0       =   0.0000
elo1       =   5.0000
elo        =   2.5000
draw_ratio =   0.9500
bias       =   0.0000
ovcor      =   1
threads    =   4
truncate   =   0
batch      =   1
elo_model  =   normalized
seed       =   1609184261

BayesElo
========
draw_elo   = 636.4258
advantage  =   0.0000
probs      =  [0.000586, 0.045994, 0.903702, 0.049052, 0.000667]

Elo         Logistic          Normalized            Bayes
===         ========          ==========            =====
Elo0         0.00000             0.00000          0.00000                   
Elo1         1.11905             5.00000         11.47029                   
Elo          0.55915             2.50000          5.73392                   

sims=3582 pass=0.498046[0.472983,0.523108] length=42403.3
sims=7134 pass=0.498038[0.480278,0.515797] length=42285.6
sims=10770 pass=0.498236[0.483782,0.512690] length=42140.3
sims=14416 pass=0.494936[0.482444,0.507429] length=42002.2
sims=18091 pass=0.495937[0.484785,0.507089] length=42015.6
sims=21758 pass=0.495680[0.485511,0.505848] length=42048.9
sims=25373 pass=0.498286[0.488869,0.507702] length=42118.0
...
```
6. The `bias` parameter is a proxy for the `RMS bias` of the opening
book.  The RMS bias is the Root Mean Square of the biases of the
openings in the book where the bias of an opening is defined as the
conversion to Elo (using the standard logistic formula) of the
expected score for white between engines of "equal
strength". Explicitly the RMS bias is the square root of the average
of the squares of the biases expressed in Elo. In the simulation we
assume that every opening has the same bias. One may show that in
first approximation this is correct.

7. The command line arguments for the simulator are in logistic Elo but to perform
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

8. We perform `dynamic overshoot correction` using `Siegmund -
Sequential Analysis - Corollary 8.33`. For a rough introduction
see http://hardy.uhasselt.be/Fishtest/dynamic_overshoot_correction.pdf.

9. Normalized Elo is discussed in http://hardy.uhasselt.be/Fishtest/normalized_elo.pdf.
However we now use a normalization factor
```
800/log(10)=347.43558552260146
```
which has the effect that for small elo differences and a perfectly balanced book the following
formula holds approximately
```
normalized_elo=logistic_elo/sqrt(1-draw_ratio)
```
Hence when the draw ratio is zero, normalized Elo and logistic Elo coincide. For other
draw ratios one has the following conversion table

| draw ratio | 0.0 | 0.3 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 |
|- |------|-----|-----|----|-|-|-|
| Normalized Elo| 5.00 | 5.00 | 5.00| 5.00 | 5.00 |5.00 |5.00|
| Logistic Elo| 5.00| 4.18 | 3.54 | 3.16 | 2.74 | 2.24 | 1.58|

A `SPRT(0,5)`  with bounds expressed in normalized Elo (the above sample case) takes about `42k` games to complete
(expected worst case) regardless of the book or the draw ratio. 

## Summary of parameters

| Parameter | Description | Default |
| --------- | ----------- |---------|
| --elo       | Actual Elo (difference) | 0 |
| --elo0      | H0 corresponds to Elo=Elo0 | 0 |
| --elo1      | H1 corresponds to Elo=Elo1 | 5 |
| --alpha     | Pass probability if H0 is true | 0.05|
| --beta     | Fail probability if H1 is true |0.05|
| --draw_ratio | Draw ratio between equal strength engines | 0.61 |
| --bias       | A proxy for `RMS bias` | 0 |
| --batch      | LLR calculation frequency |1 game pair|
| --ovcor       | Choose the discrete time correction algorithm | 1 |
| --threads       | Simultaneous runs |# of CPUs |
| --truncate      | Stop the simulation after this many runs |run forever|
| --elo_model     | logistic or normalized | logistic |
| --seed      | Seed for the random number generator| time(0) |


