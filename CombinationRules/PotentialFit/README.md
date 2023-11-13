Scripts to fit potentials for noble gases to QM dissociation curves
===================================================================
There are two important scripts in this directory. The first one to
fit nine potentials to QM curves. To run it for the first time try:
```
./fit_script -h
```
and study the flags. Without the -h (help) flag the script will fit
potentials to the CCSDT_CBS, corresponding to CCSD(T)/CBS and write
results in the RESULTS_CCSDT_CBS_20_0.1 directory, where 20 and 0.1
are the limits used for fitting.

The script will also generate all possible combination rules and evaluate
their accuracy in reproducing heterodimer dissociation curves.

Second virial coefficients
==========================
The other important script will compute second virial coefficients with
quantum corrections.
```
./calc_virials.py -h
```
again, without the help flag the script will use the default level of theory.
This script will fetch experimental data from the ../SecondVirial directory
so please make sure that you generated those numbers before running this.

Auxiliary scripts
=================
Some auxiliary scripts are available for plotting potentials (pot_plot.py), 
and the performance of combination rules (plot_comb_rules.py). 
Please have a look.