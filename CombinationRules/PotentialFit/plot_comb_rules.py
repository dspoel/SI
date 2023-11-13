#!/usr/bin/env python3

import os, sys, math
from fit_script import get_potentials
import numpy as np
import holoviews as hv
from holoviews import opts
hv.extension('matplotlib')

def read_comb(infile:str)->list:
    mylines = []
    with open(infile, "r") as inf:
        params = {}
        for line in inf:
            words = line.strip().split(",")
            if line.find("RMSE") >= 0:
                for i in range(1,len(words)):
                    params[words[i]] = i
            elif len(params)+1 == len(words):
                mylines.append(words)
            else:
                sys.exit("Inconsistency on line '%s' in %s" % ( line.strip(), infile ) )
    
    return params, mylines

def plotit2(params:list, mylines:list, outf, pot:str)->hv.HeatMap:
    print("Plotting heatmap for %s" % pot)
    newdata = []
    eqn = { "f1": 10, "f2": 11, "f3": 12, "f4": 14,
            "f5": 13, "f6": 17, "f7": 18, "f8": 16, 
            "f9": 20, "f10": 15, "f12": 19, "f14": 21 }
    for i in range(len(mylines)):
        eps = mylines[i][params["epsilon"]]
        sig = mylines[i][params["sigma"]]
        if eps in eqn and sig in eqn:
            # Making strings out of the equation to keep the
            # boxes equal in size.
            newentry = (str(eqn[eps]), str(eqn[sig]),
                        float(mylines[i][params["RMSE"]]))
            newdata.append(newentry)
    return hv.HeatMap(newdata).sort().aggregate(function=np.min)
    
if __name__ == "__main__":
    potentials = get_potentials()
    hm = {}
    cpdir = "comb_plots"
    os.makedirs(cpdir, exist_ok=True)
    lllots = { "CCSDT_CBS": "CCSD(T)/CBS", "CCSDTQ_CBS": "CCSDT(Q)/CBS" }
    for lot in lllots.keys():
        outfile = ( "comb_corr_%s.out" % lot )
        with open(outfile, "w") as outf:
            for pot in potentials.keys():
                infile = ( "RESULTS_%s_20_0.1/%s_comb.csv" % ( lot, pot ))
                params, mylines = read_comb(infile)
                if not "epsilon" in params or not "sigma" in params:
                    continue
                mytitle = lllots[lot] + " - " + pot
                hm[pot] = plotit2(params, mylines, outf, pot).opts(show_values=False, title=mytitle, xlabel="epsilon", ylabel="sigma", clabel="RMSE (kJ/mol)", colorbar=True, fontsize='x-large', clim=(0, 1.5))
                pdfname = ( "%s/%s-%s" % ( cpdir, lot, pot ))
                hv.save(hm[pot], filename=pdfname, fmt='pdf')
