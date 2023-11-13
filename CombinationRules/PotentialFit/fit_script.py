#!/usr/bin/env python3
from numpy import arange
from pandas import read_csv
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import scipy
from matplotlib import pyplot
import numpy as np 
import itertools, math, os, argparse
import pandas, sys
from lmwei2017 import Uwei, get_table1, Bohr, Hartree

debug = False
short = { "helium": "He", "neon": "Ne", "argon": "Ar", "krypton": "Kr", "xenon": "Xe", "radon": "Rn" }

#Definition of fitted potentials
###
def LJ_12_6(x, e, s):
    return (4 * e * ( ((s/x)**12) - ((s/x)**6)) )
def LJ_8_6(x, e, s):
    return ( e * (  3*((s/x)**8) - 4*((s/x)**6)) )
def LJ_14_7(x, e, s, g, d):
    if 0 == s:
        return 0
    return ( e*((((1+d)/((x/s)+d))**7) * (((1+g)/(((x/s)**7) + g   ) ) - 2 ) ))
def BHA(x, e, s, g):
    return (e*np.exp(-g*x) - s/x**6)
def WBH(x, e, s, g):
    if 0 == s:
        return 0
    return ((2 * e) / (1 - (3 / (3 + g)) )) * ((s**6)/((s**6) + (x**6))) *  (   ((3/(g + 3)) * np.exp(g * (1 - (x/s)))) - 1 )
def GBH(x, e, s, g, d):
    if 0 == g or 0 == s:
        return 0
    return ( e*((d + 2*g + 6) /(2*g)) * (1/(1 + ((x/s)**6))) * ( ((6+d)/(d + 2*g + 6))*np.exp(g*(1-(x/s))) - 1 ) - e/(1+((x/s)**d)))
def MBH(x, e, s, g):
    if 0 == g or 0 == s:
        return 0
    return ((e/(1- (6/g))) * ( (6/g)* np.exp(g*(1-(x/s)))  - ((s/x)**6) )      )
def MRS(x, e, s, g):
    return (e* (  np.exp(-2*g*(x - s))  - 2*np.exp(-g*(x - s)) ))
def TT(x:float, A:float, b:float, c6:float, c8:float, c10:float, Re:float, De:float)->float:
    params = { "A": A*De, "b": b/Re, "C6": c6*De*Re**6, "C8": c8*De*Re**8, "C10": c10*De*Re**10, "Re": Re, "De": De }
    return Uwei(x, params)

def compute(myfunction, x, params:dict)->np.array:
    y = np.copy(x)
    i = 0
    function = eval(myfunction)
    for distance in x:
        Y = "undefined"
        if function == WBH or function == MRS or function == MBH or function == BHA:
            Y = function(distance, e=params["epsilon"], s=params["sigma"], g=params["gamma"])
        elif function == GBH or function == LJ_14_7:
            Y = function(distance, params["epsilon"], params["sigma"], params["gamma"], params["delta"])
        elif function == LJ_12_6 or function == LJ_8_6:
            Y = function(distance, params["epsilon"], params["sigma"])
        elif function == TT:
            Y = function(distance, params["A"], params["b"], params["C6"], params["C8"], params["C10"], params["Re"], params["De"])
        else:
            sys.exit("function %s not recognized in RMSE calculator" % function)
        y[i] = Y
        i += 1

    return y


#Definition of combination rules
#geometric mean
def f1(a1, a2):
    return np.sqrt( a1 * a2 )
#arithmetic mean
def f2(a1, a2):
    return 0.5*( a1 + a2 )
def f3(a1, a2):
    return (2.0 * a1 * a2)/( a1 + a2 )
def f4(a1, a2):
    return (( a1 * a2 ) * ( a1 + a2 ) / ((a1**2) + (a2**2))  )
def f6(a1, a2):
    return (((a1**6) + (a2**6))/2)**(1/6)
def f7(a1, a2):
    return ((a1**3) + (a2**3) )/ ((a1**2) + (a2**2))
#Hogervorst for sigma
def f5(e, e2, g, g2, s, s2):
    tempi = abs(e * g * (s ** 6) /(g - 6 ))
    tempj = abs(e2 * g2 * (s2 ** 6) /(g2 - 6  ))
    return ((( tempi * tempj ) ** (1/2) )* abs(f2(g, g2) - 6) /(f2(g, g2) * f3(e, e2)) ) ** (1/6)
#Qi et al, proposed for epsilon
def f8(e1, e2, s1, s2):
    return (( e1 * e2 )**(1/2)) * ( 2 * (s1**3) * (s2**3) / ((s1**6) + (s2**6)) )
#"spherical" average
def f9(a, a2):
    return (((a**3) + (a2**3))/2)**(1/3)
#Kong mason, used for gamma, g
def f10(g, g2, s, s2):
    return f1(s, s2) * (0.5*((g/s)+(g2/s2)))
#used for sigma,  Matar2004, eq.  6
def f11(e1, e2, s1, s2):
    return ( (0.5**( ( (e1**(1/13))*(s1**(12/13)) + (e2**(1/13))*(s2**(12/13))  )) **(13))/(f1(e1, e2) * ((s1*s2)**3))  )**(1/6)
#matar2004, eq. 9
def f12(a1, a2):
    return 4*a1*a2/ ( (a1**(1/2) + (a2**(1/2)) )**2 )
# New stuff
def f14(a1:float, a2:float):
    denom = 1/a1**2 + 1/a2**2
    return (2/denom)**(1/2.)
def f15(a1:float, a2:float):
    denom = 1/a1**3 + 1/a2**3
    return (2/denom)**(1/3.)
def f16(a1:float, a2:float):
    denom = 1/a1**4 + 1/a2**4
    return (2/denom)**(1/4.)
def f17(a1:float, a2:float):
    denom = 1/a1**(0.5) + 1/a2**(0.5)
    return (2/denom)**(2)
###
#some unused rules. were not successful in trials
###
#def fkm_lj(e, e2, s, s2, s12):
#    #for epsilon in 12-6
#    return ( ( (e*(s**(12)))/((2**13)*(s12)**(12)) )*(1+ (e2*(s2**(12))/e*(s**12) )**(1/13) )**(13)  )
#def fmason(e, e2, g, g2, g12, s, s2, s12):
#    return (        ( (e * e2)**(1/2) ) * ( (g * g1/(g12**2) )**3 ) * ( (1-6/s12)**2 / ((1-6/s)*(1-6/s2)) )**(1/2)   )
#def fkm2(g, g2, s, s2):
#    return f7(s, s2) * (0.5*((g/s)+(g2/s2)))
#def ftest(g, g2, s, s2):
#    return (( g * g2 )**(1/2)) * ( 2 * ((s**6) + (s2**6)) / (s**3) * (s2**3) )
###

def pol_string(pol:bool):
    if pol:
        return " Electrostatics were taken into account through Coulomb interactions and polarizability."
    else:
        return ""

def latexit(name:str)->str:
    return name.replace("_", "\_")
    
#RMSE calculator
def RMSE(function, x, y, params:dict)->float:
    sumdiff = 0
    count = 0
#    print(function)
    for distance, energy in zip(x, y):
        Y = "undefined"
        if function == WBH or function == MRS or function == MBH or function == BHA:
            Y = function(distance, e=params["epsilon"], s=params["sigma"], g=params["gamma"])
        elif function == GBH or function == LJ_14_7:
            Y = function(distance, params["epsilon"], params["sigma"], params["gamma"], params["delta"])
        elif function == LJ_12_6 or function == LJ_8_6:
            Y = function(distance, params["epsilon"], params["sigma"])
        elif function == TT:
            Y = function(distance, params["A"], params["b"], params["C6"], params["C8"], params["C10"], params["Re"], params["De"])
        else:
            sys.exit("function %s not recognized in RMSE calculator" % function)
        count = count + 1
        diff2 = (Y - energy)**2
        sumdiff = sumdiff +  diff2
    RMSE = np.sqrt(sumdiff/count)
#    print(f"{function}: {RMSE}")
    return RMSE

def read_data(noble1:str, noble2:str, lot:str, upper:float, lower:float, polarizable:bool):
    dataframe = read_csv(f"{args.lot}/{noble1}#{noble2}.csv", sep=",")
    # Finding a minimum and gathering dataset based on bounds
    data = dataframe.values
    x_1, y_1 = data[:, 0], data[:, -1]

    # Subtract electrostatics if needed
    if polarizable:
        elecframe = read_csv(f"Electrostatics/{noble1}#{noble2}.csv", sep=",")
        elec      = elecframe.values
        x_p, y_p  = elec[:, 0], elec[:, -1]
        eldict    = CubicSpline(x_p, y_p)
    xmin = 0
    min1 = 100000
    x = []
    y = []
    for a, b in zip(x_1, y_1):
        if float(b) < min1:
            xmin = a
            min1 = b
    if debug:
        print(f"minimum found at r = %g E = %g for %s-%s" % ( xmin, min1, noble1, noble2 ))
    for a, b in zip(x_1, y_1):
        if polarizable:
            b -= eldict(a) 
        if b < upper and abs(b) > abs(lower*min1):
            x.append(a)
            y.append(b)
    datapoints.update({f"{noble1}{noble2}": f"{[x, y]}"})
    return x, y

def parse_args():
    parser  = argparse.ArgumentParser(description="""
    The script fits parameters of potentials of noble gas dimers, loops through combination rules to 
    reproduce parameters of heterodimer, writes and plots results, that is _fit tables with fitted 
    parameters and RMSE of fit, png files of potentials (below), comb.csv files with all combinations 
    of combination rules and their RMSE for respective pari and total, and a single best_comb file 
    with only best of them for each potential. Required architesture is: <level of theory>/(csv files 
    named <nobleatom1>#<nobleatom2>.csv), csv files has two columns, e.g. distance and energy 
    (of desired unit) separated by <,>. The plots have reference points to fit in blue dots, the 
    potential fitted throught them in blue dash, the parameters for the potential reconstructed 
    from the best combining rule in green and for comparison, also parameters reconstructed when 
    taking arithmetic (yellow) or geometric means (red) for every parameter. The latter can be found 
    in combination tables under <f1> or <f2> for all parameters, for all-geometric or all-arithmetic rules")
    """)
    
    deflot = "CCSDT_CBS"
    parser.add_argument("-lot", "--lot", help="Level of theory, default "+deflot, type=str, default=deflot)
    parser.add_argument("-up","--upper", help="upper cutoff of energy, default 20 kJ/mol. Will take only data points with energies below this value.", default=20, type=int)
    parser.add_argument("-down","--lower", help="lower cutoff. Fraction of the minimum energy, defaults to 0.1, will take only energies with absolute value greater than this percentage of energy well depth", default=0.1, type=float)
    parser.add_argument("-blink", "--blink", help="Time between frames of plotted potential, defaults to 0.25 sec. Increase to prevent epileptic seizures. If 0, plots are written only.", type=float, default=0.25)
    parser.add_argument("-v", "--verbose", help="Print more stuff", action="store_true")
    nlist = [ "He", "Ne", "Ar", "Kr", "Xe" ]
    parser.add_argument("-elem", "--elements", nargs="+", help="List of noble gas elements to consider, default all", default=nlist)
    parser.add_argument("-pol", "--polarizable", help="Subtract electrostatic energy before fitting", action="store_true")
    args = parser.parse_args()
    if args.verbose:
        debug = True
    return args

def get_Re_De(noble1:str, noble2:str, table:dict):
    elem = short[noble1] + "-" + short[noble2]
    #if not elem in table:
    #    elem = short[noble2] + "-" + short[noble1]
    Re = table[elem]["Re"]*Bohr
    De = table[elem]["De"]*Hartree
    return Re, De
    
def dofit(lot:str, upper:float, lower:float, noble1:str, noble2:str, potentials:dict, table:dict, polarizable:bool)->dict:
#Fitting of each potential. e, g, s, d for epsilon, gamma, sigma, delta
    x, y = read_data(noble1, noble2, lot, upper, lower, polarizable)
    x_fit = np.copy(x)
    rounddecimal = 3
    abs_sigma = True
    maxsteps  = 10000000 
    myfit = {}
    for pot in potentials.keys():
        params     = {}
        nparams = len(potentials[pot]["min"])
        pindex  = 0
        if 2 == nparams:
            popt, pcov = curve_fit(lambda t, e, s: ( eval(pot)(t, e, s) ), x, y, bounds=(potentials[pot]["min"], potentials[pot]["max"]), maxfev=maxsteps, absolute_sigma=abs_sigma )
        elif 3 == nparams:
            popt, pcov = curve_fit(lambda t, e, s, g: ( eval(pot)(t, e, s, g) ), x, y, bounds=(potentials[pot]["min"], potentials[pot]["max"]), maxfev=maxsteps, absolute_sigma=abs_sigma )
        elif 4 == nparams:
            popt, pcov = curve_fit(lambda t, e, s, g, d: ( eval(pot)(t, e, s, g, d) ), x, y, bounds=(potentials[pot]["min"], potentials[pot]["max"]), maxfev=maxsteps, absolute_sigma=abs_sigma )
        elif 7 == nparams:
            # The last two parameters are not real parameters, this needs a bit of a hack
            Re, De = get_Re_De(noble1, noble2, table)
            params["Re"] = Re
            params["De"] = De
            popt, pcov = curve_fit(lambda t, A, b, c6, c8, c10: ( eval(pot)(t, A, b, c6, c8, c10, Re, De) ), x, y, 
                                   bounds=(potentials[pot]["min"][:5], potentials[pot]["max"][:5]), maxfev=maxsteps, absolute_sigma=abs_sigma )
            popt = np.append(popt, [Re, De], axis=0)
        myfit[pot] = {}
        for i in range(len(potentials[pot]["label"])):
            mylabel = potentials[pot]["label"][i]
            params[mylabel] = popt[i]
            myfit[pot][mylabel] = popt[i]
        myrms     = RMSE(eval(pot), x, y, params)
        myfit[pot]["RMSE"] = myrms
    return myfit

def low_combine(myfit: dict, nn1:str, nn2:str, nn12:str, myrule:dict, pot:str)->dict:
    myparm = {}
    if "epsilon" in myfit[nn1][pot] and "sigma" in myfit[nn1][pot]:
        eps1 = myfit[nn1][pot]["epsilon"]
        eps2 = myfit[nn2][pot]["epsilon"]
        sig1 = myfit[nn1][pot]["sigma"]
        sig2 = myfit[nn2][pot]["sigma"]
        gam1 = None
        gam2 = None
        if "gamma" in myfit[nn1][pot]:
            gam1 = myfit[nn1][pot]["gamma"]
            gam2 = myfit[nn2][pot]["gamma"]
        if myrule["epsilon"] == "f8":
            myparm["epsilon"] = f8(eps1, eps2, sig1, sig2)
        else:
            myparm["epsilon"] = eval(myrule["epsilon"])(eps1, eps2)
        if myrule["sigma"] == "f11":
            myparm["sigma"] = f11(eps1, eps2, sig1, sig2)
        elif myrule["sigma"] != "f5":
            myparm["sigma"] = eval(myrule["sigma"])(sig1, sig2)
        elif gam1 and gam2:
            # Rule f5 can only be applied if we have gamma parameters
            myparm["sigma"] = f5(eps1, eps2, gam1, gam2, sig1, sig2)
        else:
            sys.exit("Inconsistency for pot %s sigma %s epsilon %s" % (pot, myrule["sigma"], myrule["epsilon"]))

        if gam1 and gam2:
            if myrule["gamma"] == "f10":
                myparm["gamma"] = f10(gam1, gam2, sig1, sig2)
            else:
                myparm["gamma"] = eval(myrule["gamma"])(gam1, gam2)
            if "delta" in myfit[nn1][pot]:
                dlt1 = myfit[nn1][pot]["delta"]
                dlt2 = myfit[nn2][pot]["delta"]
                myparm["delta"] = eval(myrule["delta"])(dlt1, dlt2)
    elif "A" in myfit[nn1][pot] and "b" in myfit[nn1][pot] and "C6" in myfit[nn1][pot]:
        myparm["A"]  = eval(myrule["A"])(myfit[nn1][pot]["A"], myfit[nn2][pot]["A"])
        myparm["b"]  = eval(myrule["b"])(myfit[nn1][pot]["b"], myfit[nn2][pot]["b"])
        myparm["C6"] = eval(myrule["C6"])(myfit[nn1][pot]["C6"], myfit[nn2][pot]["C6"])
        if "C8" in myfit[nn1][pot] and "C10" in myfit[nn1][pot]:
            myparm["C8"] = eval(myrule["C8"])(myfit[nn1][pot]["C8"], myfit[nn2][pot]["C8"])
            myparm["C10"] = eval(myrule["C10"])(myfit[nn1][pot]["C10"], myfit[nn2][pot]["C10"])
            myparm["Re"] = myfit[nn12][pot]["Re"]
            myparm["De"] = myfit[nn12][pot]["De"]
            
    return myparm

def eval_comb_rules(myfit:dict, potentials:dict, noble1:str, noble2:str, addRule:bool,
                    upper:float, lower:float, lot:str, polarizable:bool):
    #looping throught dictionaries for each parameter and combining relations. for each, homodimer parameters attained from the fit are used (global parameters with a name looping through following...
    x, y = read_data(noble1, noble2, lot, upper, lower, polarizable)
    nn1  = f"{noble1}{noble1}"
    nn2  = f"{noble2}{noble2}"
    if not nn1 in myfit or not nn2 in myfit:
        return
    nn12 = f"{noble1}{noble2}"
    for pot in potentials.keys():
        if debug:
            print("*** Evaluating combination rules for %s-%s and potential %s ***" % (noble1, noble2, pot) )
        count_combrules = 0
        allrules  = []
        if "epsilon" in myfit[nn1][pot] and "sigma" in myfit[nn1][pot]:
            epsilons = [ "f1", "f2", "f3", "f4", "f6", "f7", "f8", "f9", "f12", "f14", "f15", "f16", "f17" ] 
            gammas   = [ "f1", "f2", "f3", "f4", "f6", "f7", "f9", "f10", "f12", "f14" ]
            sigmas   = [ "f1", "f2", "f3", "f4", "f5", "f6", "f7", "f9", "f11", "f12", "f14" ]
            deltas   = [ "f1", "f2", "f3", "f4", "f6", "f7", "f9", "f12", "f14" ]
            for eps in epsilons:
                myrule = { "epsilon": eps }
                for sig in sigmas:
                    if "f5" == sig and not "gamma" in myfit[nn1][pot]:
                        continue
                    myrule["sigma"] = sig

                    if not "gamma" in myfit[nn1][pot]:
                        allrules.append(myrule.copy())
                    else:
                        # Loop over gamma as well
                        for gam in gammas:
                            myrule["gamma"] = gam
                            if not "delta" in myfit[nn1][pot]:
                                allrules.append(myrule.copy())
                            else:
                                for dlt in deltas:
                                    myrule["delta"] = dlt
                                    allrules.append(myrule.copy())
        elif "A" in myfit[nn1][pot] and "b" in myfit[nn1][pot] and "C6" in myfit[nn1][pot]:
            mysigmas   = [ "f1", "f2", "f3", "f4", "f6", "f7", "f9", "f12" ] 
            myepsilons = [ "f1", "f2", "f3", "f4", "f6", "f7", "f9", "f12" ]
            mygammas   = [ "f1", "f2", "f3", "f4", "f6", "f7", "f9", "f12" ]
            myrule     = {}
            for A in myepsilons:
                myrule["A"] = A
                for b in mygammas:
                    myrule["b"] = b
                    for C6 in mysigmas:
                        myrule["C6"] = C6
                        if not "TT" == pot:
                            allrules.append(myrule.copy())
                        else:
                            for C8 in mysigmas:
                                myrule["C8"] = C8
                                for C10 in mysigmas:
                                    myrule["C10"] = C10
                                    allrules.append(myrule.copy())
                        

        else:
            sys.exit("Incomprehensible potential "+str(potential))
        # Check rules
        if allrules[0] == allrules[1]:
            print("Rules mixed up for %s" % pot)
        for myrule in allrules:
            testparm = low_combine(myfit, nn1, nn2, nn12, myrule, pot)
            if "TT" == pot:
                # Use these two parameters from Wei
                testparm["Re"] = myfit[nn12][pot]["Re"]
                testparm["De"] = myfit[nn12][pot]["De"]
            if not nn12 in potentials[pot]["df"]:
                potentials[pot]["df"][nn12] = []
            potentials[pot]["df"][nn12].append(RMSE(eval(pot), x, y, testparm))
            if addRule:
                potentials[pot]["df"]["rule"].append(myrule)

def write_rules(workdir:str, potentials:dict):
    for pot in potentials.keys():
        pairs = []
        for k in potentials[pot]["df"].keys():
            if "rule" != k and "RMSE" != k:
                pairs.append(k)
        csv = f"{workdir}/{pot}_comb.csv"
        count = 0
        with open(csv, "w") as outf:
            labels = list(potentials[pot]["df"]["rule"][0].keys())
            for l in labels:
                outf.write(",%s" % l)
            outf.write(",RMSE")
            for pair in pairs:
                outf.write(","+pair)
            outf.write("\n")
            for i in range(len(potentials[pot]["df"]["rule"])):
                outf.write("%d" % ( i ) )
                for l in labels:
                    outf.write(",%s" % potentials[pot]["df"]["rule"][i][l] )
                outf.write(",%g" % ( potentials[pot]["df"]["RMSE"][i] ) )
                for pair in pairs:
                    outf.write(",%g" % ( potentials[pot]["df"][pair][i] ))
                outf.write("\n")

def compute_rms(rrr:list):
    r2sum = 0
    for r in rrr:
        r2sum += r**2
    if r2sum > 0:
        return math.sqrt(r2sum/len(rrr))
    else:
        return 0

def save_tex_pot(potentials:dict, myfit:dict, noble:list, workdir:str, comb_df_best:dict, polarizable:bool):
    for pot in potentials.keys():
        if not pot in comb_df_best:
            continue
        bestrule = comb_df_best[pot]["rule"]
        myrule   = potentials[pot]["df"]["rule"][bestrule]
        mycomb = "Rules: "
        for i in range(len(potentials[pot]["label"])):
            if potentials[pot]["label"][i] in myrule:
                mycomb += (" %s:\\ref{%s}" % ( potentials[pot]["texlabel"][i], myrule[potentials[pot]["label"][i]]))
        texfit = open(f"{workdir}/{pot}_fit.tex", "w")
        texfit.write("\\begin{table}[ht]\n\\centering\n\\caption{Force field parameters from fitting a %s potential to quantum chemistry data at the %s level of theory.%s Energy cut-off %g kJ/mol for repulsion and %.0f\\%% of the well-depth at the long range. %s}\n" % ( potentials[pot]["name"], thislot, pol_string(polarizable), args.upper, 100*args.lower, mycomb ))
        polstr = ""
        if polarizable:
            polstr = "_pol"
        texfit.write("\\label{fit_%s_%s%s}\n"  % (pot, args.lot, polstr))
        texfit.write("\\begin{tabular}{lc")
        for i in range(len(potentials[pot]["label"])+1):
            texfit.write("c")
        texfit.write("}\n")
        texfit.write("\\hline\n")
        texfit.write("Pair&")
        for i in range(len(potentials[pot]["label"])):
            texfit.write("& %s" % ( potentials[pot]["texlabel"][i]))
        texfit.write(f"& RMSE\\\\\n")
        for i in range(len(potentials[pot]["unit"])):
            texfit.write("& %s" % potentials[pot]["unit"][i])
        texfit.write(f"&kJ/mol\\\\\n\\hline\n")
        count    = { "fit": [], "comb": [] }
        bestrule = comb_df_best[pot]["rule"]
        myrule   = potentials[pot]["df"]["rule"][bestrule]
        for n1 in range(len(noble)):
            noble1 = noble[n1]
            for n2 in range(n1, len(noble)):
                noble2 = noble[n2]
                nn1    = f"{noble1}{noble1}"
                nn2    = f"{noble2}{noble2}"
                nn12   = f"{noble1}{noble2}"
                texfit.write("%s & %s " % ( f"{short[noble1]}-{short[noble2]}", "fit" ))
                for i in range(len(potentials[pot]["label"])):
                    mylabel = potentials[pot]["label"][i]
                    texfit.write("& %g" % myfit[nn12][pot][mylabel] )
                texfit.write("& %.4f\\\\\n" % myfit[nn12][pot]["RMSE"])
                count["fit"].append(myfit[nn12][pot]["RMSE"])
                if n2 > n1:
                    texfit.write("  & %s " % ( "comb" ))
                    myparm = low_combine(myfit, nn1, nn2, nn12, myrule, pot)
                    for i in range(len(potentials[pot]["label"])):
                        mylabel = potentials[pot]["label"][i]
                        if not ("RMSE" == mylabel):
                            texfit.write("& %g" % myparm[mylabel] )
                    myrms = potentials[pot]["df"][nn12][bestrule]
                    texfit.write("& %.4f\\\\\n" % myrms)
                    count["comb"].append(myrms)
        
        for cf in [ "fit", "comb" ]:
            if len(count[cf]) > 0:
                texfit.write("RMSE & %s" % cf)
                for i in range(len(potentials[pot]["label"])):
                    texfit.write(f" &")
                rmse = compute_rms(count[cf])
                texfit.write(" & %.4f\\\\\n" % rmse)
        texfit.write("\\hline\n\\end{tabular}\n\end{table}\n")
        texfit.close()

def save_csv_pot(potentials:dict, myfit:dict, noble:list, workdir:str, comb_df_best:dict):
    for pot in potentials.keys():
        if not pot in comb_df_best:
            continue
        bestrule = comb_df_best[pot]["rule"]
        myrule   = potentials[pot]["df"]["rule"][bestrule]
        with open(("%s/%s.csv" % ( workdir, pot)), "w") as outf:
            outf.write(",type")
            for lab in potentials[pot]["label"]:
                outf.write(",%s" % lab)
            outf.write("\n")
            for n1 in range(len(noble)):
                for n2 in range(n1, len(noble)):
                    nn1    = f"{noble[n1]}{noble[n1]}"
                    nn2    = f"{noble[n2]}{noble[n2]}"
                    nn12   = f"{noble[n1]}{noble[n2]}"
                    outf.write("%s-%s,fit" % ( short[noble[n1]], short[noble[n2]] ))
                    for ilab in range(len(potentials[pot]["label"])):
                        lab = potentials[pot]["label"][ilab]
                        if nn12 in myfit:
                            mymax = potentials[pot]["max"][ilab]
                            mymin = potentials[pot]["min"][ilab]
                            toler = (mymax-mymin)*1e-3
                            value = myfit[nn12][pot][lab]
                            if abs(value-mymin) < toler and mymin > 0:
                                print("Parameter %s %s close to minimum value %g" % ( lab, pot, mymin ))
                            if abs(value-mymax) < toler:
                                print("Parameter %s %s close to maximum value %g" % ( lab, pot, mymax ))
                            outf.write(",%g" % value)
                    outf.write("\n")
                    if n1 != n2 and nn12 in myfit:
                        outf.write("%s-%s,comb" % ( short[noble[n1]], short[noble[n2]] ))
                        myparm = low_combine(myfit, nn1, nn2, nn12, myrule, pot)
                        for i in range(len(potentials[pot]["label"])):
                            mylabel = potentials[pot]["label"][i]
                            if not ("RMSE" == mylabel):
                                outf.write(",%g" % myparm[mylabel] )
                        outf.write("\n")


def save_table2(potentials:dict, myfit:dict, noble:list, mylot:str, workdir:str, args):
    rmsepot = {}
    with open(workdir+"/table_fit.tex", "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{Root mean square error (kJ/mol) from analytical fit of potentials to quantum chemical data at the %s level of theory.%s Repulsive energy cut-off %g kJ/mol, long range cut at %d\\%% of the well-depth.}\n" % (latexit(mylot), pol_string(args.polarizable), args.upper, int(100*args.lower)))
        polstr = ""
        if args.polarizable:
            polstr = "_pol"
        outf.write("\\label{fit_%s%s}\n" % (args.lot, polstr))
        outf.write("\\begin{tabular}{l")
        for k in range(len(potentials.keys())):
            outf.write("c")
        outf.write("}\n")
        outf.write("\\hline\n")
        outf.write("Pair")
        for k in potentials.keys():
            outf.write("&%s" % potentials[k]["latex"])
        outf.write("\\\\\n\\hline\n")
        average = {}
        for pot in potentials.keys():
            average[pot] = 0
        npair   = 0
        for n1 in range(len(noble)):
            for n2 in range(n1, len(noble)):
                mycomb = f"{noble[n1]}{noble[n2]}"
                if mycomb in myfit:
                    outf.write("%s-%s" % ( short[noble[n1]], short[noble[n2]] ) )
                    for pot in potentials.keys():
                        thisrmse = myfit[mycomb][pot]["RMSE"]
                        outf.write("&%.3f" % thisrmse)
                        average[pot] += thisrmse**2
                    outf.write("\\\\\n")
                    npair += 1
        outf.write("\\hline\n")
        outf.write("RMSE")
        for pot in potentials.keys():
            if npair > 0:
                rmsepot[pot] = np.sqrt(average[pot]/npair)
            else:
                rmsepot[pot] = 0
            outf.write("&%.3f" % ( rmsepot[pot] ))
        outf.write("\\\\\n")
        outf.write("\\hline\n\\end{tabular}\n\\end{table}\n")
    return rmsepot

def save_table3(potentials:dict, myfit:dict, noble:list, mylot:str, workdir:str, args, comb_df_best:dict):
    with open(workdir+"/table_comb.tex", "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{Root mean square error (kJ/mol) from of combination-rule based interaction energies from heterodimer potentials at the %s level of theory.%s Repulsive energy cut-off %g kJ/mol, long range cut at %d\\%% of the well-depth.}\n" % (latexit(mylot), pol_string(args.polarizable), args.upper, int(100*args.lower)))
        polstr = ""
        if args.polarizable:
            polstr = "_pol"
        outf.write("\\label{comb_%s%s}\n" % (args.lot, polstr))
        width = "14mm"
        outf.write("\\begin{tabular}{l")
        for k in range(len(potentials.keys())):
            outf.write("p{%s}" % width)
        outf.write("}\n")
        outf.write("\\hline\n")
        outf.write("Pair")
        for k in potentials.keys():
            outf.write("&%s" % potentials[k]["latex"])
        outf.write("\\\\\n\n")
        nrow = 5
        outf.write("\multirow{%d}{%s}{Eq.}" % ( nrow, width ))
        for pot in potentials.keys():
            bestrule = comb_df_best[pot]["rule"]
            myrule   = potentials[pot]["df"]["rule"][bestrule]
            outf.write("&\\multirow{%d}{%s}{" % (nrow, width))
            for i in range(len(potentials[pot]["label"])):
                if potentials[pot]["label"][i] in myrule:
                    outf.write(" %s:\\ref{%s}" % ( potentials[pot]["texlabel"][i], myrule[potentials[pot]["label"][i]]))
            outf.write("}")
        outf.write("\\\\\n\\\\\n\\\\\n\\\\\n")
        outf.write("\\\\\n")
        
        outf.write("\\hline\n")
        for n1 in range(len(noble)):
            for n2 in range(n1+1, len(noble)):
                mycomb = f"{noble[n1]}{noble[n2]}"
                if mycomb in myfit:
                    outf.write("%s-%s" % ( short[noble[n1]], short[noble[n2]] ) )
                    for pot in potentials.keys():
                        bestrule = comb_df_best[pot]["rule"]
                        thisrmse = potentials[pot]["df"][mycomb][bestrule]
                        outf.write("&%.3f" % thisrmse)
                    outf.write("\\\\\n")
        outf.write("\\hline\n")
        outf.write("RMSE")
        for pot in potentials.keys():
            outf.write("&%.3f" % ( comb_df_best[pot]["RMSE"] ))
        outf.write("\\\\\n")
        outf.write("\\hline\n\\end{tabular}\n\\end{table}\n")
        
def save_table_treshold(potentials:dict, rmsepot:dict, mylot:str, workdir:str, args, comb_df_best:dict):
    with open(workdir+"/comb_best.tex", "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\caption{Fitting RMSE and best combination rules on datasets with different upper energy thresholds. Level of theory used %s.%s}\n" % (latexit(mylot), pol_string(args.polarizable)) )
        outf.write("\\label{treshold}\n")
        outf.write("\\begin{tabular}{l")
        for pot in potentials.keys():
            outf.write("c")
        outf.write("}\n")
        for pot in potentials.keys():
            outf.write("&%s" % ( potentials[pot]["latex"] ))
        outf.write("\\\\\nParam.")
        for pot in potentials.keys():
            bestrule = comb_df_best[pot]["rule"]
            myrule   = potentials[pot]["df"]["rule"][bestrule]
            outf.write(" &")
            for i in range(len(potentials[pot]["label"])):
                if potentials[pot]["label"][i] in myrule:
                    outf.write(" %s" % ( potentials[pot]["texlabel"][i] ))
        outf.write("\\\\\n\\hline\n&\multicolumn{%d}{c}{Treshold %g kJ/mol}\\\\\n" % ( len(potentials.keys()), args.upper ))
        outf.write("RMSE fit")
        for pot in potentials.keys():
            outf.write("&%.3f" % ( rmsepot[pot] ))
        outf.write("\\\\\n")
        outf.write("RMSE comb")
        for pot in potentials.keys():
            outf.write("&%.3f" % ( comb_df_best[pot]["RMSE"] ))
        outf.write("\\\\\nRule")
        for pot in potentials.keys():
            bestrule = comb_df_best[pot]["rule"]
            myrule   = potentials[pot]["df"]["rule"][bestrule]
            outf.write("&")
            for i in range(len(potentials[pot]["label"])):
                if potentials[pot]["label"][i] in myrule:
                    outf.write(" \\ref{%s}" % ( myrule[potentials[pot]["label"][i]]))
        outf.write("\\\\\n")
        outf.write("\\hline\n\\end{tabular}\n\\end{table}\n")

def get_potentials()->dict:
    epsmin = 0 # kJ/mol
    epsmax = 10 # kJ/mol
    sigmin = 2 # Angstrom
    sigmax = 20 # Angstrom
    potentials = { "LJ_12_6": { "name": "Lennard-Jones 12-6", "latex": "LJ12-6",
                                "label": [ "epsilon", "sigma" ], "unit": [ "kJ/mol", "{\AA}" ],
                                "texlabel": [ "$\epsilon$", "$\sigma$" ], 
                                "min": [epsmin, sigmin], "max": [epsmax, sigmax] },
                   "LJ_8_6":  { "name": "Lennard-Jones 8-6", "latex": "LJ8-6",
                                "label": [ "epsilon", "sigma" ], "unit": [ "kJ/mol", "{\AA}" ],
                                "texlabel": [ "$\epsilon$", "$\sigma$" ], 
                                "min": [epsmin, sigmin],          "max": [epsmax, 12*sigmax] },
                   "WBH":     { "name": "Wang-Buckingham", "latex": "WBH",
                                "label": [ "epsilon", "sigma", "gamma" ], "unit": [ "kJ/mol", "{\AA}", "-" ],
                                "texlabel": [ "$\epsilon$", "$\sigma$", "$\gamma$" ], 
                                "min": [epsmin, sigmin, 6],      "max": [epsmax, sigmax, 40] },
                   "MBH":     { "name": "modified Buckingham", "latex": "MBH",
                                "label": [ "epsilon", "sigma", "gamma" ], "unit": [ "kJ/mol", "{\AA}", "-" ],
                                "texlabel": [ "$\epsilon$", "$\sigma$", "$\gamma$" ], 
                                "min": [epsmin, sigmin, 6],      "max": [epsmax, sigmax, 30] },
                   "MRS":     { "name": "Morse", "latex": "MRS",
                                "label": [ "epsilon", "sigma", "gamma" ], "unit": [ "kJ/mol", "{\AA}", "1/{\AA}" ],
                                "texlabel": [ "$\epsilon$", "$\sigma$", "$\gamma$" ], 
                                "min": [epsmin, sigmin, 1],       "max": [epsmax, sigmax, 4]  },
                   "BHA":     { "name": "Buckingham", "latex": "BHA",
                                "label": [ "epsilon", "sigma", "gamma" ], "unit": [ "kJ/mol", "kJ/mol {\AA}$^6$", "1/{\AA}" ],
                                "texlabel": [ "A", "C$_6$", "b" ],
                                "min": [3e4, 100, 2 ],           "max": [2e6, 35000, 6 ] },
                   "GBH":     { "name": "generalized Buckingham", "latex": "GBH",
                                "label": [ "epsilon", "sigma", "gamma", "delta" ], "unit": [ "kJ/mol", "{\AA}", "-", "-" ],
                                "texlabel": [ "$\epsilon$", "$\sigma$", "$\gamma$", "$\delta$" ], 
                                "min": [epsmin, sigmin, 10,  0], "max": [epsmax, sigmax, 50, 40] },
                   "LJ_14_7": { "name": "Buffered 14-7", "latex": "LJ14-7",
                                "label": [ "epsilon", "sigma", "gamma", "delta" ], "unit": [ "kJ/mol", "{\AA}", "-", "-" ],
                                "texlabel": [ "$\epsilon$", "$\sigma$", "$\gamma$", "$\delta$" ], 
                                "min": [epsmin, sigmin, 0,  0],   "max": [epsmax, sigmax, 0.8, 2] },
                   "TT":      { "name": "Tang-Toennies", "latex": "TT",
                                "label": [ "A", "b", "C6", "C8", "C10", "Re", "De" ], 
                                "texlabel": [ "A", "b", "C$_6$", "C$_8$", "C$_{10}$", "Re", "De"  ], 
                                "unit": [ "-", "-", "-", "-", "-", "{\AA}", "kJ/mol" ],
                                "min": [2e5, 8, 0.0, 0.0, 0.0, 0, 0], "max": [1e6, 20, 5, 4, 8, 0, 0] }
                }
    return potentials

if __name__ == "__main__":
    args   = parse_args()
    noble  = []
    for elem in args.elements:
        if elem in short:
            noble.append(elem)
        else:
            for k in short.keys():
                if short[k] == elem:
                    noble.append(k)
    potentials = get_potentials()
    # Full level of theory for printing
    lots = { "CCSDT_CBS": "CCSD(T)/CBS", "CCSDTQ_CBS": "CCSDT(Q)/CBS", "SAPT_CBS_a_3.8": "SAPT/CBS", "SAPT_TZ": "SAPT/TZ", "SAPT_QZ": "SAPT/QZ", "CCSDT_TZ": "CCSDT/TZ", "CCSDT_QZ": "CCSDT/QZ" }
    thislot = args.lot
    if thislot in lots:
        thislot = lots[thislot]
    myfit = {}
    datapoints = {}
    workdir = f"RESULTS_{args.lot}_{args.upper}_{args.lower}"
    if args.polarizable:
        workdir += "_pol"
    os.makedirs(workdir, exist_ok=True)

    #cleanup of previous runs, Definition of some variables for the data storage and csv conversion later
    for potential in potentials.keys():
        for fn in [ f"{workdir}/{potential}_fit", f"{workdir}/{potential}_comb.csv", f"best_comb" ]:
            if os.path.exists(fn):
                os.remove(fn)

        potentials[potential]["df"] = {}
        potentials[potential]["df"]["rule"] = []
        potentials[potential]["df"]["RMSE"] = []

    print("Fitting curves for each potential and each atom pair")
    count = 0
    list_combinations = ["pair"]
    table  = get_table1("lmwei2017_table1.csv")
    table4 = get_table1("lmwei2017_table4.csv")
    for t4 in table4.keys():
        table[t4] = table4[t4]
    for n1 in range(len(noble)):
        for n2 in range(n1, len(noble)):
            csvfile = f"{args.lot}/{noble[n1]}#{noble[n2]}.csv"
            if os.path.exists(csvfile):
                count += 1
                mycomb = f"{noble[n1]}{noble[n2]}"
                list_combinations.append(mycomb)
                myfit[mycomb] = dofit(args.lot, args.upper, args.lower, noble[n1], noble[n2], potentials, table, args.polarizable)
            else:
                print("File %s not found" % csvfile)

    # Now dump the manuscript table 2
    rmsepot = save_table2(potentials, myfit, noble, thislot, workdir, args)

    print("...and combining fitted parameters of homodimers with a permutation of rules for each of the parameters.")
    ncombrule = 0
    for n1 in range(len(noble)):
        for n2 in range(n1+1, len(noble)):
            if os.path.exists(f"{args.lot}/{noble[n1]}#{noble[n2]}.csv"):
                eval_comb_rules(myfit, potentials, noble[n1], noble[n2],
                                (ncombrule == 0), args.upper, args.lower, args.lot, args.polarizable)
                ncombrule += 1
                
    comb_df_best = {}
    if ncombrule > 0:
        print("Writing the final RMSE and parsing combination rules to tables")
        for potential in potentials.keys():
            if debug:
                print("Finding best combination rule for %s" % potential)
            bestrule = None
            min_rmse = 1000000
            pairs = []
            for k in potentials[potential]["df"].keys():
                if "rule" != k and "RMSE" != k:
                    pairs.append(k)
            ncomb = len(potentials[potential]["df"]["rule"])
            print("Evaluating %d combination rules for %s." % ( ncomb, potentials[potential]["name"]))
            for ii in range(ncomb):
                sum2 = 0
                counter = 0
                for pair in pairs:
                    counter += 1
                    sum2 += potentials[potential]["df"][pair][ii]**2
                rule_rmse = np.sqrt(sum2/counter)
                if rule_rmse < min_rmse:
                    min_rmse = rule_rmse
                    bestrule = ii
                potentials[potential]["df"]["RMSE"].append(rule_rmse)
            comb_df_best[potential] = { "RMSE": min_rmse, "rule": bestrule }
        write_rules(workdir, potentials)
    # Now dump the manuscript table 3
    if comb_df_best:
        save_table3(potentials, myfit, noble, thislot, workdir, args, comb_df_best)
        save_table_treshold(potentials, rmsepot, thislot, workdir, args, comb_df_best)
    # Now dump the fitted parameters in a machine readable file as well as a Latex file
    save_csv_pot(potentials, myfit, noble, workdir, comb_df_best)
    save_tex_pot(potentials, myfit, noble, workdir, comb_df_best, args.polarizable)
    
    print("Now plotting the png files, taking the analytical fit, geometric, arithmetic and the individual best rule set")
    for n1 in range(len(noble)):
        for n2 in range(n1, len(noble)):
            if not os.path.exists(f"{args.lot}/{noble[n1]}#{noble[n2]}.csv"):
                continue
            x, y = read_data(noble[n1], noble[n2], args.lot, args.upper, args.lower, args.polarizable)
            x_fit = np.copy(x)
            for potential in potentials.keys():
                ax = pyplot.axes()

                comb_fit = myfit[f"{noble[n1]}{noble[n2]}"][potential]
                parm1    = myfit[f"{noble[n1]}{noble[n1]}"][potential]
                parm2    = myfit[f"{noble[n2]}{noble[n2]}"][potential]
                if not len(comb_fit) == len(parm1) or not len(comb_fit) == len(parm2):
                    sys.exit("len(comb_fit) = %d, len(parm1) = %d, len(parm2) = %d" % 
                             ( len(comb_fit), len(parm1), len(parm2) ))
                comb_geo  = {}
                comb_ari  = {}
                comb_best = {}
                if potential in comb_df_best:
                    bestrule   = comb_df_best[potential]["rule"]
                    rule       = potentials[potential]["df"]["rule"][bestrule]
                    label_best = ""
                    for i in comb_fit:
                        comb_geo[i]  = f1(parm1[i], parm2[i])
                        comb_ari[i]  = f2(parm1[i], parm2[i])
                        if "epsilon" == i:
                            if "f8" == rule[i]:
                                comb_best[i] = f8(parm1["epsilon"], parm2["epsilon"], parm1["sigma"], parm2["sigma"])
                            else:
                                comb_best[i] = eval(rule[i])(parm1[i], parm2[i])
                        elif "sigma" == i:
                            if "f5" == rule[i]:
                                comb_best[i] = f5(parm1["epsilon"], parm2["epsilon"], parm1["gamma"], parm2["gamma"], parm1["sigma"], parm2["sigma"])
                            elif "f11" == rule[i]:
                                comb_best[i] = f11(parm1["epsilon"], parm2["epsilon"], parm1["sigma"], parm2["sigma"])
                            else:
                                comb_best[i] = eval(rule[i])(parm1[i], parm2[i])
                        elif "gamma" == i:
                            if "f10" == rule[i]:
                                comb_best[i] = f10(parm1["gamma"], parm2["gamma"], parm1["sigma"], parm2["sigma"])
                            else:
                                comb_best[i] = eval(rule[i])(parm1[i], parm2[i])
                        elif "Re" == i or "De" == i or "RMSE" == i:
                            # These parameters are not fitted at all
                            comb_best[i] = comb_fit[i]
                        else:
                            comb_best[i] = eval(rule[i])(parm1[i], parm2[i])
                        if "Re" != i and "De" != i and "RMSE" != i:
                            if len(label_best) > 0:
                                label_best += " "
                            label_best += ( "%s:%s" % ( i, rule[i] ))
                    if debug:
                        print("Best rule for %s: %s" % ( potential, label_best))

                yfunc = compute(potential, x_fit, comb_fit)
                ax.plot(x_fit, yfunc, label=f'Fit of {potentials[potential]["name"]}', color="blue", linestyle='dotted', linewidth=3)
                if comb_geo:
                    ygeo  = compute(potential, x_fit, comb_geo)
                    ax.plot(x_fit, ygeo, label=f'All param. geom.', color="red", linestyle='dotted', linewidth=3)
                if comb_ari:
                    yari  = compute(potential, x_fit, comb_ari)
                    ax.plot(x_fit, yari,  label=f'All param. arit.', color="yellow", linestyle='dotted', linewidth=3)
                if comb_best:
                    ybest = compute(potential, x_fit, comb_best)
                    ax.plot(x_fit, ybest, label=f'{label_best}', color="green", linestyle='dotted', linewidth=3)
                # plotting the most successful combination rule combination for all elements and respective potential
                ax.scatter(x, y, label=f"{short[noble[n1]]}-{short[noble[n2]]}", color="blue")
                ax.set_ylabel('Energy (kJ/mol)')
                ax.set_xlabel('Distance (Å)')
                ax.legend()
                pyplot.savefig(f'RESULTS_{args.lot}_{args.upper}_{args.lower}/{potential}{noble[n1]}{noble[n2]}.pdf')
                if args.blink != "0":
                    pyplot.show(block=False)
                    pyplot.pause(args.blink)
                pyplot.close()
