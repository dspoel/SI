#!/usr/bin/env python3

import argparse, math, os, sys
from lmwei2017 import Uwei, get_table1
from constants import *
from fit_script import LJ_14_7, WBH, MRS, TT, GBH, LJ_12_6, LJ_8_6, get_potentials, pol_string
import numpy as np
from scipy.interpolate import CubicSpline

B2debug      = False

b2map = { "He-Ne": "helium#neon",
          "He-Ar": "helium#argon",
          "He-Kr": "helium#krypton",
          "He-Xe": "helium#xenon",
          "He-Rn": "helium#radon",
          "Ne-Ne": "neon#neon",
          "Ne-Ar": "neon#argon",
          "Ne-Kr": "neon#krypton",
          "Ne-Xe": "neon#xenon",
          "Ne-Rn": "neon#radon",
          "Ar-Ar": "argon#argon",
          "Ar-Kr": "argon#krypton",
          "Ar-Xe": "argon#xenon",
          "Ar-Rn": "argon#radon",
          "Kr-Kr": "krypton#krypton",
          "Kr-Xe": "krypton#xenon",
          "Kr-Rn": "krypton#radon",
          "Xe-Xe": "xenon#xenon",
          "Xe-Rn": "xenon#radon",
          "Rn-Rn": "radon#radon" }

def calculate_curve(func:str, param:dict, distance_cutoff:float):
    dx  = 0.01
    ir2 = int(distance_cutoff/dx)
    x   = np.zeros(ir2, dtype=float)
    y   = np.zeros(ir2, dtype=float)
    for ir in range(0, ir2):
        r   = dx*(ir+1)
        U   = 0
        if "Wei2017a" == func:
            U = Uwei(r/Bohr, param)*Hartree # *param["De"]
        elif func in [ "LJ_12_6", "LJ_8_6" ]:
            U = eval(func)(r, param["epsilon"], param["sigma"])
        elif func in [ "WBH", "MRS" ]:
            U = eval(func)(r, param["epsilon"], param["sigma"], param["gamma"])
        elif func in [ "LJ_14_7", "GBH" ]:
            U = eval(func)(r, param["epsilon"], param["sigma"],
                           param["gamma"], param["delta"])
        elif "TT" == func:
            U = TT(r, param["A"], param["b"], param["C6"], param["C8"], param["C10"],
                   param["Re"], param["De"])
        else:
            sys.exit("Unknown function %s" % func)
        x[ir] = r
        y[ir] = U
    return x, y

def qmcorr(cs, r:float, reverse_mass:float, Temp:float)->float:
    # Compute QM corrections to second virial according to
    # Bich et al. Molecular Physics 106, No. 6, 20 March 2008, 813-825
    # https://doi.org/10.1080/00268970801964207
    beta   = 1/(Boltz*Temp)
    # Check factor 24 or 12 (between eqn 10 and 11)
    # hbar has unit (kJ/mol) ps then lambda has (kJ/mol) ps^2 / Dalton
    # or (kJ ps^2)/(mol Dalton) -> (1e3 J 1e-24 s^2)/(g) ->
    # J = kg m^2/s^2, therefore
    # lambda has unit (1e3 * 1e3 g m^2 1e-24 / g) -> 1e-18 m^2 -> nm^2
    # Since the integrals are computed in cubic Angstrom, we have to scale
    # the final result by 100 to have square Angstrom unit.
    mylambda = 100*hbar**2*beta*reverse_mass/24
    betaV1 = beta*cs(r, 1)
    betaV2 = beta*cs(r, 2)
    betaV3 = beta*cs(r, 3)
    Bqm1 = betaV1**2
    Bqm2 = -( (6.0/5.0)*betaV2**2 + (12.0/(5.0*r**2))*betaV1**2 + 
              (4.0/(3.0*r))*betaV1**3 - (1.0/6.0)*betaV1**4)
    Bqm3 = ((36.0/35.0)*betaV3**2 + (216.0/(35.0*r**2))*betaV2**2 + (24.0/21.0)*betaV2**3 +
            (24.0/(5.0*r))*betaV1*betaV2**2 + (288.0/(315.0*r**3))*betaV1**3 -
            (6.0/5.0)*betaV1**2*betaV2**2 - (2.0/(15.0*r**2))*betaV1**4 -
            (2.0/(5.0*r))*betaV1**5 + (1.0/30.0)*betaV1**6)
    return [ mylambda*Bqm1, mylambda**2*Bqm2, mylambda**3*Bqm3 ]

def compute_mayer(func:str, param:dict, T:float, reverse_mass:float, distance_cutoff:float)->list:
    mayer = [ [], [], [], [] ]
    beta  = 1/(Boltz*T)
    x, y  = calculate_curve(func, param, distance_cutoff)
    # Do a spline interpolation to compute derivatives
    if len(x) < 4:
        sys.exit("Not enough data points")
    cs    = CubicSpline(x, y)
    for ir in range(x.size):
        r   = x[ir]
        U   = y[ir]
        mmm = [ -1, 0, 0, 0 ]
        if abs(beta*U) < 250:
            uexp = math.exp(-beta*U)
            qmc = qmcorr(cs, r, reverse_mass, T)
            mmm = [ uexp-1, uexp*qmc[0], uexp*qmc[1], uexp*qmc[2] ]

        for k in range(4):
            mayer[k].append((r, mmm[k]))
    return mayer

def qm_mayer(pair:str, qmdata:list, T:float, debug:bool, reverse_mass:float, nintermediate:int):
    mayer = [ [], [], [], [] ]
    beta  = 1/(Boltz*T)
    r0    = qmdata[0][0]
    dx    = 0.1
    ir0   = 0
    while ir0*dx < r0:
        rr = ir0*dx
        mayer[0].append( (rr, -1) )
        for k in range(1,4):
            mayer[k].append( (rr, 0) )
        ir0 += 1

    x   = np.zeros(len(qmdata), dtype=float)
    y   = np.zeros(len(qmdata), dtype=float)
    for ir in range(len(qmdata)):
        x[ir] = qmdata[ir][0]
        y[ir] = qmdata[ir][1]
        
    cs = CubicSpline(x, y)
    for ir in range(len(qmdata)):
        r = x[ir]
        U = y[ir]
        if B2debug:
            print("Data pair %s r = %g U = %g beta*U = %g" % ( pair, r, U, beta*U ))
        mmm = [ -1, 0, 0, 0 ]
        if abs(beta*U) < 250:
            uexp = math.exp(-beta*U)
            qmc  = qmcorr(cs, r, reverse_mass, T)
            mmm = [ uexp-1, uexp*qmc[0], uexp*qmc[1], uexp*qmc[2] ]
        for k in range(4):
            mayer[k].append( ( r, mmm[k] ) )
        if ir < len(qmdata)-1:
            for j in range(1,1+nintermediate):
                if x[ir+1] <= r:
                    sys.exit("Incorrect QM data for %s r[%d] = %g r[%d] = %g" % (pair, ir, r, x[ir+1], r) )
                rr   = r + j*(x[ir+1]-r)/(nintermediate+1)
                U    = cs(rr)
                print("Interpolated pair %s rr = %g UU = %g beta*UU = %g" % ( pair, rr, U, beta*U ))
                if abs(beta*U) < 250:
                    qmc  = qmcorr(cs, rr, reverse_mass, T)
                    uexp = math.exp(-beta*U)
                    mmm = [ uexp-1, uexp*qmc[0], uexp*qmc[1], uexp*qmc[2] ]
                for k in range(4):
                    mayer[k].append( ( rr, mmm[k] ) )

    return mayer

def dump_mayer(mm:list, filenm:str):
    with open(filenm, "w") as outf:
        for m in range(len(mm[0])):
            outf.write("%10g" % mm[0][m][0])
            for k in range(len(mm)):
                outf.write("  %g" % ( mm[k][m][1] ) )
            outf.write("\n")

def read_b2_data(filenm:str, minb2:float)->list:
    b2data = []
    with open(filenm, "r") as inf:
        for line in inf:
            if line.find("#") >= 0:
                continue
            words = line.strip().split(",")
            if len(words) == 2:
                try:
                    T = float(words[0])
                    B2 = float(words[1])
                    if B2 >= minb2:
                        b2data.append({ "T": T, "b2": B2 })
                except ValueError:
                    print("Cannot interpret line '%s' in %s" % line.strip(), filenm)
    return b2data

def integrate_sphere(mayer:list, filenm:str, lengthfac:float)->float:
    integ = 0
    r1    = 0
    m1    = -1
    outf  = None
    if len(filenm) > 0:
        outf = open(filenm, "w")
    for i in range(1,len(mayer)):
        ( r2, m2 ) = mayer[i]
        # Approximate trapezium by y = ax + b (a == slope)
        # then integrate that multiplied by x^2 to get
        # a/4 x^4 + b/3 x^3
        # insert old and new point
        a      = (m2-m1)/(r2-r1)
        b      = m1 - a*r1
        integ += ((a/4)*(r2*r2*r2*r2 - r1*r1*r1*r1) + 
                  (b/3)*(r2*r2*r2 - r1*r1*r1))*lengthfac
        if outf:
            outf.write(" %g  %g\n" % ( r2, integ ))
        m1 = m2
        r1 = r2
    if outf:
        outf.close()
    return 4*math.pi*integ

def read_qm(thislot:str, nlist:list, qmcut:float)->dict:
    qmdict = {}
    for e1 in range(len(nlist)):
        elem1 = nlist[e1]
        for e2 in range(e1, len(nlist)):
            elem2 = nlist[e2]
            elems = ("%s-%s" % ( elem1, elem2 ))
            if not elems in b2map:
                print("Cannot find QM data for pair %s" % elems)
                continue
            data = ("%s/%s.csv" % ( thislot, b2map[elems]))
            if not os.path.exists(data):
                print("Cannot find QM data in %s" % data)
                continue
            qmdict[elems] = []
            with open(data, "r") as inf:
                for line in inf:
                    words = line.split(",")
                    if 2 == len(words):
                        try:
                            r = float(words[0])
                            if (r <= qmcut):
                                qmdict[elems].append( ( r, float(words[1]) ) )
                            else:
                                print("Ignoring data point at r = %g since qmcut = %g" % ( r, qmcut))
                        except ValueError:
                            sys.exit("Error reading %s" % data)
    return qmdict

def latexit(name:str)->str:
    return name.replace("_", "\_")

def calc_b2(table:dict, args, resultdir:str, suffix:str):
    qmdict     = read_qm(args.lot, args.elements, args.qmcut)
    # TODO check masses
    rmse       = {}
    lots = { "CCSDT_CBS": "CCSD(T)/CBS", "CCSDTQ_CBS": "CCSDT(Q)/CBS", "SAPT_CBS_a_3.8": "SAPT/CBS", "SAPT_TZ": "SAPT/TZ", "SAPT_QZ": "SAPT/QZ", "CCSDT_TZ": "CCSDT/TZ", "CCSDT_QZ": "CCSDT/QZ" }
    if args.lot in lots:
        thislot = lots[args.lot]
    else:
        thislot = args.lot
    Tmin = {}
    Tmax = {}
    allfuncs  = [ "qm" ]
    for func in table.keys():
        allfuncs.append(func)
    for e1 in range(len(args.elements)):
        elem1 = args.elements[e1]
        for e2 in range(e1, len(args.elements)):
            elem2 = args.elements[e2]
            pair = ("%s-%s" % ( elem1, elem2 ))
            if not pair in b2map:
                print("Cannot find pair %s in name table" % pair)
                continue
            if not pair in table["WBH"]:
                print("No force field parameters for %s" % pair)
                continue
            b2data = read_b2_data("../SecondVirial/%s.csv" % b2map[pair], args.minb2)
            outf = {}
            outf["abs"] = open("%s/b2classical-%s-%s.xvg" % ( resultsdir, b2map[pair], suffix), "w")
            outf["rel"] = open("%s/b2classical-%s-%s-rel.xvg" % ( resultsdir, b2map[pair], suffix), "w")
            print("Will generate second virials with %d data points for %s" % ( len(b2data), pair ))

            relabs = [ "abs", "rel" ]
            for of in relabs:
                outf[of].write("@ title \"%s\"\n" % ( b2map[pair] ))
                outf[of].write("@ xaxis label \"T (K)\"\n")
            outf["rel"].write("@ yaxis label \"B2(T,computed)-B2(T,experiment)\"\n")
            outf["abs"].write("@ yaxis label \"B2(T)\"\n")
            setnum = { "abs": 0, "rel": 0 }
            outf["abs"].write("@ s0 legend \"Exper.\"\n")
            setnum["abs"] += 1
            for of in relabs:
                for func in allfuncs:
                    outf[of].write("@ s%d legend \"%s\"\n" % ( setnum[of], func) )
                    setnum[of] += 1
            reverse_mass = (mass[elem1]+mass[elem2])/(mass[elem1]*mass[elem2])
            if B2debug:
                print("Reverse mass for %s %s = %g" % ( elem1, elem2, reverse_mass ))
            nmse = {}
            mse  = {}
            for func in allfuncs:
                nmse[func] = 0
                mse[func]  = 0
            # Remove this hack!
            for bb in range(len(b2data)):
                Temp = b2data[bb]["T"]
                if not pair in Tmin:
                    Tmin[pair] = Temp
                else:
                    Tmin[pair] = min(Temp, Tmin[pair])
                if not pair in Tmax:
                    Tmax[pair] = Temp
                else:
                    Tmax[pair] = max(Temp, Tmax[pair])
                b2cl = {}
                b2qm = {}
                lengthfac = AVOGADRO*1e-24
                mydebug   = 0 <= bb
                for func in allfuncs:
                    if "qm" == func:
                        mayer = qm_mayer(pair, qmdict[pair], Temp, mydebug, reverse_mass, args.interpolate)
                    else:
                        mytab = table[func]
                        if not pair in mytab:
                            continue
                        mayer = compute_mayer(func, mytab[pair], Temp, reverse_mass, args.qmcut)
                    label = ""
                    if mydebug:
                        dump_mayer(mayer, ( "%s/mayer-%s-%s-%g.xvg" % ( resultdir, func, pair, Temp )))
                        label = ( "%s/integral-%s-%s-%g.xvg" % (resultdir, func, pair, Temp ))
                    b2cl[func]  = -0.5*integrate_sphere(mayer[0], label, lengthfac)
                    b2qm[func]  = []
                    for k in range(1,4):
                        b2qm[func].append(0.5*integrate_sphere(mayer[k], "", lengthfac))
                for of in relabs:
                    outf[of].write("%10g" % Temp)
                outf["abs"].write("  %10g" % ( b2data[bb]["b2"] ))
                for func in allfuncs:
                    if not func in b2cl:
                        continue
                    if B2debug:
                        print("Func %s Pair %s Temp %g B2cl %g B2qm1 %g B2qm2 %g B2qm3 %g" % ( func, pair, Temp, b2cl[func], b2qm[func][0], b2qm[func][1], b2qm[func][2]))
                    b2calc = b2cl[func] + b2qm[func][0] + b2qm[func][1] + b2qm[func][2]
                    b2diff = b2calc - b2data[bb]["b2"]
                    outf["rel"].write("  %10g" % (b2diff))
                    outf["abs"].write("  %10g" % (b2calc))
                    mse[func] += b2diff**2
                    nmse[func] += 1
                for of in relabs:
                    outf[of].write("\n")
            rmse[pair] = {}
            for func in allfuncs:
                rmse[pair][func] = 0
                if nmse[func] > 0:
                    rmse[pair][func] = math.sqrt(mse[func]/nmse[func])

            for of in relabs:
                outf[of].close()

    # For printing statistics
    write_b2_stats(resultsdir, suffix, args, rmse, thislot, allfuncs, Tmin, Tmax)

def write_b2_stats(resultsdir:str, suffix:str, args, rmse, thislot:str, allfuncs:list, Tmin:list, Tmax:list):
    potentials = get_potentials()
    exception  = ""
    if "CCSTD(Q)" in thislot:
        exception = "For the qm, the range of integration is limited by the data published elsewhere (see main text)."
    with open(("%s/b2_rmse-%s.tex" % ( resultsdir, suffix)), "w") as outf:
        outf.write("\\begin{table}[ht]\n")
        outf.write("\\caption{Root mean square error from experimental second virial coefficient for quantum chemical reference %s and %s potentials, in cm$^3$/mol.%s Only data points with experimental B(T) $\\ge$ %.0f are taken into account. For training the potentials, the repulsive energy was cut-off at %g kJ/mol and the long range tail cut at %d\\%% of the well-depth. Integration of eq.~\\ref{eqb2} was performed from zero to %g {\\AA}. %s}\n" % (thislot, suffix, pol_string(args.polarizable), args.minb2, args.upper, int(100*args.lower), args.qmcut, exception) )
        polstr = ""
        if args.polarizable:
            polstr = "_pol"
        outf.write("\\label{b2_%s_%s%s}\n" % ( suffix, args.lot, polstr ))
        outf.write("\\begin{tabular}{lcc")
        for i in range(len(allfuncs)):
            outf.write("c")
        outf.write("}\n\\hline\n")
        outf.write("Element& T (K) ")
        for func in allfuncs:
            pname = func
            if func in potentials and "latex" in potentials[func]:
                pname = potentials[func]["latex"]
            outf.write("& %s" % pname)
        outf.write("\\\\\n\\hline\n")
        funcrms = {}
        nrmse   = 0
        for pair in rmse.keys():
            include_pair = True
            outf.write("%s & %d-%d " % ( pair, int(Tmin[pair]), int(Tmax[pair])))
            for func in allfuncs:
                outf.write("& %.2f" % ( rmse[pair][func] ))
                if include_pair:
                    if not func in funcrms:
                        funcrms[func] = 0
                    funcrms[func] += rmse[pair][func]**2
            if include_pair:
                nrmse += 1
            outf.write("\\\\\n")
        outf.write("\\hline\n")
        outf.write("RMSE &")
        for func in allfuncs:
            if nrmse > 1:
                outf.write("& %.2f" % ( math.sqrt(funcrms[func]/nrmse)))
            else:
                outf.write("& -")
        outf.write("\\\\\n")
        outf.write("\\hline\n\\end{tabular}\n\\end{table}\n")
    print("Check second virial plots in %s" % resultsdir)

def read_fit_table(tablefn:str)->dict:
    if B2debug:
        print("Reading table %s" % tablefn)
    table = {}
    comb  = {}
    with open(tablefn, "r") as inf:
        header = []
        for line in inf:
            words = line.strip().split(",")
            if len(header) == 0:
                header = words
            elif len(words) == len(header):
                try:
                    if "fit" == words[1]:
                        table[words[0]] = {}
                        for k in range(2,len(header)):
                            table[words[0]][header[k]] = float(words[k])
                    else:
                        comb[words[0]] = {}
                        for k in range(2,len(header)):
                            comb[words[0]][header[k]] = float(words[k])
                except ValueError:
                    sys.exit("Wrong line '%s'" % line.strip())
    return table, comb

def parse_args():
    parser  = argparse.ArgumentParser(description="""
    This script computes the second virial coefficient for potentials from Wei2017a as well as
    those computed here. The experimental second virial is plotted as well.
    The upper and lower flags are used for finding the corresponding RESULTS directory.
    """)
    
    defqm = "CCSDT_CBS"
    parser.add_argument("-lot", "--lot", help="Level of theory for quantum calculations, default "+defqm, type=str, default=defqm)
    parser.add_argument("-up","--upper", help="Upper cutoff of energy, default 20 kJ/mole. Will take only data points with energies below this value.", default=20, type=int)
    parser.add_argument("-down","--lower", help="Lower cutoff. Fraction of the minimum energy, defaults to 0.1, will take only energies with absolute value greater than this percentage of energy well depth", default=0.1, type=float)
    nlist = [ "He", "Ne", "Ar", "Kr", "Xe" ]
    parser.add_argument("-elem", "--elements", nargs="+", help="List of noble gas elements to consider, default all", default=nlist)
    parser.add_argument("-v", "--verbose", help="Print more output, mainly for debugging", action="store_true")
    minb2 = -1000
    parser.add_argument("-minb2", "--minb2", help="Do not consider experimental data with B2 smaller than this value, default "+str(minb2), type=float, default=minb2)
    qmcut = 100.0
    parser.add_argument("-qmcut", "--qmcut", help="Cut-off distance for QM calculations, energies at distances larger than this will be disregarded. Default "+str(qmcut)+" Angstrom", type=float, default=qmcut) 
    parser.add_argument("-interpolate", "--interpolate", help="Interpolate quantum chemistry data when computing B2, with this many point between data points. Default 0", type=int, default=0)
    parser.add_argument("-pol", "--polarizable", help="Subtract electrostatic energy before fitting", action="store_true")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    if args.verbose:
        B2debug = True
    resultsdir = ( "RESULTS_%s_%d_%g" % ( args.lot, args.upper, args.lower ))
    if args.polarizable:
        resultsdir += "_pol"
    if not os.path.isdir(resultsdir):
        sys.exit("No such directory %s" % resultsdir)
    table1 = {}
    table1["Wei2017a"] = get_table1("lmwei2017_table1.csv")
    table4 = get_table1("lmwei2017_table4.csv")
    for t4 in table4.keys():
        table1["Wei2017a"][t4] = table4[t4]
    comb = table1.copy()
    pots = [ "LJ_12_6", "LJ_8_6", "WBH", "MBH", "MRS", "BHA", "GBH", "LJ_14_7", "TT" ]
    pots = [ "LJ_12_6", "MRS", "WBH", "GBH", "LJ_14_7", "TT"  ]
    for pot in pots:
        table1[pot], comb[pot]  = read_fit_table(("%s/%s.csv" % (resultsdir, pot)))
    calc_b2(table1, args, resultsdir, "fitted")
    calc_b2(comb, args, resultsdir, "combined")
