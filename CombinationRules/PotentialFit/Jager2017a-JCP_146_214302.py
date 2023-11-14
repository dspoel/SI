#!/usr/bin/env python3

import math, os
import numpy as np
from calc_virials import Boltz, Bohr, Hartree

debug   = False

def faculty(k:int)->int:
    if 0 == k:
        return 1
    else:
        return k*faculty(k-1)

def Ujager2017a(R, params:dict)->float:
    A = params["A"]
    a_1 = params["a_1"]
    a_2 = params["a_2"]
    a_min1 = params["a_-1"]
    a_min2 = params["a_-2"]
    d1 = params["d_1"]
    d2 = params["d_2"]
    d3 = params["d_3"]

    b = params["b"]

    Urep = A * np.exp(a_1 * R
                      + a_2 * R**2
                      + a_min1/R
                      + a_min2/R**2
                      )
    
    
    #x     = R/params["Re"]
    x = R
    #Astar = params["A"]/params["De"]
    #bstar = params["b"]*params["Re"]
    expbx = np.exp(-b*x)
    #Urep  = Astar*expbx
    Udisp = 0
    for n in [ 3,4,5,6,7,8 ]:
        facsum = 0
        for k in range(0, 2*n+1):
            facsum += ((b*x)**k)/faculty(k)
        if n >= 6:
            C2nmin2 = params["c_" + str(2*n-2)]
            C2nmin4 = params["c_" + str(2*n-4)]
            C2nmin6 = params["c_" + str(2*n-6)]
            C2n = C2nmin6 * (C2nmin2/C2nmin4)**3
            params["c_" + str(2*n)] = C2n
            #print("%s %f" % ("c_" + str(2*n), C2n))
        else:
            C2n = params["c_" + str(2*n)]
        Udisp += (1 - expbx*facsum)*C2n/(x**(2*n))
    if debug:
        print(R, Urep,Udisp, Urep-Udisp)
    
    return (Urep-Udisp)

def get_table(filenm:str)->dict:
    table = {}
    with open("LiteratureData/"+filenm, "r") as inf:
        for line in inf:
            words = line.strip().split(",")
            if 19 == len(words):
                try:
                    table[words[0]] = { "A": float(words[1]),
                                        "a_1": float(words[2]),
                                        "a_-1": float(words[3]),
                                        "a_2": float(words[4]),
                                        "a_-2": float(words[5]),
                                        "d_1": float(words[6]),
                                        "d_2": float(words[7]),
                                        "d_3": float(words[8]),
                                        "b": float(words[9]),
                                        "c_6": float(words[10]),
                                        "c_8": float(words[11]),
                                        "c_10": float(words[12]),
                                        "c_12": float(words[13]),
                                        "c_14": float(words[14]),
                                        "c_16": float(words[15]),
                                        "eps": float(words[16]),
                                        "R_eps": float(words[17]),
                                        "sig": float(words[18])
                                       }
                except ValueError:
                    print("Cannot interpret line '%s'" % line.strip())
            else:
                print("Incorrect number of words (%d), expected 19 in line %s" % (len(words), filenm))
    return table

def gen_table3(table:dict):
    print("Jager et al. J. Chem. Phys., 146, 214302, Table III")
    for elem in table.keys():
        A = table[elem]["A"]
        Reps = table[elem]["R_eps"]
        print("%-15s %15.12e" % ("A [K]", table[elem]["A"]) )
        print("%-15s %15.12e" % ("a_1 [Ang^-1]", table[elem]["a_1"]) )
        print("%-15s %15.12e" % ("a_2 [Ang^-2]", table[elem]["a_2"]) )
        print("%-15s %15.12e" % ("a_-1 [Ang]", table[elem]["a_-1"]) )
        print("%-15s %15.12e" % ("a_-2 [Ang**2]", table[elem]["a_-2"]) )
        print("%-15s %15.12e" % ("b [Ang^-1]", table[elem]["b"]) )
        print("%-15s %15.12e" % ("C_6 [K Ang^6]", table[elem]["c_6"]) )
        print("%-15s %15.12e" % ("C_8 [K Ang^8]", table[elem]["c_8"]) )
        print("%-15s %15.12e" % ("C_10 [K Ang^10]", table[elem]["c_10"]) )
        print("%-15s %15.12e" % ("eps/kB [K]", table[elem]["eps"]) )
        print("%-15s %15.12e" % ("R_eps [Ang]", table[elem]["R_eps"]) )

def plot_xvg_ang_kelvin(table:dict):
    for elem in table.keys():
        outfn = ("%s-ang-kelvin.xvg" % ( elem ))
        print("Generating %s" % outfn)
        with open(outfn, "w") as outf:
            outf.write("@ title \"%s\"\n" % ( elem ))
            outf.write("@ xaxis label \"r (Angstrom)\"\n")
            outf.write("@ yaxis label \"E (Kelvin)\"\n")
            for ir in range(61):
                r = (1.3 + 0.1*ir)
                U = Ujager2017a(r, table[elem])
                outf.write("%10g  %15f\n" % ( r, U ))

def plot_xvg(table:dict):
    for elem in table.keys():
        outfn = ("%s.xvg" % ( elem ))
        print("Generating %s" % outfn)
        with open(outfn, "w") as outf:
            outf.write("@ title \"%s\"\n" % ( elem ))
            outf.write("@ xaxis label \"r (Angstrom)\"\n")
            outf.write("@ yaxis label \"E (kJ/mol)\"\n")
            for ir in range(301):
                r = (1.0 + 0.2*ir)
                U = Ujager2017a(r, table[elem])
                outf.write("%10g  %15f\n" % ( r, U*Boltz ))

def gen_csv(table:dict):
    for elem in table.keys():
        outfn = ("%s.csv" % ( elem ))
        print("Generating %s" % outfn)
        with open(outfn, "w") as outf:
            for ir in range(9990):
                r = (1.0 + 0.01*ir)
                U = Ujager2017a(r, table[elem])
                outf.write("%10g,  %15e\n" % ( r, U*Boltz ))


if __name__ == "__main__":
    if debug:
        for k in range(11):
            print("%d! = %d" % ( k, faculty(k)))
    table = get_table("Jager2017a-JCP_146_214302-table3.csv")
    gen_table3(table)
    plot_xvg(table)
    plot_xvg_ang_kelvin(table)
    gen_csv(table)
