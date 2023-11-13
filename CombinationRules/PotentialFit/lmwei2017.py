#!/usr/bin/env python3

import math, os
import numpy as np

debug   = False
Bohr    = 0.529177
Hartree = 2625.5 
Boltz   = 0.008314

def faculty(k:int)->int:
    if 0 == k:
        return 1
    else:
        return k*faculty(k-1)

def Uwei(R, params:dict)->float:
    x     = R/params["Re"]
    Astar = params["A"]/params["De"]
    bstar = params["b"]*params["Re"]
    expbx = np.exp(-bstar*x)
    Urep  = Astar*expbx
    Udisp = 0
    for n in [ 3, 4, 5 ]:
        facsum = 0
        for k in range(0, 2*n+1):
            facsum += ((bstar*x)**k)/faculty(k)
        C2n = params["C" + str(2*n)]/(params["De"]*(params["Re"]**(2*n)))
        Udisp += (1 - expbx*facsum)*C2n/(x**(2*n))
    return (Urep - Udisp)*params["De"]
  
def get_table1(filenm:str)->dict:
    table = {}
    with open("LiteratureData/"+filenm, "r") as inf:
        for line in inf:
            words = line.strip().split(",")
            if 8 == len(words):
                try:
                    table[words[0]] = { "A": float(words[1]),
                                        "b": float(words[2]),
                                        "C6": float(words[3]),
                                        "C8": float(words[4]),
                                        "C10": float(words[5]),
                                        "Re": float(words[6]),
                                        "De": float(words[7]) }
                except ValueError:
                    print("Cannot interpret line '%s'" % line.strip())
            else:
                print("Incorrect number of words (%d), expected 8 in line %s" % (len(words), filenm))
    return table

def gen_table2(table:dict):
    print("L. M. Wei et al. Chem. Phys. Lett. 675 (2017) 40-45, Table 2")
    for elem in table.keys():
        De = table[elem]["De"]
        Re = table[elem]["Re"]
        print("%5s  %8.3f  %8.2f %8.3f  %8.4f  %8.4f" %
              (elem,
               table[elem]["A"]/(De*1e6),
               table[elem]["b"]*Re,
               table[elem]["C6"]/(De*Re**6),
               table[elem]["C8"]/(De*Re**8),
               table[elem]["C10"]/(De*Re**10) ))

def gen_table4(table:dict):
    print("L. M. Wei et al. Chem. Phys. Lett. 675 (2017) 40-45, Table 4")
    for elem in table.keys():
        De = table[elem]["De"]
        Re = table[elem]["Re"]
        print("%5s  %8.2f  %8.2f %8.1f  %8.2f  %8.3f" %
              (elem,
               table[elem]["C6"],
               table[elem]["C8"],
               table[elem]["C10"],
               table[elem]["A"],
               table[elem]["b"]
               ))
               
def plot_xvg(table:dict):
    for elem in table.keys():
        outfn = ("%s.xvg" % ( elem ))
        print("Generating %s" % outfn)
        with open(outfn, "w") as outf:
            outf.write("@ title \"%s\"\n" % ( elem ))
            outf.write("@ xaxis label \"r (Angstrom)\"\n")
            outf.write("@ yaxis label \"E (kJ/mol)\"\n")
            for ir in range(61):
                r = (0.8 + 0.02*ir)*table[elem]["Re"]
                U = Uwei(r, table[elem])
                outf.write("%10g  %10g\n" % ( r*Bohr, U*table[elem]["De"]*Hartree ))


if __name__ == "__main__":
    if debug:
        for k in range(11):
            print("%d! = %d" % ( k, faculty(k)))
    table1 = get_table1("lmwei2017_table1.csv")
    gen_table2(table1)
    plot_xvg(table1)
    table4 = get_table1("lmwei2017_table4.csv")
    gen_table4(table4)
    plot_xvg(table4)
