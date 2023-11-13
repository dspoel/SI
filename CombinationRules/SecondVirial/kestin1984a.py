#!/usr/bin/env python

import math, os, sys

def read_data(infile:str)->list:
    data = []
    with open(infile, "r") as inf:
        labels = []
        linenum = 0
        for line in inf:
            linenum += 1
            if line.find("#") >= 0:
                labels = line.strip().split(",")
                labels[0] = labels[0][1:]
            elif len(labels) > 0:
                words = line.strip().split(",")
                if len(words) != len(labels):
                    sys.exit("Inconsistency on line %d in %s" % ( linenum, infile ))
                else:
                    try:
                        T = float(words[0])
                        if "C" == words[1]:
                            T += 273.15
                        elif "K" != words[1]:
                            sys.exit("Incorrect temperature unit %s on line %d in %s" % ( words[1], linenum, infile ))
                        data.append( { "T": T, "B": float(words[2]),
                                       "x0": float(words[3]),
                                       "x1": float(words[4]) } )
                    except ValueError:
                        sys.exit("Incomprehensible line %d in %s" % ( linenum, infile ))
    return data
    
def plot_b12(outfile:str, dimerdata:list, mon0:dict, mon1:dict, stddev:bool, Tmin:float):
    # Get sorted temperature list
    b12  = {}
    for d in dimerdata:
        b12[d["T"]] = []
    for T in sorted(b12.keys()):
        for d in dimerdata:
            if T != d["T"]:
                continue
            if T < Tmin:
                continue
            x0 = d["x0"]
            x1 = d["x1"]
            if x0*x1 != 0:
                if T in mon0 and T in mon1:
                    bb = (d["B"] - x0**2*mon0[d["T"]] - x1**2*mon1[d["T"]])/(2*x0*x1)
                    b12[T].append(bb)
    with open(outfile, "w") as outf:
        outf.write("# Data from Kestin et al. Journal of Physical and Chemical Reference Data 13, 229 (1984)\n")
        for T in sorted(b12.keys()):
            b2aver = 0
            b2aver2 = 0
            for b in b12[T]:
                b2aver += b
                b2aver2 += b*b
            nb12 = len(b12[T])
            if 0 == nb12:
                print("Empty b12 for %s" % outfile)
                continue
            b2aver /= nb12
            if stddev:
                b2std = 0
                if nb12 > 1:
                    diff = (b2aver2/nb12) - b2aver**2
                    if diff >= 0:
                        b2std   = math.sqrt(diff/(nb12-1) )
                outf.write("%g,%g,%g\n" % ( T, b2aver, b2std ))
            else:
                outf.write("%g,%g\n" % ( T, b2aver ))

if __name__ == "__main__":
    elems = [ "helium", "neon", "argon", "krypton", "xenon" ]
    mondata = {}
    rawdata = "Kestin1984a"
    for monomer in elems:
        monomerdata = ("%s/kestin1984a_%s#%s.csv" % ( rawdata, monomer, monomer))
        mondata[monomer] = {}
        montemp = read_data(monomerdata)
        print("Found %d data points in %s" % (len(montemp), monomerdata))
        for tmp in montemp:
            mondata[monomer][tmp["T"]] = tmp["B"]

    Tmin = { "neon#xenon": 200, "argon#xenon": 200 }
    stddev = False
    for m1 in range(len(elems)):
        for m2 in range(m1+1, len(elems)):
            dimer = ( "%s#%s" % ( elems[m1], elems[m2] ))
            dimerdata = ("%s/kestin1984a_%s.csv" % ( rawdata, dimer))
            data = read_data(dimerdata)
            print("Found %d data points in %s" % (len(data), dimerdata))
            TTT = 0
            if dimer in Tmin:
                TTT = Tmin[dimer]
            plot_b12(("%s.csv" % dimer), data, mondata[elems[m1]], mondata[elems[m2]], stddev, TTT)
