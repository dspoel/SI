#!/usr/bin/env python3

import os, sys
from calc_virials import Boltz, Bohr, Hartree

CM1 = 0.0119627

def jager2017():
    infile = "LiteratureData/jager2017a-table2.txt"
    outfile = "CCSDTQ_CBS/helium#krypton.csv"
    with open(outfile, "w") as outf:
        with open(infile, "r") as inf:
            for line in inf:
                words = line.strip().split()
                if 10 == len(words):
                    try:
                        vtot = float(words[7])
                        outf.write("%s,%g\n" % ( words[0], vtot*Boltz ))
                    except ValueError:
                        sys.exit("Cannot interpret line '%s' from %s" % ( line.strip(), infile ))

def liu2022():
    infile = "LiteratureData/liu2022a-table2.txt"
    xe_he = open("CCSDTQ_CBS/helium#xenon.csv", "w")
    xe_ne = open("CCSDTQ_CBS/neon#xenon.csv", "w")
    xe_ar = open("CCSDTQ_CBS/argon#xenon.csv", "w")
    with open(infile, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if 7 == len(words):
                try:
                    r    = float(words[0])
                    xe_he.write("%s,%g\n" % ( words[0], float(words[1])*CM1 ) )
                    xe_ne.write("%s,%g\n" % ( words[0], float(words[3])*CM1 ) )
                    xe_ar.write("%s,%g\n" % ( words[0], float(words[5])*CM1 ) )
                except ValueError:
                    sys.exit("Cannot interpret line '%s' from %s" % ( line.strip(), infile ))
    xe_he.close()
    xe_ne.close()
    xe_ar.close()

def hu2022():
    infile = "LiteratureData/hu2022a-table1.txt"
    ne_kr = open("CCSDTQ_CBS/neon#krypton.csv", "w")
    ar_kr = open("CCSDTQ_CBS/argon#krypton.csv", "w")
    kr_xe = open("CCSDTQ_CBS/krypton#xenon.csv", "w")
    with open(infile, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if 7 == len(words):
                try:
                    ne_kr.write("%s,%g\n" % ( words[0], float(words[1])*CM1 ) )
                    ar_kr.write("%s,%g\n" % ( words[0], float(words[3])*CM1 ) )
                    kr_xe.write("%s,%g\n" % ( words[0], float(words[5])*CM1 ) )
                except ValueError:
                    sys.exit("Cannot interpret line '%s' from %s" % ( line.strip(), infile ))
    ne_kr.close()
    ar_kr.close()
    kr_xe.close()

def lopez2006():
    infile = "LiteratureData/lopez-cacheiro2006-table1.csv"
    he_ne = open("CCSDTQ_CBS/helium#neon.csv", "w")
    he_ar = open("CCSDTQ_CBS/helium#argon.csv", "w")
    ne_ar = open("CCSDTQ_CBS/neon#argon.csv", "w")
    mH = Hartree/1e6
    with open(infile, "r") as inf:
        for line in inf:
            words = line.strip().split(",")
            if line.find("R") >= 0:
                continue
            if 6 == len(words):
                for w in range(len(words)):
                    words[w] = words[w].replace('"','')
                try:
                    he_ne.write("%g,%g\n" % ( float(words[0])*Bohr, float(words[1])*mH ) )
                    if len(words[2]) > 2:
                        he_ar.write("%g,%g\n" % ( float(words[2])*Bohr, float(words[3])*mH ) )
                    if len(words[4]) > 2:
                        ne_ar.write("%g,%g\n" % ( float(words[4])*Bohr, float(words[5])*mH ) )
                except ValueError:
                    sys.exit("Cannot interpret line '%s' from %s" % ( line.strip(), infile ))
    he_ne.close()
    he_ar.close()
    ne_ar.close()
    
def hellmann2007a():
    infile = "LiteratureData/hellmann2007a_table5.txt"
    he_he = open("CCSDTQ_CBS/helium#helium.csv", "w")
    with open(infile, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if 2 == len(words):
                try:
                    he_he.write("%g,%g\n" % ( float(words[0])*Bohr, float(words[1])*Boltz ) )
                except ValueError:
                    sys.exit("Cannot interpret line '%s' from %s" % ( line.strip(), infile ))
    he_he.close()
    
if __name__ == "__main__":
    jager2017()
    liu2022()
    hu2022()
    lopez2006()
    hellmann2007a()
