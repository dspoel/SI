#!/usr/bin/env python3

import math, os, sys
from calc_b2   import b2map, Bohr, Hartree
from lmwei2017 import faculty

def sheng(r:float, Re:float, De:float, Za:int, Zb:int)->float:
    a1X    =  10.34329 #  10.343
    a2X    = -23.68667 # -23.687
    a3X    =  12.34334 #  12.343
    alphaX =  22.76733 #  22.767
    Ax     =  45.612 * 1e5
    gammaX =   1.531
    betaX  =  15.485
    bX     =  13.950
    params = { 6: 1.3499, 8: 0.4147, 10: 0.1716 }
    
    x      = r/Re
    
    expaX  = math.exp(-alphaX*x)
    # Make sure that the terms add up till 0 for x = 1.
    a0X    = 1 - a1X - a2X - a3X
    Ushort = (1.0/x) * (a0X + a1X*x + a2X*x**2 + a3X*x**3) * expaX
    
    expbx  = math.exp(-bX*x)
    Uvdw   = Ax * (x**gammaX) * math.exp(-betaX*x)
    for n in [ 3, 4, 5 ]:
        facsum = 0
        for k in range(2*n+1):
            facsum += ((bX*x)**k)/faculty(k)
        C2n   = params[2*n]
        Uvdw -= (1 - expbx*facsum)*C2n/(x**(2*n))
    Ulong = (1 - expaX) * Uvdw
    
    return (Za*Zb/Re) * Ushort, De * Ulong

def get_table2()->dict:
    tab2  = {}
    stab2 = "LiteratureData/Sheng2020a_table2.txt"
    with open(stab2, "r") as inf:
        for line in inf:
            words = line.strip().split()
            if 3 == len(words):
                try:
                    tab2[words[0]] = { "Re": float(words[1]), "De": float(words[2])*1e-5 }
                except ValueError:
                    sys.exit("Cannot read line '%s' in %s" % ( line.strip(), stab2 ))
    return tab2

atomicnum = { "He": 2, "Ne": 10, "Ar": 18, "Kr": 36, "Xe": 54, "Rn": 86 }
if __name__ == "__main__":
    tab2 = get_table2()
    sdir = "Sheng2020a"
    os.makedirs(sdir, exist_ok=True)
    for pair in tab2.keys():
        ppp = pair.split("-")
        Za = atomicnum[ppp[0]]
        Zb = atomicnum[ppp[1]]
        Ushort, Ulong = sheng(tab2[pair]["Re"], tab2[pair]["Re"],
                              tab2[pair]["De"], Za, Zb )
        print("%s  Re  %g  De  %g  Ushort  %g  Ulong %g sum %g" % 
              ( pair, tab2[pair]["Re"], tab2[pair]["De"], Ushort, Ulong, 
                ( Ushort + Ulong )*Hartree ))
        if pair in b2map:
            fn = ("%s/%s.csv" % ( sdir, b2map[pair]) )
        else:
            fn = ("%s/helium#helium.csv" % ( sdir ) )
        with open(fn, "w") as outf:
            # R is in Bohr (a.u.)
            r = 1
            while r <= 100:
                Ushort, Ulong = sheng(r/Bohr, tab2[pair]["Re"], tab2[pair]["De"], Za, Zb )
                # Radius used for calculation is interpreted as Bohr so we
                # must scale back to get Angstrom.
                # Also, energy returned is in Hartree so we have to multiply by
                # that factor to get kJ/mole.
                outf.write("%g,%e\n" % ( r, Hartree*(Ushort+Ulong) ) )
                r += 0.02

