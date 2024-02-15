#!/usr/bin/env python3

import argparse, math
from calc_virials import  read_fit_table
from constants import mass, longname

def parse_args():
    parser  = argparse.ArgumentParser(description="""
    This script will generate an OpenMM compatible force field for
    a user specified level of theory and potential. Note however, that
    not all potentials are supported, just LJ_12_6, LJ_8_6, LJ_14_7,
    WBH and GBH.
    """)
    
    defqm = "CCSDT_CBS"
    parser.add_argument("-lot", "--lot", help="Level of theory for quantum calculations, default "+defqm, type=str, default=defqm)
    defpot = "LJ_14_7"
    parser.add_argument("-pot", "--pot", help="Potential to make force field file for", type=str, default=defpot)
    parser.add_argument("-up","--upper", help="Upper cutoff of energy, default 20 kJ/mole. Will take only data points with energies below this value.", default=20, type=int)
    parser.add_argument("-down","--lower", help="Lower cutoff. Fraction of the minimum energy, defaults to 0.1, will take only energies with absolute value greater than this percentage of energy well depth", default=0.1, type=float)
    nlist = [ "He", "Ne", "Ar", "Kr", "Xe" ]
    args = parser.parse_args()
    return args

def print_ff(table:dict, outfn:str, pot:str):
    with open(outfn, "w") as outf:
        outf.write("""<?xml version="1.0" encoding="ISO-8859-1"?>
<!DOCTYPE ForceField.dtd PUBLIC "ForceField.dtd" "ForceField.dtd">
<ForceField>
  <AtomTypes>
""")
        elems = {}
        for pair in table.keys():
            ppp = pair.split("-")
            if ppp[0] == ppp[1]:
                elems[ppp[0]] = pair
        for elem in elems:
            outf.write("    <Type name=\"%s\" class=\"%s\" element=\"%s\" mass=\"%g\"/>\n" %
                       ( elem, elem, elem, mass[elem]) )
        outf.write("  </AtomTypes>\n")
        outf.write("  <Residues>\n")
        for elem in elems.keys():
            outf.write("    <Residue name=\"%s\">\n" % longname[elem])
            outf.write("      <Atom name=\"%s\" type=\"%s\" charge=\"0\"/>\n" % (elem, elem))
            outf.write("    </Residue>\n")
        outf.write("  </Residues>\n")
        outf.write("  <CustomNonbondedForce energy=\"0\" bondCutoff=\"3\">\n")
        outf.write("    <UseAttributeFromResidue name=\"charge\"/>\n")
        sigeps = [ "sigma", "epsilon" ]
        for se in sigeps:
            if se in table[elems["Xe"]]:
                outf.write("      <PerParticleParameter name=\"%s\"/>\n" % se)
        for param in table[elems["Xe"]]:
            if not param in sigeps:
                outf.write("      <PerParticleParameter name=\"%s\"/>\n" % param)
        outf.write("      <PerParticleParameter name=\"charge\"/>\n")
        outf.write("      <PerParticleParameter name=\"zeta\"/>\n")
        for elem in elems.keys():
            outf.write("      <Atom type=\"%s\"" % elem)
            pair = elems[elem]
            for param in sigeps:
                fac = 1
                if "sigma" == param:
                    fac = 0.1
                if param in table[pair]:
                    outf.write(" %s=\"%g\"" % ( param, table[pair][param]*fac ))
            for param in table[pair]:
                if not param in sigeps:
                    outf.write(" %s=\"%g\"" % ( param, table[pair][param] ))
            outf.write(" zeta=\"13\"/>\n")
        outf.write("  </CustomNonbondedForce>\n")
        outf.write("  <NonbondedForce coulomb14scale=\"1\" lj14scale=\"1\" energy=\"0\" bondCutoff=\"3\">\n")
        outf.write("    <UseAttributeFromResidue name=\"charge\"/>\n")
        for elem in elems.keys():
            pair = elems[elem]
            epsilon = table[pair]["epsilon"]
            sigma   = table[pair]["sigma"]
            if pot in [ "GBH", "WBH", "LJ_14_7", "LJ_8_6" ]:
                sigma /= 2**(1.0/6.0)
            outf.write("      <Atom type=\"%s\" sigma=\"%g\" epsilon=\"%g\"/>\n" % ( elem, sigma*0.1, epsilon ))
        outf.write("  </NonbondedForce>\n")
        outf.write("</ForceField>\n")

if __name__ == "__main__":
    args       = parse_args()
    resultsdir = ( "RESULTS_%s_%d_%g" % ( args.lot, args.upper, args.lower ))
    table, comb = read_fit_table(("%s/%s.csv" % (resultsdir, args.pot)))
    pot2OpenMM  = { "LJ_12_6": "LJ12_6", "LJ_8_6": "LJ8_6", "LJ_14_7": "LJ14_7",
                    "WBH": "WBH", "GBH": "GBHAM" }
    outfn       = ( "%s/ff_%s.xml" % ( resultsdir, pot2OpenMM[args.pot] ) ) 
    print_ff(table, outfn, args.pot)
    print("Generated force field file %s" % outfn)
