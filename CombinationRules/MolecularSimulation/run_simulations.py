#!/usr/bin/env python3

import argparse, os

def parse_args():
    parser  = argparse.ArgumentParser(
    description="""This script will run OpenMM simulations to study melting of
    noble gas crystals. A directory tree will be created with the structure of
    potential/element/temperature.
    Note that only the following potentials are supported: 
    LJ12_6, LJ8_6, LJ14_7, WBH and GBH and that force field files must be created using
    the generate_ff.py script in the ../PotentialFit directory.
    By default, it will be assumed that you run this script in a cluster and that jobs
    can be spawned using the sbatch command. If not, you can supply the -submit python3 argument
    in which case the simulations will be run sequentially on your computer.
    """)
    
    defpot = [ "LJ14_7" ]
    parser.add_argument("-pot", "--pot", nargs="+", help="Potential(s) to run simulations with, default "+str(defpot), default=defpot)
    nlist = [ "He", "Ne", "Ar", "Kr", "Xe" ]
    parser.add_argument("-elem", "--elements", nargs="+", help="List of noble gas elements to consider, default all", default=nlist)
    mysubmit = "sbatch"
    parser.add_argument("-submit", "--submit", help="Command to run the simulations, default "+mysubmit, type=str, default=mysubmit)
    args = parser.parse_args()
    return args

Trange = { 
    "He": { "min":   0, "max":  10, "delta": 1 },
    "Ne": { "min":   0, "max":  30, "delta": 1 },
    "Ar": { "min":  70, "max": 100, "delta": 5 },
    "Kr": { "min": 100, "max": 150, "delta": 5 },
    "Xe": { "min": 140, "max": 210, "delta": 5 } 
}

if __name__ == "__main__":
    args        = parse_args()
    root        = os.getcwd()
    ppath       = "PYTHONPATH"
    if ppath in os.environ:
        os.environ[ppath] = ("%s:%s" % ( root, ppath ))
    else:
        os.environ[ppath] = root

    simdatfile  = root + "/Parameters/meltsim.dat"
    for pot in args.pot:
        forcefield  = root + ( "/ff_%s.xml" % pot )
        os.makedirs(pot, exist_ok=True)
        os.chdir(pot)
        for elem in args.elements:
            structure   = root + ( "/Structures/%sbiphase.pdb" % elem )
            os.makedirs(elem, exist_ok=True)
            os.chdir(elem)
            mydatfile = "meltsim.dat"
            for Temp in range(Trange[elem]["min"], Trange[elem]["max"], Trange[elem]["delta"]):
                print("Running %s/%s/%s" % ( pot, elem, Temp))
                os.makedirs(str(Temp), exist_ok=True)
                os.chdir(str(Temp))
                os.system("cp %s %s" % (simdatfile, mydatfile))
                with open(mydatfile, "a") as outf:
                    outf.write("temperature_c = %g\n" % Temp)
                    outf.write("vanderwaals   = %s\n" % pot)
                jobpy = ( "%s-%s-%g.py" % ( pot, elem, Temp ))
                with open(jobpy, "w") as outf:
                    outf.write("#!/usr/bin/env python3\n")
                    outf.write("#SBATCH -c 4\n")
                    outf.write("#SBATCH -t 24:00:00\n")
                    outf.write("from act_openmm import ActOpenMMSim\n")
                    outf.write("sim = ActOpenMMSim(pdbfile=\"%s\", datfile=\"%s\", xmlfile=\"%s\")\n" %
                               ( structure, mydatfile, forcefield))
                    outf.write("sim.run()\n")
                os.system("%s %s\n" % (args.submit, jobpy))
                os.chdir("..")
            os.chdir("..")
        os.chdir("..")
