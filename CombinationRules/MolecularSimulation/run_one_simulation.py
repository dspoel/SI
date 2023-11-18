#!/usr/bin/env python3

import os, shutil

root        = os.getcwd()

simdatfile  = root + "/Parameters/meltsim.dat"
for pot in [ "LJ14_7" ]:
    forcefield  = root + ( "/ff_%s.xml" % pot )
    os.makedirs(pot, exist_ok=True)
    os.chdir(pot)
    for elem in [ "Xe" ]:
        structure   = root + ( "/Structures/%sbiphase.pdb" % elem )
        os.makedirs(elem, exist_ok=True)
        os.chdir(elem)
        mydatfile = "meltsim.dat"
        os.system("cp %s %s" % (simdatfile, mydatfile))
        with open(mydatfile, "a") as outf:
            outf.write("temperature_c = 300\n")
            outf.write("vanderwaals   = %s\n" % pot)
        jobpy = "job.py"
        shutil.copy(("%s/act_openmm.py" % root), ".")
        with open(jobpy, "w") as outf:
            outf.write("#!/usr/bin/env python3\n")
            outf.write("from act_openmm import ActOpenMMSim\n")
            outf.write("sim = ActOpenMMSim(pdbfile=\"%s\", datfile=\"%s\", xmlfile=\"%s\")\n" %
                       ( structure, mydatfile, forcefield))
            outf.write("sim.run()\n")
        os.system("python3 job.py\n")
        os.chdir("..")
    os.chdir("..")
