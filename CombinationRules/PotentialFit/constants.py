
import math
        
Bohr         = 0.529177
Hartree      = 2625.5 
Boltz        = 8.314462618e-3
SpeedOfLight = 2.99792458e5 # nm/ps
AVOGADRO     = 6.02214076e23
PLANCK1      = 6.62606957e-34 # J s, NIST 2010 CODATA
PLANCK       = PLANCK1*AVOGADRO*1e9 # (kJ/mol) ps
hbar         = PLANCK/(2*math.pi)
mass         = { "He": 4.002602, "Ne": 20.1797, "Ar": 39.948, "Kr": 83.798, "Xe": 131.293 }
longname     = { "He": "helium", "Ne": "neon", "Ar": "argon", "Kr": "krypton", "Xe": "xenon" }
