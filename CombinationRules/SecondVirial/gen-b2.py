#!/usr/bin/env python3

def header(outf):
    outf.write("#T,B2(T)\n")

def helium():
    with open("helium#helium.csv", "w") as outf:
        header(outf)
        for T in range(2,36,1):
            outf.write("%g,%g\n" % ( T, 1.5943*10 - 3.4601e2/T - 5.9545e2/T**2 + 1.9929e3/T**3 - 2.2269e3/T**4))
        for T in range(40,1471,10):
            outf.write("%g,%g\n" % ( T, 9.2479 + 1.0876e3/T - 1.0880e5/T**2 + 2.3869e6/T**3 ))

def helium_xenon():
    with open("helium#xenon.csv", "w") as outf:
        outf.write("# Data taken from Hurly et al. Int. J. Thermophysics. Vol. 18, 1997, p. 579\n")
        for T in range(225,400,5):
            val = 1e6*(2.756941e-5 + 4.617880e-3/T - 1.741538/T**2 + 1.364453e2/T**3)
            outf.write("%g,%g\n" % ( T, val ))

def helium_neon_lb():
    with open("helium#neon.csv", "w") as outf:
        outf.write("# Data taken from LB\n")
        for T in range(15,350,5):
            val = 1.3558e1 - 2.5654e2/T - 2.8502e4/T**2 + 2.2060e5/T**3
            outf.write("%g,%g\n" % ( T, val ))

def helium_argon_lb():
    with open("helium#argon.csv", "w") as outf:
        outf.write("# Data taken from LB\n")
        for T in range(100,750,10):
            val = 1.4965e1 + 6.8188e3/T - 2.6724e6/T**2 + 3.2953e8/T**3 - 1.3616e10/T**4
            outf.write("%g,%g\n" % ( T, val ))

            
def helium_krypton_lb():
    with open("helium#krypton.csv", "w") as outf:
        outf.write("# Data taken from LBi, except first value Kestin\n")
        outf.write("50,-29.28\n")
        for T in range(90,350,5):
            val = 2.2445e1 + 1.2894e3/T - 5.1960e5/T**2 + 2.1073e7/T**3
            outf.write("%g,%g\n" % ( T, val ))

def helium_xenon_lb():
    with open("helium#xenon.csv", "w") as outf:
        outf.write("# Data taken from LB\n")
        for T in range(120,350,5):
            val = 2.2759e1 + 3.5574e3/T - 5.8539e5/T**2 - 4.5983e6/T**3
            outf.write("%g,%g\n" % ( T, val ))

def neon_krypton_lb():
    with open("neon#krypton.csv", "w") as outf:
        outf.write("# Data taken from LB\n")
        for T in range(100,500,10):
            val = 3.0219e1 - 6.1482e3/T - 6.7100e4/T**2
            outf.write("%g,%g\n" % ( T, val ))
    
def neon_argon_lb():
    with open("neon#argon.csv", "w") as outf:
        outf.write("# Data taken from LB\n")
        for T in range(100,475,5):
            val = 2.9286 + 1.2477e4/T - 4.1259e6/T**2 + 3.7687e8/T**3 - 1.2338e10/T**4
            outf.write("%g,%g\n" % ( T, val ))

def neon_xenon_lb():
    with open("neon#xenon.csv", "w") as outf:
        outf.write("# Data taken from LB\n")
        for T in range(160, 500, 5):
            val = 1.4791e1 + 4.7169e3/T - 1.1208e6/T**2 - 8.5658e7/T**3
            outf.write("%g,%g\n" % ( T, val ))

def argon_krypton_lb():
    with open("argon#krypton.csv", "w") as outf:
        outf.write("# Data taken from LB\n")
        for T in range(110,700,10):
            val = 3.8147e1 - 1.5052e4/T - 1.5687e6/T**2
            outf.write("%g,%g\n" % ( T, val ))

def argon_xenon():
    with open("argon#xenon.csv", "w") as outf:
        data = []
        for src in [ "Brewer1967", "Rentschler1977", "Schramm1977a" ]:
            with open(("Ar-Xe/%s.dat" % src), "r") as inf:
                for line in inf:
                    words = line.strip().split()
                    if len(words) == 2:
                        data.append( ( float(words[0]), float(words[1]) ) )
        for d in sorted(data):
            outf.write("%s,%s\n" % ( d[0], d[1] ))
                    
def krypton_xenon_lb():
    with open("krypton#xenon.csv", "w") as outf:
        outf.write("# Data taken from LB\n")
        for T in range(150,750,10):
            val = 4.6804e1 - 2.5419e4/T - 4.1375e6/T**2 + 2.6799e8/T**3 - 3.5667e10/T**4
            outf.write("%g,%g\n" % ( T, val ))

def neon():
    with open("neon#neon.csv", "w") as outf:
        header(outf)
        for T in range(50,871,4):
            outf.write("%g,%g\n" % ( T, 1.5894*10 - 9.9406e2/T - 1.2641e5/T**2 + 2.2721e6/T**3))

def argon():
    with open("argon#argon.csv", "w") as outf:
        header(outf)
        for T in range(76,1001,4):
            outf.write("%g,%g\n" % ( T, 3.4162*10 - 1.2087e4/T - 7.6702e5/T**2 - 1.9600e7/T**3 ))

def krypton():
    with open("krypton#krypton.csv", "w") as outf:
        header(outf)
        for T in range(110,871,10):
            outf.write("%g,%g\n" % ( T, 3.8030*10 - 1.9985e4/T - 1.4767e6/T**2 - 1.3450e8/T**3) )

def xenon():
    with open("xenon#xenon.csv", "w") as outf:
        header(outf)
        for T in range(165,971,10):
            outf.write("%g,%g\n" % ( T, 6.7836*10 - 5.3169e4/T + 2.2239e6/T**2 - 1.3905e9/T**3 + 5.8510e10/T**4 ) )
            
if __name__ == "__main__":
    print("# All equations taken from Landolt-Bornstein, New Series IV/21A")
    print("# https://materials.springer.com/googlecdn/assets/sm_lbs/025/sm_lbs_978-3-540-70732-5_3/sm_lbs_978-3-540-70732-5_3.pdf?trackRequired=true&originUrl=/lb/docs/sm_lbs_978-3-540-70732-5_3&componentId=Download%20Chapter")
    helium()
    # This formula is more limited than the Data from Kestin
    #    helium_xenon()
    neon()
    argon()
    krypton()
    xenon()
    helium_neon_lb()
    helium_argon_lb()
    helium_krypton_lb()
    helium_xenon_lb()
    neon_argon_lb()
    neon_krypton_lb()
    neon_xenon_lb()
    argon_krypton_lb()
    argon_xenon()
    krypton_xenon_lb()
