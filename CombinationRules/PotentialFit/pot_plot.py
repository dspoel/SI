#!/usr/bin/env python3

import os
from gen_plots import make_xvg

def get_bounds()->dict:
    bounds = {}
    csv = "plot_bounds.csv"
    with open(csv, "r") as inf:
        for line in inf:
            words = line.strip().split(",")
            if len(words) != 3:
                print("Incorrect line '%s' in %s" % ( line.strip(), csv ))
                continue
            try:
                bounds[words[0]] = { "xmin": float(words[2]), "ymin": float(words[1]) }
            except ValueError:
                print("Strange line '%s' in %s" % (line.strip(), csv ))
    return bounds
    
def plot(bounds:dict):
    with open("potentials.tex", "w") as outf:
        plots = "Plots"
        os.makedirs(plots, exist_ok=True)
        os.chdir(plots)
        outf.write("\\begin{figure}[ht]\n")
        outf.write("\\centering\n")
        outf.write("\\caption{Dimer energies (kJ/mol) as a function of distance (\\AA) at different levels of theory.}\n")
        outf.write("\\label{dimer_plots}\n")
        for pair in bounds.keys():
            mypairs = []
            mylabel = []
            label = { "CCSDTQ_CBS": "CCSDT(Q)/CBS",
                      "CCSDT_CBS" : "CCSD(T)/CBS",
                      "CCSDT_QZ" : "CCSD(T)/QZ",
                      "CCSDT_TZ" : "CCSD(T)/TZ",
                      "SAPT_QZ" : "SAPT/QZ",
                      "SAPT_TZ" : "SAPT/TZ" }
                      #,"Sheng2020a": "Sheng2020a" }
            for src in label.keys():
                myxvg = ( "../%s/%s.xvg" % ( src, pair ))
                mycsv = myxvg[:-3] + "csv"
                if os.path.exists(mycsv):
                    mypairs.append(mycsv)
                    mylabel.append("'"+label[src]+"'")
            outpdf = ("potential-%s.pdf" % pair).replace("#", "-")
            cmd = ("viewxvg -f")
            for src in mypairs:
                cmd += ( " %s" % src )
            cmd += " -label"
            for src in mylabel:
                cmd += ( " %s" % src )
            cmd += (" -ymin %f -ymax %f -xmin %f -xmax %f -pdf %s -noshow -title %s -alfs 14 -lfs 14 -legend_x 0.4 -legend_y 0.95 -tickfs 14 -xframe 8 -yframe 6 -tfs 14" % 
                    ( bounds[pair]["ymin"], -bounds[pair]["ymin"],
                      bounds[pair]["xmin"], 2*bounds[pair]["xmin"],
                      outpdf, pair ))
            os.system(cmd)
            outf.write("\includegraphics[width=75mm]{%s/%s}\n" % ( plots, outpdf))
        outf.write("\\end{figure}\n")

if __name__ == "__main__":
    bounds = get_bounds()
    print("There are %d pairs to plot" % len(bounds.keys()))
    plot(bounds)
