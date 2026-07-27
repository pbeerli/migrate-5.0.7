#!/usr/bin/env python
"""Overdispersed-start test: does the M chain's answer depend on where it starts?

Both conditionals p(G|theta,M) and p(M|G) are verified exact, so if starts at
low and high M fail to meet, the chain is not converging (trapping); if they all
collapse to the same value, the stationary distribution itself is wrong.
"""
import subprocess, sys, tempfile, numpy as np
from pathlib import Path

MIGRATE = str(Path.home()/"src/migrate-5.0.7/src/migrate-n")
THETA=0.05; SITES=10; NPOP=2; NPERPOP=4; SCALE=20000.0
SAMPLE=int(sys.argv[1]) if len(sys.argv)>1 else 20000

PARM="""################################################################################
# Parmfile for Migrate
menu=NO
datatype=SequenceData
datamodel=JC69
freqs-from-data=NO:0.25,0.25,0.25,0.25
infile=infile
progress=NO
print-data=NO
pdf-outfile=NO
heating=NO
replicate=NO
long-chains=1
long-inc=20
long-sample={sample}
burn-in=0
auto-tune=NO
use-M=YES
custom-migration={{**}}
updatefreq= tree:0.70 parameter:0.30 scaler:0.00
theta=OWN:{theta}
migration=OWN:{mstart}
bayes-priors= THETA * * UNIFORMPRIOR: {tlo:.7f} {thi:.7f} 0.0000005
bayes-priors= MIG * * EXPPRIOR: 0.000000 {scale:.6f} {scale:.6f}
random-seed=OWN:{seed}
outfile=outfile
bayes-allfile=YES:1:bayesall
"""
COLS=("Steps Locus Replicate lnPost lnDataL lnPrbGParam lnPrior treeintervals "
      "treelength Theta_1 Theta_2 M_2_1 M_1_2").split()

def run(mstart,seed):
    with tempfile.TemporaryDirectory() as wd:
        wd=Path(wd)
        with open(wd/"infile","w") as f:
            f.write("   %i 1 . allmiss\n%i\n"%(NPOP,SITES))
            for p in range(1,NPOP+1):
                f.write("   %i Pop%i\n"%(NPERPOP,p))
                for i in range(NPERPOP):
                    f.write("%-10s%s\n"%("p%si%i"%(p,i),"?"*SITES))
        (wd/"parmfile").write_text(PARM.format(theta=THETA,tlo=THETA*0.99998,
            thi=THETA*1.00002,scale=SCALE,seed=seed,mstart=mstart,sample=SAMPLE))
        subprocess.run([MIGRATE,"parmfile","-nomenu"],cwd=wd,
                       stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True)
        rows=[]
        for line in (wd/"bayesall").read_text().splitlines():
            if line.startswith("#") or line.startswith("Steps"): continue
            p=line.split("\t")
            if len(p)>=len(COLS): rows.append([float(x) for x in p[:len(COLS)]])
    A=np.array(rows); d={c:A[:,i] for i,c in enumerate(COLS)}
    return 0.5*(d["M_2_1"]+d["M_1_2"]), d["treeintervals"]-2*NPOP*NPERPOP

if __name__=="__main__":
    tmean=(0+SCALE)  # analytic below
    z=1-np.exp(-1.0); target=(SCALE-(SCALE+SCALE)*np.exp(-1.0))/z
    print("prior truncExp(scale %g) on [0,%g]: target mean %.1f"%(SCALE,SCALE,target))
    print("recorded samples per run: %d   (burn-in disabled)\n"%SAMPLE)
    print(" start M   seed |  fifths of the chain: mean M in each 20%% block          | 2nd-half mean  n_mig(last 20%%)")
    print("-"*118)
    for mstart in (50.0, 5000.0, 19500.0):
        for seed in (131,231,331):
            M,nm=run(mstart,seed)
            k=len(M)//5
            blocks=" ".join("%8.0f"%M[i*k:(i+1)*k].mean() for i in range(5))
            print("%8.0f %6d | %s |   %9.0f      %7.0f"%(
                  mstart,seed,blocks,M[len(M)//2:].mean(),nm[-k:].mean()))
