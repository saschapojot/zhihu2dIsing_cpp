import numpy as np

#this script generates MCMC computing scripts

TAll=np.linspace(0.05,30,100)

counter=0
for T in TAll:
    bashContents = []
    bashContents.append("#!/bin/bash\n")
    bashContents.append("#SBATCH -n 12\n")
    bashContents.append("#SBATCH -N 1\n")
    bashContents.append("#SBATCH -t 0-40:00\n")
    bashContents.append("#SBATCH -p CLUSTER\n")
    bashContents.append("#SBATCH --mem=5GB\n")
    bashContents.append("#SBATCH -o outising2d" + str(counter) + ".out\n")
    bashContents.append("#SBATCH -e outising2d" + str(counter) + ".err\n")
    bashContents.append("cd /home/cywanag/liuxi/Documents/cppCode/zhihu2dIsing_cpp\n")
    command="./ising "+str(T)
    bashContents.append(command)
    bsFileName="./isingBash/ising"+str(counter)+".sh"
    fptrTmp=open(bsFileName,"w+")
    for oneline in bashContents:
        fptrTmp.write(oneline)
    fptrTmp.close()
    counter+=1