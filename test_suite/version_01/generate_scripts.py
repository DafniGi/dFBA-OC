#!/home/users/mgotsmy/.conda/envs/2302test3.10/bin/python

import os
import itertools as it

nFE_list  = [5,10,15,20]
tmax_list = [10,15,20]
phi_list  = [1e0]

base_name = "base_v01"
version = base_name.split("_")[1]
with open(f"{base_name}.jl","r") as file:
    base = file.read()

run_file = "#!/bin/bash \n"
executable = "/home/users/mgotsmy/julia/julia-1.8.5/bin/julia -Cnative -J/home/users/mgotsmy/julia/julia-1.8.5/lib/julia/sys.so -g1 --color=yes"
for n,(nFE, tmax, phi) in enumerate(it.product(nFE_list,tmax_list,phi_list)):
    tmp = base
    tmp = tmp.replace("##nFE##",f" {nFE} ")
    tmp = tmp.replace("##TMAX##",f" {tmax} ")
    tmp = tmp.replace("##PHI##",f" {phi} ")
    tmp_name = f"{version}.{n+1:02d}.jl"
    with open(tmp_name,"w") as file:
        file.write(tmp)
    run_file += f"{executable} {tmp_name} \n"
    print(nFE,tmax,phi,tmp_name)

with open("run.sh","w") as file:
    file.write(run_file)
    
os.system("chmod +x run.sh")
    
print("DONE")