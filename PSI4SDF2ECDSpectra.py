#!/usr/bin/env python
from rdkit import Chem
import psi4
import numpy as np
import pandas as pd
from psi4.driver.procrouting.response.scf_response import tdscf_excitations
import sys
import string
import os

###
def insert_first_line(file,value):
    f = open(file, "r")
    contents = f.readlines()
    f.close()
    contents.insert(0, value)
    f = open(file, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()
###
#
# Author: Gaokeng Xiao
# Email: info@molcalx.com
# @ Phone: 020-38261356
# Organization: Guangzhou Molcalx.Information & Technology Ltd
# Home page: http://www.molcalx.com.cn.
# Date: 2022-09-06
# Version: 0.1
#
# This script is used to calculate ECD spectra with PSI4 at B3LYP-D3BJ/6-311+g(d,p) level.
#
if len(sys.argv) < 4:
    print("")
    print('This script is used to calculate ECD spectra with PSI4 at B3LYP/6-311+g(d,p) level')
    print("")
    print('There are something wrong with the input parameters:')
    print('No input file. Please specify a input SDF file, output file and solvent.')
    print("")
    print("usage:  ")
    print(sys.argv[0]+" <SDF file>  <output file> <Solvent>")
    print("")
    print("For example:")
    print(sys.argv[0]+" CONF_1.sdf CONF_1_tddft_output.dat Methanol")
    print('')
    print('The available solvent: https://pcmsolver.readthedocs.io/en/stable/users/input.html#available-solvents')
    print('')
    print("Any question, Please feel free to contact me. info@molcalx.com")
    sys.exit()

# input SDF file
isdf = sys.argv[1]
if not os.path.exists(isdf):
   print("Sorry, I cannot find the %s file" % isdf)
   sys.exit()


# output PSI4 file
output = sys.argv[2]
if os.path.exists(output):
   print("the %s file already exist, please use another one." % output)
   sys.exit()


# solvent
solvent = sys.argv[3]

# Read first molecule from input SDF file
#
suppl = Chem.SDMolSupplier(isdf,removeHs=False)
mol = suppl[0]
title = mol.GetProp('_Name')

if mol.HasProp('QM_ENERGY'):
    heat = mol.GetProp('QM_ENERGY')
else:
    print("the input SDF file didn't include tag of QM_ENERGY, exit...")
    sys.exit()

n = mol.GetNumAtoms()
qm_mol=[]
formalcharge = str(Chem.rdmolops.GetFormalCharge(mol))
charge_and_multiplicity=str(formalcharge)+' 1'
qm_mol.append(charge_and_multiplicity)

for i in range(n):
    pos = mol.GetConformer().GetAtomPosition(i)
    element = mol.GetAtoms()[i].GetSymbol()
    x = str(round(pos.x,4))
    y = str(round(pos.y,4))
    z = str(round(pos.z,4))
    line = element+" "+x+" "+y+" "+z
    qm_mol.append(line)

qm_molblock = '\n'.join(qm_mol[:])

# Output file path
# arg 0: Path of output file
# arg 1: Whether overwrite or not
psi4.core.set_output_file(output,True)

# set threads and memory
psi4.core.set_num_threads(24)
psi4.set_memory(int(24e9))

# Scratch directory path
psi4.core.IOManager.shared_object().set_default_path("./tmp")


# Set solvent model
pcm_string = """
    Units = Angstrom
    Medium {
    SolverType = IEFPCM
    Solvent = %s
    Nonequilibrium = True
    }
    Cavity {
    Type = GePol
    Area = 1.0
    }
""" %(solvent)

# molecule
CONF = psi4.geometry(qm_molblock, name=title)

psi4.set_options({
    "save_jk": True,
    "pcm": True,
    "pcm__input": pcm_string
})

e, wfn = psi4.energy("B3LYP-D3BJ/6-311+g(d,p)", return_wfn=True, molecule=CONF)
res = tdscf_excitations(wfn, states=100)


# get poles and residues to plot OPA and ECD spectra
poles = [r["EXCITATION ENERGY"] for r in res]
opa_residues = [np.linalg.norm(r["ELECTRIC DIPOLE TRANSITION MOMENT (LEN)"])**2 for r in res]
ecd_residues = [r["ROTATORY STRENGTH (LEN)"] for r in res]
print(poles)
print(opa_residues)
print(ecd_residues)
data = {'poles':poles,'opa':opa_residues,'ecd':ecd_residues}
df = pd.DataFrame(data)

#（1）光量子方程E=h·ν，h为普朗克常数，ν为频率；
#（2）光速方程c=λ·ν，λ为波长，ν为频率。
# 运算两个方程得到:
#                 E=hc/λ
# 而h=4.13×10-15 eV·s，c=3×1017 nm·s，故:
# hc=1240 eV·nm。
# 因此：1240/λ(nm)=(能量) eV
#
# 其中poles列的单位是au,需要转化波长
# 1 au = 27.21139664 eV
# wave length (nM) = 1240/(27.21139664* 原子能量)
df['wavelength']=df['poles'].map(lambda x:45.56914209163503/x)
output_uv_ecd = title+'_output.csv'
df.to_csv(output_uv_ecd,index=False)

# save heat file
# convert hatree to Kcal/mol
heat = float(heat)*627.5
file_heat = open(title+'.heat', "w")
file_heat.write("HEAT OF FORMATION  =  %.11f kcal/mol" %(heat))
file_heat.close()
# save UV spectra for Specdis
uv_file=title+'.uv.bil'
uv_first_line = 'SpecDis bil-file (length formalism) UV\n'
uv = df[['wavelength','opa']]
uv.to_csv(uv_file,index=False,sep='\t',header=None)
insert_first_line(uv_file,uv_first_line)

# Save ECD spectra for Specdis
ecd_file=title+'.cd.bil'
ecd_first_line = 'SpecDis bil-file (length formalism) ECD\n'
ecd = df[['wavelength','ecd']]
ecd.to_csv(ecd_file,index=False,sep='\t',header=None)
insert_first_line(ecd_file,ecd_first_line)

