#Python script to compile and run the fortran code. Will include the initial condition generation (for now in python) and any plotting, using as much fortran as is practicable

import os

import numpy as np
from python import pfield_3d
import shutil
import time


######################################
#SYSTEM PARAMETERS

G = 5              #Grid resolution. Corresponds to (2**G)^3 cells per segment. Integers from 3 to 6 are tested.
rmax = 2.5         #Maximum radius in solar radii. Only tested at rmax = 2.5.
nprocs = 4         #Number of processors. Must be 1,2,4,5,10,20,40 or 80. NOT TESTED PROPERLY FOR MORE THAN 20: Runs fine but diagnostics might go wrong.
run_id = 0         #Run id, keeping output data in separate folders etc.

data_directory = '/extra/tmp/trcn27/icos_data/%02d/' % run_id   #Data directory location. IF THIS CHANGES - MAKE SURE FORTRAN KNOWS TOO AND HAS THE RIGHT CHARACTER LENGTH! (in src/import_export.f90)

if True:
    if os.path.exists(data_directory):
        for filename in os.listdir(data_directory):
            os.remove(data_directory + filename)
    else:
        os.mkdir(data_directory)

ndiags = 1000        #Number of diagnostic outputs (open flux etc.)
nplots = 10          #Number of full magnetic field saves (which can be rather large)

# CODE PARAMETERS - in code units unless otherwise specified (days/solar radii)

rsun = 1.0                 #Radius of the sun in code units
hamilton_flag = 0           #Set to 1 if using hamilton (outputs to /nobackup directory). Also need to change output directories etc.
constant_friction = 1       #Set to 1 for constant friction, 0 for friction varying with latitude and altitude (recommend reducing cfl if this is the case)

time_unit = 1.0                #Code time unit in days
nu0 = 0.05           #Magnetofriction relaxation rate
eta = 5e-3*nu0     #Coronal diffusion rate
eta0 = 10.0*nu0     #Supergraunular diffusion rate

crot = 2160                #Carrington rotation of initial condition, or specified lower boundary condition, set in python/pfield_3d.py. To use function 1, use crot = 1001 etc.

shear_frame = 2   #1 = background stars, 2 = carrington frame
voutfact = 50.0*86400/696000   #Outflow speed
shearfact = 1.0                #Shearing rate (set to 1.0 for the correct rate in code units)
tmax = 1.0/time_unit          #Maximum time (in days)

#different cfls are allowed depending on whether friction or diffusion dominates. Can be higher for nu0, for instance.
cfl = 1.0                      #cfl condition factor (set to 1.0 seems to be OK)
delta = 1e-6                     #Magnetofrictional smoothing term

print('nu0 = ', nu0)
print('eta0 = ', eta0)
print('eta = ', eta)

def save_parameters():
    paras = -1.0*np.ones((20))
    paras[0] = G; paras[1] = rmax; paras[2] = nprocs; paras[3] = tmax
    paras[4] = eta; paras[5] = nu0; paras[6] = eta0; paras[7] = voutfact; paras[8] = shearfact
    paras[9] = cfl; paras[10] = delta
    paras[11] = crot; paras[12] = ndiags; paras[13] = nplots
    paras[14] = run_id; paras[15] = time_unit
    paras[16] = shear_frame
    paras[17] = rsun
    paras[18] = hamilton_flag
    paras[19] = constant_friction
    np.savetxt('parameters%02d.txt' % run_id, paras, delimiter = ',')

save_parameters()  #Save out parameters, to be read in by python

init_filename = './inits/pfield%d_%d_%d.nc' % (G, 0, crot)   #filename of the initial condition

print('Initial Condition filename', init_filename)

#Create data directory
if os.path.exists(data_directory):
    shutil.rmtree(data_directory)
    os.mkdir(data_directory)
else:
    os.mkdir(data_directory)

#Generate initial condition (if one does not already exist) and output as netcdf into the 'inits' folder.
if not os.path.exists(init_filename):
    print('Calculating initial potential field using multigrid in python')
    pfield_3d.pfss(G,rmax, import_data=True,lim = 1e-12, crot = crot)        #lim is the l2 norm that the potential field needs to converge to
    print('Initial condition calculated and saved to file', init_filename)
else:
    print('Initial condition already exists, filename', init_filename)

#Compile and run the code
print('Compiling Fortran')
os.system('make')


print('Running code with', nprocs, 'processes')
print('')
print('Output directory acording to python:', data_directory)
time.sleep(1.0)

if nprocs <= 4:
    os.system('/usr/lib64/openmpi/bin/mpiexec -n %d ./bin/mf3d %d' % (nprocs, run_id))
else:
    os.system('/usr/lib64/openmpi/bin/mpiexec -n %d --oversubscribe ./bin/mf3d %d' % (nprocs, run_id))





