# icos
A new magnetofrictional code based on an icosahedral grid. Details of its development and uses to follow in my thesis.

Initialise the code using the python file 'run.py', in which you input the model parameters etc. It should run without modification except:
1. Change the location of your Fortran 90 compiler in the makefile (near the top). You will also require netcdf to be installed to output the data and transfer from Fortran to Python. The installation of netcdf is fraught with danger and frustration.
2. Change the location of the data directory (which is usually not in the same folder as the code as the files are large). This needs to be done in run.py, plot.py and /src/main.f90. I recommend not changing the length of the filename I have provided or Fortran will be unhappy with you.

If you're reading this you probably are familiar with the equations and parameters, but more information may follow in due course. The main limitation currently is the requirement of having a reasonably high ratio of coronal diffusion to magnetofrictional relaxation rate. Hopefully this can be fixed in due course with improvements to the numerical scheme. Maybe.

The python wrapper will check to see if a suitable initial condition exists, and if not it will calculate one. This can take a while for high resolutions - an hour or so at G = 6. Sorry about this.
It will then compile and run the code, and will print some diagnostic data as it runs.
Once the code is running, you can use 'plot.py' to plot diagnostics and the magnetic field, using a field line tracer (hidden in /fltrace/). 

More to follow, one hopes...

