#Python script to plot data outputted by the icosahedral code
#Can plot diagnostics, pyvista 3d plot (in any orientation), surface plot of closed field lines or a lonitudinal plane plot for axisymmetric ones
#Calls the fortran field line tracer in the folder /fltrace/, which contains many of the functions used in the main code (including the gridpoint averaging), but only runs on one core.
#Will by default produce animations every 10 plots

import os

from python import grid_lite as grid_data
import pyvista as pv

import numpy as np
from scipy.io import netcdf_file
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt, colors
from matplotlib import cm
import matplotlib
import shutil
import time
from matplotlib.collections import PolyCollection

run_id = 0    #To identify the run that you need to plot.

#Plotting options
vista = True        #Plots 3D pyvista
surface = False      #Plots the closed field lines on a Mollweide projection. A bit slow...
plane = False       #Plots in the longitudinal plane (good for the axisymetric ones)
diags = True       #Plots some diagnostics

plot_directory = './plots/'
this_directory = os.path.abspath(os.getcwd())
print(this_directory)
data_directory = '/extra/tmp/trcn27/icos_data/'

moving_frame = False #set to true to make the pyvista plot move in the carrington frame
nlines = 2500   #Number of field lines to be plotted from the surface
nopen = 180     #Number of open field lines to be plotted

dt_fact = 0.02  #'timestep' in the field line tracer - as a proportion of the size of each grid cell
rmax = 2.5
export_length = 500   #max number of points in each field line)
start = 0     #Plot number to begin with

init_view_angle = 0 #Angle of pyvista plot (in degrees from facing the Earth (I think))

parameters = np.loadtxt('parameters%02d.txt' % run_id)


if True and start == 0:
    if os.path.exists(plot_directory):
        for filename in os.listdir(plot_directory):
            os.remove(plot_directory + filename)

    if os.path.exists('./temp/'):
        for filename in os.listdir('./temp/'):
            os.remove('./temp/' + filename)
    else:
        os.mkdir('./temp/')

G = int(parameters[0])
rmax = parameters[1]
tmax = parameters[3]
nplots = int(parameters[13])
time_unit = parameters[15]
crot = int(parameters[11])

end = nplots-1


print('Running')
def initialise_grid(G, pts_init = [[0]]):
    '''
    Sets up the required data using the already-calculated tweaked grids (done up to 7), or from the specified pts_init.
    '''
    n = 2**G
    pts, triangs, nbrs, corners = grid_data.icosahedron()
    if len(pts_init) == 1:  #not specified, load tweaked points
        pts_init = np.loadtxt('./tweaked_points/tweakpoints%d.txt' % G)   #tweaked points, just the first segment

    all_pts = grid_data.rotate_pts(pts_init, pts)  #rotate the points in such a manner as to cover the whole sphere
    grids = [grid_data.Subdomain(pts[triangs], nbrs, corners, n, all_pts[seg], rmax, seg) for seg in range(20)]   #all the grids
    #Bit of a mess doing it here, but grid_whole only does one segment so this can't be paralellised.
    pole_map = np.zeros((12, 2, 5), dtype='int')
    for i in range(12):
        pole_map[i] = np.array(np.where(triangs == i))
    return grids, pole_map, triangs

grid, _, _ = initialise_grid(G)  #highest level grid calculated, using tweaked points.

def transform(xyz):
    #transforms proper coordinates onto a 2D plane
    if abs(1.0- xyz[2]) < 1e-4:  #on top
        return 0.0, np.pi/2
    if abs(-1.0- xyz[2]) < 1e-4:  #on top
        return 0.0, -np.pi/2
    u = np.arctan2(xyz[1],xyz[0])
    v = np.arctan2(xyz[2],np.sqrt(xyz[0]**2+xyz[1]**2))
    return u, v

verts = []
allcentres = []

for seg in range(20):
    for k in range(grid[0].ntri):
        vert = []
        centre = grid[seg].centres[k]
        uc, vc = transform(centre)  #should give quadrant of centre
        allcentres.append([uc, vc])
        #Make sure the triangle is in the same hemisphere as the centre of the triangle
        if uc >= 0 and vc >= 0:
            quad = 0; base = np.pi/4
        if uc < 0  and vc >= 0:
            quad = 1; base = 3*np.pi/4
        if uc <  0 and vc < 0:
            quad = 2; base = 5*np.pi/4
        if uc >= 0 and vc < 0:
            quad = 3; base = 7*np.pi/4
        for i in range(3):
            pt = grid[seg].pts[grid[seg].triangs[k,i]]
            u, v = transform(pt)
            if abs(u-uc) > 2*np.pi/5:
                if u > uc:
                    u = u - 2*np.pi
                else:
                    u = u + 2*np.pi
                if abs(u-uc) > 2*np.pi/5:
                    u = 0.0
            vert.append([u,v])
        verts.append(vert)

print('Plotting grid established in python (only done once for a series of plots)')

def plot_surface(grid, br, plot_num, vmax):
    #plots the surface grid and the radial magnetic field distribution
    #vmax = np.max(np.abs(br))
    vmin = -vmax
    cmap = cm.seismic
    norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
    map2 = cm.ScalarMappable(norm = norm, cmap = 'seismic')

    count = 0
    zdata = []
    for seg in range(20):
        for k in range(grid[0].ntri):
            zdata.append(br[seg,k,0])
            count = count + 1
    if G < 5:
        lw = 0.1
    if G == 5:
        lw = 0.05
    if G >= 6:
        lw = 0.02
    coll = PolyCollection(verts, array = zdata, cmap=cmap, edgecolors = 'black', linewidth = lw, clim = [vmin,vmax])
    plt.gca().add_collection(coll)
    plt.plot([-np.pi,np.pi,np.pi,-np.pi,-np.pi],[np.pi/2,np.pi/2,-np.pi/2,-np.pi/2,np.pi/2],c='black', linewidth = 4)
    fig.colorbar(map2,ax = plt.gca(),label = 'Radial Magnetic Field Strength')
    plt.xticks([])
    plt.title('Time = %.1f days' % (t*time_unit))
    plt.savefig('plots/a%d.png' % plot_num,dpi=200)
    print('Saved heatmap figure', plot_num)
    #plt.show()
    plt.close()

def plot_fieldlines(fieldlines):
    for i, line in enumerate(fieldlines[:]):
        #transform to the correct coordinates and plot
        doplot = True
        line = line[:-1]
        if np.sum(np.abs(line[0])) < 1e-10:
            doplot = False
        if np.sqrt(np.sum(np.array(line[-1])**2)) > 2.4:
            doplot = False
        if np.max(np.sqrt(np.sum(np.array(line)**2, axis = 1))) < 1.025:
            doplot = False

        if doplot:
            xs=  []; ys = []
            dist = 0
            for j, pt in enumerate(line):
                if np.sum(np.abs(pt)) < 1e-10 or j == len(line) - 1:  #This is the end of the line
                    if np.sum(line[j-1]**2) < 1e-10:
                        openclosed = -1
                    elif np.sum(line[j-1]**2) < 1.1:
                        openclosed = 0
                    else:
                        openclosed = 1
                    break

                u,v = transform(pt)
                if len(xs) > 0:
                    if abs(u-xs[-1]) > 0.1 or abs(v-ys[-1]) > 0.1:  #line has passed into the other hemisphere
                        if len(xs) > 1:
                            plt.plot(xs[:-1],ys[:-1],c='black',linewidth = 0.25)
                            xs = [u]; ys = [v]
                            dist = 0
                    else:
                        xs.append(u); ys.append(v)
                        dist = dist + abs(xs[-1]-xs[-2])
                else:
                    xs.append(u); ys.append(v)
                    dist = 0

            if len(xs) > 1 and openclosed == 0:
                plt.plot(xs,ys,c='black',linewidth = 0.25)

def plot_plane(fieldlines):
    fig = plt.figure(figsize = (7,7))

    norm = colors.Normalize(vmin=0, vmax=np.pi/2)
    cmap = plt.cm.Reds

    maxwind = 0
    for li, line in enumerate(fieldlines):
        ud = line[-1][0]
        line = line[:-1]
        xs = []; ys = []
        wind = 0  #distance around the sun longitudinally. Need to be clever with the breaks
        for i, pt in enumerate(line):
            if np.sum(pt**2) < 1e-10:
                break
            else:
                xs.append(np.sqrt(pt[0]**2 + pt[1]**2)); ys.append(pt[2])
                if i > 0:
                    a1 = np.arctan2(line[i-1][1],line[i-1][0])
                    a2 = np.arctan2(line[i][1],line[i][0])
                    wind = wind + (a2 - a1)
        if len(xs) > 1:
            if np.sqrt(np.sum(xs[0]**2 + ys[0]**2)) < 1.1:  #the ones starting from the bottom
                plt.plot(xs, ys, c = cmap(norm(abs(wind))), linewidth = 0.5)
                maxwind = max(wind,abs(maxwind))

    plt.plot(np.sin(np.linspace(0,np.pi,1000)), np.cos(np.linspace(0,np.pi,1000)), c= 'black', linewidth = 1.0)
    #plt.plot(2.5*np.sin(np.linspace(0,np.pi,1000)), 2.5*np.cos(np.linspace(0,np.pi,1000)), c= 'white', linewidth = 1.0)
    plt.fill_between(np.sin(np.linspace(0,np.pi/2,1000)), -np.cos(np.linspace(0, np.pi/2, 1000)),np.cos(np.linspace(0, np.pi/2, 1000)), color = 'darkorange')

    plt.gca().set_facecolor('grey')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    fig.colorbar(sm, ax=plt.gca())
    plt.axis('equal')
    plt.xlim(xmin=-0.0, xmax = 2.5)
    plt.title('Time = %.1f days' % (t*time_unit))
    plt.savefig('plots/b%d.png' % plot_num,dpi=200)
    print('Saved plane figure', plot_num)
    plt.close()

def plot_vista(grid, br, plot_num, camera, vmax):
    off_screen = True
    p = pv.Plotter(off_screen=off_screen)
    p.background_color = "black"
    for seg in range(20):
        faces_local = np.array([[3, t[0], t[1], t[2]]
                                for t in grid[seg].triangs[:grid[seg].ntri]])
        surf_local = pv.PolyData(grid[seg].pts[:grid[seg].npts], faces_local)
        p.add_mesh(surf_local, scalars=br[seg][:,0], show_edges=False, cmap='hot', clim = [-vmax, vmax])

    data = netcdf_file('flines.nc', 'r', mmap=False)   #plotted from the surface
    fieldlines = data.variables['lines'][:]
    data.close()
    for line in range(0,len(fieldlines)):
        pts = []
        for pt in fieldlines[line]:
            if np.sum(np.abs(pt)) > 1e-10:
                pts.append(pt)
            else:
                break
        if len(pts) > 2:
            if np.sqrt(np.sum(pts[-1]**2)) < 2.0:
            #if True:
                p.add_mesh(pv.Spline(pts, len(pts)),color='white',line_width=2.0)

    if nopen > 0:  #also ones from the top
        data = netcdf_file('flines_plane.nc', 'r', mmap=False)   #plotted from the top
        fieldlines = data.variables['lines'][:]
        data.close()

        for line in range(0,len(fieldlines)):
            pts = []
            for pt in fieldlines[line][:-1]:
                if np.sum(np.abs(pt)) > 1e-10:
                    pts.append(pt)
                else:
                    break
            if len(pts) > 2:
                if np.sqrt(np.sum(pts[-1]**2)) < 2.0:
                    p.add_mesh(pv.Spline(pts, len(pts)),color='white',line_width=2.0)

    theta = 4*np.pi/8
    r = 12.
    phi = view_angle*np.pi/180.
    p.remove_scalar_bar()
    p.camera.position = (r*np.sin(theta)*np.sin(phi), r*np.sin(theta)*np.cos(phi), r*np.cos(theta))
    p.camera.focal_point = (0,0,0)
    p.add_title('%d days' % t, font='times', color='white', font_size=40)
    if off_screen:
        #p.export_html('pv.html')
        p.show(screenshot='plots/plot%d.png' % plot_num,window_size=[3840, 2160])
    else:
        p.show()
    return camera

def plot_diags():
    fname = this_directory + '/diagnostics/run0.nc'

    fname = ('./diagnostics/run%02d.nc' % run_id)

    try:
        data = netcdf_file(fname, 'r', mmap=False)   #plotted from the surface
        data.close()
    except:
        print('Diagnostic file not found')
        return False

    data = netcdf_file(fname, 'r', mmap=False)   #plotted from the surface

    time = data.variables['time'][:]*time_unit
    oflux = data.variables['oflux'][:]
    lflux = data.variables['lflux'][:]
    energy = data.variables['energy'][:]
    sumj = data.variables['current squared'][:]
    v1 = data.variables['frictional velocity'][:]

    data.close()

    end = len(time)
    for i in range(1,len(time)):
        if time[i] < 1e-12:
            end = i-1
            break
    fig, ax = plt.subplots(2,3, figsize=(10, 5))

    ax[0,0].plot(time[:end],oflux[:end])
    ax[0,0].set_ylabel('Open Flux')
    ax[0,0].set_xlabel('Time (days)')

    ax[0,1].plot(time[:end],energy[:end])
    ax[0,1].set_ylabel('Energy')
    ax[0,1].set_xlabel('Time (days)')

    ax[1,0].plot(time[:end],lflux[:end])
    ax[1,0].set_ylabel('Total Flux')
    ax[1,0].set_xlabel('Time (days)')

    ax[1,1].plot(time[:end],sumj[:end])
    ax[1,1].set_ylabel('Total Current')
    ax[1,1].set_xlabel('Time (days)')
    #ax[1,1].set_yscale('log')


    ax[0,2].plot(time[:end],v1[:end])
    ax[0,2].set_ylabel('MF Velocity')
    ax[0,2].set_xlabel('Time (days)')
    #ax[0,2].set_yscale('log')

    plt.suptitle('Diagnostics')
    plt.tight_layout()
    #plt.show()
    plt.savefig('diagplots/diags%02d.png' % run_id)
    plt.close()

for plot_num in range(start,end + 1):
    print('_______________________________')
    print('Plotting snapshot', plot_num)
    print('_______________________________')

    if plot_num == start:
        camera = [] #For the plotter
        if os.path.exists('flparameters.txt'):
            os.remove('flparameters.txt')

    filename = data_directory + '%02d' % run_id + '/%04d.nc' % (plot_num)
    if not os.path.exists(filename):
        while not os.path.exists(filename):
            time.sleep(10.0)

    time.sleep(2.0)  #Let netcdf catch up

    data = netcdf_file(filename, 'r', mmap=False)

    t = plot_num*(tmax/float(nplots-1))

    br = data.variables['br'][:]
    bh = data.variables['bh'][:]

    data.close()

    if diags:
        print('Plotting diagnostics')
        plot_diags()
        print('Diagnostics plotted')

    if plot_num == start:
        vmax = np.max(abs(br[:,:,0]))

    if moving_frame:
        view_angle = init_view_angle + t*360.0/27.3   #Rotates the plot at the Carrington rotation speed
    else:
        view_angle = init_view_angle   #Angle at which to plot any open field lines from (in degrees)

    flparas = np.zeros(20)
    flparas[0] = plot_num
    flparas[1] = G
    flparas[2] = nlines
    flparas[3] = dt_fact
    flparas[4] = crot
    flparas[5] = rmax
    flparas[6] = export_length
    flparas[7] = view_angle
    flparas[8] = 0
    flparas[9] = run_id
    flparas[10] = nopen

    fl_directory = this_directory + '/fltrace/'

    while os.path.exists('flparameters.txt'):  #don't want to overwrite things
        time.sleep(0.1)

    np.savetxt('flparameters.txt', np.array(flparas), delimiter = ',')

    if plot_num == 0:
        os.system('make -C' + fl_directory)

    os.system('/usr/lib64/openmpi/bin/mpiexec -n 1 ' + fl_directory + 'bin/fltrace')

    if os.path.exists('flparameters.txt'):
        os.remove('flparameters.txt')

    time.sleep(1.0)

    data = netcdf_file('flines.nc', 'r', mmap=False)
    fieldlines = data.variables['lines'][:]

    data.close()

    vmax = np.max(np.abs(br[:,:,0]))


    if surface:
        fig= plt.figure(figsize = (10,6))
        plt.subplot(111,projection="mollweide")
        print('Field lines calculated, plotting...')
        plot_fieldlines(fieldlines)
        print('Field lines plotted')
        plot_surface(grid, br, plot_num, vmax)
        print('Surface plotted')
        if plot_num%(10) == 0:
            os.system('ffmpeg -y -framerate 20 -i ./plots/a%d.png -b:v 10M  mov_a.mp4')

    if plane:
        print('Plotting in-plane magnetic field')
        data = netcdf_file('flines_plane.nc', 'r', mmap=False)
        fieldlines = data.variables['lines'][:]
        plot_plane(fieldlines)
        data.close()
        print('Finished plotting')
        if plot_num%(10) == 0:
            os.system('ffmpeg -y -framerate 20 -i ./plots/b%d.png -b:v 10M  mov_b.mp4')
    if vista:
        print('Plotting pyvista')

        camera = plot_vista(grid, br, plot_num, camera, vmax)

        print('Pyvista plotted')
        if plot_num%10 == 0 and plot_num > 1:
            os.system('ffmpeg -y -framerate 20 -i ./plots/plot%d.png -b:v 10M  mov_c.mp4')

    os.remove('flines.nc')
    os.remove('flines_plane.nc')


