# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:57:24 2022

#Script for calculating potential fields, as initial conditions for the icosahadral code
#Initialised using run.py, and uses either boundary data from HMI (downloaded if necessary), or a prescribed boundary function (set in init)
#Uses the multigrid scheme defined in my thesis

@author: eleph
"""

import numpy as np
from python import grid_pfield as grid_whole
np.set_printoptions(precision=16)
import time
import matplotlib.pyplot as plt

from scipy.interpolate import RegularGridInterpolator

import os
import drms
from astropy.io import fits
from scipy.io import netcdf_file
import time

class pfss:
    '''
    Class for all the multigrid functions, designed to solve Laplace's equation in 2D, on the surface of the sphere at grid cetnres
    '''
    def __init__(self, G_max, rmax, rhs_fn = 0, import_data = False, lim = 1e-10, crot = 2148):
        '''
        Requires the maximum grid level and the right hand side function (arguments to be determined)
        '''
        #INPUT LOWER BOUNDARY FUNCTIONS HERE IF NECESSARY, USE crot = 1000, 1001, 1002 etc. to use these
        def fn0(pt, seg):
            return 0

        def fn1(pt, seg):
            r = np.sqrt(np.sum(pt**2))
            s = pt[2]/r
            d1 = s - np.cos(0.4*np.pi)
            return s**7# + 100*d1*np.exp(-100*d1**2) 
    
        def fn2(pt, seg):
            r = np.sqrt(np.sum(pt**2))
            s = pt[2]/r
            theta = np.arccos(s)
            return (np.cos(4*theta - np.pi/2)) + s**7

        def fn3(pt, seg):
            r = np.sqrt(np.sum(pt**2))
            s = pt[2]/r
            theta = np.arccos(s)
            d1 = s - np.cos(0.35*np.pi)
            return s**7 + 5*d1*np.exp(-10*d1**2)# + 100*d2*np.exp(-100*d2**2)
    
        def fn4(pt, seg):
            r = np.sqrt(np.sum(pt**2))
            s = pt[2]/r
            theta = np.arccos(s)
            d1 = s - np.cos(0.3*np.pi)
            d2 = s - np.cos(0.55*np.pi)
            return s**7 + 100*d1*np.exp(-100*d1**2) + 100*d2*np.exp(-100*d2**2)

        def fn5(pt, seg):   #for testing rotation speed
            nstripes = 10
            phi = np.arctan2(pt[1],pt[0])

            return 2*(int(10*phi/(2*np.pi))%2)-1

        G = G_max
        self.G_max = G_max; self.rmax = rmax
        self.all_grids = self.generate_all_grids()  

        self.omegas = [1.0,1.41,1.65,1.75,1.49,1.35,1.28,1.25]   #SOR omega parameters

        self.rhs = self.set_scalar_faces(fn0, G)[:,:,:self.all_grids[G-1][0].ntri]
        if crot > 1500:
            print('Using downloaded boundary data, Carrington Rotation', crot)
            self.lower_bound = self.readmap(G, crot, smooth = 0.01/(2**G))
        else:
            bound_fn = crot-1000
            print('Using boundary function', bound_fn)
            if bound_fn == 0:
                self.lower_bound = self.set_lower_bound(fn0, G)
            if bound_fn == 1:
                self.lower_bound = self.set_lower_bound(fn1, G)
            if bound_fn == 2:
                self.lower_bound = self.set_lower_bound(fn2, G)
            if bound_fn == 3:
                self.lower_bound = self.set_lower_bound(fn3, G)
            if bound_fn == 4:
                self.lower_bound = self.set_lower_bound(fn4, G)
            if bound_fn == 5:
                self.lower_bound = self.set_lower_bound(fn5, G)

        self.lower_bound = self.correct_flux_balance(G, self.lower_bound)
        self.u = self.MGV(G, lim)
        self.br, self.bh = self.find_b_field(G, self.u)
        
        if not os.path.exists('./inits/'):
            os.mkdir('./inits/')
            
        grid = self.all_grids[0][0]
        
        do_netcdf = True
        save_id = 0
        if do_netcdf:
            netcdf_filename = './inits/pfield%d_%d_%d.nc' % (G_max, save_id, crot)
            """
            Saves out the magnetic fields. Only two arrays needed (coordinates should just happen anyways and will complicate things).
            All done in the fortran numbering scheme to save messing around when it gets there.
            """
            grid = self.all_grids[G_max-1][0]
            nr = grid.nr
            ntri = grid.ntri
            
            fid = netcdf_file(netcdf_filename, 'w')
            
            fid.createDimension('rs', nr+1)
            fid.createDimension('rc', nr)

            fid.createDimension('tri', ntri)   #has dimension but doesn't need a variable
            fid.createDimension('nsides', 3)    #three coordinate directions. There's probably a proper way of doing this but I don't know what it is.
            fid.createDimension('nsegs', 20)    #three coordinate directions. There's probably a proper way of doing this but I don't know what it is.

            vid = fid.createVariable('rs', 'd', ('rs',))
            vid[:] = grid.rs[1:-1]
            vid = fid.createVariable('rc', 'd', ('rc',))
            vid[:] = grid.rc[1:-1]
            
            
            vid = fid.createVariable('br', 'd', ('nsegs','tri','rs'))
            br_export = np.swapaxes(np.swapaxes(self.br[:,1:-1,:ntri], 1, 2), 0, 0)
            
            vid[:] = br_export

            vid = fid.createVariable('bh', 'd', ('nsegs','tri','nsides','rc'))
            
            bh_export = np.swapaxes(np.swapaxes(self.bh[:,1:-1,:ntri], 1, 2), 2, 3)
            vid[:] = bh_export

        
            fid.close()
            print('Wrote B to netcdf file '+ netcdf_filename)
            
         
        if not do_netcdf:
            print(np.shape(self.br))
            print(np.shape(self.bh))
            np.save('./inits/%dbr%d.npy' % (crot, G_max), self.br)
            np.save('./inits/%dbh%d.npy' % (crot, G_max), self.bh)

    
    def readmap(self, G, rot, smooth=0):
        """
            Reads the synoptic map for Carrington rotation rot.
            Also reads in the neighbouring maps, and puts them together for smoothing.
            Interpolates onto the icosaheadral grid and corrects for flux.
            
            ARGUMENTS:
                rot is the number of the required Carrington rotation (e.g. 2190)
                ns and nph define the required grid (e.g. 180 and 360)
                smooth [optional] controls the strength of smoothing (default 0 is no smoothing)
            
            [Important: the output map is corrected for flux balance, but the same correction is applied to the
            two neighbouring maps (for continuity), so they are not balanced.]
        """

        if not os.path.exists('./hmi/'):
            os.mkdir('./hmi/')
        if not os.path.exists('./hmi/%d.npy' % rot):
            #data not already downloaded. Download it...
            if True:   #Automatic thingy isn't working
                c = drms.Client()
                seg = c.query(('hmi.synoptic_mr_polfil_720s[%4.4i]' % rot), seg='Mr_polfil')
                segr = c.query(('hmi.synoptic_mr_polfil_720s[%4.4i]' % (rot-1)), seg='Mr_polfil')
                segl = c.query(('hmi.synoptic_mr_polfil_720s[%4.4i]' % (rot+1)), seg='Mr_polfil')

                with fits.open('http://jsoc.stanford.edu' + seg.Mr_polfil[0]) as fid:
                    brm = fid[1].data
                with fits.open('http://jsoc.stanford.edu' + segl.Mr_polfil[0]) as fid:
                    brm_l = fid[1].data
                with fits.open('http://jsoc.stanford.edu' + segr.Mr_polfil[0]) as fid:
                    brm_r = fid[1].data
            else:
                with fits.open('./hmi/CR%04d.fits' % rot) as fid:
                    brm = fid[0].data
                with fits.open('./hmi/CR%04d.fits' % (rot-1)) as fid:
                    brm_l = fid[0].data
                with fits.open('./hmi/CR%04d.fits' % (rot+1)) as fid:
                    brm_r = fid[0].data
                brm[np.isnan(brm)] = 0.0
                brm_l[np.isnan(brm_l)] = 0.0
                brm_r[np.isnan(brm_r)] = 0.0

            # Stitch together:
            brm3 = np.concatenate((brm_l, brm, brm_r), axis=1)

            del(brm, brm_l, brm_r)
            np.save('./hmi/%d.npy' % rot, brm3)
        else:
            brm3 = np.load('./hmi/%d.npy' % rot)
        # Remove NaNs:
        #brm3 = np.nan_to_num(brm3)
        nsm = np.size(brm3, axis=0)
        npm = np.size(brm3, axis=1)
        dsm = 2.0/nsm
        dpm = 2*np.pi/npm
        scm = np.linspace(-1 + 0.5*dsm, 1 - 0.5*dsm, nsm)  
        pcm = np.linspace(0.5*dpm, 2*np.pi - 0.5*dpm, npm)  
        
        def plgndr(m,x,lmax):
            """
                Evaluate associated Legendre polynomials P_lm(x) for given (positive)
                m, from l=0,lmax, with spherical harmonic normalization included.
                Only elements l=m:lmax are non-zero.
                
                Similar to scipy.special.lpmv except that function only works for
                small l due to overflow, because it doesn't include the normalization.
            """
            nx = np.size(x)
            plm = np.zeros((nx, lmax+1))
            pmm = 1
            if (m > 0):
                somx2 = (1-x)*(1+x)
                fact = 1.0
                for i in range(1,m+1):
                    pmm *= somx2*fact/(fact+1)
                    fact += 2
            
            pmm = np.sqrt((m + 0.5)*pmm)
            pmm *= (-1)**m
            plm[:,m] = pmm
            if (m < lmax):
                pmmp1 = x*np.sqrt(2*m + 3)*pmm
                plm[:,m+1] = pmmp1
                if (m < lmax-1):
                    for l in range(m+2,lmax+1):
                        fact1 = np.sqrt(((l-1.0)**2 - m**2)/(4.0*(l-1.0)**2-1.0))
                        fact = np.sqrt((4.0*l**2-1.0)/(l**2-m**2))
                        pll = (x*pmmp1 - pmm*fact1)*fact
                        pmm = pmmp1
                        pmmp1 = pll
                        plm[:,l] = pll
            return plm
            
        print('Data imported. Smoothing...')
        # (2) SMOOTH COMBINED MAP WITH SPHERICAL HARMONIC FILTER
        # ------------------------------------------------------
        if (smooth > 0):
            # Azimuthal dependence by FFT:
            brm3 = np.fft.fft(brm3, axis=1)
    
            # Compute Legendre polynomials on equal (s, ph) grid,
            # with spherical harmonic normalisation:
            lmax = 2*int((3*2**G-1)/2)  # note - already lower resolution
            nm = 2*lmax+1  # only need to compute this many values
            plm = np.zeros((nsm, nm, lmax+1))
            for m in range(lmax+1):
                plm[:,m,:] = plgndr(m, scm, lmax)
            plm[:,nm-1:(nm-lmax-1):-1,:] = plm[:,1:lmax+1,:]
            
            # Compute spherical harmonic coefficients:
            blm = np.zeros((nm,lmax+1), dtype='complex')
            for l in range(lmax+1):
                blm[:lmax+1,l] = np.sum(plm[:,:lmax+1,l]*brm3[:,:lmax+1]*dsm, axis=0)
                blm[lmax+1:,l] = np.sum(plm[:,lmax+1:,l]*brm3[:,-lmax:]*dsm, axis=0)
                # Apply smoothing filter:
                blm[:,l] *= np.exp(-smooth*l*(l+1))
    
            # Invert transform:
            brm3[:,:] = 0.0
            for j in range(nsm):
                brm3[j,:lmax+1] = np.sum(blm[:lmax+1,:]*plm[j,:lmax+1,:], axis=1)
                brm3[j,-lmax:] = np.sum(blm[lmax+1:,:]*plm[j,lmax+1:,:], axis=1)
    
            brm3 = np.real(np.fft.ifft(brm3, axis=1))

        def smooth_rotations(brm3):   #uses data from adjacent maps to get rid of discontinuities. Need to generalise this for all times.
            pn = len(brm3[0])//3
            p1 = pn; p2 = 2*pn
            bl = brm3[:,:p1]; bc = brm3[:,p1:p2]; bu = brm3[:,p2:]
            b = bc*0
            for p in range(pn//2):  #first half
                b[:,p] = bc[:,p]*(p + pn//2)/pn - bu[:,p]*(p - pn//2)/pn
            for p in range(pn//2, pn):  #second half
                b[:,p] = bc[:,p]*(3*pn//2 - p)/pn + bl[:,p]*(p - pn//2)/pn
            return b
        
        brm3 = smooth_rotations(brm3)
        nsm = np.size(brm3, axis=0)
        npm = np.size(brm3, axis=1)
        dsm = 2.0/nsm
        dpm = 2*np.pi/npm
        scm = np.linspace(-1 + 0.5*dsm, 1 - 0.5*dsm, nsm)  
        pcm = np.linspace(0.5*dpm, 2*np.pi - 0.5*dpm, npm)       
        #plt.pcolormesh(brm3)
        #plt.show()
        lower_bound = self.interpolate_to_grid(G, scm, pcm, brm3)
        print('Correcting flux balance')
        lower_bound = self.correct_flux_balance(G, lower_bound)

        return lower_bound
        
    def initialise_grid(self, G, pts_init = [[0]]):   
        '''
        Sets up the required data using the already-calculated tweaked grids (done up to 7), or from the specified pts_init.
        '''
        n = 2**G
        pts, triangs, nbrs, corners = grid_whole.icosahedron()
        if len(pts_init) == 1:  #not specified, load tweaked points
            pts_init = np.loadtxt('./tweaked_points/tweakpoints%d.txt' % G)   #tweaked points, just the first segment
        
        all_pts = grid_whole.rotate_pts(pts_init, pts)  #rotate the points in such a manner as to cover the whole sphere
        grids = [grid_whole.Subdomain(pts[triangs], nbrs, corners, n, all_pts[seg], self.rmax, seg) for seg in range(20)]   #all the grids
        #Bit of a mess doing it here, but grid_whole only does one segment so this can't be paralellised.
        pole_map = np.zeros((12, 2, 5), dtype='int')   
        for i in range(12):
            pole_map[i] = np.array(np.where(triangs == i))
        return grids, pole_map, triangs

        
    def generate_all_grids(self): 
        '''
        Imports the tweaked grid at G_max level. Then removes various points to produce the reduced grids, down to G=1.
        This is necessary for doing multgrid things, and it only needs to be done once. 
        '''
        print('Generating all grids. This can take some time.')
        grid, _, _ = self.initialise_grid(self.G_max)  #highest level grid calculated, using tweaked points.
        all_grids = []   #each of the grids at various levels
        all_grids.append(grid)
        print('Grid', self.G_max, 'done.')
        G = self.G_max
        while G > 1:  #remove points until this gets down to zero. Need to set up the correct pts_init each time.
            G = G - 1
            grid_upper = all_grids[-1]   #all_grids has the number of cells in descending order. Flip at the end.
            n_l = 2**G
            npts_l = (n_l+1)*(n_l+2)//2
            pts_new = np.zeros((npts_l, 3))
            for u in range(0,n_l+1):              #run through points in the first segment
                for v in range(u,n_l+1):
                    k_upper = self.findk(2*u, 2*v, 2**(G+1))
                    k_lower = self.findk(u, v, n_l)
                    pts_new[k_lower] = grid_upper[0].pts[k_upper]
            grid, _, _ = self.initialise_grid(G, pts_init = pts_new)

            all_grids.append(grid)
            print('Grid', G, 'done.')
        all_grids.reverse()   #flip into ascending order
        print('All grids generated.')
        return all_grids
        
    def set_scalar_faces(self, fn, G):    
        #for a fn of the coordinates (triple), outputs a scalar array on the grid faces (not points)
        u = np.zeros((20, len(self.all_grids[G-1][0].rc), len(self.all_grids[G-1][0].centres)))
        for seg, grid in enumerate(self.all_grids[G-1]):
            for r in range(len(grid.rc)):
                for k in range(len(grid.centres)):
                    u[seg,r,k] = fn(grid.all_centres[r,k], seg)
        return np.array(u)            

    def set_lower_bound(self, fn, G):    
        #for a fn of the coordinates (triple), outputs a scalar array on the grid faces (not points)
        lb = np.zeros((20, self.all_grids[G-1][0].ntri))
        for seg, grid in enumerate(self.all_grids[G-1]):
            for k in range(grid.ntri):
                lb[seg,k] = fn(grid.centres[k], seg)
        return np.array(lb) 
    
    def findk(self, u, v, n):
        return int(((2*n + 3)*u - u**2)//2 + v - u)

    def findt(self, u,v,n,up):
        if up:
            return 2*(v-1) + 2*(n-1)*u - u**2
        else:
            return 2*v-1 + 2*(n-1)*u - u**2

    def transfer_faces(self, u, G, seg_transfer = np.arange(0,20,dtype='int16')):
        """
        Transfers ghost point data between adjacent segments. Radial ghost points are unaffected
        U can either be with ghosts or without, and a ghosted array is produced.
        """
        u = np.array(u)
        grids = self.all_grids[G-1]
        u = u[:,:,:grids[0].ntri]   #remove any existing ghost points
        u = np.append(u, np.zeros((20, grids[0].nr + 2, len(grids[0].centres) - grids[0].ntri)), axis = 2)
        for seg in seg_transfer:   #segments to transfer information to
            grid = grids[seg]
            for g in range(grid.ntri, len(grids[seg].centres)):
                ghost_map = grid.tri_transfer[g - grid.ntri]
                u[seg,:,g] = u[ghost_map[0], :, ghost_map[1]]
        return np.array(u)

    def restrict_faces(self, f, G):
        """
        Input data f on grid size G. 
        Restricts this to grid size G-1 just by taking the average of the four faces.
        The reverse of this is probably adequate for interpolation, but will check with Anthony.
        """
        n_u = 2**G
        n_l = 2**(G-1)  #n of the restricted segment
        nr = self.all_grids[G-2][0].nr
        new_f = np.zeros((20, nr+2, self.all_grids[G-2][0].ntri))
        def findt(u,v,n,updown):   #up is 1, down is 0
            if updown:
                return 2*(v-1) + 2*(n-1)*u - u**2
            else:
                return 2*v-1 + 2*(n-1)*u - u**2
        for seg in range(20):
            for r in range(1, nr+1):
                r_uppers = [2*r-1, 2*r]
                for u in range(0,n_l):
                    for v in range(u + 1, n_l + 1):   #there is an upward and downward triangle (on lower grid)
                        #upward triangle
                        k_lower = findt(u,v,n_l,1)
                        k_uppers = [findt(2*u,2*v-1,n_u,1),findt(2*u,2*v-1,n_u,0),findt(2*u,2*v,n_u,1),findt(2*u+1,2*v,n_u,1)]
                        new_f[seg][r, k_lower] = np.sum(f[seg][r_uppers[0], k_uppers])/8 + np.sum(f[seg][r_uppers[1], k_uppers])/8
                        if v < n_l:   #downward triangle
                            k_lower = findt(u,v,n_l,0)
                            k_uppers = [findt(2*u,2*v,n_u,0),findt(2*u+1,2*v+1,n_u,1),findt(2*u+1,2*v+1,n_u,0),findt(2*u+1,2*v,n_u,0)]
                            new_f[seg][r, k_lower] = np.sum(f[seg][r_uppers[0], k_uppers])/8 + np.sum(f[seg][r_uppers[1], k_uppers])/8
        return new_f
    
    def extend_faces(self, f, G):
        """
        Input data f on grid size G-1. 
        Converts this onto a grid size G. 
        No interpolation or averaging as I'm not sure how that could help without more information.
        """
        n_u = 2**(G+1)
        n_l = 2**G  
        nr = self.all_grids[G-1][0].nr
        temp_f = np.zeros((20, nr*2 + 2, len(self.all_grids[G][0].centres)))
        new_f = np.zeros((20, nr*2 + 2, len(self.all_grids[G][0].centres)))
        def findt(u,v,n,updown):   #up is 1, down is 0
            if updown:
                return 2*(v-1) + 2*(n-1)*u - u**2
            else:
                return 2*v-1 + 2*(n-1)*u - u**2
        for seg in range(20):
            for r in range(1, nr+1):
                r_uppers = [2*r-1, 2*r]
                for u in range(0,n_l):
                    for v in range(u + 1, n_l + 1):   #there is an upward and downward triangle (on lower grid)
                        #upward triangle
                        k_lower = findt(u,v,n_l,1)
                        k_uppers = [findt(2*u,2*v-1,n_u,1),findt(2*u,2*v-1,n_u,0),findt(2*u,2*v,n_u,1),findt(2*u+1,2*v,n_u,1)]
                        temp_f[seg][r_uppers[0], k_uppers] = f[seg][r, k_lower]
                        temp_f[seg][r_uppers[1], k_uppers] = f[seg][r, k_lower]
                        if v < n_l:   #downward triangle as well (usually both)
                            k_lower = findt(u,v,n_l,0)
                            k_uppers = [findt(2*u,2*v,n_u,0),findt(2*u+1,2*v+1,n_u,1),findt(2*u+1,2*v+1,n_u,0),findt(2*u+1,2*v,n_u,0)]
                            temp_f[seg][r_uppers[0], k_uppers] = f[seg][r, k_lower]
                            temp_f[seg][r_uppers[1], k_uppers] = f[seg][r, k_lower]
        #Populate ghost points, for averaging.
        temp_f = self.transfer_faces(temp_f, G+1)
        #Smooth the new function, taking a simple mean of the surrounding cells
        for seg in range(20):
            for r in range(1, 2*nr+1):
                for kt in range(self.all_grids[G][0].ntri): #triangles in higher resolution
                    new_f[seg][r][kt] += (np.sum(temp_f[seg][r][self.all_grids[G][0].nbrfaces[kt]]) + temp_f[seg][r][kt])/4
        return new_f

    def l2_norm(self, G, scalar):
        #Returns l2 norm of a scalar field on the grid, using the Heikes formula and the volumes
        lsum = 0
        nr = 2**G
        for seg in range(20):
            grid = self.all_grids[G-1][seg]
            lsum += np.sum(grid.cell_volumes[1:nr+1,:grid.ntri]*scalar[seg][1:nr+1,:grid.ntri]**2)
        lsum = lsum**0.5
        lsum = lsum/(20*np.sum(grid.cell_volumes[1:nr+1,:grid.ntri]))
        return lsum

    def lap_faces(self, u, G, lower_bound, segs = np.arange(0,20,dtype='int16')):
        """
        Compute discrete laplacian of array u at grid faces.
        """
        u = self.transfer_faces(u, G)
        for seg in range(20):
            grid = self.all_grids[G-1][seg]   #the grids for the required resolution
            for k in range(grid.ntri):
                u[seg][0,k] = u[seg][1,k] - lower_bound[seg][k]*(grid.rc[1] - grid.rc[0])
                u[seg][grid.nr+1,k] = -(grid.rc[grid.nr+1]/grid.rc[grid.nr])*u[seg][grid.nr,k]

        lap = np.zeros((20, self.all_grids[G-1][0].nr+2, self.all_grids[G-1][0].ntri))
        for seg in segs:   #G-1 as this starts from zero. I think this makes more intuitive sense.
            grid = self.all_grids[G-1][seg]
            for r in range(1, self.all_grids[G-1][0].nr+1):
                for k in range(grid.ntri):
                    for j in range(3):   #Vertical faces first
                        dj = grid.tri_d_all[r,k,j]
                        Aj = grid.v_areas[r,k,j]
                        lap[seg][r,k] += Aj*(u[seg][r,grid.nbrfaces[k,j]] - u[seg][r,k])/dj
                    #lower one
                    dj = grid.rc[r] - grid.rc[r-1]
                    Aj = grid.h_areas[r,k]
                    lap[seg][r,k] += Aj*(u[seg][r-1,k] - u[seg][r,k])/dj
                    #upper one
                    dj = grid.rc[r+1] - grid.rc[r]
                    Aj = grid.h_areas[r+1,k]
                    lap[seg][r,k] += Aj*(u[seg][r+1,k] - u[seg][r,k])/dj
                    #print(r,k,Aj*(u[seg][r+1,k] - u[seg][r,k])/dj)
                    lap[seg][r,k] /= grid.cell_volumes[r, k]
        return np.array(lap)

    def lap_mat(self, u, G):
        """
        Does the non-boundary bit of the Laplacian matrix
        """
        u = self.transfer_faces(u, G)
        Au = np.zeros((20, self.all_grids[G-1][0].nr+2, self.all_grids[G-1][0].ntri))
        for seg in range(20):   #G-1 as this starts from zero. I think this makes more intuitive sense.
            grid = self.all_grids[G-1][seg]
            for r in range(1, self.all_grids[G-1][0].nr+1):
                for k in range(grid.ntri):
                    for j in range(3):   #Vertical faces first
                        dj = grid.tri_d_all[r,k,j]
                        Aj = grid.v_areas[r,k,j]
                        Au[seg][r,k] += Aj*(u[seg][r,grid.nbrfaces[k,j]] - u[seg][r,k])/dj

                    if r in range(2,grid.nr+1):
                        #lower one
                        dj = grid.rc[r] - grid.rc[r-1]
                        Aj = grid.h_areas[r,k]
                        Au[seg][r,k] += Aj*(u[seg][r-1,k] - u[seg][r,k])/dj

                    if r in range(1,grid.nr):
                        #upper one
                        dj = grid.rc[r+1] - grid.rc[r]
                        Aj = grid.h_areas[r+1,k]
                        Au[seg][r,k] += Aj*(u[seg][r+1,k] - u[seg][r,k])/dj
                    else:
                        dj = grid.rc[r+1] - grid.rc[r]
                        Aj = grid.h_areas[r+1,k]
                        Au[seg][r,k] += Aj*(-(1.0 + grid.rc[r+1]/grid.rc[r])*u[seg][r,k])/dj
                    Au[seg][r,k] /= grid.cell_volumes[r, k]
        return np.array(Au)

    def b_vector(self,lower_bound,G):
        """
        Generates the rhs vector using the lower boundary condition
        """
        b = np.zeros((20, self.all_grids[G-1][0].nr+2, self.all_grids[G-1][0].ntri))
        for seg in range(20):
            grid = self.all_grids[G-1][seg]
            #r = 1 here - only affects the bottom boundary
            r = 1
            for k in range(grid.ntri):
                dj = grid.rc[r] - grid.rc[r-1]
                Aj = grid.h_areas[r,k]
                b[seg][r,k] -= Aj*(lower_bound[seg][k]*(grid.rc[1] - grid.rc[0]))/dj
                b[seg][r,k] /= grid.cell_volumes[r, k]
        return b

    def MGV(self, G_max, lim = 1e-10):
        #Full multigrid function converging to the l2 norm limit (or could test l inf I suppose)
        lower_bound_init = self.lower_bound
        error = 1e10
        rhs = -self.b_vector(lower_bound_init,G_max)
        t0 = time.time()
        self.id_nits = [[0,0],[0,0],[9,16],[10,18],[2,11],[2,5],[2,5],[2,5]]
        f = 0.*rhs
        rem = rhs - self.lap_mat(f,G_max)
        lap_init = self.l2_norm(G_max,self.lap_faces(f, G_max, lower_bound_init))
        t0 = time.time()
        while True:
            if False:
                eps = self.f_cycle(G_max, 0*f, rem)
            else:
                print('Starting V-cycle')
                eps = self.v_cycle(G_max, 0.*f, rem, just_sor = False)
            f = f + eps
            rem = rhs - self.lap_mat(f,G_max)
            frem = self.lap_faces(f,G_max,lower_bound_init)
            print('Time = ', time.time() - t0)
            lap_after = self.l2_norm(G_max,self.lap_faces(f, G_max, lower_bound_init))
            t1 = time.time()
            print('Lalpacian Norm', lap_after)
            alpha = np.log(lap_after/lap_init)/(t1-t0)
            t_converge = t1 + (1./alpha)*np.log(lim/lap_after)
            print('Predicted remaining convergence time', t_converge-t1, 'seconds')
            print('Sorry this is so slow... Fortran will be better but hard to code.')
            if lap_after < lim:
                break
        print('Solution found, L2 error', lap_after, 'in', time.time() - t0, 'seconds')
        return f

    def f_cycle(self, G_max, f, rhs):
        '''
        Try to do this a bit more neatly. Each F-cycle using multgrid goodness
        '''
        print('F cycle level', G_max)
        #print(G_max, 'Initial error', self.l2_norm(G_max, rhs - self.lap_mat(0*f,G_max)))

        G = G_max
        if G == 1:
            f = self.v_cycle(G,f, rhs)
        else:
            #Restrict and run F-cycle at the lower resolution, then do a v-cycle at this resolution
            f = self.smooth_faces(G, f, rhs, nits = self.id_nits[G][0])[:,:,:(2**(G))**2]  #pre-condition
            rem = rhs - self.lap_mat(f,G)
            rem = self.restrict_faces(rem, G)
            G = G - 1
            eps = self.f_cycle(G,0*rem,rem)  #call itself because why not?
            eps = self.extend_faces(eps,G)[:,:,:(2**(G+1))**2]
            G = G + 1
            f = f + eps

            f = self.smooth_faces(G, f, rhs, nits = self.id_nits[G][0])[:,:,:(2**(G))**2]
            rem = rhs - self.lap_mat(f,G)
            rem = self.restrict_faces(rem, G)
            G = G - 1
            eps = self.v_cycle(G,0*rem,rem)  #call itself because why not?
            eps = self.extend_faces(eps,G)[:,:,:(2**(G+1))**2]
            G = G + 1

            f = f + eps
            f = self.smooth_faces(G, f, rhs, nits = self.id_nits[G][1])[:,:,:(2**(G))**2]

        print(G, 'Solution found, final error', self.l2_norm(G_max, rhs - self.lap_mat(f,G_max)))

        return f[:,:,:(2**G_max)**2]


    def v_cycle(self, G_top, f_init, rhs, nits_1 = 5, nits_2 = 5, just_sor = False):
        '''
        The v-cycle with G_top as its top level, going down to level 1 and back up again.
        '''
        f = f_init.copy()
        G = G_top
        if G == 1 or just_sor:   #just smooth directly using SOR
            f = self.smooth_faces(G, f_init, rhs, nits = 0)
        else:   #resrict to coarser grid and go from there
            f = self.smooth_faces(G, f_init, rhs, nits = self.id_nits[G][0])  #pre-condition
            rem = rhs - self.lap_mat(f, G)
            rem = self.restrict_faces(rem, G)
            G = G - 1
            eps = 0.*rem
            eps = self.v_cycle(G, eps, rem)
            eps = self.extend_faces(eps, G)
            G = G + 1
            f = f + eps
            f = self.smooth_faces(G, f, rhs, nits = self.id_nits[G][1])  #post-condition
        return f[:,:,:(2**G_top)**2]


    def smooth_faces(self, G, u, rhs, nits = 0, lim = 1e-13, check_freq = 10, omega = 1.65, iterate_order = [0]):
        '''
        nits iterations of the Gauss-Seidel process on faces
        if nits is set to zero then it will converge (almost) exactly
        '''
        if nits != 0:
            for i in range(nits):
                u = self.iterate_faces(G, u, rhs, omega = self.omegas[G])
        else:     #keep going until converged
            count = 0
            while True:
                u = self.iterate_faces(G, u, rhs, omega = self.omegas[G])
                error = self.l2_norm(G,self.lap_mat(u,G) - rhs)
                if error < lim:
                    break
                count += 1
        return u

    def iterate_faces(self, G, u, rhs, omega = 1.45, count = 0):  #one iteration of the smoothing. SOR method.
        #Needs u to be extended beforehand, or segments won't talk to each other.
        u = self.transfer_faces(u, G)   #obtain points from other segments that have just been iterated.
        grid = self.all_grids[G-1][0]   #the grids for the required resolution
        D_update = u[0]*0
        for j in range(3):
            D_update[1:grid.nr+1,:grid.ntri] -= grid.v_areas[1:grid.nr+1,:grid.ntri,j]/grid.tri_d_all[1:grid.nr+1,:grid.ntri,j]
        D_update[2:grid.nr+1,:grid.ntri] -= grid.h_areas[2:grid.nr+1,:grid.ntri]/(grid.rc[2:grid.nr+1,np.newaxis] - grid.rc[1:grid.nr,np.newaxis])
        D_update[1:grid.nr,:grid.ntri] -= grid.h_areas[2:grid.nr+1,:grid.ntri]/(grid.rc[2:grid.nr+1,np.newaxis] - grid.rc[1:grid.nr,np.newaxis])
        D_update[grid.nr,:grid.ntri] -= (1.0 + grid.rc[grid.nr+1]/grid.rc[grid.nr])*grid.h_areas[grid.nr+1,:grid.ntri]/(grid.rc[grid.nr+1,np.newaxis] - grid.rc[grid.nr,np.newaxis])

        for seg in range(20):  #run through the points and update them so the laplacian is zero at that point (momentarily)
            grid = self.all_grids[G-1][seg]   #the grids for the required resolution
            u = self.transfer_faces(u, G, [seg])  #obtain information from adjacent segments
            for r in range(1,grid.nr+1):   #interior points
                for k in range(grid.ntri):
                    #Contribution from boundary
                    b_update = rhs[seg][r,k]*grid.cell_volumes[r, k]
                    #Contribution from other cells
                    A_update = 0.0
                    for j in range(3):   #Vertical faces first
                        A_update += grid.v_areas[r,k,j]*(u[seg][r,grid.nbrfaces[k,j]])/grid.tri_d_all[r,k,j]
                    if r in range(2,grid.nr+1):
                        #lower one
                        A_update+= grid.h_areas[r,k]*(u[seg][r-1,k])/(grid.rc[r] - grid.rc[r-1])
                    if r in range(1,grid.nr):
                        #upper one
                        A_update += grid.h_areas[r+1,k]*(u[seg][r+1,k])/(grid.rc[r+1] - grid.rc[r])
                    #Contribution from THIS cell
                    u[seg][r,k] = omega*(b_update - A_update)/D_update[r,k] + (1-omega)*u[seg][r,k]
        return u

    def plotsnap(self, G, seg_plot = -1, u = [0], plot_num = -1, show=True, vmax = 0):   #plots grid level G with various options
        if show:
            p = pv.Plotter(off_screen = False)
        else:
            p = pv.Plotter(off_screen = True)
        if vmax == 0:
            vmax = np.max(np.abs(u))
            vmin = -np.max(np.abs(u))
        else:
            vmin = -vmax
        if seg_plot == -1:  #plot all segments
            for seg, grid in enumerate(self.all_grids[G-1]):
                faces_local = np.array([[3, t[0], t[1], t[2]] for t in grid.triangs[:grid.ntri]])
                surf_local = pv.PolyData(grid.pts[:grid.npts], faces_local)
                if len(u) == 1:
                    p.add_mesh(surf_local, show_edges=False, clim = [-1, 1])
                else:
                    if len(u[seg]) == len(self.all_grids[G-1][seg].triangs) or len(u[seg]) == self.all_grids[G-1][seg].ntri:  #Plot scalars at grid centres
                        p.add_mesh(surf_local, scalars = u[seg][:self.all_grids[G-1][seg].ntri], show_edges=False, clim = [vmin, vmax])
                    elif len(u[seg]) == len(self.all_grids[G-1][seg].pts) or len(u[seg]) == self.all_grids[G-1][seg].npts:  #Plot scalars at grid points
                        p.add_mesh(surf_local, scalars = u[seg][:self.all_grids[G-1][seg].npts], show_edges=False, clim = [vmin, vmax])
                    else:
                        raise Exception('Plot failed. Scalars do not fit on grid.')        
        else:
            seg = seg_plot; grid = self.all_grids[G-1][seg_plot]
            faces_local = np.array([[3, t[0], t[1], t[2]] for t in grid.triangs[:grid.ntri]])
            surf_local = pv.PolyData(grid.pts[:grid.npts], faces_local)
            if len(u) == 1:
                p.add_mesh(surf_local, show_edges=False, clim = [-1, 1])
            else:
                if len(u[seg]) == len(self.all_grids[G-1][seg].triangs) or len(u[seg]) == self.all_grids[G-1][seg].ntri:  #Plot scalars at grid centres
                    p.add_mesh(surf_local, scalars = u[seg][:self.all_grids[G-1][seg].ntri], show_edges=False, clim = [vmin, vmax])
                elif len(u[seg]) == len(self.all_grids[G-1][seg].pts) or len(u[seg]) == self.all_grids[G-1][seg].npts:  #Plot scalars at grid points
                    p.add_mesh(surf_local, scalars = u[seg][:self.all_grids[G-1][seg].npts], show_edges=False, clim = [vmin, vmax])
                else:
                    raise Exception('Plot failed. Scalars do not fit on grid.')
        if plot_num >= 0:
            p.show(screenshot='snaps/snap%d' % plot_num)
        else:
            p.show()
                
    def interpolate_to_grid(self, G, scm, pcm, brm3):
        #Reads in cartesian grid data and interpolatees onto the triangular grid
        print('Interpolating read-in data onto triangular mesh. May take a while as the new scipy function is horrendously slow for some reason...')
        lb = np.zeros((20, self.all_grids[G-1][0].ntri))
        fn = RegularGridInterpolator((pcm, scm), brm3.T, method='linear', bounds_error=False,fill_value=None)
        #interpolation function
        for seg, grid in enumerate(self.all_grids[G-1]):
            pts = []
            #print(seg, '/20 done. Apologies for the wait.')
            for k in range(grid.ntri):
                #FLIP THE COORDINATE HERE or it reads in incorrectly
                p = -np.arctan2(grid.centres[k,0], grid.centres[k,1]) + np.pi
                s = grid.centres[k,2]   #the z coordinate
                pts.append([p,s])
            lb[seg][:] = fn(pts)
        return lb

    def correct_flux_balance(self, G, f):
        """
        Correct the flux balance in the map f.
        Cells do not have equal area so this isn't trivial. Need to use face areas and adjust using them.
        """
        # Compute net flux through the sphere
        # Compute positive and negative fluxes:
        
        fluxp = 0
        fluxn = 0
        for seg, grid in enumerate(self.all_grids[G-1]):
            ipos = f[seg] > 0
            ineg = f[seg] < 0
            fluxp += np.sum(f[seg][ipos]*grid.areas[:grid.ntri][ipos])
            fluxn += np.sum(f[seg][ineg]*grid.areas[:grid.ntri][ineg])
        a = (fluxp - fluxn)/(2*fluxp); b = -(fluxp - fluxn)/(2*fluxn)
        for seg, grid in enumerate(self.all_grids[G-1]):
            ipos = f[seg] > 0
            ineg = f[seg] < 0
            f[seg][ipos] *= a
            f[seg][ineg] *= b  
        return f
        
    def find_b_field(self, G, u):  
        #Calculates the magnetic FIELD STRENGTH through each of the faces.
        for seg in range(20):
            grid = self.all_grids[G-1][seg]   #the grids for the required resolution
            #Apply boundary condition to the field u
            for k in range(grid.ntri):
                u[seg][0,k] = u[seg][1,k] - self.lower_bound[seg][k]*(grid.rc[1] - grid.rc[0])
                u[seg][grid.nr+1,k] = -(grid.rc[grid.nr+1]/grid.rc[grid.nr])*u[seg][grid.nr,k]

        br = np.zeros((20, self.all_grids[G-1][0].nr + 3, self.all_grids[G-1][0].ntri))
        bh = np.zeros((20, self.all_grids[G-1][0].nr + 2, self.all_grids[G-1][0].ntri, 3))
        u = self.transfer_faces(u, G)
        ntri = self.all_grids[G-1][0].ntri
        rc = self.all_grids[G-1][0].rc
        for seg, grid in enumerate(self.all_grids[G-1][:]):
            for r in range(1, grid.nr + 2):
                br[seg][r] = (u[seg,r,:ntri] - u[seg,r-1,:ntri])/(rc[r]-rc[r-1])
            for r in range(0, grid.nr + 2):
                for k in range(grid.ntri):
                    for j in range(3):
                        bh[seg][r, k, j] = (u[seg,r,grid.nbrfaces[k, j]] - u[seg,r,k])/(grid.rc[r]*grid.tri_d[k, j])

        print('Total magfield', np.sum(br[:,1,:ntri]))
        magflux = 0
        for seg in range(20):
            magflux = magflux + np.sum(br[seg,1,:ntri]*grid.h_areas[1,:ntri])
        print('Total magflux', magflux)


        return br, bh

    def testplot(self, G, u, r_indices = [1], seg_indices = np.arange(0,20,1)):
        p = pv.Plotter(off_screen = False)
        vmin = np.min(u)
        vmax = np.max(u)
        for seg, grid in enumerate(self.all_grids[G-1][:]):
            for r in [1]:
                faces_local = np.array([[3, t[0], t[1], t[2]] for t in grid.triangs[:grid.ntri]])
                surf_local = pv.PolyData(grid.rs[r]*grid.pts[:grid.npts], faces_local)
                p.add_mesh(surf_local, scalars = u[seg][r,:grid.ntri], show_edges=False, cmap = 'hot', clim = [vmin, vmax])
        for seg_test in seg_indices:
            grid = self.all_grids[G-1][seg_test]; seg= seg_test
            for r in r_indices:
                faces_local = np.array([[3, t[0], t[1], t[2]] for t in grid.triangs[:grid.ntri]])
                surf_local = pv.PolyData(grid.rs[r]*grid.pts[:grid.npts], faces_local)
                p.add_mesh(surf_local, scalars = u[seg][r,:grid.ntri], show_edges=False, cmap = 'hot', clim = [vmin, vmax])
        p.show()
        
    def boundplot(self, G, bound):
        p = pv.Plotter(off_screen = False)
        for seg, grid in enumerate(self.all_grids[G-1][:]):
            faces_local = np.array([[3, t[0], t[1], t[2]] for t in grid.triangs[:grid.ntri]])
            surf_local = pv.PolyData(grid.pts[:grid.npts], faces_local)
            p.add_mesh(surf_local, scalars = bound[seg][:grid.ntri], show_edges=False, cmap = 'hot')
        p.show()
  
def fn1(pt, seg):
    r = np.sqrt(np.sum(pt**2))
    s = pt[2]/r
    return s**7
        

#mgrid = pfss(7,2.5, import_data=True,lim = 1e-10, crot = 2130, bound_fn = 0)



