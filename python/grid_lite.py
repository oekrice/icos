"""
    Routines for constructing icosahedral-hexagonal-pentagonal grid on 2d spherical surface.
    AR Yeates, Jan-2022.
"""
import numpy as np

def icosahedron():
    """
    Generate icosahedron inscribed in unit sphere.
    Returns array of points and array of triangles (with points ordered anti-clockwise).
    Also returns neighbours of each triangle and their corresponding orientations.
    """
    s, c = 2/np.sqrt(5), 1/np.sqrt(5)
    topPts = [(0,0,1)] + [(s*np.cos(i*2*np.pi/5.-np.pi/5), s*np.sin(i*2*np.pi/5.-np.pi/5), c) for i in range(5)]
    bottomPts = [(0,0,-1)] + [(s*np.cos(-i*2*np.pi/5.+4*np.pi/5), s*np.sin(-i*2*np.pi/5.+4*np.pi/5), -c) for i in range(5)]
    icoPts = topPts + bottomPts
    icoTriangs = [(0,i+1,(i+1)%5+1) for i in range(5)] +\
                 [(6,i+7,(i+1)%5+7) for i in range(5)] +\
                 [(i+1,(i+1)%5+1,(7-i)%5+7) for i in range(5)] +\
                 [(i+1,(7-i)%5+7,(8-i)%5+7) for i in range(5)]

    icoPts = np.array(icoPts)
    icoTriangs = np.array(icoTriangs)

    # Ensure that points in each triangle are ordered anti-clockwise:
    for triang in icoTriangs:
        p1 = icoPts[triang[1]]
        v10 = icoPts[triang[0]] - icoPts[triang[1]]
        v12 = icoPts[triang[2]] - icoPts[triang[1]]
        orx = p1[0]*(v12[1]*v10[2] - v12[2]*v10[1])
        ory = p1[1]*(v12[2]*v10[0] - v12[0]*v10[2])
        orz = p1[2]*(v12[0]*v10[1] - v12[1]*v10[0])
        orientn = orx + ory + orz
        if (orientn < 0):
            triang[0], triang[1] = triang[1], triang[0]
            
    # Identify the three neighbours of each triangle (which share an edge) and also which of the neighbour's edges it is:
    icoNbrs = np.zeros((20, 3, 2), dtype='int16')   #also has opposite triangles now as they're important
    # Triangles share an edge if they have two vertices in common.
    for k, triang in enumerate(icoTriangs):
        # Consider each edge in turn:
        for i in range(3):
            p0 = triang[i]
            p1 = triang[(i+1)%3]
            for j in range(20):
                otherTriang = icoTriangs[j]
                if ((j != k) & np.isin(p0, otherTriang) & np.isin(p1, otherTriang)):
                    icoNbrs[k,i,0] = j  # id of neighbour
                    icoNbrs[k,i,1] = np.where(otherTriang == p1)[0][0]  # which edge of nbr
                    
    icoCorners = -1*np.ones((20, 3, 4, 2), dtype = 'int16')  #the four other segments which share the segment points. Anticlockwise.
    for k, triang in enumerate(icoTriangs):
        for i in range(3):   #corner to populate
            p0 = triang[i]
            p1 = triang[(i+2)%3]
            for j in range(4):
                for t_check in range(20):
                    otherTriang = icoTriangs[t_check]
                    if (not np.isin(t_check, icoCorners[k,i,:,0])) and t_check != k:
                        if np.isin(p0, otherTriang) and np.isin(p1, otherTriang):
                            icoCorners[k,i,j, 0] = t_check
                            icoCorners[k,i,j,1] = np.where(otherTriang == p0)[0][0] 
                            p1 = otherTriang[(icoCorners[k,i,j,1] + 2)%3]
                            break   
    return icoPts, icoTriangs, icoNbrs, icoCorners

# ==================================================================

def pts_create(triang, nbrs, n):  #the initial distribution of raw grid points (including subdivisions)
    corners = triang
    npts = (n+1)*(n+2)//2
    pts = np.zeros((npts, 3))  # (global) cartesian coordinates
    # - equally spaced points in (u,v) on original triangle:
    k = 0
    for u in range(n+1):
        uh = u/float(n)
        for v in range(u, n+1):
            vh = v/float(n)
            pts[k,:] = corners[0] + uh*(corners[1]-corners[0]) + (vh-uh)*(corners[2]-corners[0])
            k += 1
    # - project onto the unit sphere:
    r = np.sqrt(pts[:,0]**2 + pts[:,1]**2 + pts[:,2]**2)
    pts[:,0] /= r
    pts[:,1] /= r
    pts[:,2] /= r
    return pts

def rotate_on_sphere(pts, axis, angle):
    """
    Rotate the (npts, 3) array pts (in cartesians) around the axis given by a
    position vector (3) on sphere, through given anti-clockwise angle in radians.
    - use Rodrigues' rotation formula
    """
    cosa = np.cos(angle)
    sina = np.sin(angle)
    newpts = pts.copy()*0
    for k, pt in enumerate(pts):
        dot = axis[0]*pt[0] + axis[1]*pt[1] + axis[2]*pt[2]
        newpts[k,0] = pt[0]*cosa+ (axis[1]*pt[2] - axis[2]*pt[1])*sina + axis[0]*dot*(1 - cosa)
        newpts[k,1] = pt[1]*cosa+ (axis[2]*pt[0] - axis[0]*pt[2])*sina + axis[1]*dot*(1 - cosa)
        newpts[k,2] = pt[2]*cosa+ (axis[0]*pt[1] - axis[1]*pt[0])*sina + axis[2]*dot*(1 - cosa)
    return newpts

def rotate_pts(seg_pts, icoPts):   #this provides the point data for the other segments, by a series of complicated rotations
    all_pts = np.zeros((20, len(seg_pts), 3))  #all of the points, in the whole sphere
    all_pts[0] = seg_pts
    all_pts[1] = rotate_on_sphere(seg_pts, [0,0,1], 2*np.pi/5)
    all_pts[2] = rotate_on_sphere(seg_pts, [0,0,1], 4*np.pi/5)
    all_pts[3] = rotate_on_sphere(seg_pts, [0,0,1], 6*np.pi/5)
    all_pts[4] = rotate_on_sphere(seg_pts, [0,0,1], 8*np.pi/5)
    
    all_pts[10] = rotate_on_sphere(all_pts[0], icoPts[1], 8*np.pi/5)
    all_pts[11] = rotate_on_sphere(all_pts[1], icoPts[2], 8*np.pi/5)
    all_pts[12] = rotate_on_sphere(all_pts[2], icoPts[3], 8*np.pi/5)
    all_pts[13] = rotate_on_sphere(all_pts[3], icoPts[4], 8*np.pi/5)
    all_pts[14] = rotate_on_sphere(all_pts[4], icoPts[5], 8*np.pi/5)
    
    all_pts[15] = rotate_on_sphere(all_pts[10], icoPts[1], 8*np.pi/5)
    all_pts[16] = rotate_on_sphere(all_pts[11], icoPts[2], 8*np.pi/5)
    all_pts[17] = rotate_on_sphere(all_pts[12], icoPts[3], 8*np.pi/5)
    all_pts[18] = rotate_on_sphere(all_pts[13], icoPts[4], 8*np.pi/5)
    all_pts[19] = rotate_on_sphere(all_pts[14], icoPts[5], 8*np.pi/5)

    all_pts[5] = rotate_on_sphere(all_pts[17], icoPts[8], 8*np.pi/5)
    all_pts[6] = rotate_on_sphere(all_pts[16], icoPts[9], 8*np.pi/5)
    all_pts[7] = rotate_on_sphere(all_pts[15], icoPts[10], 8*np.pi/5)
    all_pts[8] = rotate_on_sphere(all_pts[19], icoPts[11], 8*np.pi/5)
    all_pts[9] = rotate_on_sphere(all_pts[18], icoPts[7], 8*np.pi/5)
    return all_pts  

def initialise_grid(G, rmax, pts_init = [[0]]):   
    '''
    Sets up the required data using the already-calculated tweaked grids (done up to 7), or from the specified pts_init.
    '''
    n = 2**G
    pts, triangs, nbrs, corners = icosahedron()
    
    #pts_init  = pts_create(pts[triangs[0]], nbrs, n)

    if len(pts_init) == 1:  #not specified, load tweaked points
        pts_init = np.load('tweaking/tweaked_points/tweakpoints%d.npy' % G)   #tweaked points, just the first segment
    #pts_init = pts_untweaked
    
    all_pts = rotate_pts(pts_init, pts)  #rotate the points in such a manner as to cover the whole sphere
    
    
    grids = [Subdomain(pts[triangs], nbrs, corners, n, all_pts[seg], rmax, seg) for seg in range(20)]   #all the grids
    #Bit of a mess doing it here, but grid_whole only does one segment so this can't be paralellised.
    pole_map = np.zeros((12, 2, 5), dtype='int')   
    for i in range(12):
        pole_map[i] = np.array(np.where(triangs == i))
        
    return grids, pole_map, triangs

# ==================================================================
    
class Subdomain:
    """
    A single spherical triangular subdomain of the global mesh.
    - divided into triangular "raw" cells and hexagonal/pentagonal dual cells.
    """
    
    def __init__(self, triang, nbrs, corners, n, pts, rmax, seg):
        """
        Initialise this subdomain.
        """
        def findk(u, v, n):
            return int(((2*n + 3)*u - u**2)//2 + v - u)
        self.corners = corners    
        self.nbrs = nbrs   # 3x2 array giving neighbour ids for each edge and which edge of neighbour they match with.
        self.n = n   # number of cells within subgrid in u and v geodesic coordinates.
        self.G = int(np.log(n)/np.log(2))
        self.rmax = rmax
        # Compute raw grid pts
        # --------------------
        npts = (n+1)*(n+2)//2
        self.npts = int(npts)
        self.nintpts = int((n-1)*(n-2)//2)
        self.pts = pts
        # Specify points in radial direction
        self.r_fact =  max(1,int((np.exp(np.log(rmax)/(2**self.G))-1)//(1.107/(2**self.G))))  
        self.nr = (2**self.G)*self.r_fact
        self.dp = np.log(self.rmax)/self.nr
        self.rs = np.exp((np.linspace(-self.dp,np.log(self.rmax)+self.dp, self.nr+3)))
        self.nr = len(self.rs)-3
        self.rc = np.exp(0.5*(np.log(self.rs[:-1]) + np.log(self.rs[1:])))
        
        # Specify raw grid triangles
        # --------------------------
        ntri = n**2
        self.ntri = ntri
        self.triangs = np.zeros((ntri, 3), dtype='int16')  # indices to self.pts for vertices of each triangle
        
        k = 0
        kt = 0
        for u in range(n+1):
            for v in range(u, n+1):
                # upward pointing triangle:
                if (u < n) & (v < n):
                    self.triangs[kt,:] = np.array([k, k+n-u+1, k+1])
                    kt += 1
                    if (v < n-1):
                        # downward pointing triangle:
                        self.triangs[kt,:] = np.array([k+1, k+n-u+1, k+n-u+2])
                        kt += 1
                k += 1
        # Add global ghost (raw) grid points
        # ---------------------------
        if n%2 != 0:
            raise Exception('Needs to be even for now')
        self.tri_centre = self.circumcentre(pts[n], pts[0], pts[npts-1])
        
        nghost = 3*(n+1)
        self.nghost = nghost
        ghst_pts = np.zeros((nghost, 3))
        rots = np.zeros((3,3))  #rotation points
        # Ghost points along each edge are produced by rotation of appropriate row around the corner point. Can remove this step if not close to the edge (How to tell this...)
        row0 = [n+1 + i for i in range(n)]
        rots[0] = pts[n//2]; rots[1] = pts[findk(n/2, n/2, n)]; rots[2] = pts[findk(n/2, n, n)]  #midpoints of each side, to rotate around
                
        row1 = []
        for u in range(n-1, -1, -1):
            row1.append(((2*n+3)*u - u**2)//2 + 1)
        row2 = []
        for u in range(0, n):
            row2.append(((2*n+3)*u - u**2)//2 + n - u - 1)
            
            
        #sides
        ghst_pts[:n] = self.rotate_on_sphere(pts[row0,:], rots[0], np.pi)
        ghst_pts[n+1:2*n+1] = self.rotate_on_sphere(pts[row1,:], rots[1], np.pi)
        ghst_pts[2*n+2:3*n+2] = self.rotate_on_sphere(pts[row2,:], rots[2], np.pi)
        #corners
        ghst_pts[n,:] = self.rotate_on_sphere(pts[1:2,:], pts[0,:], 4*np.pi/5)
        ghst_pts[2*n+1,:] = self.rotate_on_sphere(pts[npts-2:npts-1,:], pts[npts-1,:], -4*np.pi/5)
        
        ghst_pts[3*n+2,:] = self.rotate_on_sphere(pts[n-1:n,:], pts[n,:], -4*np.pi/5)
        # - append ghost points to end of array of grid points:
        self.pts = np.append(self.pts, ghst_pts, axis=0)
            

        # Add ghost triangles
        # -------------------
        nghosttri = 6*n+3
        self.nghosttri = nghosttri
        ghst_tri = np.zeros((nghosttri, 3), dtype='int16')
        k = 0
        for v in range(n, 0, -1):
            ghst_tri[k,:] = np.array([npts + n-v, v-1, v])
            k += 1
            ghst_tri[k,:] = np.array([npts+n-v+1, v-1, npts+n-v])
            k += 1
        ghst_tri[k-1,:] = np.array([0, npts + n- 1, npts+n])
        ghst_tri[k,:] = np.array([0, npts+n, npts+n+1])
        k += 1
        for u in range(n):
            i = ((2*n + 3)*u - u**2)//2
            i1 = ((2*n + 3)*(u+1) - (u+1)**2)//2
            ghst_tri[k,:] = np.array([i, npts+n+1+u, i1])
            k += 1
            ghst_tri[k,:] = np.array([npts+n+1+u, npts+n+2+u, i1])
            k += 1
        ghst_tri[k-1,:] = np.array([i1, npts+2*n, npts+2*n+1])
        ghst_tri[k,:] = np.array([npts-1, npts+2*n+1, npts+2*n+2])
        k += 1
        for u in range(n,0,-1):
            i = ((2*n + 3)*(u+1) - (u+1)**2)//2 - 1
            i1 = ((2*n + 3)*u - u**2)//2 - 1
            ghst_tri[k,:] = np.array([i1, i, npts + 2*n+2+n-u])
            k += 1
            ghst_tri[k,:] = np.array([i1, npts + 2*n+2+n-u, npts + 2*n+3+n-u])
            k += 1
        ghst_tri[k,:] = np.array([n, npts+3*n+2, npts])
        # - append ghost triangles to end of array of triangles:
        self.triangs = np.append(self.triangs, ghst_tri, axis=0)
        
        # Compute centre of each triangle (i.e. dual grid points)
        # ------------------------------------------------------
        self.centres = np.zeros((ntri + nghosttri, 3))
        for kt, triang in enumerate(self.triangs):
            self.centres[kt,:] = self.circumcentre(self.pts[triang[0]], self.pts[triang[1]], self.pts[triang[2]])

        # Construct the dual grid faces as arrays of hexagons/pentagons
        # -------------------------------------------------------------
        # [index in same order as grid pts]
        self.faces = np.zeros((npts, 6), dtype='int16')
        self.faces[:,:] = -1
        
        # - do faces around interior points first (whose neighbours are all within this grid)
        # - all of these are hexagons
        for u in range(1,n-1):
            for v in range(u+1, n):
                k = ((2*n+3)*u - u**2)//2 + v-u
                f1 = (u-1)*(2*n-u+1) + 2*(v-u) + 1
                f2 = (u-1)*(2*n-u+1) + 2*(v-u)
                f3 = (u-1)*(2*n-u+1) + 2*(v-u) - 1
                f4 = u*(2*n-u) + 2*(v-u) - 2
                f5 = u*(2*n-u) + 2*(v-u) - 1
                f6 = u*(2*n-u) + 2*(v-u)
                self.faces[k,:] = np.array([f3, f4, f5, f6, f1, f2])
        # - now add hexagons around boundary points [except corners]
        for v in range(n-1, 0, -1):
            k = v
            f1 = 2*v - 2
            f2 = 2*v - 1
            f3 = 2*v
            f4 = ntri + 2*(n-1-v)
            f5 = ntri + 2*(n-1-v) + 1
            f6 = ntri + 2*(n-1-v) + 2
            self.faces[k,:] = np.array([f6, f1, f2, f3, f4, f5])
        for u in range(1, n):
            k = ((2*n + 3)*u - u**2)//2
            f1 = u*(2*n-u)
            f2 = (u-1)*(2*n-u+1) + 1
            f3 = (u-1)*(2*n-u+1)
            f4 = ntri + 2*n+1 + 2*(u-1)
            f5 = ntri + 2*n+1 + 2*(u-1) + 1
            f6 = ntri + 2*n+1 + 2*(u-1) + 2
            self.faces[k,:] = np.array([f4, f5, f6, f1, f2, f3])
        for u in range(n-1, 0, -1):
            k = ((2*n+3)*(u+1) - (u+1)**2)//2 - 1
            f1 = (u-1)*(2*n-u+1) + 2*(n-u)
            f2 = (u-1)*(2*n-u+1) + 2*(n-u) - 1
            f3 = u*(2*n-u) + 2*(n-u) - 2
            f4 = ntri + 4*n+2 + 2*(n-1-u)
            f5 = ntri + 4*n+2 + 2*(n-1-u) + 1
            f6 = ntri + 4*n+2 + 2*(n-1-u) + 2
            self.faces[k,:] = np.array([f2, f3, f4, f5, f6, f1])
        # - now add pentagons around corner points [for these, repeat the first point]
        self.faces[0,:] = np.array([0, ntri+2*(n-1), ntri+2*n-1, ntri+2*n, ntri+2*n+1, 0])
        self.faces[npts-1,:] = np.array([ntri-1, ntri+4*n-1, ntri+4*n, ntri+4*n+1, ntri+4*n+2, ntri-1])
        self.faces[n,:] = np.array([2*(n-1), ntri+6*n, ntri+6*n+1, ntri+6*n+2, ntri, 2*(n-1)])
            
        # Record identity of neighbouring (raw) grid pts [not ghost pts]:
        # ---------------------------------------------------------------
        # (order corresponds to ordering of dual cell)
        # - number of neighbours:
        self.nnbrs = np.zeros(npts, dtype='int16')
        self.nnbrs[:] = 6
        self.nnbrs[0], self.nnbrs[npts-1], self.nnbrs[n] = 5, 5, 5
        self.nbrpts = np.zeros((npts, 6), dtype='int16')
        # interior points:
        for u in range(1,n-1):
            for v in range(u+1, n):
                k = ((2*n+3)*u - u**2)//2 + v-u
                self.nbrpts[k,:] = np.array([k-n+u-2, k-1, k+n-u, k+n-u+1, k+1, k-n+u-1])
        # boundary points except corners:
        for v in range(n-1, 0, -1):
            k = v
            self.nbrpts[k,:] = np.array([npts + n-v, k-1, k+n, k+n+1, k+1, npts + n-1-v])
        for u in range(1, n):
            k = ((2*n + 3)*u - u**2)//2
            self.nbrpts[k,:] = np.array([k-n+u-2, npts + n+u, npts + n+u+1, k+n-u+1, k+1, k-n+u-1])
        for u in range(n-1, 0, -1):
            k = ((2*n+3)*(u+1) - (u+1)**2)//2 - 1
            self.nbrpts[k,:] = np.array([k-n+u-2, k-1, k+n-u, npts + 2*n+1 + n-u, npts + 2*n+2+n-u, k-n+u-1])
        # corners points [only 5 neighbours]:
        self.nbrpts[0,:] = np.array([ n+1, 1, npts+n-1, npts+n, npts+n+1, -1])
        #self.nbrpts[npts-1,:] = np.array([npts-3, npts+2*n, npts+2*n+1, npts+2*n+2, npts-2, -1])
        #self.nbrpts[n,:] = np.array([ npts,n-1,  2*n, npts+3*n+1, npts+3*n+2, -1])
        self.nbrpts[npts-1,:] = np.array([npts-2, npts-3, npts+2*n, npts+2*n+1, npts+2*n+2, -1])
        self.nbrpts[n,:] = np.array([n-1,  2*n, npts+3*n+1, npts+3*n+2,npts, -1])
        
        
    def findk(self, u, v, n):
        return int(((2*n + 3)*u - u**2)//2 + v - u)

    def interior_index(self, side, index):   #finds the index of a required interior point. Numbered anticlockwise according to the ghost point location
        if side == 0:
            return int(self.findk(self.n - index, self.n-index + 1, self.n))
        if side == 1:
            return int(self.findk(index - 1, self.n-1, self.n))
        if side == 2:
            return int(self.findk(1, index, self.n))
        
    def edge_index(self, side, index):   #index runs from 1 to n-1
        if side == 0:
            return int(self.findk(self.n - index, self.n-index, self.n))
        if side == 1:
            return int(self.findk(index, self.n, self.n))
        if side == 2:
            return int(self.findk(0, index, self.n))
        
        
    def rotate_on_sphere(self, pts, axis, angle):
        """
        Rotate the (npts, 3) array pts (in cartesians) around the axis given by a
        position vector (3) on sphere, through given anti-clockwise angle in radians.
        - use Rodrigues' rotation formula
        """
        cosa = np.cos(angle)
        sina = np.sin(angle)
        newpts = pts.copy()*0
        for k, pt in enumerate(pts):
            dot = axis[0]*pt[0] + axis[1]*pt[1] + axis[2]*pt[2]
            newpts[k,0] = pt[0]*cosa+ (axis[1]*pt[2] - axis[2]*pt[1])*sina + axis[0]*dot*(1 - cosa)
            newpts[k,1] = pt[1]*cosa+ (axis[2]*pt[0] - axis[0]*pt[2])*sina + axis[1]*dot*(1 - cosa)
            newpts[k,2] = pt[2]*cosa+ (axis[0]*pt[1] - axis[1]*pt[0])*sina + axis[2]*dot*(1 - cosa)
        return newpts
        
    def ghost_transform_point(self, pt0, side): 
        """ 
        Rotates the point pt according to the side in question.
        """
        if side == 0:   #flip in second coordinate and rotate around the pole +2pi/5
            pt = pt0.copy()
            pt[1] = -pt[1] #flipped
            pt = self.rotate_on_sphere(pt, [0,0,1], 2*np.pi/5 )
        if side == 1:   #flip in second coordinate and rotate around the segment tri_centre +2pi/3
            pt = pt0.copy()
            pt[1] = -pt[1] #flipped
            pt = self.rotate_on_sphere(pt, self.tri_centre, 4*np.pi/3 )
        if side == 2:   #Just flip
            pt = pt0.copy()
            pt[1] = -pt[1] 
        return pt
    
    def distance(self, p0, p1):
        """
        Calculate great-circle distance between points p0 and p1 on surface of (unit) sphere.
        - it's just the angle in radians between the position vectors
        """
        
        return np.arccos(min(1.0,p0[0]*p1[0] + p0[1]*p1[1] + p0[2]*p1[2]))
        
    def circumcentre(self, p0, p1, p2):
        """
        Compute circumcentre of triangle with vertices p0, p1, p2. Project on sphere.
        """
        v01, v02 = p1 - p0, p2 - p0
        cx = v01[1]*v02[2] - v01[2]*v02[1]
        cy = v01[2]*v02[0] - v01[0]*v02[2]
        cz = v01[0]*v02[1] - v01[1]*v02[0]
        mag = np.sqrt(cx**2 + cy**2 + cz**2)
        return np.array([cx/mag, cy/mag, cz/mag])
        
    def triangle_area(self, p0, p1, p2):
        """
        Calculate area of spherical triangle with vertices p0, p1, p2 (3) in anti-clockwise order.
        For formula see Eriksson: https://www.jstor.org/stable/2691141?seq=1#metadata_info_tab_contents
        """
        stp = np.abs(p0[0]*(p1[1]*p2[2] - p1[2]*p2[1]) + p0[1]*(p1[2]*p2[0] - p1[0]*p2[2]) + p0[2]*(p1[0]*p2[1] - p1[1]*p2[0]))
        dot01 = p0[0]*p1[0] + p0[1]*p1[1] + p0[2]*p1[2]
        dot02 = p0[0]*p2[0] + p0[1]*p2[1] + p0[2]*p2[2]
        dot12 = p1[0]*p2[0] + p1[1]*p2[1] + p1[2]*p2[2]
        return 2*np.arctan(stp/(1 + dot01 + dot02 + dot12))
    
    def set_u(self, u_exact):
        """
        Set array of u at grid points given function u_exact(s, phi).
        Including ghost points.
        No area weighting.
        """
        self.u = self.pts[:,0]*0
        for k, pt in enumerate(self.pts):  # include ghost pts
            # Convert point to spherical coordinates:
            s = pt[2]
            ph = (np.arctan2(pt[1], pt[0]) + 2*np.pi) % (2*np.pi)
            # Interpolate value of u at this point:
            self.u[k] = u_exact(s, ph)
        
    def laplacian_u(self):
        """
        Compute discrete laplacian of array self.u at grid points.
        """
        lap = np.zeros(self.npts)
        for k in range(self.npts): # omit ghost points
            for j in range(self.nnbrs[k]):
                lap[k] += self.nbr_lens[k,j]*(self.u[self.nbrpts[k,j]] - self.u[k])/self.nbr_dist[k,j]
            lap[k] /= self.dual_areas[k]
            
        return lap

    def midpoint(self, p0, p1): 
        """
        Finds the midpoint of two vectors on the spherical surface
        """
        mid = (p0 + p1)/(np.linalg.norm(p0 + p1, ord = 2))
        return mid  
'''    
    def plot(self, scalars, fname, show=True, dual=False):
        if show:
            p = pv.Plotter()
        else:
            p = pv.Plotter(off_screen=True)
        if not dual:
            faces_local = np.array([[3, t[0], t[1], t[2]] for t in self.triangs[:self.ntri]])
            surf_local = pv.PolyData(self.pts[:self.npts], faces_local)
            p.add_mesh(surf_local, scalars=scalars, show_edges=True)
            p.show(screenshot=fname)
        else:
            dual_faces_local = np.array([[6, t[0], t[1], t[2], t[3], t[4], t[5]] for t in self.faces])
            dual_surf_local = pv.PolyData(self.centres, dual_faces_local)
            scalars = [np.sum((self.lambdas[k,:self.nnbrs[k]]/self.d_lens[k,:self.nnbrs[k]])**4) for k in range(self.npts)]
            p.add_mesh(dual_surf_local, scalars=scalars, show_edges=True)
            p.show(off_screen=True)
'''           
