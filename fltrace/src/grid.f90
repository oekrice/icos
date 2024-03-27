!*******************************************************************************
MODULE grid
!*******************************************************************************
! Generates all the shared grid data, as grid3d does in the python code. Maybe additionally put shared grid in arrays in a different file - not sure yet. Most (if not all) of the grid arrays should not depend on the process, but I guess we'll see about that.
!*******************************************************************************
    USE shared_data
    USE gridpoint_tools

    IMPLICIT NONE

    contains

    subroutine establish_grid()
        !Finds initial point locations, for ALL the segments (not just the ones dealt with by the individual process)
        call icosahedron()

        n = 2**G

        if (use_tweaked_pts) then
            call load_tweaked_pts()   !Loads in the tweaked points from the text file
        else
            call initial_pts()   !Establishes initial (untweaked) points.
        end if
        !The initial points have been found in the main process. These now need to be rotated and sent to the other processes. May as well just keep as one large array here rather than being clever
        call find_all_points()

        call populate_global_grids()

        call extend_to_3d()

        call calculate_grid_extras()

    end subroutine establish_grid

    subroutine calculate_grid_extras()
        implicit none

        real(num), dimension(:):: t0(0:2), n0(0:2), t_side(0:2)
        real(num), dimension(:):: x(0:2)
        integer:: i, k, ri, seg
        ! Calculate the cell normals, not normals, coefficients etc.
        ! This allows for averaging to gridpoints, but isn't necessary otherwise. Is quite slow...

        !Trisides and trinorms are the (unit) cartesian vectors along the triangle sides and triangle normals, respectively.
        !Order coordinate, side, k, segment
        allocate(trisides(0:2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(trinorms(0:2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))

        allocate(alphas(0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(betas(0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(gammas(0:ntri+nghosttri-1,seg_min:seg_max))

        allocate(jbases(0:2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(bbases(0:2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))

        allocate(jcoeffs(0:2,0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(bcoeffs(0:2,0:ntri+nghosttri-1,seg_min:seg_max))
        !do seg = seg_min, seg_max
        do seg = seg_min, seg_max
        !Find triangle side vectors and normals
        do k = 0, ntri+nghosttri - 1
            do i = 0, 2
                t0 = -all_base_pts(:,triangs(i,k),seg) + all_base_pts(:,triangs(mod(i+1,3),k),seg)
                trisides(:,i,k,seg) = t0/sqrt(sum(t0**2))
                t_side = 0.5*(all_base_pts(:,triangs(i,k),seg) + all_base_pts(:,triangs(mod(i+1,3),k),seg))
                n0 = cross1(t0, t_side)
                trinorms(:,i,k,seg) = n0/sqrt(sum(n0**2))
            end do
        end do

        ! Find alpha, gamma coordinates of the centres of each triangle

        do k = 0, ntri+nghosttri-1
        x = coeffs(centres(:,k,seg), all_base_pts(:,triangs(0,k),seg),&
        all_base_pts(:,triangs(1,k),seg),  all_base_pts(:,triangs(2,k),seg))
        x = x/sum(x)
        alphas(k,seg) = x(0)
        betas(k,seg) = x(1)
        gammas(k,seg) = x(2)
        end do

        !Estalish basis functions using the above data
        do k = 0, ntri+nghosttri-1
            jbases(:,0,k,seg) = alphas(k,seg)*trinorms(:,2,k,seg)*tri_l_faces(1,2,k,seg) - &
            betas(k,seg)*trinorms(:,1,k,seg)*tri_l_faces(1,1,k,seg)
            jbases(:,1,k,seg) = betas(k,seg)*trinorms(:,0,k,seg)*tri_l_faces(1,0,k,seg) - &
            gammas(k,seg)*trinorms(:,2,k,seg)*tri_l_faces(1,2,k,seg)
            jbases(:,2,k,seg) = gammas(k,seg)*trinorms(:,1,k,seg)*tri_l_faces(1,1,k,seg) - &
            alphas(k,seg)*trinorms(:,0,k,seg)*tri_l_faces(1,0,k,seg)

            bbases(:,0,k,seg) = alphas(k,seg)*trisides(:,2,k,seg)*tri_l_faces(1,2,k,seg) - &
            betas(k,seg)*trisides(:,1,k,seg)*tri_l_faces(1,1,k,seg)
            bbases(:,1,k,seg) = betas(k,seg)*trisides(:,0,k,seg)*tri_l_faces(1,0,k,seg) - &
            gammas(k,seg)*trisides(:,2,k,seg)*tri_l_faces(1,2,k,seg)
            bbases(:,2,k,seg) = gammas(k,seg)*trisides(:,1,k,seg)*tri_l_faces(1,1,k,seg) - &
            alphas(k,seg)*trisides(:,0,k,seg)*tri_l_faces(1,0,k,seg)
        end do
        !Then jcoeffs and bcoeffs
        do k = 0, ntri+nghosttri-1
        do i = 0, 2
        jcoeffs(i,k,seg) = 1.0_num/(sum(trinorms(:,mod(i+2,3),k,seg)*trisides(:,i,k,seg))*tri_l_faces(1,mod(i+2,3),k,seg))
        bcoeffs(i,k,seg) = 1.0_num/(sum(trisides(:,mod(i+2,3),k,seg)*trinorms(:,i,k,seg))*tri_l_faces(1,mod(i+2,3),k,seg))
        !CHECK THIS IS OK!
        end do
        end do

        end do

        ! Extend some of these matrices into 3D to allow for faster multiplication in the code (removing one do loop)
        ! Essentially the same as using np.newaxis in python
        allocate(trisides_extend(0:2,0:nr+3,0:2,0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(centres_extend(0:2,0:nr+3,0:ntri+nghosttri-1,0:19))

        allocate(jbases_extend(0:2,0:nr+3,0:2,0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(bbases_extend(0:2,0:nr+3,0:2,0:ntri+nghosttri-1,seg_min:seg_max))

        allocate(jcoeffs_extend(0:nr+3,0:2,0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(bcoeffs_extend(0:nr+3,0:2,0:ntri+nghosttri-1,seg_min:seg_max))

        do ri = 0, nr+3
        trisides_extend(0:2,ri,:,:,:) = trisides(:,:,:,:)
        centres_extend(0:2,ri,:,:) = centres(:,:,:)

        jbases_extend(0:2,ri,:,:,:) = jbases(:,:,:,:)
        bbases_extend(0:2,ri,:,:,:) = bbases(:,:,:,:)

        jcoeffs_extend(ri,0:2,:,:) = jcoeffs(:,:,:)
        bcoeffs_extend(ri,0:2,:,:) = bcoeffs(:,:,:)
        end do

        return

    end subroutine calculate_grid_extras

    subroutine load_tweaked_pts()
        implicit none

        character(len=31):: filename
        real(num), allocatable:: load_pts(:)
        integer:: i
        npts = (n+1)*(n+2)/2

        allocate(load_pts(0:3*npts-1))
        allocate(base_pts(0:2, 0:npts-1))
        write (filename, "(A26,I1,A4)") "tweaked_points/tweakpoints", G, ".txt"
        open(1,file= filename )
        read(1, *) load_pts
        close(1)

        do i = 0, npts-1
            base_pts(:,i) = load_pts(3*i:3*i+2)
        end do

        deallocate(load_pts)
        return
    end subroutine load_tweaked_pts

    subroutine extend_to_3d()
        implicit none
        integer:: ri
        allocate(all_pts(0:nr +2,0:2,0:npts+nghosts-1,0:19))
        allocate(h_areas(0:nr +2,0:ntri+nghosttri-1,0:19))
        allocate(tri_l_faces(0:nr +2,0:2,0:ntri+nghosttri-1,0:19))

        allocate(all_centres(0:nr +1,0:2,0:ntri+nghosttri-1,0:19))
        allocate(all_ds(0:nr +1,0:2,0:ntri-1,0:19))
        allocate(all_point_ds(0:nr +1,0:5,0:npts-1,0:19))
        allocate(tri_d_all(0:nr +1,0:2,0:ntri-1,0:19))
        allocate(all_duals(0:nr +1,0:npts-1,0:19))
        allocate(tri_l_mids(0:nr +1,0:2,0:ntri+nghosttri-1,0:19))

        allocate(v_areas(0:nr +1,0:2,0:ntri+nghosttri-1,0:19))
        allocate(cell_volumes(0:nr +1,0:ntri+nghosttri-1,0:19))

        allocate(v_duals(1:nr +1,0:2,0:ntri-1,0:19))

        do ri = 0, nr +2  !The arrays saved on radial grid points. Only the local section is required
            all_pts(ri,:,:,:) = rs_global(ri + seg_layer*nr)*all_base_pts(:,:,:)
            h_areas(ri,:,:) = rs_global(ri + seg_layer*nr)**2*areas(:,:)
            tri_l_faces(ri,:,:,:) = rs_global(ri + seg_layer*nr)*tri_l(:,:,:)

        end do

        do ri = 0, nr +1  !The arrays saved on radial centres
            all_centres(ri,:,:,:) = rc_global(ri + seg_layer*nr)*centres(:,:,:)
            all_ds(ri,:,:,:) = rc_global(ri + seg_layer*nr)*tri_d(:,:,:)
            all_point_ds(ri,:,:,:) = rc_global(ri + seg_layer*nr)*d_lens(:,:,:)
            tri_d_all(ri,:,:,:) = rc_global(ri + seg_layer*nr)*tri_d(:,:,:)
            all_duals(ri,:,:) = rc_global(ri + seg_layer*nr)**2*dual_areas(:,:)
            tri_l_mids(ri,:,:,:) = rc_global(ri + seg_layer*nr)*tri_l(:,:,:)
        end do
        !Arrays based on the difference between grid points (vertical areas and volumes)
        do ri = 0, nr +1  !The arrays saved on radial centres
            v_areas(ri,:,:,:) = 0.5_d*(rs_global(ri + seg_layer*nr+1)**2 - rs_global(ri + seg_layer*nr)**2)*tri_l(:,:,:)
            cell_volumes(ri,:,:) = (1.0_d/3.0_d)*(rs_global(ri + seg_layer*nr+1)**3 - rs_global(ri + seg_layer*nr)**3)*areas(:,:)
        end do

        do ri = 1, nr +1
            !Don't have to do the exception to the rule!! Careful to pass this correctly though, numbering starts at 1
            v_duals(ri,:,:,:) = 0.5_d*(rc_global(ri + seg_layer*nr)**2 - rc_global(ri + seg_layer*nr - 1)**2)*tri_d(:,:,:)

        end do


        return
    end subroutine

    subroutine icosahedron()
        ! Generate icosahedron inscribed in unit sphere.
        ! Returns array of points and array of segments (with points ordered anti-clockwise).
        ! Also returns neighbours of each segment and their corresponding orientations.
        implicit none
        real(num):: s, c, orx, ory, orz
        real(num), dimension(0:2, 0:5):: top_pts, bottom_pts
        real(num), dimension(0:2):: p0, p1, v10, v12
        integer:: i, j, k, i1, i2, k1
        integer, dimension(0:2):: tri_temp, tri_check
        integer:: t0, t1
        integer:: loc1, loc2
        logical:: flag1

        s = (2.0_d/sqrt(5.0_d)); c = (1.0_d/sqrt(5.0_d))
        top_pts(:,0) = [0.0_d,0.0_d,1.0_d]; bottom_pts(:,0) = [0.0_d,0.0_d,-1.0_d]

        do i = 0, 4
            top_pts(:,i+1) = [s*cos(i*2.0_d*pi/5.0_d - pi/5.0_d), s*sin(i*2.0_d*pi/5.0_d - pi/5.0_d), c]
            bottom_pts(:,i+1) = [s*cos(-i*2.0_d*pi/5.0_d + 4.0_d*pi/5.0_d), s*sin(-i*2.0_d*pi/5.0_d + 4.0_d*pi/5.0_d), -c]
        end do

        ico_pts(:,0:5) = top_pts(:,:)
        ico_pts(:,6:11) = bottom_pts(:,:)

        do i = 0, 4   !Determine the points at the vertices of each triangle
            ico_triangs(i,:) =    [0, i+1, mod(i+1,5)+1]
            ico_triangs(i+5,:) =  [6, i+7, mod(i+1,5)+7]
            ico_triangs(i+10,:) = [i+1, mod(i+1,5)+1, mod(7-i,5) + 7]
            ico_triangs(i+15,:) = [i+1, mod(7-i,5)+7, mod(8-i,5) + 7]
        end do

        do i = 0, 19  !Ensure that the triangles are numbered in the correct direction (anti-clockwise), by means of a cross product
            p0 = ico_pts(:,ico_triangs(i,0))
            p1 = ico_pts(:,ico_triangs(i,1))
            v10 = ico_pts(:,ico_triangs(i,0)) - ico_pts(:,ico_triangs(i,1))
            v12 = ico_pts(:,ico_triangs(i,2)) - ico_pts(:,ico_triangs(i,1))

            orx = p1(0)*(v12(1)*v10(2) - v12(2)*v10(1))
            ory = p1(1)*(v12(2)*v10(0) - v12(0)*v10(2))
            orz = p1(2)*(v12(0)*v10(1) - v12(1)*v10(0))

            if (orx + ory + orz < 0) then  !Swap triangle sides
                tri_temp = ico_triangs(i,:)
                ico_triangs(i,0) = tri_temp(1)
                ico_triangs(i,1) = tri_temp(0)
            end if
        end do
        !Iterate through segments to check location of segment neighbours. Segments share an edge if they have two vertices in common
        do k = 0, 19
            do i  = 0, 2   !Check side of each segment
                t0 = ico_triangs(k, i)
                t1 = ico_triangs(k, mod(i+1, 3))  !Indices of the points on each side of the segment
                do j = 0, 19
                    tri_check = ico_triangs(j,:)
                    loc1 = -1; loc2 = -1  !Flags to replace np.isin
                    if (j .ne. k) then
                        do i1 = 0, 2
                            if (t0 == tri_check(i1)) loc1 = i1
                            if (t1 == tri_check(i1)) loc2 = i1
                        end do
                        if (loc1 > -1 .and. loc2 > -1) then   !Does indeed share a side
                            ico_nbrs(k, i, 0) = j
                            ico_nbrs(k, i, 1) = loc2
                        end if
                    end if
                end do
            end do
        end do

        !Determine the segments which share corner points. Each segment has three corners, which are shared with four other segments, so the dimension here is (20,3,4,2), where the last index gives the 'position' of the corner in the other segment
        ico_corners = -1   !Set to -1 so we know that they've been populated
        do k = 0, 19
            do i = 0, 2
                t0 = ico_triangs(k, i)
                t1 = ico_triangs(k, mod(i+2, 3))  !Indices of the points on each side of the segment
                do j = 0, 4  !Find the four other segments sharing this icosahedron point. Anticlockwise
                    do k1 = 0, 19  !Checking other segments
                        tri_check = ico_triangs(k1, :)
                        flag1 = .false.
                        do i2 = 0, 4   !check not already in list
                            if (k1 == ico_corners(k, i, i2, 0)) flag1 = .true.
                        end do
                        if (.not. flag1 .and. k1 .ne. k) then  !This is a valid point sharer.
                            loc1 = -1; loc2 = -1  !Flags to replace np.isin
                            do i1 = 0, 2
                                if (t0 == tri_check(i1)) loc1 = i1
                                if (t1 == tri_check(i1)) loc2 = i1
                            end do
                            if (loc1 > -1 .and. loc2 > -1) then   !Does indeed share a side
                                ico_corners(k, i, j, 0) = k1
                                ico_corners(k, i, j, 1) = loc1
                                !New t1 to check (this makes sense now!)
                                t1 = tri_check(mod(ico_corners(k, i, j, 1) + 2, 3))
                                exit
                            end if
                        end if
                    end do
                end do
            end do
        end do

        pole_map = -1
        do i = 0, 11
            j = 0  !Count index for each pole (should go to 5)
            do k = 0, 19  !Check segments
                do i1 = 0, 2
                    if (ico_triangs(k, i1) == i) then
                        pole_map(i, 0, j) = k
                        pole_map(i, 1, j) = i1
                        j = j + 1
                    end if
                end do
            end do
        end do
        return
    end subroutine icosahedron

    subroutine initial_pts()
    ! Given the icosahedral corners of a segment, produces the `untweaked' point grid, for a single segment (zero). Other segments are to be ascertained b rotations
        implicit none
        INTEGER:: k, i
        INTEGER:: u, v
        real(num):: uh, vh, r
        real(num), dimension(0:2, 0:2):: corners

        npts = (n+1)*(n+2)/2
        allocate(base_pts(0:2, 0:npts-1))

        do i = 0, 2
            corners(:,i) = ico_pts(:,ico_triangs(0,i))  !Coordinates of corner points
        end do

        k = 0
        do u = 0, n
            uh = dble(u)/dble(n)
            do v = u, n
                vh = dble(v)/dble(n)
                base_pts(:,k) = corners(:,0) + uh*(corners(:,1) - corners(:,0)) + (vh - uh)*(corners(:,2) - corners(:,0))
                k = k + 1
            end do
        end do
        do k = 0, npts-1
            r = sqrt(base_pts(0,k)**2 + base_pts(1,k)**2 + base_pts(2,k)**2)
            base_pts(:,k) = base_pts(:,k)/r
        end do
        return

    end subroutine initial_pts

    subroutine populate_global_grids()
        ! Do most of what goes on in the __init__ section of the subdomain class
        implicit none
        real(num), dimension(0:n+2):: ps
        real(num), dimension(0:n+1):: pc
        real(num), dimension(0:2,0:2, 0:19):: rots !Ghost point rotation axes
        real(num), dimension(0:2,0:n-1,0:19):: row0_pts, row1_pts, row2_pts !Actual row point locations
        integer, dimension(0:2):: triang
        integer:: i, i1, j, k, kt, tk, u, v, g1
        integer:: f1,f2,f3,f4,f5,f6
        integer:: n0, n1, n2, n3, n4, n5
        integer:: ghostk, seg
        real(num), dimension(0:2):: pt0, pt1
        integer, dimension(0:2):: corner_faces
        integer, dimension(0:npts-1):: tri_c  !Counter for the trimap

        nr_global = n
        nr = nr_global/1
        dp = log(rmax)/nr_global

        allocate(rs_global(0:nr_global +2))
        allocate(rc_global(0:nr_global +1))

        allocate(rs(0:nr +2))
        allocate(rc(0:nr +1))

        do j = 0, nr_global +2
            ps(j) = -dp + j*(log(rmax) + 2*dp)/(nr_global +2)
        end do
        pc(:) = 0.5*(ps(0:nr_global +1)+ps(1:nr_global +2))
        rs_global = exp(ps)
        rc_global = exp(pc)

        rs = rs_global(seg_layer*nr:seg_layer*nr+nr+2)
        rc = rc_global(seg_layer*nr:seg_layer*nr+nr+1)

        ntri = n**2   !number of non-ghost triangles
        nghosttri = 6*n+3

        allocate(triangs(0:2,0:ntri+nghosttri-1))
        triangs = -1
        k = 0; kt = 0
        do u = 0, n
            do v = u, n
                if (u < n .and. v < n) then
                    !Upward pointing triangle
                    triangs(:,kt) = [k, k +n-u+1, k+1]
                    kt = kt + 1
                    if (v < n-1) then
                    !Also downward-pointing triangle
                        triangs(:,kt) = [k+1, k +n-u+1, k+n-u+2]
                        kt = kt + 1
                    end if
                end if
                k = k + 1
            end do
        end do

        if (mod(n, 2) .ne. 0) call mpi_abort(comm, ierr)  !If odd number of side points (not usually a problem...) then these ghost point transformations will not work

        !Add ghost (raw) grid points, obtained through rotations (rather than from other segments, which would perhaps be neater...)
        do seg = 0 ,19
            tri_centre(:,seg) = circumcentre(all_base_pts(:,n,seg), all_base_pts(:,0,seg), all_base_pts(:,npts-1,seg))
        end do
        nghosts = 3*(n+1)

        do seg = 0, 19
            do u = 0, n-1
                row0_pts(:,u,seg) = all_base_pts(:,n + 1 + u,seg)
                row2_pts(:,u,seg) = all_base_pts(:,(int(((2*n+3)*u - u**2)/2) + n - u - 1),seg)
            end do
            i = 0
            do u = n-1, 0, -1
                row1_pts(:,i,seg) = all_base_pts(:,(int(((2*n+3)*u - u**2)/2) + 1),seg)
                i = i + 1
            end do
        end do

        !Add ghost points to the base_pts array
        do seg = 0, 19
            rots(:,0,seg) = all_base_pts(:,n/2, seg)   !Ghost point rotation axes
            rots(:,1,seg) = all_base_pts(:,findk(n/2,n/2,n), seg)
            rots(:,2,seg) = all_base_pts(:,findk(n/2,n,n), seg)
            !Sides
            all_base_pts(:,npts:npts + n - 1, seg) = rotate_pts(row0_pts(:,:,seg), rots(:,0,seg), pi, n)
            all_base_pts(:,npts + n + 1:npts + 2*n , seg) = rotate_pts(row1_pts(:,:,seg), rots(:,1,seg), pi, n)
            all_base_pts(:,npts + 2*n + 2:npts + 3*n+1, seg) = rotate_pts(row2_pts(:,:,seg), rots(:,2,seg), pi, n)
            !Corners (which have to be sourced elsewhere)
            all_base_pts(:,npts+n:npts+n, seg) = &
            rotate_pts(all_base_pts(:,1, seg), all_base_pts(:,0, seg), 4.0_d*pi/5.0_d, 1)
            all_base_pts(:,npts+2*n+1:npts+2*n+1, seg) = &
            rotate_pts(all_base_pts(:,npts-2, seg), all_base_pts(:,npts-1, seg), -4.0_d*pi/5.0_d, 1)
            all_base_pts(:,npts+3*n+2:npts+3*n+2, seg) = &
            rotate_pts(all_base_pts(:,n-1, seg), all_base_pts(:,n, seg), -4.0_d*pi/5.0_d, 1)

        end do

        !Add ghost triangles. These are the indices of the points in each triangle, rather than the coordinates so don't need to split by segment.
        k = 0
        do v = n, 1, -1
            triangs(:, k + ntri) = [npts + n - v, v-1, v]
            k = k + 1
            triangs(:, k + ntri) = [npts + n - v + 1, v-1, npts + n - v]
            k = k + 1
        end do
        triangs(:, ntri + k - 1) = [0, npts + n -1, npts + n]
        triangs(:, ntri + k) = [0, npts + n, npts + n + 1]
        k = k + 1
        do u = 0, n-1
            i = ((2*n + 3)*u - u**2)/2
            i1 = ((2*n + 3)*(u+1) - (u+1)**2)/2
            triangs(:, k + ntri) = [i, npts+n+1+u, i1]
            k = k + 1
            triangs(:, k + ntri) = [npts+n+1+u, npts+n+2+u, i1]
            k = k + 1
        end do
        i1 = ((2*n + 3)*(n-1+1) - (n-1+1)**2)/2
        triangs(:, ntri + k - 1) = [i1, npts+2*n, npts+2*n+1]
        triangs(:, ntri + k) = [npts-1, npts+2*n+1, npts+2*n+2]
        k = k + 1
        do u = n, 1, -1
            i = ((2*n + 3)*(u+1) - (u+1)**2)/2 - 1
            i1 = ((2*n + 3)*u - u**2)/2 - 1
            triangs(:, k + ntri) = [i1, i, npts + 2*n+2+n-u]
            k = k + 1
            triangs(:, k + ntri) = [i1, npts + 2*n+2+n-u, npts + 2*n+3+n-u]
            k = k + 1
        end do
        triangs(:, k + ntri) = [n, npts+3*n+2, npts]

        !Calculate centre of each triangles (dual grid points)
        allocate(centres(0:2, 0:ntri+nghosttri-1, 0:19))


        do seg = 0, 19
            do kt = 0, ntri + nghosttri - 1
                triang(:) = triangs(:,kt)
                centres(:,kt,seg) = &
                circumcentre(all_base_pts(:,triang(0),seg),all_base_pts(:,triang(1),seg), all_base_pts(:,triang(2),seg))
            end do

        end do

        tri_c = 0
        allocate(trimap(0:1,0:5,0:npts-1))
        trimap = -1
        !Construct inverse map of points to triangles. For the interpolation code (for now)
        do k = 0, ntri+nghosttri-1  !Can't do ghost points but need ghost tri
        do i = 0, 2
        if (triangs(i,k) < npts) then
            trimap(0,tri_c(triangs(i,k)),triangs(i,k)) = k
            trimap(1,tri_c(triangs(i,k)),triangs(i,k)) = i
            tri_c(triangs(i,k)) = tri_c(triangs(i,k)) + 1
        end if
        end do
        end do

        !Construct dual grid faces as arrays of hexagons/pentagons
        !Indexed in the same order as grid points
        allocate(faces(0:5, 0:npts-1))  !Not ghost points here, for reasons I used to understand

        faces = -1.0_d
        !All the interior hexagons. Ordering definitely does matter here
        do u = 1, n-2
            do v = u+1, n-1
                k = ((2*n+3)*u - u**2)/2 + v-u
                f1 = (u-1)*(2*n-u+1) + 2*(v-u) + 1
                f2 = (u-1)*(2*n-u+1) + 2*(v-u)
                f3 = (u-1)*(2*n-u+1) + 2*(v-u) - 1
                f4 = u*(2*n-u) + 2*(v-u) - 2
                f5 = u*(2*n-u) + 2*(v-u) - 1
                f6 = u*(2*n-u) + 2*(v-u)
                faces(:,k) = [f3, f4, f5, f6, f1, f2]
            end do
        end do
        !Hexagons around the boundary points (except the corners)
        do v = n-1, 1, -1
            k = v
            f1 = 2*v - 2
            f2 = 2*v - 1
            f3 = 2*v
            f4 = ntri + 2*(n-1-v)
            f5 = ntri + 2*(n-1-v) + 1
            f6 = ntri + 2*(n-1-v) + 2
            faces(:,k) = [f6, f1, f2, f3, f4, f5]
        end do
        do u = 1, n-1
            k = ((2*n + 3)*u - u**2)/2
            f1 = u*(2*n-u)
            f2 = (u-1)*(2*n-u+1) + 1
            f3 = (u-1)*(2*n-u+1)
            f4 = ntri + 2*n+1 + 2*(u-1)
            f5 = ntri + 2*n+1 + 2*(u-1) + 1
            f6 = ntri + 2*n+1 + 2*(u-1) + 2
            faces(:,k) = [f4, f5, f6, f1, f2, f3]
        end do
        do u = n-1, 1, -1
            k = ((2*n+3)*(u+1) - (u+1)**2)/2 - 1
            f1 = (u-1)*(2*n-u+1) + 2*(n-u)
            f2 = (u-1)*(2*n-u+1) + 2*(n-u) - 1
            f3 = u*(2*n-u) + 2*(n-u) - 2
            f4 = ntri + 4*n+2 + 2*(n-1-u)
            f5 = ntri + 4*n+2 + 2*(n-1-u) + 1
            f6 = ntri + 4*n+2 + 2*(n-1-u) + 2
            faces(:,k) = [f2, f3, f4, f5, f6, f1]
        end do
        !Add pentagons around the corner points, repeating the first point in the sixth index
        faces(:,0) = [0, ntri+2*(n-1), ntri+2*n-1, ntri+2*n, ntri+2*n+1, 0]
        faces(:,npts-1) = [ntri-1, ntri+4*n-1, ntri+4*n, ntri+4*n+1, ntri+4*n+2, ntri-1]
        faces(:,n) = [2*(n-1), ntri+6*n, ntri+6*n+1, ntri+6*n+2, ntri, 2*(n-1)]

        !Record identity of neighbouring raw grid points (not ghost points)
        !The ordering of these corresponds to the ordering of the respective dual cell
        allocate(nnbrs(0:npts-1))
        nnbrs = 6
        nnbrs(0) = 5; nnbrs(npts-1) = 5; nnbrs(n) = 5
        allocate(nbrpts(0:5, 0:npts-1))

        nbrpts = -1
        do u = 1, n-2
            do v = u+1, n-1
                k = ((2*n+3)*u - u**2)/2 + v-u
                nbrpts(:,k) = [k-n+u-2, k-1, k+n-u, k+n-u+1, k+1, k-n+u-1]
            end do
        end do
        !boundary points except corners:
        do v = n-1, 1, -1
            k = v
            nbrpts(:,k) = [npts + n-v, k-1, k+n, k+n+1, k+1, npts + n-1-v]
        end do
        do u = 1, n-1
            k = ((2*n + 3)*u - u**2)/2
            nbrpts(:,k) = [k-n+u-2, npts + n+u, npts + n+u+1, k+n-u+1, k+1, k-n+u-1]
        end do
        do u = n-1, 1, -1
            k = ((2*n+3)*(u+1) - (u+1)**2)/2 - 1
            nbrpts(:,k) = [k-n+u-2, k-1, k+n-u, npts + 2*n+1 + n-u, npts + 2*n+2+n-u, k-n+u-1]
        end do

        nbrpts(:,0) = [ n+1, 1, npts+n-1, npts+n, npts+n+1, -1]

        nbrpts(:,npts-1) = [npts-2, npts-3, npts+2*n, npts+2*n+1, npts+2*n+2, -1]
        nbrpts(:,n) = [n-1,  2*n, npts+3*n+1, npts+3*n+2,npts, -1]


        !Record identity of neighbouring grid FACES (not ghosts)
        !These are ordered anticlockwise, and do have a coherent ordering system beyond that
        allocate(updown(0:ntri-1))
        allocate(nbrfaces(0:2, 0:ntri-1))
        !Populate neighbours that are not ghosts. Ghost points are set as -1 for now.
        do u = 0, n-1
            do v = u + 1, n
                tk = findt(u, v, n, 1)  !Upward triangle
                n1 = -1; n3 =  -1; n5 =  -1
                if (v > u + 1)  n1 = findt(u, v-1, n, 0)
                if (v < n) n3 = findt(u, v, n, 0)
                if (u > 0)  n5 = findt(u-1, v-1, n, 0)

                nbrfaces(:, tk) = [n1,n3,n5]
                updown(tk) = 1

                if (v < n) then   !there is also a downward triangle. Always on the segment interior
                    tk = findt(u, v, n, 0)   !downward triangle
                    n0 = -1; n2 =  -1; n4 =  -1
                    n0 = findt(u, v, n, 1)
                    n2 = findt(u+1, v+1, n, 1)
                    n4 = findt(u, v+1, n, 1)
                    nbrfaces(:, tk) = [n0,n2,n4]
                    updown(tk) = 0
                end if
            end do
        end do

        do i = 0, n-1
            ghostk = n**2 + 2*i
            nbrfaces(2, findt(0, n-i, n, 1)) = ghostk
            ghostk = n**2 + 2*n + 1 + 2*i
            nbrfaces(0, findt(i, i+1, n, 1)) = ghostk
            ghostk = n**2 + 4*n + 2 + 2*i
            nbrfaces(1, findt(n-i-1, n, n, 1)) = ghostk
        end do

    allocate(nbr_dist(0:5,0:npts-1,0:19))
    allocate(nbr_lens(0:5,0:npts-1,0:19))
    allocate(areas(0:ntri+nghosttri-1,0:19))
    allocate(dual_areas(0:npts-1,0:19))

    do seg = 0, 19

        !Calculate distances from neighbouring raw grid points
        do k = 0, npts-1
            pt0 = all_base_pts(:,k,seg)
            do j = 0, nnbrs(k)-1
                pt1 = all_base_pts(:,nbrpts(j,k),seg)
                nbr_dist(j,k,seg) = distance(pt0, pt1)
            end do
        end do

        !Calculate lengths of dual edges around each (raw) grid point
        do k = 0, npts-1
            do j = 0, nnbrs(k)-1
                pt0 = centres(:,faces(j,k),seg)
                pt1 = centres(:,faces(mod(j+1,6),k),seg)
                nbr_lens(j,k,seg) = distance(pt0, pt1)
            end do
        end do

        !Calculate raw grid cell areas
        do k = 0, npts-1
            do j = 0, nnbrs(k)-1
                pt0 = centres(:,faces(j,k),seg)
                pt1 = centres(:,faces(mod(j+1,6),k),seg)
                nbr_lens(j,k,seg) = distance(pt0, pt1)
            end do
        end do
        ! Compute areas of raw grid cells (triangles)
        do kt = 0, ntri + nghosttri - 1
            triang = triangs(:,kt)
            areas(kt, seg) = &
            triangle_area(all_base_pts(:,triang(0),seg),all_base_pts(:,triang(1),seg),all_base_pts(:,triang(2),seg))
        end do
        dual_areas(:,seg) = 0.0_d
        ! Compute areas of dual grid cells (hexagons or pentagons)
        do k = 0, npts-1
            do j = 0, nnbrs(k) - 1
                pt0 = centres(:,faces(j,k),seg)
                pt1 = centres(:,faces(mod(j+1,6),k),seg)
                dual_areas(k, seg) = dual_areas(k, seg)  + triangle_area(all_base_pts(:,k,seg),pt0,pt1)
            end do

        end do
    end do

    !Weights for dual grids in global integration
    !Not sure this is entirely necessary in the main code, but it is worthwhile checking it works
    allocate(dual_weights(0:npts-1))
    dual_weights = 1.0_d
    !Edge points
    do v = n-1, 1, -1
        dual_weights(v) = 0.5_d
    end do
    do u = 1, n-1
        dual_weights(((2*n + 3)*u - u**2)/2) = 0.5_d
    end do
    do u = n-1, 1, -1
        dual_weights(((2*n+3)*(u+1) - (u+1)**2)/2 - 1) = 0.5_d
    end do
    !Corner points
    dual_weights(0) = 0.2_d
    dual_weights(npts-1) = 0.2_d
    dual_weights(n) = 0.2_d

    if (.false.) then  !Print area errors. Nice check at this point I think
        do seg = 0, seg-1
            print*, 'segment', seg, 'area error' , sum(dual_weights(:)*dual_areas(:,seg)) - pi/5
        end do
    end if

    allocate(d_lens(0:5,0:npts-1,0:19))
    allocate(tri_d(0:2,0:ntri-1,0:19))
    allocate(tri_l(0:2,0:ntri+nghosttri-1,0:19))

    d_lens = 0.0
    do seg = 0, seg-1
        do k = 0, npts-1
            do j = 0, nnbrs(k) -1
                pt0 = centres(:, faces(mod(j,nnbrs(k)), k), seg)
                pt1 = centres(:, faces(mod(j+1,nnbrs(k)), k), seg)
                d_lens(j,k,seg) = distance(pt0, pt1)
            end do
        end do

        do kt  = 0, ntri - 1
            do j = 0, 2
                tri_d(j,kt,seg) = distance(centres(:,kt,seg), centres(:,nbrfaces(j, kt),seg))
            end do
        end do
        do kt  = 0, ntri + nghosttri - 1
            do j = 0, 2
                tri_l(j,kt,seg) = distance(all_base_pts(:,triangs(j,kt),seg), all_base_pts(:,triangs(mod(j+1,3),kt),seg))
            end do
        end do

    end do

    allocate(edge_transfer(0:2,0:(n-1)*3-1,0:19))
    allocate(tri_transfer(0:2,0:nghosttri-1,0:19))

    do seg  = 0, 19
        do i = 1, n-1
            edge_transfer(0, i-1, seg) = findk(i,i,n)
            edge_transfer(1, i-1, seg) = ico_nbrs(seg, 0, 0)
            edge_transfer(2, i-1, seg) = edge_index(ico_nbrs(seg, 0, 1), i, n)

            edge_transfer(0, i+n-2, seg) = findk(n-i, n, n)
            edge_transfer(1, i+n-2, seg) = ico_nbrs(seg, 1, 0)
            edge_transfer(2, i+n-2, seg) = edge_index(ico_nbrs(seg, 1, 1), i, n)

            edge_transfer(0, i + 2*n-3, seg) = findk(0, n-i, n)
            edge_transfer(1, i + 2*n-3, seg) = ico_nbrs(seg, 2, 0)
            edge_transfer(2, i + 2*n-3, seg) = edge_index(ico_nbrs(seg, 2, 1), i, n)
        end do
        !call printmat(dble(edge_transfer(:,:,seg)), 3, (n-1)*3)
    end do

    corner_faces = [0, n**2 - 1, 2*(n-1)]
    do seg  = 0, 19

        i = 0
        do g1 = n**2 + 2*n - 1, n**2 + 2*n
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_corners(seg,0,i+1,0), corner_faces(ico_corners(seg,0,i+1,1)), ico_corners(seg,0,i+1, 1)]
            i = i + 1
        end do
        i = 0
        do g1 = n**2 + 4*n, n**2 + 4*n + 1
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_corners(seg,1,i+1,0), corner_faces(ico_corners(seg,1,i+1,1)), ico_corners(seg,1,i+1, 1)]
            i = i + 1
        end do
        i = 0
        do g1 = n**2 + 6*n + 1, n**2 + 6*n + 2
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_corners(seg,2,i+1,0), corner_faces(ico_corners(seg,2,i+1,1)), ico_corners(seg,2,i+1, 1)]
            i = i + 1
        end do
        i = 0
        do g1 = n**2 + 2*n + 1, n**2 + 4*n -1  !side 0
            if (ico_nbrs(seg,0,1) == 0) then
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_nbrs(seg,0,0), side_tri(ico_nbrs(seg,0,1),i,n), -1]
            end if
            if (ico_nbrs(seg,0,1) == 1) then
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_nbrs(seg,0,0), side_tri(ico_nbrs(seg,0,1),i,n), 2*(mod(g1,2))]
            end if
            if (ico_nbrs(seg,0,1) == 2) then
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_nbrs(seg,0,0), side_tri(ico_nbrs(seg,0,1),i,n), mod(g1-1,2)]
            end if
            i = i + 1
        end do

        i = 0
        do g1 = n**2 + 4*n + 2, n**2 + 6*n   !side 1
            if (ico_nbrs(seg,1,1) == 0) then
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_nbrs(seg,1,0), side_tri(ico_nbrs(seg,1,1),i,n), mod(g1,2)]
            end if
            if (ico_nbrs(seg,1,1) == 1) then
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_nbrs(seg,1,0), side_tri(ico_nbrs(seg,1,1),i,n), -1]
            end if
            if (ico_nbrs(seg,1,1) == 2) then
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_nbrs(seg,1,0), side_tri(ico_nbrs(seg,1,1),i,n), 2*(mod(g1-1,2))]
            end if
            i = i + 1
        end do

        i = 0
        do g1 = n**2, n**2 + 2*n - 2   !side 2
            if (ico_nbrs(seg,2,1) == 0) then
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_nbrs(seg,2,0), side_tri(ico_nbrs(seg,2,1),i,n), 2*(mod(g1-1,2))]
            end if
            if (ico_nbrs(seg,2,1) == 1) then
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_nbrs(seg,2,0), side_tri(ico_nbrs(seg,2,1),i,n), mod(g1,2)]
            end if
            if (ico_nbrs(seg,2,1) == 2) then
            tri_transfer(:,g1 - ntri,seg) = &
            [ico_nbrs(seg,2,0), side_tri(ico_nbrs(seg,2,1),i,n), -1]
            end if
            i = i + 1
        end do

    end do

    return

    end subroutine populate_global_grids

    subroutine find_all_points()
        implicit none

        !Do the series of rotations
        nghosts = 3*(n+1)

        allocate(all_base_pts(0:2,0:npts + nghosts-1,0:19))

        all_base_pts(:,0:npts-1,0) = base_pts(:,:)
        all_base_pts(:,0:npts-1,1) = rotate_pts(all_base_pts(:,0:npts-1,0), dble([0.,0.,1.]), 2.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,2) = rotate_pts(all_base_pts(:,0:npts-1,0), dble([0.,0.,1.]), 4.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,3) = rotate_pts(all_base_pts(:,0:npts-1,0), dble([0.,0.,1.]), 6.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,4) = rotate_pts(all_base_pts(:,0:npts-1,0), dble([0.,0.,1.]), 8.0_d*pi/5.0_d, npts)

        all_base_pts(:,0:npts-1,10) = rotate_pts(all_base_pts(:,0:npts-1,0), ico_pts(:,1), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,11) = rotate_pts(all_base_pts(:,0:npts-1,1), ico_pts(:,2), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,12) = rotate_pts(all_base_pts(:,0:npts-1,2), ico_pts(:,3), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,13) = rotate_pts(all_base_pts(:,0:npts-1,3), ico_pts(:,4), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,14) = rotate_pts(all_base_pts(:,0:npts-1,4), ico_pts(:,5), 8.0_d*pi/5.0_d, npts)

        all_base_pts(:,0:npts-1,15) = rotate_pts(all_base_pts(:,0:npts-1,10), ico_pts(:,1), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,16) = rotate_pts(all_base_pts(:,0:npts-1,11), ico_pts(:,2), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,17) = rotate_pts(all_base_pts(:,0:npts-1,12), ico_pts(:,3), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,18) = rotate_pts(all_base_pts(:,0:npts-1,13), ico_pts(:,4), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,19) = rotate_pts(all_base_pts(:,0:npts-1,14), ico_pts(:,5), 8.0_d*pi/5.0_d, npts)

        all_base_pts(:,0:npts-1,5) = rotate_pts(all_base_pts(:,0:npts-1,17), ico_pts(:,8), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,6) = rotate_pts(all_base_pts(:,0:npts-1,16), ico_pts(:,9), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,7) = rotate_pts(all_base_pts(:,0:npts-1,15), ico_pts(:,10), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,8) = rotate_pts(all_base_pts(:,0:npts-1,19), ico_pts(:,11), 8.0_d*pi/5.0_d, npts)
        all_base_pts(:,0:npts-1,9) = rotate_pts(all_base_pts(:,0:npts-1,18), ico_pts(:,7), 8.0_d*pi/5.0_d, npts)

        return
    end subroutine


    function rotate_pts(old_pts, axis, angle, nrotate)
        !Rotate the (npts, 3) array pts (in cartesians) around the axis given by a
        !position vector (3) on sphere, through given anti-clockwise angle in radians.
        !- use Rodrigues' rotation formula
        implicit none
        real(num):: angle, dot
        real(num), dimension(0:2):: axis, pt
        real(num), dimension(0:2, 0:nrotate-1):: old_pts, rotate_pts
        real(num):: cosa, sina
        integer:: k, nrotate
        cosa = cos(angle)
        sina = sin(angle)
        do k = 0, nrotate - 1
            pt = old_pts(:,k)
            dot = axis(0)*pt(0) + axis(1)*pt(1) + axis(2)*pt(2)
            rotate_pts(0,k) = pt(0)*cosa+ (axis(1)*pt(2) - axis(2)*pt(1))*sina + axis(0)*dot*(1 - cosa)
            rotate_pts(1,k) = pt(1)*cosa+ (axis(2)*pt(0) - axis(0)*pt(2))*sina + axis(1)*dot*(1 - cosa)
            rotate_pts(2,k) = pt(2)*cosa+ (axis(0)*pt(1) - axis(1)*pt(0))*sina + axis(2)*dot*(1 - cosa)
        end do

    end function rotate_pts

    function findk(u, v, n)
        integer:: u, v, n, findk
        findk = int(((2*n + 3)*u - u**2)/2) + v - u
    end function findk

    function edge_index(side, i, n)
        integer:: side, i, n, edge_index

        if (side == 0) then
            edge_index = findk(n -i, n- i, n)
        else if (side == 1) then
            edge_index = findk(i, n, n)
        else
            edge_index = findk(0, i, n)
        end if
    end function edge_index

    function side_tri(side, i, n)
        integer:: side, i, i2, n, side_tri
        !returns index of triangles along a certain side.
        !Ordered CLOCKWISE to match up with the other segment, where it will go in the other direction.
        if (side == 0) then
            i2 = 2*n - i - 2
            side_tri = 2*n*(i2/2) - (i2/2)**2 + mod(i2, 2)
        else if (side == 1) then
            side_tri = (2*n-2)*((i/2)+1) - (i/2)**2 - mod(i, 2)
        else
            side_tri = i
        end if

    end function side_tri

    function findt(u, v, n, up)
        integer:: u, v, n, findt
        integer:: up
        if (up == 1) then
            findt = 2*(v-1) + 2*(n-1)*u - u**2
        else
            findt = 2*v-1 + 2*(n-1)*u - u**2
        end if
    end function


    function circumcentre(p0, p1, p2)
        !Compute circumcentre of triangle with vertices p0, p1, p2. Project on sphere.
        real(num), dimension(0:2):: circumcentre
        real(num), dimension(0:2):: p0, p1, p2, v01, v02
        real(num):: cx, cy, cz, mag
        v01 =  p1 - p0; v02 = p2 - p0
        cx = v01(1)*v02(2) - v01(2)*v02(1)
        cy = v01(2)*v02(0) - v01(0)*v02(2)
        cz = v01(0)*v02(1) - v01(1)*v02(0)
        mag = sqrt(cx**2 + cy**2 + cz**2)
        circumcentre = [cx/mag, cy/mag, cz/mag]

    end function circumcentre

    function distance(p0, p1)
        ! Calculate great-circle distance between points p0 and p1 on surface of (unit) sphere.
        real(num):: distance
        real(num), dimension(0:2):: p0, p1
        distance = acos(min(1.0_d, p0(0)*p1(0) + p0(1)*p1(1) + p0(2)*p1(2) ))

    end function distance

    function triangle_area(p0, p1, p2)
        ! Calculate area of spherical triangle with vertices p0, p1, p2 (3) in anti-clockwise order.
        !For formula see Eriksson: https://www.jstor.org/stable/2691141?seq=1#metadata_info_tab_contents
        real(num):: triangle_area
        real(num), dimension(0:2):: p0, p1, p2
        real(num):: stp, dot01, dot02, dot12
        stp = abs(p0(0)*(p1(1)*p2(2) - p1(2)*p2(1)) + p0(1)*(p1(2)*p2(0) - p1(0)*p2(2)) + p0(2)*(p1(0)*p2(1) - p1(1)*p2(0)))
        dot01 = p0(0)*p1(0) + p0(1)*p1(1) + p0(2)*p1(2)
        dot02 = p0(0)*p2(0) + p0(1)*p2(1) + p0(2)*p2(2)
        dot12 = p1(0)*p2(0) + p1(1)*p2(1) + p1(2)*p2(2)
        triangle_area = 2.0_d*atan(stp/(1.0_d + dot01 + dot02 + dot12))

    end function triangle_area

    subroutine printmat(array,nx,ny)  !prints a matrix with dimensions nx, ny the right way around
        implicit none
        real(num):: array(nx,ny)
        integer:: nx, ny, i, j
        if (proc_num == 0) then
            do j = 1, ny, 1
                write(*,10)(array(i,j), i=1,nx)
            end do
            print*, ''
            10 format(100f10.4)  !Maximum matrix size 100
        end if
        return
    end subroutine

END MODULE grid
