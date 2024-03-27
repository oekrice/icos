!*******************************************************************************
PROGRAM fltrace
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code. Theoretically. Want to run with intel debuggers on, otherwise definitely not going to pick up on all the mistakes
!*******************************************************************************
    USE grid

    USE import_export
    USE shared_data

    IMPLICIT NONE

    !Import parameters from text file (saved out by python in the fltrace directory)
    open(1,file= "flparameters.txt")
    read(1, *) flparameters
    close(1)

    use_tweaked_pts = .true.

    ! Put some of the major variables in here - things that can be changed occasionally but not in a series of runs
    plot_num = int(flparameters(0))
    G = int(flparameters(1))
    nlines = int(flparameters(2))
    ds_fact = flparameters(3)
    crot = int(flparameters(4))
    rmax = flparameters(5)
    export_length = int(flparameters(6))
    view_angle = flparameters(7)
    hamilton_flag = int(flparameters(8))
    run_id = int(flparameters(9))
    nopen = int(flparameters(10))

    seg_min = 0; seg_max = 19

    call establish_grid()  !Establish the grid

    call read_bfield()  !Read in the appropriate magnetic field

    call bfield_ghosts_1()  !uses MPI
    call bfield_bound()
    call bfield_ghosts_2()   !This is the same as bfield_ghosts_1 but only transfers the cells at the top and bottom (those affected by bfield_bound)

    call find_starts()

    call b_to_gridpts()  !Averages B to gridpoints using the RMS approach
    print*, 'Magnetic field averaged to gridpoints'

    call integrate_fieldlines()

    call export_fieldlines()

    if (nopen > 0) then
        ! Also calculate the fieldlines necessary for a side-on approach.
        ! Uses the same tracing concepts but different conditions and export filename
        ntop = nopen; nbottom = nopen  !number of lines to trace from the top and bottom
        call find_plane_starts()

        call integrate_plane_fieldlines()

        call export_plane_fieldlines()

    end if
    contains

    subroutine b_to_gridpts()
        implicit none
        real(num), dimension(:,:,:,:):: a_vec(0:2,0:11)  !Transposed normal vectors at each point
        real(num), dimension(:,:,:,:):: a_mat(0:2,0:2)  !3x3 matrix at each point
        real(num), dimension(:):: rhs(0:2)
        integer:: seg, k, j, i1, i2, ri

        allocate(b0(0:2,1:nr+1,0:npts-1,0:19))  !The magnetic field that has been interpolated, in cartesian coordinates
        !Establish matrix of normal vectors for each point. Require 12 of these (in general)
        do seg = 0, 19
        do k = 0, npts-1
        !Find horizontal normal vectors
        a_vec = 0.0_num
        a_mat = 0.0_num

        do j = 0, nnbrs(k)-1  !Some of them have 5 unfortunately so need to be careful
        a_vec(:,j) = trinorms(:,trimap(1,j,k),trimap(0,j,k),seg)
        end do
        do j = 0, nnbrs(k)-1
        a_vec(:,j+nnbrs(k)) = centres(:,trimap(0,j,k),seg)
        end do

        do i1 = 0,2  !Components of the matrix obtained with the ATA of a_vec^T
        do i2 = 0,2
        a_mat(i1,i2) = sum(a_vec(i1,0:2*nnbrs(k)-1)*a_vec(i2,0:2*nnbrs(k)-1))
        end do
        !print*, a_mat(i1,:)
        end do
        !print*, ''
        !!!!!! ALL OF THE ABOVE ONLY NEEDS CALCULATING INITIALLY
        !Multiply by magnetic field to obtain the rhs that needs to be solved
        do ri = 1, nr+1  !Is now radial variation so need to do this
        rhs = 0.0_num

        do j = 0, nnbrs(k)-1

            rhs(:) = rhs(:) + 0.5*a_vec(:,j)*(bh(ri-1,trimap(1,j,k),trimap(0,j,k),seg) + &
            bh(ri,trimap(1,j,k),trimap(0,j,k),seg)) !Horizontal vectors

            rhs(:) = rhs(:) + a_vec(:,j+nnbrs(k))*br(ri,trimap(0,j,k),seg)  !Vertical vectors


        end do
        !Solve for this system at all the grid points
        b0(0:2,ri,k,seg) = solve3d(a_mat, rhs)

        end do

        end do
        end do

        return

    end subroutine b_to_gridpts


    subroutine find_plane_starts()
        !Find start points for evenly distributed lines on the top and bottom
        implicit none
        integer:: line, k, j, seg_test, count1
        real(num):: dtheta, theta, phi
        real(num), dimension(:):: startpt(0:2)
        integer, dimension(:,:):: start_ind(0:2)

        deallocate(starts); deallocate(start_inds)
        allocate(starts(0:2,0:ntop+nbottom-1))
        allocate(start_inds(0:2,0:ntop+nbottom-1))

        count1 = 0
        do j = 0, 1
            do line = 0, ntop/2-1
                !Tops
                dtheta = pi/float(ntop/2)
                theta = dtheta*line + dtheta/2
                phi = 0.01_num + pi*j - view_angle*(pi/180.)
                startpt = 0.99*rmax*(/sin(theta)*cos(phi), sin(theta)*sin(phi), cos(dtheta*line+dtheta/2)/)
                !Find initial cell indices of these
                do seg_test  = 0, 19
                do k = 0, ntri-1
                if (isin_triangle(startpt, (/nr,k,seg_test/))) then
                    start_ind = (/nr,k,seg_test/)
                    starts(0:2,count1) = startpt
                    start_inds(0:2,count1) = start_ind
                    count1 = count1 + 1
                end if
                end do

                end do
            end do
        end do
        do line = 0, nbottom-1
            !Bottoms
            dtheta = pi/float(nbottom)
            theta = dtheta*line + dtheta/2
            phi = 0.01_num
            startpt = 1.01*(/sin(theta)*cos(phi), sin(theta)*sin(phi), cos(dtheta*line+dtheta/2)/)
            !Find initial cell indices of these
            do seg_test  = 0, 19
            do k = 0, ntri-1
            if (isin_triangle(startpt, (/1,k,seg_test/))) then
                start_ind = (/1,k,seg_test/)
                starts(0:2,count1) = startpt
                start_inds(0:2,count1) = start_ind
                count1 = count1 + 1
                continue
            end if
            end do
            end do

        end do


        return

    end subroutine find_plane_starts

    subroutine integrate_plane_fieldlines()
        implicit none
        integer:: line
        !real(num), dimension(:):: b_direction(0:2)
        real(num), dimension(:):: start_grad(0:2)
        integer:: ud, i, skip

        ds = ds_fact*minval(d_lens(0,:,:))

        max_length = min(100000, int(50.0/ds))

        deallocate(new_line); deallocate(export_lines)
        allocate(new_line(0:2,0:max_length-1))
        allocate(export_lines(0:2,0:export_length-1,0:ntop+nbottom-1))

        export_lines = 0.0_num

        do line = 0, ntop+nbottom-1
            !print*, starts(:, line), start_inds(:,line), ud

            start_grad = b_direction_new(starts(:,line), start_inds(:,line))
            if (sum(abs(start_grad)) < 1e-10) then
                ud = 0
            else
                ud = int(abs(sum(start_grad*starts(:,line)))/sum(start_grad*starts(:,line)))
            end if

            if (ud .ne. 0) then
                if (start_inds(0,line) == 1) then
                    call integrate_line(starts(:,line),start_inds(:,line),ud,max_length)
                else
                    call integrate_line(starts(:,line),start_inds(:,line),-ud,max_length)
                end if
            end if

            skip = max(1,int((current_length-1)/(export_length-2)))

            if (current_length > 0) then
             do i = 0, min(export_length, current_length) - 2
                export_lines(:,i,line) = new_line(:,i*skip)
             end do
            export_lines(:,min(export_length, current_length)-2,line) = new_line(:,current_length-1)
            end if
            export_lines(:,export_length-1,line) = float(ud)   !Add flag saying which way the line was integrated
            !print*, 'Line', line, 'integrated'

        end do
        print*, 'Open fieldlines integrated'

        return

    end subroutine

    subroutine find_starts()
        !Find the start points for randomly distributed points over the whole surface
        implicit none
        integer:: seg, k, startcount, newcount, skip
        real(num), dimension(:,:):: allstarts(0:2,0:ntri*20 - 1)
        integer, dimension(:,:):: all_inds(0:2,0:ntri*20-1)

        allocate(starts(0:2,0:nlines-1))
        allocate(start_inds(0:2,0:nlines-1))

        startcount = 0
        do seg = 0, 19
        do k = 0, ntri-1
        allstarts(:,startcount) = 1.001*centres(:,k,seg)
        all_inds(:,startcount) = (/1,k,seg/)
        startcount = startcount + 1
        end do
        end do
        if (nlines > ntri*20) nlines = ntri*20
        !Thin these out based on nlines
        skip = int((ntri*20)/nlines) + 1
        newcount = 0; startcount = skip-1

        do while (newcount < nlines)
            starts(:,newcount) = allstarts(:,mod(startcount,ntri*20))
            start_inds(:,newcount) = all_inds(:,mod(startcount,ntri*20))
            newcount = newcount + 1
            startcount = startcount + skip
        end do
        return

    end subroutine find_starts


    subroutine integrate_fieldlines()
        implicit none
        integer:: line, smoothfact, sf_real
        real(num), dimension(:):: start_grad(0:2)
        integer:: ud, i, skip

        ds = ds_fact*minval(d_lens(0,:,:))

        max_length = min(1000000, int(1000.0/ds))

        allocate(new_line(0:2,0:max_length-1))
        allocate(export_lines(0:2,0:export_length-1,0:nlines-1))
        allocate(smooth_line(0:2,0:max_length-1))

        export_lines = 0.0_num

        do line = 0, nlines-1
            start_grad = b_direction_new(starts(:,line), start_inds(:,line))
            if (sum(abs(start_grad)) < 1e-10) then
                ud = 0
            else
                ud = int(abs(sum(start_grad*starts(:,line)))/sum(start_grad*starts(:,line)))
            end if

            if (ud .ne. 0) then
                call integrate_line(starts(:,line),start_inds(:,line),ud,max_length)  !Actually integrates the line. Not trivial...
            end if

            skip = max(1,int((current_length-1)/(export_length-1)))

            smoothfact = 0!int(1/ds_fact)  !How much the lines should be smoothed
            if (current_length > 0) then
            do i = 0, current_length - 1
                sf_real = 1 + min(smoothfact, i) + min(current_length - 1 - i, smoothfact)
                smooth_line(:,i) = sum(new_line(:,max(0,i-smoothfact):min(current_length-1,i+smoothfact)),2)/sf_real
             end do

             do i = 0, min(export_length, current_length) - 2
                export_lines(:,i,line) = smooth_line(:,i*skip)  !May be better to take an average?
             end do
            export_lines(:,min(export_length, current_length)-1,line) = new_line(:,current_length-1)
            !print*, 'Line', line, 'integrated'
            end if
        end do
        print*, 'Closed fieldlines integrated'
        return

    end subroutine

    subroutine integrate_line(pt,ind,ud,max_length)
        implicit none
        real(num), dimension(:):: pt(0:2)
        real(num), dimension(:,:):: line_grads(0:2,0:max_length-1)
        real(num):: r

        integer:: max_length, lc, ud, fc
        integer, dimension(:):: ind(0:2)
        integer:: ri, k, seg, i, i1, i2, kt, kt1, kt2, seg_test
        logical:: keep_going, tri_found

        ri = ind(0); k = ind(1); seg = ind(2)

        lc = 0; fc = 0
        new_line = 0.0_num
        new_line(:,lc) = pt
        !line_grads(0:2,lc) = b_direction(pt, ind)
        line_grads(0:2,lc) = b_direction_new(pt, ind)

        keep_going = .true.
        do while (keep_going)
            !print*, lc, pt, line_grads(0:2,lc), seg,k,ri
            lc = lc + 1
            !Add to line
            if (lc == 1) then
                pt = pt + ud*ds*line_grads(0:2,lc-1)
            else
                pt = pt + ud*1.5*ds*line_grads(0:2,lc-1) - ud*0.5*ds*line_grads(0:2,lc-2)
            end if
            !Find new cell indices
            r = sqrt(sum(pt**2))
            if (r >= rs(ri) .and. r <= rs(ri+1)) then
                ri = ri
            else if (r >= rs(ri+1) .and. r <= rs(ri+2)) then
                ri = ri + 1
            else if (r >= rs(ri-1) .and. r <= rs(ri)) then
                ri = ri - 1
            else
                keep_going = .false.  !A whole segment has been skipped
                print*, 'Entire radial segment skipped, suggest reduce integration step'
            end if
            !Find triangle. If not in the same cell, check the neghbours and then neighbours of neighbours (which should be rare)
            tri_found = .false.
            if (isin_triangle(pt, ind)) then
                !Triangle still correct, don't move
                k = k
                tri_found = .true.
            end if
            if (.not. tri_found) then
                do i = 0, 2
                    kt = nbrfaces(i,k)
                    if (isin_triangle(pt, (/ri,kt,seg/))) then
                        !Neighbour triangle is the correct one
                        tri_found = .true.

                        if (kt < ntri) then  !Not a ghost point
                            k = kt
                        else !Is a ghost point, switch segment
                            k = tri_transfer(1,kt-ntri,seg)
                            seg = tri_transfer(0,kt-ntri,seg)
                        end if
                    end if
                end do
            end if
            if (.not. tri_found) then  !Check neighbours of neighbours. Will miss a few around the poles but meh. Reduce dt if this is a problem
                do i1 = 0, 2
                    kt1 = nbrfaces(i1,k)
                    if (kt1 < ntri) then
                        do i2 = 0, 2
                            kt2 = nbrfaces(i2, kt1)
                            if (isin_triangle(pt, (/ri,kt2,seg/))) then
                                !Neighbour triangle is the correct one
                                tri_found = .true.
                                if (kt2 < ntri) then  !Not a ghost point
                                    k = kt2
                                else !Is a ghost point, switch segment
                                    k = tri_transfer(1,kt2-ntri,seg)
                                    seg = tri_transfer(0,kt2-ntri,seg)
                                end if
                            end if
                        end do
                    end if
                end do
            end if
            if (fc < 100 .and. (.not. tri_found)) then  !Check all triangles. Can turn off for speed...
                fc = fc + 1
                do seg_test = 0, 19
                do kt = 0, ntri-1
                if (isin_triangle(pt, (/ri,kt,seg_test/))) then
                    k = kt
                    seg = seg_test
                    tri_found = .true.
                end if
                end do
                end do

            end if

            if (.not. tri_found) then
                keep_going = .false.
                new_line = 0.0_num
                !Failed to find new triangle. Discard this line
                lc = 0
            end if

            if (lc + 1 == max_length) then
                keep_going = .false.
                print*, 'Infinite rope may be present'
                !Line is too long
            end if

            if ((r > 2.5) .or. (r < 1.0)) then
                keep_going = .false.
                !Reached the boundary
            end if

            if (keep_going) then
                ind = (/ri, k, seg/)
                new_line(:,lc) = pt
                line_grads(:,lc) = b_direction_new(pt, ind)
            end if
        end do
        current_length = lc
        !print*, sqrt(sum(new_line(:,0)**2)), sqrt(sum(new_line(:,lc-1)**2)), lc
        return
    end subroutine

    function b_direction_new(pt, ind)
    real(num), dimension(:):: pt(0:2), b_direction_new(0:2), e_r(0:2)
    real(num), dimension(:):: b_g(0:2), x(0:2)
    real(num), dimension(:,:):: tripts(0:2,0:2)
    integer, dimension(:):: ind(0:2)
    real(num):: r, prop_r, alpha, beta, gamm, mag
    integer:: seg, k, ri, i

    ri = ind(0); k = ind(1); seg = ind(2)

    r = sqrt(sum(pt**2))
    e_r = pt/r

    b_direction_new = (/0.0, 0.0, 0.0/)

    do i = 0, 2
        tripts(:,i) = all_base_pts(:,triangs(i,k),seg)
    end do

    prop_r = (r - rs_global(ri)) / ((rs_global(ri+1) - rs_global(ri))) !Distance up this cell

    x = coeffs(e_r, tripts(:,0), tripts(:,1), tripts(:,2))
    x = x/sum(x)  !Project into plane

    alpha = x(0); beta = x(1); gamm = 1 - x(0) - x(1)

    b_g = 0.0

    do i = 0,2
        b_g = b_g + 0.5*x(i)*((1.0-prop_r)*b0(0:2,ri,triangs(i,k),seg) + prop_r*b0(0:2,ri+1,triangs(i,k),seg))


    end do

    mag = sqrt(sum(b_g**2))

    if (mag < 1e-8) then
        b_direction_new = (/0.0, 0.0, 0.0/)
    else
        b_direction_new = b_g/mag
    end if

    end function b_direction_new


    subroutine bfield_ghosts_1()
        !As for bfield_ghosts_1 but only transfers the top and bottom layers, for speed purposes
        implicit none
        integer:: gi, i
        integer:: seg, seg_source
        !Transfer ghosts within the process
        do seg = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg)
                br(0:2, ntri+gi, seg) = br(0:2, tri_transfer(1,gi,seg), seg_source)
                br(nr:nr+2, ntri+gi, seg) = br(nr:nr+2, tri_transfer(1,gi,seg), seg_source)

                do i = 0, 2
                    bh(nr-1:nr+1,i,ntri+gi,seg) = &
                    bh(nr-1:nr+1,mod(i+tri_transfer(2,gi,seg), 3),  tri_transfer(1,gi,seg), seg_source)

                    bh(0:2,i,ntri+gi,seg) = &
                    bh(0:2,mod(i+tri_transfer(2,gi,seg), 3),  tri_transfer(1,gi,seg), seg_source)
                end do
            end do
        end do
        return
    end subroutine bfield_ghosts_1

    subroutine bfield_ghosts_2()
        !Populates magnetic field ghost points, obtaining data from adjacent segments.
        implicit none
        integer:: gi, i
        integer:: seg, seg_source
        !Transfer ghosts within the process
        do seg = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg)
                br(0:nr+2, ntri+gi, seg) = br(0:nr+2, tri_transfer(1,gi,seg), seg_source)
                do i = 0, 2
                    bh(0:nr+1,i,ntri+gi,seg) = &
                    bh(0:nr+1,mod(i+tri_transfer(2,gi,seg), 3),  tri_transfer(1,gi,seg), seg_source)
                end do
            end do
        end do
        return
    end subroutine bfield_ghosts_2

    subroutine bfield_bound()
        !Imposes the boundary conditions on the upper and lower boundaries.
        !Does some MPI between the layers as necessary.
        !Will eventually have the option to choose radial or current-free condition on the top. But will pick current-free for now

        implicit none
        real(num), dimension(:,:):: bh0_field(0:2,0:ntri-1,seg_min:seg_max)   !Lower boundary
        real(num), dimension(:,:):: bh1_field(0:2,0:ntri-1,seg_min:seg_max)   !Upper boundary
        real(num):: diff   !Difference in current. Keep within a segment
        integer:: i,k,seg
        !TOP BOUNDARY
        if (seg_layer == nlayers - 1) then
        !On the actual top of the domain, so apply the real boundary condition.
        do seg  = seg_min, seg_max
            do i = 0, 2
                do k = 0, ntri-1
                    diff = 0.0_num
                    diff = diff - bh(nr,i,k,seg)*all_ds(nr,i,k,seg)
                    diff = diff + br(nr+1,k,seg)*(rc(nr+1) - rc(nr))
                    diff = diff - br(nr+1,nbrfaces(i,k),seg)*(rc(nr+1) - rc(nr))
                    if (.true.) then  !Current-free condition
                    bh1_field(i,k,seg) = -diff/all_ds(nr+1,i,k,seg)
                    else  !Radial condition
                    bh1_field(i,k,seg) = -bh(nr,i,k,seg)
                    end if
                    bh(nr+1,i,k,seg) = bh1_field(i,k,seg)
                end do
            end do
        end do
        end if

        !BOTTOM BOUNDARY
        if (seg_layer == 0) then
        !On the actual bottom of the domain, so apply the real boundary condition.
        do seg  = seg_min, seg_max
            do i = 0, 2
                do k = 0, ntri-1
                    diff = 0.0_num
                    diff = diff - bh(1,i,k,seg)*all_ds(1,i,k,seg)
                    diff = diff - br(1,k,seg)*(rc(1) - rc(0))
                    diff = diff + br(1,nbrfaces(i,k),seg)*(rc(1) - rc(0))
                    bh0_field(i,k,seg) = -diff/all_ds(0,i,k,seg)
                    bh(0,i,k,seg) = bh0_field(i,k,seg)
                end do
            end do
        end do
        end if

        if (seg_layer == 0) then
        !Radial field ghosts are determined by the solenoidal condition (I missed this, thanks Anthony)
        do seg = seg_min, seg_max
            do k = 0, ntri-1
                diff = 0
                do i = 0, 2
                    diff = diff - bh(0,i,k,seg)*v_areas(0,i,k,seg)
                end do
                diff = diff - br(1,k,seg)*h_areas(1,k,seg)
                br(0,k,seg) = -diff/h_areas(0,k,seg)
            end do
        end do
        end if

        if (seg_layer == nlayers - 1) then
        !Radial field ghosts are determined by the solenoidal condition (I missed this, thanks Anthony)
        do seg = seg_min, seg_max
            do k = 0, ntri-1
                diff = 0
                do i = 0, 2
                    diff = diff - bh(nr+1,i,k,seg)*v_areas(nr+1,i,k,seg)
                end do
                diff = diff + br(nr+1,k,seg)*h_areas(nr+1,k,seg)
                br(nr+2,k,seg) = diff/h_areas(nr+2,k,seg)
            end do
        end do
        end if

        return
    end subroutine bfield_bound

    function isin_triangle(pt, ind)

    logical:: isin_triangle
    integer, dimension(:):: ind(0:2)
    integer:: i, ri, k, seg
    real(num), dimension(:):: pt(0:2), x(0:2), e_r(0:2)
    real(num), dimension(:,:):: tripts(0:2,0:2)

    ri = ind(0); k = ind(1); seg = ind(2)

    do i = 0, 2
        tripts(:,i) = all_base_pts(:,triangs(i,k),seg)
    end do
    e_r = pt/sqrt(sum(pt**2))
    x = coeffs(e_r, tripts(:,0), tripts(:,1), tripts(:,2))

    if (minval(x) >= 0.0) then  !Linear combination of all triangle points is strictly positive
        isin_triangle = .true.
    else
        isin_triangle = .false.
    end if

    end function isin_triangle

END PROGRAM fltrace

























