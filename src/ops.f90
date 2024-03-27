!*******************************************************************************
MODULE ops
!*******************************************************************************
! Module containing all the operations (or at least some of them...). The actual timestep will be in main.f90 but try to put everything else in here. Should be fully paralell.
!*******************************************************************************
    USE shared_data
    USE grid

    IMPLICIT NONE

    contains

    subroutine allocate_variable_arrays()
        ! Allocates arrays. No need to do this more than necessary
        implicit none
        allocate(br(0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(bh(0:nr+1,0:2,0:ntri+nghosttri-1,seg_min:seg_max))
        br = 0.0_num; bh = 0.0_num

        allocate(jr(0:nr+1,0:npts+nghosts-1,seg_min:seg_max))
        allocate(jh(0:nr+2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))
        jr = 0.0_num; jh = 0.0_num

        allocate(diff_r(0:nr+1,0:npts-1,seg_min:seg_max))
        allocate(diff_h(0:nr+2,0:2,0:ntri-1,seg_min:seg_max))
        diff_r = 0.0_num; diff_h = 0.0_num

        allocate(er(0:nr+1,0:npts+nghosts-1,seg_min:seg_max))
        allocate(eh(0:nr+2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))

        allocate(e0(0:2,0:nr+2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))   !Directions, then radial coordinate then triangle side.
        allocate(j0(0:2,0:nr+2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))   !Directions, then radial coordinate then triangle side.


        allocate(b0(0:2,0:nr+2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))   !Directions, then radial coordinate then triangle side.
        allocate(v0(0:2,0:nr+2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))   !Directions, then radial coordinate then triangle side.
        allocate(b_up(0:2,0:nr+2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))   !As above, but saved on radial grid faces

        allocate(j1(0:2,0:nr+2,0:npts-1,seg_min:seg_max))
        allocate(b1(0:2,0:nr+2,0:npts-1,seg_min:seg_max))  !Saved on radial grid points
        allocate(v1(0:2,0:nr+2,0:npts-1,seg_min:seg_max))
        allocate(bu(0:2,0:nr+2,0:npts-1,seg_min:seg_max))

        allocate(jh1(0:2,0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max))
        allocate(eh1(0:2,0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max))

        allocate(e1(0:2,0:nr+2,0:npts-1,seg_min:seg_max))
        allocate(e1_vout(0:2,0:nr+2,0:npts-1,seg_min:seg_max))


        allocate(vout_field(0:nr+2))
        allocate(fric_field(0:nr+2,0:npts-1,seg_min:seg_max))


        allocate(b2(0:nr+2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))

        f_orders(0,:) = (/2,0,2,2,2/)
        f_orders(1,:) = (/1,2,2,2,1/)
        f_orders(2,:) = (/0,1,2,2,0/)
        corner_pts_1 = (/0,n,npts-1/)
        corner_pts_2 = (/0,npts-1,n/)
        return

    end subroutine

    subroutine find_dt()
        implicit none

        real(num):: dx
        if (proc_num == 0) print*, ''
        if (proc_num == 0) print*, 'Determining timestep...'
        dx = minval(tri_l)*rsun

        dt = 1e6
        if (nu0 > 0.0) then
            if (constant_friction) then
               dt = min(dt, 0.2*dx**2/nu0)
            else
               dt = min(dt, 0.2*dx**2/(nu0*rmax))
            end if
            if (proc_num == 0) print*, 'nu0 dt', 0.2*dx**2/nu0
        end if
        if (eta > 0.0) then
            dt = min(dt, 0.075*dx**2/eta)
            if (proc_num == 0) print*, 'eta dt', 0.075*dx**2/eta
        end if
        if (voutfact > 0.0) then
            dt = min(dt, 2.0*dx/voutfact)  !Bodging this doesn't appear to be harmful in general and speeds it up
            if (proc_num == 0) print*, 'outflow dt', 2.0*dx/voutfact
        end if
        if (shearfact > 0.0) then
            dt = min(dt, 0.2*dx/shearfact)
            if (proc_num == 0) print*, 'shear dt', 0.2*dx/shearfact
        end if
        if (eta0 > 0.0) then
            dt = min(dt, 0.075*dx**2/eta0)
            if (proc_num == 0) print*, 'eta0 dt', 0.075*dx**2/eta0
        end if

        if (proc_num == 0) print*, 'nocfldt' ,dt
        dt = dt*cfl
        if (proc_num == 0) print*, 'Required timestep', dt
        nt = int(tmax/dt) + 1
        dt = tmax/nt
        if (proc_num == 0) print*, 'Actual timestep  ', dt
        if (proc_num == 0) print*, ''

        return

    end subroutine find_dt

    function vout_fn(r)
    ! The outflow velocity function
    real(num):: r_crit, vout_fn, r

    r_crit = 10.0
    if (parker) then
        vout_fn = voutfact*((rmax**2*exp(-2.0*r_crit/r))/(r**2*exp(-2.0*r_crit/rmax)))
    else
        vout_fn = voutfact*((r)/(rmax))**11.5
    end if

    end function vout_fn

    function shear_fn(lat)
    ! The shearing velocity as a function of latitude (radians from the equator). In radians/day which is the same as solar radii/day at the equator. Setting a time unit as 1 day.
    real(num):: lat, shear_fn

    if (shear_frame == 1) then  !Frame 1, rotates relative to background stars, not the carrington frame
        shear_fn = shearfact*(0.2507-0.04009*sin(lat)**2 - 0.02834*sin(lat)**4)*cos(lat)
    else if (shear_frame == 2) then  !Frame 2, carrington frame (subtract 0.2304 rads/day, corresponding to rotation time of 27.26 days)
        shear_fn = shearfact*(0.2507-0.2304-0.04009*sin(lat)**2 - 0.02834*sin(lat)**4)*cos(lat)
    else if (shear_frame == 3) then  !Testing frame
        shear_fn = shearfact*(0.2507-0.2200-0.04009*sin(lat)**2 - 0.02834*sin(lat)**4)*cos(lat)
    else  !not real shear. Sort of like the 2D ones.
        shear_fn = (shearfact*(1.0-1.1455*sin(lat)**2 - 0.8544*sin(lat)**4)*cos(lat))
    end if

    end function shear_fn

    subroutine average_edges()
        !As certain bits of information are stored twice (like ah), this checks that everything actually matches up
        ! It is probably possible to avoid storing things twice by being smart, but I'm not so there
        ! This is also a useful check that there aren't any bugs and everything is stored fine

        !Does an a ghost transfer and then does the averaging. I think that should be fairly foolproof
        implicit none
        integer:: ri, k, i, gi, nbr, proc_source, seg_source,seg_target
        integer:: seg_check, proc_target, id
        real(num):: test1, test2
        real(num), dimension(:):: ah_data(1:nr+1)
        !Transfer ghosts within the process.
        do seg = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg)
                proc_source = seg_source/seg_groups + seg_layer*20
                if (proc_source == proc_num) then !Just transfer data over if it's in the process
                    do i = 0, 2
                        ah(1:nr+1, i,ntri+gi, seg) = &
                        ah(1:nr+1,mod(i+tri_transfer(2,gi,seg), 3),  tri_transfer(1,gi,seg), seg_source)
                    end do
                end if
            end do
        end do
        !MPI sends
        do seg_source = seg_min, seg_max
        !Seg_source is the segment that is sending the information.
        !This is from process proc_num
            do seg_check = 0, 19
            !Run through all (other) segments and check if they need any information from here
            !Check ghost points of all other segments (only order n, luckily). A bit messy but meh.
            proc_target = seg_check/seg_groups + seg_layer*20
            if (proc_target .ne. proc_num) then !Do need MPI here
            do gi = 0, nghosttri-1
                if (tri_transfer(0,gi,seg_check) == seg_source) then
                    do i = 0, 2
                        ah_data = ah(1:nr+1, mod(tri_transfer(2,gi,seg_check) + i,3), tri_transfer(1,gi,seg_check), seg_source)
                        id = gi + 1000*seg_check + 100000*(i+1)
                        call mpi_send(ah_data, nr+1, MPI_DOUBLE_PRECISION, proc_target, id, comm, ierr)
                    end do
                end if
            end do
            end if
            end do
        end do
        call MPI_barrier(comm, ierr)
        !MPI receives
        do seg_target = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg_target)
                proc_source = seg_source/seg_groups + seg_layer*20
                if (proc_source .ne. proc_num) then
                do i = 0, 2
                    id = gi + 1000*seg_target + 100000*(i+1)
                    call mpi_recv(ah_data, nr+1, MPI_DOUBLE_PRECISION, proc_source, id, comm, MPI_STATUS_IGNORE, ierr)
                    ah(1:nr+1, i, ntri+gi, seg_target) = ah_data
                end do
                end if
            end do
        end do
        call MPI_barrier(comm, ierr)
        !Ghosts populated. Now do the averaging
        !Interior points first - no MPI required (yay!)
        test1 = -1; test2 = -1
        do seg = seg_min, seg_max
            do ri = 1, nr+1 !Radial layer
            do i = 0 ,2  !Triangle side
                do k = 0, ntri-1  !Triangle
                    test1 = ah(ri,i,k,seg)
                    nbr = nbrfaces(i,k)
                    if (nbr < ntri) then  !This side is not on an edge
                    if (updown(k) == 1) then
                        test2 = ah(ri,mod(i+2,3),nbr,seg)  !test1 and test2 should be equal and opposite
                    else
                        test2 = ah(ri,mod(i+1,3),nbr,seg)
                    end if
                    if (test1+test2 > 1d-14 .and. seg == 0) then
                        if (updown(k) == 1) then
                            print*, 'Interior edges not matching', test1, test2, seg, ri, k, i
                            call mpi_abort(comm, ierr)
                        else
                            print*, 'Interior edges not matching', test1, test2, seg, ri, k, i
                            call mpi_abort(comm, ierr)
                        end if
                    end if
                    ah(ri,i,k,seg) = sign(0.5_num*(abs(test1) + abs(test2)),ah(ri,i,k,seg))
                    if (updown(k) == 1) then
                        ah(ri,mod(i+2,3),nbr,seg) = -test1
                    else
                        ah(ri,mod(i+1,3),nbr,seg) = -test1
                    end if
                    end if
                end do
            end do
            end do
        end do

        return
    end subroutine


    subroutine calculate_magnetic()
        ! In a local process, calculates the magnetic field from the vector potential.
        ! Only strictly internal points here
        ! This outputs the magnetic field strength, NOT fluxes
        implicit none

        integer:: i, k
        br = 0.0_num; bh = 0.0_num

        do seg = seg_min, seg_max
            !Calculate fluxes (Stokes' Theorem)
            do k = 0, ntri-1
                do i = 0, 2
                    br(:,k,seg) = br(:,k,seg) + ah(:,i,k,seg)

                    bh(1:nr,i,k,seg) = bh(1:nr ,i,k,seg) - ah(2:nr+1,i,k,seg)
                    bh(1:nr,i,k,seg) = bh(1:nr ,i,k,seg) + ah(1:nr  ,i,k,seg)

                    bh(1:nr,i,k,seg) = bh(1:nr ,i,k,seg) - ar(1:nr,triangs(i,k),seg)
                    bh(1:nr,i,k,seg) = bh(1:nr ,i,k,seg) + ar(1:nr,triangs(mod(i+1,3),k),seg)
                end do
            end do
            ! Divide by areas to give field strength (NOT flux)
            br(:,0:ntri-1,seg)   = br(:,0:ntri-1,seg)  /h_areas(:,0:ntri-1,seg)
            bh(:,:,0:ntri-1,seg) = bh(:,:,0:ntri-1,seg)/v_areas(:,:,0:ntri-1,seg)
        end do
        return

    end subroutine calculate_magnetic

    subroutine bfield_ghosts_1()
        !As for bfield_ghosts_1 but only transfers the top and bottom layers, for speed purposes
        implicit none
        integer:: gi, i, seg_source, proc_source, seg_target
        integer:: seg_check, proc_target, id
        real(num), dimension(:):: br_data(0:nr+2), bh_data(0:nr+1)
        !Transfer ghosts within the process
        do seg = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg)
                proc_source = seg_source/seg_groups + seg_layer*20
                if (proc_source == proc_num) then !Just transfer data over if it's in the process
                    br(0:nr+2, ntri+gi, seg) = br(0:nr+2, tri_transfer(1,gi,seg), seg_source)
                    do i = 0, 2
                        bh(0:nr+1,i,ntri+gi,seg) = &
                        bh(0:nr+1,mod(i+tri_transfer(2,gi,seg), 3),  tri_transfer(1,gi,seg), seg_source)
                    end do
                end if
            end do
        end do

        !MPI sends
        do seg_source = seg_min, seg_max
        !Seg_source is the segment that is sending the information.
        !This is from process proc_num
            do seg_check = 0, 19
            !Run through all (other) segments and check if they need any information from here
            !Check ghost points of all other segments (only order n, luckily). A bit messy but meh.
            proc_target = seg_check/seg_groups + seg_layer*20
            if (proc_target .ne. proc_num) then !Do need MPI here
            do gi = 0, nghosttri-1
                if (tri_transfer(0,gi,seg_check) == seg_source) then
                    br_data = br(0:nr+2, tri_transfer(1,gi,seg_check), seg_source)
                    id = gi + 1000*seg_check
                    call mpi_send(br_data, nr+3, MPI_DOUBLE_PRECISION, proc_target, id, comm, ierr)
                    do i = 0, 2
                        bh_data = bh(0:nr+1, mod(tri_transfer(2,gi,seg_check) + i,3), tri_transfer(1,gi,seg_check), seg_source)
                        id = gi + 1000*seg_check + 100000*(i+1)
                        call mpi_send(bh_data, nr+2, MPI_DOUBLE_PRECISION, proc_target, id, comm, ierr)
                    end do
                end if
            end do
            end if
            end do
        end do
        call MPI_barrier(comm, ierr)
        !MPI receives
        do seg_target = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg_target)
                proc_source = seg_source/seg_groups + seg_layer*20
                if (proc_source .ne. proc_num) then
                id = gi + 1000*seg_target
                call mpi_recv(br_data, nr+3, MPI_DOUBLE_PRECISION, proc_source, id, comm, MPI_STATUS_IGNORE, ierr)
                br(0:nr+2, ntri + gi, seg_target) = br_data
                do i = 0, 2
                    id = gi + 1000*seg_target + 100000*(i+1)
                    call mpi_recv(bh_data, nr+2, MPI_DOUBLE_PRECISION, proc_source, id, comm, MPI_STATUS_IGNORE, ierr)
                    bh(0:nr+1, i, ntri+gi, seg_target) = bh_data
                end do
                end if
            end do
        end do
        call MPI_barrier(comm, ierr)
        return
    end subroutine bfield_ghosts_1

    subroutine bfield_ghosts_2()
        !Populates magnetic field ghost points, obtaining data from adjacent segments.
        implicit none
        integer:: gi, i, seg_source, proc_source, seg_target
        integer:: seg_check, proc_target, id
        real(num), dimension(:):: br_data(0:2), bh_data(0:2)
        !Transfer ghosts within the process
        do seg = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg)
                proc_source = seg_source/seg_groups + seg_layer*20
                if (proc_source == proc_num) then !Just transfer data over if it's in the process
                    br(0:2, ntri+gi, seg) = br(0:2, tri_transfer(1,gi,seg), seg_source)
                    br(nr:nr+2, ntri+gi, seg) = br(nr:nr+2, tri_transfer(1,gi,seg), seg_source)

                    do i = 0, 2
                        bh(nr-1:nr+1,i,ntri+gi,seg) = &
                        bh(nr-1:nr+1,mod(i+tri_transfer(2,gi,seg), 3),  tri_transfer(1,gi,seg), seg_source)

                        bh(0:2,i,ntri+gi,seg) = &
                        bh(0:2,mod(i+tri_transfer(2,gi,seg), 3),  tri_transfer(1,gi,seg), seg_source)
                    end do
                end if
            end do
        end do
        !MPI sends
        do seg_source = seg_min, seg_max
        !Seg_source is the segment that is sending the information.
        !This is from process proc_num
            do seg_check = 0, 19
            !Run through all (other) segments and check if they need any information from here
            !Check ghost points of all other segments (only order n, luckily). A bit messy but meh.
            proc_target = seg_check/seg_groups + seg_layer*20
            if (proc_target .ne. proc_num) then !Do need MPI here
            do gi = 0, nghosttri-1
                if (tri_transfer(0,gi,seg_check) == seg_source) then
                    br_data = br(nr:nr+2, tri_transfer(1,gi,seg_check), seg_source)
                    id = gi + 1000*seg_check
                    call mpi_send(br_data, 3, MPI_DOUBLE_PRECISION, proc_target, id, comm, ierr)
                    do i = 0, 2
                        bh_data = bh(nr-1:nr+1, mod(tri_transfer(2,gi,seg_check) + i,3), tri_transfer(1,gi,seg_check), seg_source)
                        id = gi + 1000*seg_check + 100000*(i+1)
                        call mpi_send(bh_data, 3, MPI_DOUBLE_PRECISION, proc_target, id, comm, ierr)
                    end do
                end if
            end do
            end if
            end do
        end do
        call MPI_barrier(comm, ierr)
        !MPI receives
        do seg_target = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg_target)
                proc_source = seg_source/seg_groups + seg_layer*20
                if (proc_source .ne. proc_num) then
                id = gi + 1000*seg_target
                call mpi_recv(br_data, 3, MPI_DOUBLE_PRECISION, proc_source, id, comm, MPI_STATUS_IGNORE, ierr)
                br(nr:nr+2, ntri + gi, seg_target) = br_data
                do i = 0, 2
                    id = gi + 1000*seg_target + 100000*(i+1)
                    call mpi_recv(bh_data, 3, MPI_DOUBLE_PRECISION, proc_source, id, comm, MPI_STATUS_IGNORE, ierr)
                    bh(nr-1:nr+1, i, ntri+gi, seg_target) = bh_data
                end do
                end if
            end do
        end do
        call MPI_barrier(comm, ierr)

        !MPI sends
        do seg_source = seg_min, seg_max
        !Seg_source is the segment that is sending the information.
        !This is from process proc_num
            do seg_check = 0, 19
            !Run through all (other) segments and check if they need any information from here
            !Check ghost points of all other segments (only order n, luckily). A bit messy but meh.
            proc_target = seg_check/seg_groups + seg_layer*20
            if (proc_target .ne. proc_num) then !Do need MPI here
            do gi = 0, nghosttri-1
                if (tri_transfer(0,gi,seg_check) == seg_source) then
                    br_data = br(0:2, tri_transfer(1,gi,seg_check), seg_source)
                    id = gi + 1000*seg_check
                    call mpi_send(br_data, 3, MPI_DOUBLE_PRECISION, proc_target, id, comm, ierr)
                    do i = 0, 2
                        bh_data = bh(0:2, mod(tri_transfer(2,gi,seg_check) + i,3), tri_transfer(1,gi,seg_check), seg_source)
                        id = gi + 1000*seg_check + 100000*(i+1)
                        call mpi_send(bh_data, 3, MPI_DOUBLE_PRECISION, proc_target, id, comm, ierr)
                    end do
                end if
            end do
            end if
            end do
        end do
        call MPI_barrier(comm, ierr)
        !MPI receives
        do seg_target = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg_target)
                proc_source = seg_source/seg_groups + seg_layer*20
                if (proc_source .ne. proc_num) then
                id = gi + 1000*seg_target
                call mpi_recv(br_data, 3, MPI_DOUBLE_PRECISION, proc_source, id, comm, MPI_STATUS_IGNORE, ierr)
                br(0:2, ntri + gi, seg_target) = br_data
                do i = 0, 2
                    id = gi + 1000*seg_target + 100000*(i+1)
                    call mpi_recv(bh_data, 3, MPI_DOUBLE_PRECISION, proc_source, id, comm, MPI_STATUS_IGNORE, ierr)
                    bh(0:2, i, ntri+gi, seg_target) = bh_data
                end do
                end if
            end do
        end do
        call MPI_barrier(comm, ierr)

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
        integer:: i,k
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

        if (seg_layer > 0) then !Send MPI information downwards
            do seg = seg_min, seg_max
                do k = 0, ntri-1 !These do loops seem to be necessary because MPI is SILLY
                    do i = 0, 2
                        call mpi_send(bh(1,i,k,seg), 1, MPI_DOUBLE_PRECISION, proc_num-20, 10*k+i, comm, ierr)
                    end do
                    call mpi_send(br(2,k,seg), 1, MPI_DOUBLE_PRECISION, proc_num-20, 10*k+4, comm, ierr)
                    call mpi_send(br(1,k,seg), 1, MPI_DOUBLE_PRECISION, proc_num-20, 10*k+5, comm, ierr)
                end do
            end do
        end if

        call mpi_barrier(comm, ierr)
        if (seg_layer < nlayers - 1) then !Receive MPI information from above
            do seg = seg_min, seg_max
                do k = 0, ntri-1
                    do i = 0, 2
                    call mpi_recv(bh(nr+1,i,k,seg), 1, &
                    MPI_DOUBLE_PRECISION, proc_num+20,10*k+i, comm, MPI_STATUS_IGNORE,ierr)
                    end do
                call mpi_recv(br(nr+2,k,seg),1, &
                MPI_DOUBLE_PRECISION, proc_num+20, 10*k+4, comm, MPI_STATUS_IGNORE,ierr)
                call mpi_recv(br(nr+1,k,seg),1, &
                MPI_DOUBLE_PRECISION, proc_num+20, 10*k+5, comm, MPI_STATUS_IGNORE,ierr)
                end do
            end do

        end if
        call mpi_barrier(comm, ierr)

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
        if (seg_layer < nlayers - 1) then !Send MPI information upwards
            do seg = seg_min, seg_max
                do k = 0, ntri-1
                    do i = 0, 2
                        call mpi_send(bh(nr,i,k,seg), 1, MPI_DOUBLE_PRECISION, proc_num+20, 10*k+i, comm, ierr)
                    end do
                    call mpi_send(br(nr,k,seg), 1, MPI_DOUBLE_PRECISION, proc_num+20, 10*k+4, comm, ierr)
                end do
            end do
        end if
        call mpi_barrier(comm, ierr)
        if (seg_layer > 0) then !Receive MPI information from below
            do seg = seg_min, seg_max
                do k = 0, ntri-1
                    do i = 0, 2
                    call mpi_recv(bh(0,i,k,seg), 1, &
                    MPI_DOUBLE_PRECISION, proc_num-20,10*k+i, comm, MPI_STATUS_IGNORE,ierr)
                    end do
                call mpi_recv(br(0,k,seg),1, &
                MPI_DOUBLE_PRECISION, proc_num-20, 10*k+4, comm, MPI_STATUS_IGNORE,ierr)
                end do
            end do
        end if
        call mpi_barrier(comm, ierr)

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

    subroutine calculate_current()
        !Calculates the current on the segment interior, to be added on directly to the potential A for the heat equation test
        implicit none
        integer:: ccount !Corner count
        integer:: k, f
        real(num), dimension(:):: jh_data(1:nr+1)
        !Define the ordering of the triangles around the corner points
        !This doesn't have a reasonable formula, unfortunately
        do seg = seg_min, seg_max
            !Not corner points
            jr(:,:,seg) = 0.0_num
            do k = 0, npts-1
                do f = 0, 5
                jr(0:nr+1,k,seg) = jr(0:nr+1,k,seg) + bh(0:nr+1,mod(f/2 + 1,3),faces(f,k),seg)*all_point_ds(0:nr+1,f,k,seg)
                end do
            end do
            jh(:,:,:,seg) = 0.0_num
            do f = 0, 5
                do k = 0, npts-1
                jh_data(:) = 0.0_num
                jh_data(:) = bh(1:nr+1,mod(f/2 + 1,3),faces(f,k),seg)*all_point_ds(1:nr+1,f,k,seg)
                jh_data(:) = jh_data(:) - bh(0:nr,mod(f/2 + 1,3),faces(f,k),seg)*all_point_ds(0:nr,f,k,seg)
                jh_data(:) = jh_data(:) + br(1:nr+1,faces(f,k),seg)*(rc(1:nr+1) - rc(0:nr))
                jh_data(:) = jh_data(:) - br(1:nr+1,faces(mod(f+1,6),k),seg)*(rc(1:nr+1) - rc(0:nr))
                jh(1:nr+1,mod(f/2 + 1,3),faces(f,k),seg) = jh_data(:)
                end do
            end do
            do ccount = 0, 2
                k = corner_pts_1(ccount)
                jr(:,k,seg) = 0.0_num
                do f = 0, 4
                jr(0:nr+1,k,seg) = jr(0:nr+1,k,seg) + bh(0:nr+1,f_orders(ccount,f),faces(f,k),seg)*all_point_ds(0:nr+1,f,k,seg)

                end do
            end do
            do ccount = 0, 2
                k = corner_pts_1(ccount)
                do f = 0, 4
                if (faces(f,k) < ntri) then
                jh_data(:) = bh(1:nr+1,f_orders(ccount,f),faces(f,k),seg)*all_point_ds(1:nr+1,f,k,seg)
                jh_data(:) = jh_data(:) - bh(0:nr,f_orders(ccount,f),faces(f,k),seg)*all_point_ds(0:nr,f,k,seg)
                jh_data(:) = jh_data(:) + br(1:nr+1,faces(f,k),seg)*(rc(1:nr+1) - rc(0:nr))
                jh_data(:) = jh_data(:) - br(1:nr+1,faces(mod(f+1,5),k),seg)*(rc(1:nr+1) - rc(0:nr))
                jh(1:nr+1,f_orders(ccount,f),faces(f,k),seg) = jh_data(:)
                end if
                end do
            end do
            jr(0:nr+1,0:npts-1,seg) = jr(0:nr+1,0:npts-1,seg)/all_duals(0:nr+1,0:npts-1,seg)
            jh(1:nr+1,:,0:ntri-1,seg) = jh(1:nr+1,:,0:ntri-1,seg)/v_duals(1:nr+1,:,0:ntri-1,seg)
        end do
        return

    end subroutine calculate_current

    subroutine jfield_ghosts
        !Populates the horizontal current field ghost points, for use in the new averaging algorithm
        implicit none
        integer:: gi, i, seg_source, proc_source, seg_target
        integer:: seg_check, proc_target, id
        real(num), dimension(:):: jh_data(1:nr+1)
        !Transfer ghosts within the process.
        do seg = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg)
                proc_source = seg_source/seg_groups + seg_layer*20
                if (proc_source == proc_num) then !Just transfer data over if it's in the process
                    do i = 0, 2
                        jh(1:nr+1, i,ntri+gi, seg) = &
                        jh(1:nr+1,mod(i+tri_transfer(2,gi,seg), 3),  tri_transfer(1,gi,seg), seg_source)
                    end do
                end if
            end do
        end do
        !MPI sends
        do seg_source = seg_min, seg_max
        !Seg_source is the segment that is sending the information.
        !This is from process proc_num
            do seg_check = 0, 19
            !Run through all (other) segments and check if they need any information from here
            !Check ghost points of all other segments (only order n, luckily). A bit messy but meh.
            proc_target = seg_check/seg_groups + seg_layer*20
            if (proc_target .ne. proc_num) then !Do need MPI here
            do gi = 0, nghosttri-1
                if (tri_transfer(0,gi,seg_check) == seg_source) then
                    do i = 0, 2
                        jh_data = jh(1:nr+1, mod(tri_transfer(2,gi,seg_check) + i,3), tri_transfer(1,gi,seg_check), seg_source)
                        id = gi + 1000*seg_check + 100000*(i+1)
                        call mpi_send(jh_data, nr+1, MPI_DOUBLE_PRECISION, proc_target, id, comm, ierr)
                    end do
                end if
            end do
            end if
            end do
        end do
        call MPI_barrier(comm, ierr)
        !MPI receives
        do seg_target = seg_min, seg_max
            do gi = 0, nghosttri-1
                seg_source = tri_transfer(0,gi,seg_target)
                proc_source = seg_source/seg_groups + seg_layer*20
                if (proc_source .ne. proc_num) then
                do i = 0, 2
                    id = gi + 1000*seg_target + 100000*(i+1)
                    call mpi_recv(jh_data, nr+1, MPI_DOUBLE_PRECISION, proc_source, id, comm, MPI_STATUS_IGNORE, ierr)
                    jh(1:nr+1, i, ntri+gi, seg_target) = jh_data
                end do
                end if
            end do
        end do
        call MPI_barrier(comm, ierr)

        if (proc_num == -1) then
        print*, ''
        print*, 'jmaxes'
        print*, maxval(abs(jr(1:nr+1,0:npts-1,:)))
        print*, maxval(abs(jh(1:nr,:,0:ntri+nghosttri-1,:)))
        end if

        return
    end subroutine jfield_ghosts

    subroutine j_to_gridpts()
        !Averages j to gridpoints using the new approach. Fingers crossed...

        implicit none
        real(num), dimension(:):: rhs(0:2)
        integer:: seg, k, j, ri

        j1 = 0.0_num
        !Establish matrix of normal vectors for each point. Require 12 of these (in general)
        do seg = seg_min,seg_max
        do k = 0, npts-1
        !Find horizontal tangent vectors
        !Multiply by magnetic field to obtain the rhs that needs to be solved
        do ri = 1, nr+1  !Is now radial variation so need to do this
        rhs = 0.0_num
        do j = 0, nnbrs(k)-1
            rhs(:) = rhs(:) + jvec(:,j,k,seg)*(jh(ri,trimap(1,j,k),trimap(0,j,k),seg))!Horizontal vectors
        end do
        rhs(:) = rhs(:) + 0.5*jvec(:,nnbrs(k),k,seg)*(jr(ri-1,k,seg) + jr(ri,k,seg))  !Vertical vectors
        j1(0:2,ri,k,seg) = solve3d(jmat(:,:,k,seg), rhs)
        end do
        end do
        end do

        return

    end subroutine j_to_gridpts

    subroutine b_to_gridpts()
        !Does the main b field but also the horizontal upwinded one for the outflow term
        implicit none
        real(num), dimension(:):: rhs(0:2)
        integer:: seg, k, j, ri

        b1 = 0.0_num
        !Establish matrix of normal vectors for each point. Require 12 of these (in general)
        do seg = seg_min,seg_max
        do k = 0, npts-1
        !Multiply by magnetic field to obtain the rhs that needs to be solved
        do ri = 1, nr+1  !Is now radial variation so need to do this
        rhs = 0.0_num
        do j = 0, nnbrs(k)-1
            rhs(:) = rhs(:) + 0.5*bvec(:,j,k,seg)*(bh(ri-1,trimap(1,j,k),trimap(0,j,k),seg) + &
            bh(ri,trimap(1,j,k),trimap(0,j,k),seg)) !Horizontal vectors
            rhs(:) = rhs(:) + bvec(:,j+nnbrs(k),k,seg)*br(ri,trimap(0,j,k),seg)  !Vertical vectors
        end do
        !Solve for this system at all the grid points
        b1(0:2,ri,k,seg) = solve3d(bmat(:,:,k,seg), rhs)
        end do

        end do
        end do

        return

    end subroutine b_to_gridpts

    subroutine electric_integral()
        !Converts the electric integral (averaged to gridpoints) into line integrals
        !Currently set up just to be an integral of the current to test if it works
        implicit none

        integer:: k, i, j
        real(num):: eh_temp(1:nr+1)

        er = 0.0; eh = 0.0
        do seg = seg_min, seg_max
            !Radial direction first - electric field saved on grid points so this is simple
            do k = 0, npts-1
                do i = 0, 2
                er(1:nr,k,seg) = er(1:nr,k,seg) + &
                0.5*(e1(i,1:nr,k,seg)+e1(i,2:nr+1,k,seg))*(rs(2:nr+1)-rs(1:nr))*all_base_pts(i,k,seg)
                end do
            end do
            !Non-radial direction. No vertical averaging needed. Remove radial component
            do k = 0, ntri-1
            do j =  0, 2 !Triangle side

            do i = 0, 2 !Coordinate direction
            !Average field
            !print*, tri_l_faces(1:nr+1,i,k,seg)
            eh_temp = 0.5*(e1(i,1:nr+1,triangs(j,k),seg)+e1(i,1:nr+1,triangs(mod(j+1,3),k),seg))
            eh(1:nr+1,j,k,seg) = eh(1:nr+1,j,k,seg) + eh_temp*trisides(i,j,k,seg)*tri_l_faces(1:nr+1,j,k,seg)
            end do
            end do
            end do
        end do
        !er = 0.0_num
        return

    end subroutine electric_integral

    subroutine add_diffusion()
        implicit none

        integer:: k
        do seg = seg_min, seg_max
        do k = 0, ntri-1
            eh(1:nr+1,:,k,seg) = eh(1:nr+1,:,k,seg) + 1.0*eta*jh(1:nr+1,:,k,seg)*tri_l_faces(1:nr+1,:,k,seg)
        end do
        do k = 0, npts-1
        er(0:nr+1,k,seg) = er(0:nr+1,k,seg) + 1.0*eta*jr(0:nr+1,k,seg)*(rs(1:nr+2)-rs(0:nr+1))
        end do
        end do

        return
    end subroutine add_diffusion

    subroutine calculate_electric()
        !This is the base electric field - the friction crossed with the magnetic field
        !To be modified to use an upwinded magnetic field instead, as that may be more interesting
        e1 = 0.0*eta*j1
        e1 = e1 - cross2(v1, b1)

        return

    end subroutine

    subroutine calculate_friction()
        !Calculates the magnetofrictional term with a couple of cross products. Not too tricky
        implicit none
        real(num), dimension(:,:,:,:):: lf(0:2,0:nr+2,0:npts-1,seg_min:seg_max)
        real(num), dimension(:,:,:)  :: soft(0:nr+2,0:npts-1,seg_min:seg_max)
        real(num), dimension(:,:,:)  :: b2(0:nr+2,0:npts-1,seg_min:seg_max)
        integer:: seg, k, i
        real(num):: eps_rel

        lf = cross2(j1, b1)  !Top of the magnetofrictional term (The Lorentz Force)

        if (seg_layer == 0) lf(:,0:1,:,seg_min:seg_max) = 0.0_num ! Line-tied at the surface

        do seg = seg_min, seg_max
        b2(:,:,seg) = 0.0_num
        do k = 0, npts-1
        b2(1:nr+1,k,seg) = sum(b1(:,1:nr+1,k,seg)**2, 1)
        end do
        end do
        eps_rel = eps*maxb_init**2

        soft = b2 + eps_rel*exp(-b2/eps_rel) !Softening term
        do i = 0, 2
        v1(i,:,:,:) = fric_field*lf(i,:,:,:)/soft(:,:,:)
        end do

        if (seg_layer == 0) v1(:,0:1,:,seg_min:seg_max) = 0.0_num ! Line-tied at the surface. To be sure.

        return

    end subroutine calculate_friction


    subroutine add_shearing()
        !Adds the shearing on the lower boundary. Add to e1 vector I think, no need to be smart
        implicit none
        real(num), dimension(:,:,:):: shear_v(0:2,0:npts-1,seg_min:seg_max)
        integer:: k, seg
        real(num):: lat
        real(num), dimension(:):: rot_dir(0:2), loc(0:2)  !Rotation direction


        if (seg_layer == 0) then  !Only do this if it's on the bottom

        do seg = seg_min, seg_max
        do k = 0, npts-1
            loc = all_base_pts(:,k,seg)
            if (abs(loc(2)) > 1.0-1e-6) then
                rot_dir = (/0.0_num,0.0_num,0.0_num/)
            else
                lat = atan(loc(2)/sqrt((loc(0)**2 + loc(1)**2)))
                rot_dir = cross1((/0.0_num,0.0_num,1.0_num/), (/loc(0), loc(1), loc(2)/))
                rot_dir = rot_dir/sqrt(sum(rot_dir**2))
            end if
            !Latidude (in radians from the equator)
            shear_v(:,k,seg) = cos(lat)*shear_fn(lat)*rot_dir

        end do
        end do
        v1(:,0,0:npts-1,seg_min:seg_max) = shear_v
        v1(:,1,0:npts-1,seg_min:seg_max) = shear_v
        end if
        return
    end subroutine add_shearing

    subroutine make_vout_field()
        !Calculates the matrix to be multiplied with the magnetic field in the outflow term.
        !Only do this once (at the start). Variable array allocated in the usual place
        integer:: ri
        vout_field = 0.0_num
        do ri = 1, nr+1
        vout_field(ri) = vout_fn(rs(ri))
        end do

        if (seg_layer == 0) then
            vout_field(0:2) = 0.0_num
        end if

    end subroutine make_vout_field

    subroutine make_fric_field()
        !Necessary if the magnetofriction relaxation rate depends on latitude or radius. If not will just be a constant matrix
        integer:: seg, ri, k
        real(num):: lat, s
        real(num), dimension(:):: loc(0:2)
        fric_field = 0.0_num
        do seg = seg_min, seg_max
        do ri = 1, nr+1
        do k = 0, npts-1
        loc = all_base_pts(:,k,seg)
        lat = atan(loc(2)/sqrt((loc(0)**2 + loc(1)**2)))
        s = sin(pi/2-abs(lat))
        if (constant_friction) then
            fric_field(ri,k,seg) = nu0
        else
            fric_field(ri,k,seg) = nu0*s**2*rs(ri)**2
        end if
        end do
        end do
        end do

    end subroutine make_fric_field

    subroutine add_outflow()
        !Needs upwinding from the magnetic field, so need a special cross product with the averaged magnetic field.
        !As the vout field is vertical, only the horizontal component of the magnetic field is important. It is advantageous to use from half a cell downwards as a result.
        !This should also be fine with MPI, but not certain about that
        implicit none
        integer:: j,k
        real(num):: eh_temp(1:nr+1)

        do seg = seg_min, seg_max
        do k = 0, ntri-1
        do j = 0,2
        eh_temp = -bh(0:nr,j,k,seg)*vout_field(1:nr+1)
        eh(1:nr+1,j,k,seg) = eh(1:nr+1,j,k,seg) + eh_temp*tri_l_faces(1:nr+1,j,k,seg)
        end do
        end do
        end do

        return

    end subroutine add_outflow

    subroutine add_surface_diffusion()
        !Adds the surface diffusion by creating a fake current on the surface
        implicit none
        real(num), dimension(:,:,:):: diff(0:2,0:ntri-1,seg_min:seg_max)
        real(num), dimension(:,:,:):: jh_surf(0:2,0:ntri+nghosttri-1,seg_min:seg_max)
        real(num):: jh_data

        integer:: i, f, k, seg, ccount
        if (seg_layer == 0) then  !Only do this if it's on the bottom

        do seg = seg_min, seg_max
        do f = 0, 5
        do k = 0, npts-1
        !Do corner points (incorrectly) and then overwrite them if necessary
        jh_data = br(1,faces(f,k),seg)*(rc(1) - rc(0))
        jh_data = jh_data - br(1,faces(mod(f+1,6),k),seg)*(rc(1) - rc(0))
        jh_surf(mod(f/2+1,3),faces(f,k),seg) = jh_data
        end do
        end do

        do ccount = 0, 2
            k = corner_pts_1(ccount)
            do f = 0, 4
            jh_data = br(1,faces(f,k),seg)*(rc(1) - rc(0))
            jh_data = jh_data - br(1,faces(mod(f+1,5),k),seg)*(rc(1) - rc(0))
            jh_surf(f_orders(ccount,f),faces(f,k),seg) = jh_data
            end do
        end do

        jh_surf(:,0:ntri-1,seg) = jh_surf(:,0:ntri-1,seg)/v_duals(1,:,0:ntri-1,seg)

        do i = 0, 2  !Integrate to get into line integral form (like the electric field)
            diff(i,:,seg) = jh_surf(i,0:ntri-1,seg)*tri_l_faces(1,i,0:ntri-1,seg)
        end do

        if (seg_layer == 0) then
            eh(1,:,0:ntri-1,seg) = eh(1,:,0:ntri-1,seg) + eta0*diff(:,0:ntri-1,seg)
        end if
        end do
        end if
        return


    end subroutine add_surface_diffusion

    subroutine check_divfree_interior()
        !Check locally divergeence-free. Don't need to do this every time I suppose but it doesn't take long.
        !If cfl is too high this tends to detect it fairly rapidly
        implicit none
        integer:: i, k
        real(num), dimension(:,:,:):: divergence(1:nr ,0:ntri-1,0:19)  !Divergence of internal raw cells.
        real(num), dimension(:,:,:):: br_flux(0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max)
        real(num), dimension(:,:,:):: bh_flux(0:nr+1,0:2,0:ntri+nghosttri-1,seg_min:seg_max)
        divergence = 0.0_d
        do seg = seg_min, seg_max
            br_flux(:,:,seg) = br(:,:,seg)*h_areas(:,:,seg)
            bh_flux(:,:,:,seg) = bh(:,:,:,seg)*v_areas(:,:,:,seg)
            !print*,'lbound flux', sum(br_flux(1,:,:))
            divergence(:,:,seg) = divergence(:,:,seg) - br_flux(1:nr ,0:ntri-1,seg)  !In at the bottom
            divergence(:,:,seg) = divergence(:,:,seg) + br_flux(2:nr +1,0:ntri-1,seg)  !out at the top
            do k = 0, ntri-1
                do i = 0, 2
                    divergence(1:nr ,k,seg) = divergence(1:nr ,k,seg) + bh_flux(1:nr,i,k,seg)  !Through the sides
                end do
            end do
            if (maxval(abs(divergence(:,:,seg))) > 1d-14) print*, 'some cells not divergence free'
            if (maxval(abs(divergence(:,:,seg))) > 1d-14) print*, maxval(abs(divergence(:,:,seg)))
            if (maxval(abs(divergence(:,:,seg))) > 1d-14) call mpi_abort(comm, ierr)
        end do

        return

    end subroutine check_divfree_interior


END MODULE ops
