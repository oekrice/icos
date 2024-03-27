!*******************************************************************************
MODULE import_export
    !*******************************************************************************
    ! Imports a potential (or otherwise) field, integrates to find the vector potential A and then distributes this as appropriate to all the processes
    ! Calculates some basic diagnostics (very crudely) and exports these to /diagnostics
    ! Collects data from individual processess and then saves to netcdf file as appropriate
    ! TO CHANGE THE LENGTH OF DATA DIRECTORY DO SO IN 'export_with_checks'
    !*******************************************************************************
    USE shared_data
    USE netcdf

    IMPLICIT NONE

    INTEGER, PARAMETER:: nfreal = NF90_DOUBLE  ! type of output variables

    contains
    subroutine import_initial_condition()
        implicit none
        if (proc_num == 0) then   !Do all the imports in one process
            call read_initial_file()  !Read in initial condition (as magnetic fluxes)
            print*, 'Initial condition read in...'
            call check_magfield()   !Check the imported field is divergence-free,
            print*, 'Magnetic field is reasonably divergence-free...'
            call calculate_vector_potential()

            call check_mag_current() !Check the vector potential produces the desired initial condition and prints the error
        end if
        allocate(ar(0:nr +1,0:npts+nghosts-1,seg_min:seg_max))   !Local variables for the vector potential. Fortran is nice so can segment the end bits in a nice way
        allocate(ah(0:nr +2,0:2,0:ntri+nghosttri-1,seg_min:seg_max))

        call distribute_potential()  !Distribute the calculated potential field among the processes

        return
    end subroutine

    subroutine reset_afield()
        implicit none
        if (proc_num == 0) then   !Do all the imports in one process
            call read_bfield()  !Read in initial condition (as magnetic fluxes)
            call check_magfield()   !Check the imported field is divergence-free,
            call calculate_vector_potential()

            call check_mag_current() !Check the vector potential produces the desired initial condition and prints the error
        end if

        call distribute_potential()  !Distribute the calculated potential field among the processes

        return
    end subroutine reset_afield


    subroutine try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
        call mpi_abort(comm, ierr)
    end if

    end subroutine try

    subroutine distribute_potential()
        !Once the vector potential has been calculated, distribute it among the processes. This is really quite non-trivial...
        implicit none

        real(num), allocatable:: sum_test(:)
        real(num):: sum_global
        integer:: i
        !real(num):: test2

        ar = 0.0_num; ah = 0.0_num
        allocate(sum_test(0:nprocs-1))
        if (proc_num > 0) then
            allocate(ar_global(0:nr_global +1,0:npts-1,0:19))   !Just interior points. These should be deallocated once the information is passed to individual processes
            allocate(ah_global(0:nr_global +2,0:2,0:ntri-1,0:19))
        end if
        call MPI_Bcast(ar_global, (nr_global+2)*npts*20, MPI_DOUBLE_PRECISION, 0, comm, ierr)
        call MPI_Bcast(ah_global, (nr_global+3)*ntri*20*3, MPI_DOUBLE_PRECISION, 0, comm, ierr)

        sum_global = 0.0_num
        if (proc_num == 0) then
            sum_global = sum(abs(ah_global(2:nr_global,:,:,:))) + 0.5*sum(abs(ah_global(1,:,:,:))) + &
            0.5*sum(abs(ah_global(nr_global+1,:,:,:)))
        end if

        if (nprocs < 2) then
            ar(0:nr_global +1,0:npts-1,0:19) = ar_global
            ah(0:nr_global +2,0:2,0:ntri-1,0:19) = ah_global
            sum_test(proc_num) =  sum(abs(ah(2:nr,:,:,:))) + 0.5*sum(abs(ah(1,:,:,:))) + 0.5*sum(abs(ah(nr+1,:,:,:)))

        else if (nprocs < 21) then  !one layer, shared segments
            ar(0:nr_global +1,0:npts-1,seg_min:seg_max) = ar_global(0:nr_global +1,0:npts-1,seg_min:seg_max)
            ah(0:nr_global +2,0:2,0:ntri-1,seg_min:seg_max) = ah_global(0:nr_global +2,0:2,0:ntri-1,seg_min:seg_max)
            sum_test(proc_num) = sum(abs(ah(2:nr,:,:,:))) + 0.5*sum(abs(ah(1,:,:,:))) + 0.5*sum(abs(ah(nr+1,:,:,:)))

        else   !Multiple layers, single segment per process. Careful about nr vs nr_global here
            ar(0:nr +1,0:npts-1,seg_min:seg_max) = ar_global(seg_layer*nr:seg_layer*nr+nr+1,0:npts-1,seg_min:seg_max)
            ah(0:nr +2,:,0:ntri-1,seg_min:seg_max) = ah_global(seg_layer*nr:seg_layer*nr+nr+2,:,0:ntri-1,seg_min:seg_max)

            sum_test(proc_num) = sum(abs(ah(2:nr,:,:,:))) + 0.5*sum(abs(ah(1,:,:,:))) + 0.5*sum(abs(ah(nr+1,:,:,:)))
        end if

        if (proc_num > 0) call mpi_send(sum_test(proc_num), 1, MPI_DOUBLE_PRECISION, 0, 4, comm, ierr)
        call MPI_Barrier(comm,ierr)  !Wait for t to be broadcast everywhere.


        if (proc_num == 0)  then
            if (nprocs > 1) then
                do i = 1, nprocs-1
                    call mpi_recv(sum_test(i), 1, MPI_DOUBLE_PRECISION, i, 4, comm, MPI_STATUS_IGNORE, ierr)
                end do
            end if
            if (abs((sum(sum_test) - sum_global)/sum_global) > 1d-10) then
                print*, 'Initial information not distributed properly. Check grid is acceptably divided'
                print*, 'Absolute error:', abs(sum(sum_test) - sum_global), sum_global, sum(sum_test)
                call mpi_abort(comm, ierr)
            end if
            if (step == 0)  print*, 'Initial condition distributed to processes seemingly without problems'
            print*, '__________________________________________________________'

        end if

        do i = 1, nprocs-1
            if (proc_num == 0) call mpi_send(maxb_init, 1, MPI_DOUBLE_PRECISION,i,5,comm,ierr)
        end do
        call MPI_Barrier(comm,ierr)  !Wait for t to be broadcast everywhere.
        do i = 1, nprocs-1
            if (proc_num == i) call mpi_recv(maxb_init, 1, MPI_DOUBLE_PRECISION,0,5,comm,MPI_STATUS_IGNORE,ierr)
        end do

        !Clean up
        deallocate(ah_global); deallocate(ar_global)

        return
    end subroutine distribute_potential

    subroutine calculate_vector_potential()
        !Given the horizontal and vertical magnetic fluxes, find an appropriate vector potential.
        !This only works if the imported field is roughly divergence-free
        !If the initial field is not quite divergence-free, then the resulting current will be altered slightly and any symmetry will be lost
        implicit none
        real(num):: diff
        integer, dimension(:,:,:):: fixed(0:2,0:ntri-1,0:19)  !Records the fixed integrals
        integer, dimension(:,:,:):: fixed_sums(0:ntri-1,0:19)
        integer, dimension(:):: seg_list(0:19), othertri(0:1), loc(0:2)
        integer:: seg_count, seg, i, i2, k, ri
        integer:: dof, nbrind, seg_side
        logical:: carry_on
        allocate(ar_global(0:nr_global +1,0:npts-1,0:19))   !Just interior points. These should be deallocated once the information is passed to individual processes
        allocate(ah_global(0:nr_global +2,0:2,0:ntri-1,0:19))


        seg_list = (/0,1,2,3,4,10,11,12,13,14,15,16,17,18,19,5,6,7,8,9/)
        fixed = 0; fixed_sums = 0; carry_on = .true.; nbrind = -1
        ar_global = 0.0; ah_global = 0.0

        if (step == 0) print*, 'Integrating horizontally at r = 1...'
        do seg_count = 0, 19
            seg = seg_list(seg_count)
            carry_on = .true.; fixed_sums = 0
            do i = 0, 2
                fixed_sums(:,seg) = fixed_sums(:,seg) + fixed(i,:,seg)
            end do

            do while (carry_on)
                dof = 0
                do k = 0, ntri-1
                    if (fixed_sums(k, seg) < 3 .and. fixed_sums(k, seg) > dof) dof = fixed_sums(k, seg)
                end do
                do k = 0, ntri-1
                    !call printmat(dble(fixed(:,:,seg)), 3, ntri)
                    !Do the triangles where the degree of freedom is minimum
                    if (fixed_sums(k,seg) == dof) then
                        diff = br_global(1,k,seg)*h_areas_global(1,k,seg) - sum(ah_global(1,:,k,seg))
                        do i = 0, 2
                            if (fixed(i,k,seg) == 0) then  !Can change this side, do so.

                                ah_global(1,i,k,seg) = diff
                                fixed(:,k,seg) = 1

                                if (updown(k) == 1) then  !Upwards triangle. Determine neighbour side that needs setting
                                    othertri = (/nbrfaces(i,k), mod(i+2, 3)/)   !The neigbouring face and respective side
                                    if (othertri(0) < ntri) then  !Not a ghost point (within the segment)
                                        ah_global(1,othertri(1),othertri(0),seg) = -diff
                                    else  !Is a ghost point. Need to find the location in the next segment
                                        loc = tri_transfer(:,othertri(0) - ntri,seg)  !Ghost point location
                                        !Find the index of the neighbour
                                        do i2 = 0, 2
                                            if (ico_nbrs(seg,i2,0) == loc(0)) nbrind = i2
                                        end do
                                        seg_side = ico_nbrs(seg, nbrind, 1)
                                        ah_global(1,seg_side,loc(1),loc(0)) = -diff
                                    end if
                                else   !Downwards-facing triangle (shouldn't have a ghost point...)
                                    othertri = (/nbrfaces(i,k), mod(i+1, 3)/)
                                    ah_global(1,othertri(1),othertri(0),seg) = -diff
                                end if
                            exit
                            end if
                        end do
                        do i = 0, 2  !Fix ALL affected edges here
                            if (updown(k) == 1) then  !Upwards triangle. Determine neighbour side that needs setting
                                othertri = (/nbrfaces(i,k), mod(i+2, 3)/)   !The neigbouring face and respective side
                                if (othertri(0) < ntri) then  !Not a ghost point (within the segment)
                                    fixed(othertri(1), othertri(0), seg) = 1  !might need to do this elsewhere...
                                else  !Is a ghost point. Need to find the location in the next segment
                                    loc = tri_transfer(:,othertri(0) - ntri,seg)  !Ghost point location
                                    do i2 = 0, 2
                                        if (ico_nbrs(seg,i2,0) == loc(0)) nbrind = i2
                                    end do
                                    seg_side = ico_nbrs(seg, nbrind, 1)
                                    fixed(seg_side,loc(1),loc(0)) = 1
                                end if
                            else   !Downwards-facing triangle (shouldn't have a ghost point...)
                                othertri = (/nbrfaces(i,k), mod(i+1, 3)/)
                                fixed(othertri(1), othertri(0), seg) = 1  !might need to do this elsewhere...
                            end if
                        end do
                    end if
                    !At this point the segment in question should have been sorted, so check the flux
                end do
                fixed_sums = 0
                do i = 0, 2
                    fixed_sums(:,seg) = fixed_sums(:,seg) + fixed(i,:,seg)
                end do
                if (minval(fixed_sums(:,seg)) == 3) carry_on = .false.
            end do
        end do
        do seg = 0, 19
            do k = 0, ntri-1
                if (br_global(1,k,seg)*h_areas_global(1,k,seg) - sum(ah_global(1,:,k,seg)) > 1e-10) &
                print*, 'Integration problem. Perhaps not div-free'
            end do
        end do
        if (step == 0) print*, 'Integrating radially...'
        !Vertical ar is zero (initially) in the De Vore/Coulomb gauge. Update other horizontal cells based on those below and horizontal fluxes.
        do seg = 0, 19
            do ri = 2, nr_global  + 1
                do k = 0, ntri-1
                    do i = 0, 2
                    ah_global(ri, i, k, seg) = -bh_global(ri-1, i, k, seg)*v_areas_global(ri-1,i,k,seg) + ah_global(ri-1, i, k, seg)
                    end do
                end do
            end do
        end do
        if (step == 0) print*, 'Vector potential calculated in main process.'

        return

    end subroutine calculate_vector_potential

    subroutine check_mag_current()
        !Additional version of the magnetic field calculator. Prints the error compared to the initial condition once it has been cleaned of divergenceness.
        implicit none

        integer:: seg, k, i
        real(num), allocatable:: br_test(:,:,:), bh_test(:,:,:,:)
        real(num):: error
        !real(num), allocatable:: jr_test(:,:,:,:), jh_test(:,:,:)
        !Check that the magnetic field resulting from this is adequately similar to the inputted one. Doesn't take long and can reuse the code (one hopes)
        allocate(br_test(0:nr_global +2,0:ntri+nghosttri-1,0:19))   !Just interior points. These should be deallocated once the information is passed to individual processes
        allocate(bh_test(0:nr_global +1,0:2,0:ntri+nghosttri-1,0:19))

        br_test = 0.0; bh_test = 0.0
        do seg = 0, 19

            do k = 0, ntri-1
                do i = 0, 2
                    br_test(:,k,seg) = br_test(:,k,seg) + ah_global(:,i,k,seg)

                    bh_test(1:nr_global ,i,k,seg) = bh_test(1:nr_global ,i,k,seg) - ah_global(2:nr_global +1,i,k,seg)
                    bh_test(1:nr_global ,i,k,seg) = bh_test(1:nr_global ,i,k,seg) + ah_global(1:nr_global   ,i,k,seg)

                    bh_test(1:nr_global ,i,k,seg) = bh_test(1:nr_global ,i,k,seg) - ar_global(1:nr_global,triangs(i,k),seg)
                    bh_test(1:nr_global ,i,k,seg) = bh_test(1:nr_global ,i,k,seg) + ar_global(1:nr_global,triangs(mod(i+1,3),k),seg)
                end do
            end do

        end do
        error = max(maxval(abs(br_global(1:nr_global +1,0:ntri-1,:)*h_areas_global(1:nr_global+1,0:ntri-1,:) - &
        br_test(1:nr_global +1,0:ntri-1,:))), maxval(abs(bh_global(1:nr_global ,:,0:ntri-1,:)*&
        v_areas_global(1:nr_global ,:,0:ntri-1,:) - bh_test(1:nr_global ,:,0:ntri-1,:))))
        print*, 'Max absolute difference from specified initial condition = ', error

        deallocate(bh_test)
        deallocate(br_test)

        return

    end subroutine check_mag_current

    subroutine check_magfield()
    ! Check the magnetic field is sufficiently divergence-free at this point in time
    ! Checks the overall flux through each layer as well as all the cells
        implicit none
        integer:: ri, i, k
        real(num), dimension(:,:,:):: divergence(1:nr_global ,0:ntri-1,0:19)  !Divergence of internal raw cells
        do ri = 1, nr_global +1     !Running through radial layers
            if (sum(br_global(ri,0:ntri-1,:)*h_areas_global(ri,0:ntri-1, 0:19)) > 1d-10) print*, &
            'layer', ri, 'not divergence-free. Check lower boundary condition?'
            !print*, ri, sum(br_global(ri,0:ntri-1,:))
        end do
        divergence = 0.0_d
        divergence = divergence - br_global(1:nr_global ,0:ntri-1,:)*h_areas_global(1:nr_global,0:ntri-1, 0:19)  !In at the bottom
        divergence = divergence + br_global(2:nr_global +1,0:ntri-1,:)*h_areas_global(2:nr_global+1,0:ntri-1, 0:19)  !out at the top

        do k = 0, ntri-1
            do i = 0, 2
                divergence(1:nr_global ,k,:) = divergence(1:nr_global ,k,:) + &
                bh_global(1:nr_global ,i,k,:)*v_areas_global(1:nr_global,i,k, 0:19)  !Through the sides
            end do
        end do

        if (maxval(abs(divergence)) > 1d-10) print*, 'some cells not divergence free', maxval(abs(divergence))
        if (maxval(abs(divergence)) > 1d-10) call mpi_abort(comm, ierr)
        return

    end subroutine check_magfield

    subroutine read_initial_file()
        !Read in the initial condition from a netcdf file
        implicit none
        integer:: ncid
        integer:: vid
        character(len=23):: filename
        allocate(br_import(1:nr_global +1,1:ntri,1:20))   !No extras here - can add on later on if necessary. These are AVERAGE FIELD
        allocate(bh_import(1:nr_global ,1:3,1:ntri,1:20))

        filename = init_filename

        call try(nf90_open(filename, nf90_nowrite, ncid))  !Have to be VERY careful with the axes here, as they are stupid and don't make sense
        call try(nf90_inq_varid(ncid, 'br', vid))
        call try(nf90_get_var(ncid, vid, br_import))

        call try(nf90_inq_varid(ncid, 'bh', vid))
        call try(nf90_get_var(ncid, vid, bh_import))

        call try(nf90_close(ncid))

        allocate(br_global(0:nr_global +2,0:ntri+nghosttri-1,0:19))   !No extras here - can add on later on if necessary
        allocate(bh_global(0:nr_global +1,0:2,0:ntri+nghosttri-1,0:19))

        br_global(1:nr_global +1,0:ntri-1, 0:19) = br_import
        bh_global(1:nr_global ,:,0:ntri-1, 0:19) = bh_import

        maxb_init = sum(abs(br_global(1:nr+1,:,:)))/size(br_global(1:nr+1,:,:))

        deallocate(bh_import)
        deallocate(br_import)
        return

    end subroutine read_initial_file

    subroutine read_bfield()
        !Read in an already-saved file
        implicit none
        integer:: ncid
        integer:: vid
        character(len=100):: filename
        allocate(br_import(1:nr_global +1,1:ntri,1:20))   !No extras here - can add on later on if necessary. These are the field strength, not fluxes
        allocate(bh_import(1:nr_global ,1:3,1:ntri,1:20))

        filename = trim(output_filename)

        call try(nf90_open(filename, nf90_nowrite, ncid))  !Have to be VERY careful with the axes here, as they are stupid and don't make sense
        call try(nf90_inq_varid(ncid, 'br', vid))
        call try(nf90_get_var(ncid, vid, br_import))

        call try(nf90_inq_varid(ncid, 'bh', vid))
        call try(nf90_get_var(ncid, vid, bh_import))

        call try(nf90_close(ncid))

        br_global(1:nr_global +1,0:ntri-1, 0:19) = br_import
        bh_global(1:nr_global ,:,0:ntri-1, 0:19) = bh_import

        !print*, 'Top boundary flux (in read_bfield)', sum(br_global(nr_global,0:ntri-1,0:19)*h_areas(1,0:ntri-1, 0:19))
        if (abs(sum(br_global(nr_global,0:ntri-1,0:19)*h_areas_global(nr_global,0:ntri-1, 0:19))) > 1e-10) then
        print*, 'Top boundary flux non-zero', sum(br_global(nr_global,0:ntri-1,0:19)*h_areas(nr_global,0:ntri-1, 0:19))
        call mpi_abort(comm,ierr)
        end if

        print*, 'Global Open Flux', sum(abs(br_global(nr_global,0:ntri-1,0:19))*h_areas_global(nr_global,0:ntri-1, 0:19))

        maxb_init = sum(abs(br_global(1:nr+1,:,:)))/size(br_global(1:nr+1,:,:))

        deallocate(bh_import)
        deallocate(br_import)
        return

    end subroutine read_bfield


    subroutine export_with_checks()
        !Exports the magnetic field, but as this sometimes can go screwy due to a netcdf bug, tries again until it is filename
        implicit none
        logical:: export_is_fine
        integer:: ncid, vid

        real(num), dimension(:,:,:):: br_test(1:nr_global +1,1:ntri,1:20)
        real(num), dimension(:,:,:,:):: bh_test(1:nr_global ,1:3,1:ntri,1:20)

        if (hflag < 0.5) then
            if (plot_num < 10) then
                write (output_filename, "(A31,A3,I1,A3)") trim(data_directory), "000", plot_num, ".nc"
            else if (plot_num < 100) then
                write (output_filename, "(A31,A2,I2,A3)") trim(data_directory), "00", plot_num, ".nc"
            else if (plot_num < 1000) then
                write (output_filename, "(A31,A1,I3,A3)") trim(data_directory), "0", plot_num, ".nc"
            else if (plot_num < 10000) then
                write (output_filename, "(A31,I4,A3)") trim(data_directory), plot_num, ".nc"
            end if
        else
            if (plot_num < 10) then
                write (output_filename, "(A30,A3,I1,A3)") trim(data_directory), "000", plot_num, ".nc"
            else if (plot_num < 100) then
                write (output_filename, "(A30,A2,I2,A3)") trim(data_directory), "00", plot_num, ".nc"
            else if (plot_num < 1000) then
                write (output_filename, "(A30,A1,I3,A3)") trim(data_directory), "0", plot_num, ".nc"
            else if (plot_num < 10000) then
                write (output_filename, "(A30,I4,A3)") trim(data_directory), plot_num, ".nc"
             end if
        end if

        export_is_fine = .false.

        do while (.not. export_is_fine)
            !Try to export
            call export_magnetic_field
            !Read back in to check numbers are resonable

            !Read back in to check
            call try(nf90_open(output_filename, nf90_nowrite, ncid))  !Have to be VERY careful with the axes here, as they are stupid and don't make sense
            call try(nf90_inq_varid(ncid, 'br', vid))
            call try(nf90_get_var(ncid, vid, br_test))

            call try(nf90_inq_varid(ncid, 'bh', vid))
            call try(nf90_get_var(ncid, vid, bh_test))

            call try(nf90_close(ncid))
            call mpi_barrier(comm, ierr)

            if ( (maxval(abs(br_test)) < 1e10) .and. (maxval(abs(bh_test)) < 1e10) ) then
                export_is_fine = .true.
            else
                if (proc_num == 0) print*, 'Export number', plot_num, 'failed, trying again'
            end if
            call mpi_barrier(comm, ierr)

        end do

        !if (proc_num == 0) print*, 'Exported magnetic field number', plot_num, 'at time', (step-1)*dt
        return

    end subroutine export_with_checks

    subroutine export_magnetic_field()
        !Exports the magnetic field at this plot_num to an appropriate netcdf file
        implicit none
        integer:: ncid
        integer:: rs_id, rc_id, tri_id, side_id, seg_id
        integer:: br_id, bh_id
        integer:: proc_test

        if (proc_num == 0) then
        call try(nf90_create(trim(output_filename), nf90_clobber, ncid))
        !Define variables (in main process)
        call try(nf90_def_dim(ncid, 'rs', nr_global+1, rs_id))  !Make up fake dimensions here
        call try(nf90_def_dim(ncid, 'rc', nr_global, rc_id))  !Make up fake dimensions here
        call try(nf90_def_dim(ncid, 'tri', ntri, tri_id))  !Make up fake dimensions here

        call try(nf90_def_dim(ncid, 'nsides', 3, side_id))  !Make up fake dimensions here
        call try(nf90_def_dim(ncid, 'nsegs', 20, seg_id))  !Make up fake dimensions here

        call try(nf90_def_var(ncid, 'br', nf90_double, (/rs_id,tri_id,seg_id/), br_id))
        call try(nf90_def_var(ncid, 'bh', nf90_double, (/rc_id,side_id,tri_id,seg_id/), bh_id))
        call try(nf90_enddef(ncid))
        call try(nf90_close(ncid))


        end if
        call mpi_barrier(comm, ierr)

        !Each process writes data to the file in turn

        do proc_test = 0, nprocs-1
        call mpi_barrier(comm, ierr)

        if (proc_num == proc_test) then
            call try(nf90_open(trim(output_filename), nf90_write, ncid))
            call try(nf90_inq_varid(ncid, 'br', br_id))
            call try(nf90_inq_varid(ncid, 'bh', bh_id))

            if (seg_layer == 0) then
            
            call try(nf90_put_var(ncid, br_id, br(1:nr+1,0:ntri-1,seg_min:seg_max), &
            start = (/seg_layer*nr+1,1,seg_min+1/), count = (/nr+1,ntri,seg_groups/)))

            call try(nf90_put_var(ncid, bh_id, bh(1:nr,0:2,0:ntri-1,seg_min:seg_max), &
            start = (/seg_layer*nr+1,1,1,seg_min+1/), count = (/nr,3,ntri,seg_groups/)))

            else

            call try(nf90_put_var(ncid, br_id, br(2:nr+1,0:ntri-1,seg_min:seg_max), &
            start = (/seg_layer*nr+2,1,seg_min+1/), count = (/nr,ntri,seg_groups/)))

            call try(nf90_put_var(ncid, bh_id, bh(1:nr,0:2,0:ntri-1,seg_min:seg_max), &
            start = (/seg_layer*nr+1,1,1,seg_min+1/), count = (/nr,3,ntri,seg_groups/)))

            end if
            
            call try(nf90_close(ncid))
        end if

        call mpi_barrier(comm, ierr)

        end do

        return
    end subroutine export_magnetic_field

    subroutine diagnostics()
        !Calculates the time, open flux, total flux and magnetic energy, outputs into an appropriate netcdf file.
        !Need to deal with MPI - although that's a good check it's working I suppose
        implicit none
        real(num), dimension(:,:):: time(0:nprocs-1), oflux(0:nprocs-1)
        real(num), dimension(:,:):: energy(0:nprocs-1), lflux(0:nprocs-1)
        real(num), dimension(:,:):: sumj(0:nprocs-1), v_glob(0:nprocs-1)
        character(len=100):: filename

        integer:: i, id_1, id_2, id_3, id_4, id_5, id_6, ncid, nd_id
        if (diag_num == 0) then
            allocate(diag_time(0:ndiags-1)); allocate(diag_oflux(0:ndiags-1))
            allocate(diag_energy(0:ndiags-1)); allocate(diag_lflux(0:ndiags-1))
            allocate(diag_sumj(0:ndiags-1)); allocate(diag_v(0:ndiags-1))
            diag_time = 0.0_num; diag_oflux = 0.0_num; diag_energy = 0.0_num; diag_lflux = 0.0_num
            diag_sumj = 0.0_num; diag_v = 0.0_num
        end if

        !Do diagnostics locally
        time = 0.0_num; oflux = 0.0_num; energy = 0.0_num; lflux = 0.0_num; sumj = 0.0_num; v_glob = 0.0_num

        time(proc_num) = (step-1)*dt
        oflux(proc_num) = sum(abs(br(nr+1,0:ntri-1,seg_min:seg_max))*h_areas(nr+1,0:ntri-1,seg_min:seg_max))
        lflux(proc_num) = sum(abs(br(1,0:ntri-1,seg_min:seg_max))*h_areas(1,0:ntri-1,seg_min:seg_max))

        do i = 0,2
        energy(proc_num) = energy(proc_num) + 0.25_num*sum(dual_volumes(1:nr,0:npts-1,seg_min:seg_max)*&
        (b1(i,1:nr,0:npts-1,seg_min:seg_max)**2+b1(i,2:nr+1,0:npts-1,seg_min:seg_max)**2))
        end do
        !CHECKING various other things instead
        do i = 0,2
        sumj(proc_num) = sumj(proc_num) + 0.5_num*sum(dual_volumes(1:nr,0:npts-1,seg_min:seg_max)*&
        (j1(i,1:nr,0:npts-1,seg_min:seg_max)**2+j1(i,2:nr+1,0:npts-1,seg_min:seg_max)**2))
        end do

        do i = 0,2
        v_glob(proc_num) = v_glob(proc_num) + 0.5_num*sum(dual_volumes(1:nr,0:npts-1,seg_min:seg_max)*&
        (v1(i,1:nr,0:npts-1,seg_min:seg_max)**2+v1(i,2:nr+1,0:npts-1,seg_min:seg_max)**2))
        end do

        !Send to initial process
        if (proc_num > 0) call mpi_send(time(proc_num), 1, MPI_DOUBLE_PRECISION, 0, 11, comm, ierr)
        if (proc_num > 0) call mpi_send(oflux(proc_num), 1, MPI_DOUBLE_PRECISION, 0, 12, comm, ierr)
        if (proc_num > 0) call mpi_send(lflux(proc_num), 1, MPI_DOUBLE_PRECISION, 0, 13, comm, ierr)
        if (proc_num > 0) call mpi_send(energy(proc_num), 1, MPI_DOUBLE_PRECISION, 0, 14, comm, ierr)
        if (proc_num > 0) call mpi_send(sumj(proc_num), 1, MPI_DOUBLE_PRECISION, 0, 15, comm, ierr)
        if (proc_num > 0) call mpi_send(v_glob(proc_num), 1, MPI_DOUBLE_PRECISION, 0, 15, comm, ierr)

        call MPI_Barrier(comm,ierr)

        if (proc_num == 0) then
            do i = 1, nprocs-1
                call mpi_recv(time(i), 1, MPI_DOUBLE_PRECISION, i, 11, comm, MPI_STATUS_IGNORE, ierr)
                call mpi_recv(oflux(i), 1, MPI_DOUBLE_PRECISION, i, 12, comm, MPI_STATUS_IGNORE, ierr)
                call mpi_recv(lflux(i), 1, MPI_DOUBLE_PRECISION, i, 13, comm, MPI_STATUS_IGNORE, ierr)
                call mpi_recv(energy(i), 1, MPI_DOUBLE_PRECISION, i, 14, comm, MPI_STATUS_IGNORE, ierr)
                call mpi_recv(sumj(i), 1, MPI_DOUBLE_PRECISION, i, 15, comm, MPI_STATUS_IGNORE, ierr)
                call mpi_recv(v_glob(i), 1, MPI_DOUBLE_PRECISION, i, 15, comm, MPI_STATUS_IGNORE, ierr)

            end do
        end if

        call MPI_Barrier(comm,ierr)

        !Add up diagnostics and export
        if (proc_num == 0 .and. nprocs .le. 20) then
            diag_time(diag_num) = time(0)
            diag_oflux(diag_num) = sum(oflux)
            diag_energy(diag_num) = sum(energy)
            diag_lflux(diag_num) = sum(lflux)
            diag_sumj(diag_num) = sqrt(sum(sumj))
            diag_v(diag_num) = sqrt(sum(v_glob))

        else if (proc_num == 0) then
            diag_time(diag_num) = time(0)
            diag_oflux(diag_num) = sum(oflux(nprocs-20:nprocs-1))
            diag_energy(diag_num) = sum(energy)
            diag_lflux(diag_num) = sum(lflux(0:19))
            diag_sumj(diag_num) = sqrt(sum(sumj))
            diag_v(diag_num) = sqrt(sum(v_glob))
        end if

        !Export diagnostics
        if (run_id < 10) then
            write (filename, "(A18,I1,A3)") "./diagnostics/run0", run_id, ".nc"
        else
            write (filename, "(A17,I2,A3)") "./diagnostics/run", run_id, ".nc"
        end if


        if (proc_num == 0) then
        call try(nf90_create(trim(filename), nf90_clobber, ncid))
        !Define variables (in main process)
        call try(nf90_def_dim(ncid, 'ndiags', ndiags, nd_id))  !Make up fake dimensions here

        call try(nf90_def_var(ncid, 'time', nf90_double, (/nd_id/), id_1))
        call try(nf90_def_var(ncid, 'oflux', nf90_double, (/nd_id/), id_2))
        call try(nf90_def_var(ncid, 'lflux', nf90_double, (/nd_id/), id_3))
        call try(nf90_def_var(ncid, 'energy', nf90_double, (/nd_id/), id_4))
        call try(nf90_def_var(ncid, 'current squared', nf90_double, (/nd_id/), id_5))
        call try(nf90_def_var(ncid, 'frictional velocity', nf90_double, (/nd_id/), id_6))

        call try(nf90_enddef(ncid))

        call try(nf90_put_var(ncid, id_1, diag_time))
        call try(nf90_put_var(ncid, id_2, diag_oflux))
        call try(nf90_put_var(ncid, id_3, diag_lflux))
        call try(nf90_put_var(ncid, id_4, diag_energy))
        call try(nf90_put_var(ncid, id_5, diag_sumj))
        call try(nf90_put_var(ncid, id_6, diag_v))

        call try(nf90_close(ncid))

        if (diag_num > 0) then
        print*, 'Time = ', step*dt
        print*, 'Avg. current l_2', diag_sumj(diag_num)
        print*, 'Avg. fric l_2', diag_v(diag_num)
        print*, 'Energy', diag_energy(diag_num)
        print*, 'Open Flux', diag_oflux(diag_num)
        print*, '__________________________________________________________'

        end if
        end if

        return

    end subroutine diagnostics

END MODULE import_export
