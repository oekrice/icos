!*******************************************************************************
MODULE mpi_tools
!*******************************************************************************
! Contains some basic things done at the start used to initialise the MPI. Seems to work.
!*******************************************************************************
    USE shared_data

    IMPLICIT NONE

    contains

    subroutine start_mpi()
        ! - Have already established the global rank.
        ! - This establishes the position within the grid depending on the number of processes used
        ! - seg_layer is the radial distribution of the processes - only nonzero if there are more than 20
        ! - seg_min and seg_max are the maximum and minimum segments that each process deals with. If there are more than 20 processes there will only be one segment each

        call mpi_init(ierr)  !Tells it to start using MPI

        call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr) !Number of processes globally.
        call mpi_comm_rank(MPI_COMM_WORLD, proc_num, ierr) !Returns the rank of current process

        if (nprocs < 20) then
            if (mod(20, nprocs) .ne. 0) then
                if (proc_num == 0) then
                    print*, 'Unacceptable number of processors. Number must be a multiple or factor of 20'
                    call mpi_abort(comm, 0, ierr)  !Assuming this stops the whole thing.
                end if
            else
                nlayers = 1; seg_layer = 0; seg_groups = 20/nprocs
                seg_min = proc_num*seg_groups
                seg_max = (proc_num + 1)*seg_groups - 1
            end if
            proc_up = -1; proc_down = -1   !No processors above or below - just one layer

        else if (nprocs > 20) then !More than 20 processes
            if (mod(nprocs, 20) .ne. 0) then
                if (proc_num == 0) then
                    print*, 'Unacceptable number of processors. Number must be a multiple or factor of 20'
                    call mpi_abort(comm, 0, ierr)  !Assuming this stops the whole thing.
                end if
            !Check layers are not too small
            else if (2**G/(nprocs/20) < 4) then
                if (proc_num == 0) then
                    print*, 'Unacceptable number of processors. Layers are too thin to make any sense - reduce nprocs'
                    call mpi_abort(comm, 0, ierr)  !Assuming this stops the whole thing.
                end if

            else
                nlayers = nprocs/20; seg_layer = proc_num/20; ; seg_groups = 1

                if (mod(2**G, nlayers) .ne. 0) then
                    if (proc_num == 0) then
                        print*, 'Unacceptable number of processors. Radial grid cells do not divide evenly.'
                        call mpi_abort(comm, 0, ierr)  !Assuming this stops the whole thing.
                    end if
                end if

                seg_min = mod(proc_num*seg_groups, 20)
                seg_max = mod((proc_num + 1)*seg_groups - 1, 20)
            end if

            proc_up = proc_num + 20
            proc_down = proc_num - 20
            if (proc_num < 20) then
                proc_down = -1
            end if
            if (proc_num >= nprocs - 20) then
                proc_up = -1
            end if
        else
            nlayers = 1; seg_layer = 0; seg_groups = 1
            seg_min = proc_num*seg_groups
            seg_max = proc_num*seg_groups
            proc_up = -1; proc_down = -1 !No processors above or below - just one layer

        end if

        return
    end subroutine start_mpi


END MODULE mpi_tools





















