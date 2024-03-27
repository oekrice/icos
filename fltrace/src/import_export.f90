!*******************************************************************************
MODULE import_export
!*******************************************************************************
! Imports a potential (or otherwise) field, integrates to find the vector potential A and then distributes this as appropriate to all the processes
!*******************************************************************************
    USE shared_data
    USE netcdf

    IMPLICIT NONE

    INTEGER, PARAMETER:: nfreal = NF90_DOUBLE  ! type of output variables

    contains

    subroutine try(status)
    ! Catch error in reading netcdf fild.
    INTEGER, INTENT(IN):: status

    if (status /= NF90_noerr) THEN
        PRINT*,TRIM(ADJUSTL(NF90_STRERROR(status)))
    end if

    end subroutine try

    subroutine read_bfield()
        !Read in magnetic field from a netcdf file. This is FIELD STRENGTH, not fluxes
        implicit none
        integer:: ncid
        integer:: vid
        character(len=100):: filename
        character(len=100):: local_directory


        if (hamilton_flag > 0) then  !Read field lines from hamilton data instead  (copied into the plots folder)
            if (plot_num < 10) then
                write (filename, "(A54,I1,A3)") "/home/grads/trcn27/Documents/Thesis/icos/code/temp/000", plot_num, ".nc"
            else if (plot_num < 100) then
                write (filename, "(A53,I2,A3)") "/home/grads/trcn27/Documents/Thesis/icos/code/temp/00", plot_num, ".nc"
            else if (plot_num < 1000) then
                write (filename, "(A52,I3,A3)") "/home/grads/trcn27/Documents/Thesis/icos/code/temp/0", plot_num, ".nc"
            else if (plot_num < 10000) then
                write (filename, "(A51,I4,A3)") "/home/grads/trcn27/Documents/Thesis/icos/code/temp/", plot_num, ".nc"
            end if
        else
            if (run_id < 10) then
            write (local_directory, "(A29,I1,A1)") "/extra/tmp/trcn27/icos_data/0", run_id, "/"
            else
            write (local_directory, "(A28,I2,A1)") "/extra/tmp/trcn27/icos_data/", run_id, "/"
            end if
            if (plot_num < 10) then  !Read file from local machine
                write (filename, "(A31,A3,I1,A3)") local_directory, '000', plot_num, ".nc"
            else if (plot_num < 100) then
                write (filename, "(A31,A2,I2,A3)") local_directory, '00', plot_num, ".nc"
            else if (plot_num < 1000) then
                write (filename, "(A31,A1,I3,A3)") local_directory, '0', plot_num, ".nc"
            else if (plot_num < 10000) then
                write (filename, "(A31,I4,A3)") local_directory, plot_num, ".nc"
            end if

        end if

        if (.false.) then   !Read field lines from potential field plotter
            write (filename, "(A58,I1,A3,I4,A3)") "/home/grads/trcn27/Documents/Thesis/icos/pfss/inits/pfield",G,"_0_",crot,".nc"
        end if

        allocate(br_import(1:nr_global +1,1:ntri,1:20))   !No extras here - can add on later on if necessary. These are FIELD STRENGTH
        allocate(bh_import(1:nr_global ,1:3,1:ntri,1:20))

        call try(nf90_open(trim(filename), nf90_nowrite, ncid))  !Have to be VERY careful with the axes here, as they are stupid and don't make sense
        call try(nf90_inq_varid(ncid, 'br', vid))
        call try(nf90_get_var(ncid, vid, br_import))

        call try(nf90_inq_varid(ncid, 'bh', vid))
        call try(nf90_get_var(ncid, vid, bh_import))

        call try(nf90_close(ncid))

        allocate(br(0:nr +2,0:ntri+nghosttri-1,0:19))   !No extras here - can add on later on if necessary
        allocate(bh(0:nr +1,0:2,0:ntri+nghosttri-1,0:19))

        br = 0.0
        bh = 0.0

        br(1:nr_global +1,0:ntri-1, 0:19) = br_import
        bh(1:nr_global ,:,0:ntri-1, 0:19) = bh_import

        deallocate(bh_import)
        deallocate(br_import)

        return

    end subroutine read_bfield

    subroutine export_fieldlines()
        implicit none
        character(len=9):: filename
        integer:: aid, bid, cid, vid, ncid

        filename = "flines.nc"

        call try(nf90_create(trim(filename), nf90_clobber, ncid))

        !Define variables
        call try(nf90_def_dim(ncid, 'a', 3, aid))  !Make up fake dimensions here
        call try(nf90_def_dim(ncid, 'b', export_length, bid))  !Make up fake dimensions here
        call try(nf90_def_dim(ncid, 'c', nlines, cid))  !Make up fake dimensions here

        call try(nf90_def_var(ncid, 'lines', nf90_double, (/aid,bid,cid/), vid))
        call try(nf90_enddef(ncid))
        !Write variables

        call try(nf90_put_var(ncid, vid, export_lines))
        call try(nf90_close(ncid))

        return

    end subroutine export_fieldlines

    subroutine export_plane_fieldlines()
        implicit none
        character(len=15):: filename
        integer:: aid, bid, cid, vid, ncid

        filename = "flines_plane.nc"

        call try(nf90_create(trim(filename), nf90_clobber, ncid))

        !Define variables
        call try(nf90_def_dim(ncid, 'a', 3, aid))  !Make up fake dimensions here
        call try(nf90_def_dim(ncid, 'b', export_length, bid))  !Make up fake dimensions here
        call try(nf90_def_dim(ncid, 'c', ntop+nbottom, cid))  !Make up fake dimensions here

        call try(nf90_def_var(ncid, 'lines', nf90_double, (/aid,bid,cid/), vid))
        call try(nf90_enddef(ncid))
        !Write variables

        call try(nf90_put_var(ncid, vid, export_lines))
        call try(nf90_close(ncid))

        return

    end subroutine export_plane_fieldlines

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

END MODULE import_export
