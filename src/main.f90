!*******************************************************************************
PROGRAM main
!*******************************************************************************
! Main fortran file. Use this to initialise parameters, grid, produce/read in initial conditions and then run the code.
! Most subroutines to calculate things are in ops.f90
! Subroutines to do other things are in import_export.f90, or gridpoint_tools for one or two things
!*******************************************************************************
    USE shared_data
    USE mpi_tools
    USE grid
    USE import_export
    USE ops

    IMPLICIT NONE

    character(len=2):: var
    character(len=16)::paras_fname

    use_tweaked_pts = .true.
    parker = .true.    !Set to true for the profile used in the outflow code stuff, set to false for the power law used by DumFric previously

    !Import parameters from text file (saved out by python)
    call get_command_argument(1, var)
    read(unit=var,fmt=*) run_id

    !Find the filename of the parameters to be imported
    if (run_id < 10) then
        write (paras_fname, "(A11,I1,A4)") "parameters0", run_id, ".txt"
    else
        write (paras_fname, "(A10,I2,A4)") "parameters", run_id, ".txt"
    end if
    !Import parameters from text file (saved out by python)
    open(1,file= paras_fname)
    read(1, *) parameters
    close(1)

    ! Read in parameters that have been saved out by 'run.py'
    G = int(parameters(0))
    rmax = parameters(1)
    tmax = parameters(3)

    eta = parameters(4)
    nu0 = parameters(5)
    eta0 = parameters(6)
    voutfact = parameters(7)
    shearfact = parameters(8)

    cfl = parameters(9)
    eps = parameters(10)
    crot = int(parameters(11))  !Carrington rotation of initial condition, or index of initial condition if not using hmi data
    ndiags = int(parameters(12))  !Number of diagnostic outputs (open flux etc.)
    nplots = int(parameters(13))  !Number of full magnetic field saves (which could be rather large)

    run_id = int(parameters(14))
    shear_frame = int(parameters(16))
    rsun = parameters(17)
    hflag = int(parameters(18))
    cfric = int(parameters(19))


    ! Do some admin
    if (hflag < 0.5) then
        if (run_id < 10) then
            write (data_directory, "(A29,I1,A1)") "/extra/tmp/trcn27/icos_data/0" , run_id, "/"
        else
            write (data_directory, "(A28,I2,A1)") "/extra/tmp/trcn27/icos_data/" , run_id, "/"
        end if

    else
        if (run_id < 10) then
            write (data_directory, "(A28,I1,A1)") "/nobackup/trcn27/icos_data/0", run_id, "/"  !For hamilton
        else
            write (data_directory, "(A27,I2,A1)") "/nobackup/trcn27/icos_data/", run_id, "/"  !For hamilton
        end if
    end if

    write (init_filename, "(A12,I1,A1,I1,A1,I4,A3)") "inits/pfield", G, "_", 0, "_", crot, ".nc"

    if (cfric > 0.5_num) then
        constant_friction = .true.
    else
        constant_friction = .false.
    end if

    !Actually run the code

    !Initialise everything:
    call start_mpi()  !Establish process position within the grid

    if (proc_num == 0) then
        print*, 'Parameters read in:'
        print*, 'nu0:', nu0, 'eta0:', eta0, 'eta:', eta
        print*, 'outflow:', voutfact, 'shearing:', shearfact
        print*, '__________________________________________________________'
    end if

    if (proc_num == 0) print*, 'Number of processors', nprocs
    if (proc_num == 0) print*, 'Output directory according to fortran: ', trim(data_directory)
    if (proc_num == 0) print*, 'Initial condition filename = ', init_filename

    call establish_grid()     !Import the parameters and set up the grid


    if (proc_num == 0) print*, 'Grid established, G = ', G

    call import_initial_condition()  !This loads in the initial potential field from the netcdf file. One process then integrates to find the vector potential, which is then distributed to the individual processes
    !Run everything:
    call evolve()
    !Finish everything:
    call mpi_finalize(ierr)  !Let's MPI know that it's finished

    contains

    subroutine evolve()

        implicit none


        call allocate_variable_arrays()
        call find_dt()  !Calculates the timestep. This function to be in ops
        call make_vout_field()
        call make_fric_field()

        if (proc_num == 0) print*, 'Running code'
        if (proc_num == 0) print*, '__________________________________________________________'

        diag_num = 0; plot_num = 0
        do step = 1,  nt
            !Calculate magnetic field from the vector potential
            call calculate_magnetic()
            !Check everything is still divergence-free. Will terminate the code if not.
            call check_divfree_interior()
            !Populate magnetic field ghost points
            call bfield_ghosts_1()  !uses MPI
            call bfield_bound()
            call bfield_ghosts_2()   !This is the same as bfield_ghosts_1 but only transfers the cells at the top and bottom (those affected by bfield_bound)
            !Calculate current
            call calculate_current()
            !Find some of the current ghost points (don't need all of them)
            call jfield_ghosts()  !Populates the horizontal ghost triangles for the current field. No need to do radial current.
            !Average the current and magnetic field to the grid points
            call j_to_gridpts()
            call b_to_gridpts()
            !Calculate magnetofrictional velocity
            call calculate_friction()
            !Add shearing
            if (seg_layer == 0) then
                call add_shearing()
            end if
            !Calculate the electric field (averaged to gridpoints)
            call calculate_electric()
            !Integrate the electric field to be in the same place as teh vector potential
            call electric_integral()
            !Add outflow
            call add_outflow()
            !Add coronal diffusion directly to eh, er
            call add_diffusion()
            !Add supergranular diffusion
            if (seg_layer == 0) then
                call add_surface_diffusion()
            end if

            ! Output magnetic field for plotting

            if ((step-1)*dt .ge. plot_num*(tmax/float(nplots-1))) then
                call export_with_checks
                if (proc_num == 0) print*, 'time = ', step*dt
                call reset_afield

                plot_num = plot_num + 1
            end if

            ! Calculate and export diagnostics
            if ((step-1)*dt .ge. diag_num*(tmax/float(ndiags-1))) then
                call diagnostics
                diag_num = diag_num + 1
                if (proc_num == 0) then
                end if

            end if

            !Do the timestep
            ar(1:nr,0:npts-1,seg_min:seg_max) = ar(1:nr,0:npts-1,seg_min:seg_max) - dt*er(1:nr,0:npts-1,seg_min:seg_max)
            ah(1:nr +1,0:2,0:ntri-1,seg_min:seg_max) = ah(1:nr +1,0:2,0:ntri-1,seg_min:seg_max) - &
            dt*eh(1:nr +1,0:2,0:ntri-1,seg_min:seg_max)

            call MPI_barrier(comm, ierr)

        end do

        !At the end, export one last time
        call export_with_checks
        step = int(1 + tmax/dt)
        call diagnostics
        diag_num = diag_num + 1

        if (proc_num == 0) print*, 'Code terminated sucessfully'

    end subroutine evolve

END PROGRAM main

























