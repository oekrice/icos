!*******************************************************************************
MODULE shared_data
!*******************************************************************************
! Establish the data that needs to be shared among all the bits of the program
!*******************************************************************************

    IMPLICIT NONE

!*******************************************************************************

    INTEGER, PARAMETER:: d = KIND(0.d0) ! precision for floats
    REAL(d), PARAMETER:: PI = 4.0_d * atan(1.0_d)   ! pi
    INTEGER, PARAMETER :: num = KIND(1.D0)

    ! MPI numbers (total number of processors, etc.)
    INTEGER:: ierr, nprocs, proc_num, comm

    ! Process position information
    INTEGER:: seg_min, seg_max, seg_groups    !If one segment per process, just use SEG. Otherwise need a loop within
    INTEGER:: nlayers, seg_layer, proc_up, proc_down, hamilton_flag

    REAL(num):: send_test

    !Control logical operators
    logical:: use_tweaked_pts

    !Parameters
    real(num), dimension(:,:):: parameters(0:20)
    real(num), dimension(:):: flparameters(0:19)
    real(num):: eta, eta0, nu0, eps, tmax, voutfact, shearfact, view_angle
    real(num):: dt, cfl
    integer:: nt, step, crot, run_id
    integer:: diag_num, plot_num, ndiags, nplots, nopen
    character(len=23):: init_filename

    ! Grid information
    real(num), dimension(0:2, 0:11):: ico_pts  !This contains actual coordinates so the first dimension should be 3

    integer, dimension(0:19,0:2):: ico_triangs
    integer, dimension(0:19,0:2,0:1):: ico_nbrs
    integer, dimension(0:19,0:2,0:4,0:1):: ico_corners
    integer, dimension(0:11,0:1,0:4):: pole_map

    integer:: G, n, npts, nr, nr_global, ntri, nghosttri, nghosts
    real(num):: rmax
    real(num), allocatable:: base_pts(:,:), all_base_pts(:,:,:)
    real(num), dimension(0:2,0:19):: tri_centre
    real(num), allocatable:: centres(:,:,:)

    !Grid population things
    real(num):: dp
    real(num), allocatable:: rs_global(:), rc_global(:)
    real(num), allocatable:: rs(:), rc(:)
    integer, allocatable:: triangs(:,:) , faces(:,:)  !

    integer, allocatable:: nnbrs(:), nbrpts(:,:)
    integer, allocatable:: updown(:), nbrfaces(:,:)

    real(num), allocatable:: nbr_dist(:,:,:), nbr_lens(:,:,:)
    real(num), allocatable:: areas(:,:), dual_areas(:,:), dual_weights(:)

    real(num), allocatable:: d_lens(:,:,:), tri_d(:,:,:), tri_l(:,:,:)

    integer, allocatable:: edge_transfer(:,:,:), tri_transfer(:,:,:)   !These are different for each segment - the only integer ones that are
    integer, allocatable:: trimap(:,:,:) !Triangles adjoining a given point, with orientation data
    integer, dimension(:,:):: f_orders(0:2,0:4), corner_pts_1(0:2), corner_pts_2(0:2)
    real(num), allocatable:: all_pts(:,:,:,:), h_areas(:,:,:), tri_l_faces(:,:,:,:)
    real(num), allocatable:: all_centres(:,:,:,:), all_ds(:,:,:,:), all_point_ds(:,:,:,:)
    real(num), allocatable:: tri_d_all(:,:,:,:), all_duals(:,:,:), tri_l_mids(:,:,:,:)

    real(num), allocatable:: v_areas(:,:,:,:), cell_volumes(:,:,:), v_duals(:,:,:,:)

    !Data import arrays

    real(num), allocatable:: br_import(:,:,:), bh_import(:,:,:,:)
    real(num), allocatable:: br_global(:,:,:), bh_global(:,:,:,:)  !Stored as FLUXES
    real(num), allocatable:: ar_global(:,:,:), ah_global(:,:,:,:)

    !Local variable arrays
    real(num), allocatable:: ar(:,:,:), ah(:,:,:,:)  !Stored as LINE INTEGRALS
    real(num), allocatable:: er(:,:,:), eh(:,:,:,:)  !Stored as LINE INTEGRALS

    real(num), allocatable:: br(:,:,:), bh(:,:,:,:)  !Stored as FIELD STRENGTH
    real(num), allocatable:: jr(:,:,:), jh(:,:,:,:)  !Stored as FIELD STRENGTH

    real(num), allocatable:: j1(:,:,:,:), b1(:,:,:,:), bu(:,:,:,:), e1(:,:,:,:), v1(:,:,:,:)

    real(num), allocatable:: b0(:,:,:,:)


    real(num), allocatable:: vout_field(:,:,:,:)

    real(num), allocatable:: diff_r(:,:,:), diff_h(:,:,:,:)  !Stored as FIELD STRENGTH
    real(num), allocatable::trisides(:,:,:,:), trinorms(:,:,:,:)
    real(num), allocatable:: alphas(:,:), betas(:,:), gammas(:,:)
    real(num), allocatable:: jcoeffs(:,:,:), bcoeffs(:,:,:)
    real(num), allocatable:: jbases(:,:,:,:), bbases(:,:,:,:)

    real(num), allocatable:: trisides_extend(:,:,:,:,:), centres_extend(:,:,:,:)
    real(num), allocatable:: jcoeffs_extend(:,:,:,:), bcoeffs_extend(:,:,:,:)
    real(num), allocatable:: jbases_extend(:,:,:,:,:), bbases_extend(:,:,:,:,:)

    !Field line things
    real(num):: ds_fact, t, ds
    integer:: nlines
    real(num), allocatable:: starts(:,:)
    integer, allocatable:: start_inds(:,:)

    real(num), allocatable:: new_line(:,:), smooth_line(:,:)
    real(num), allocatable:: all_lines(:,:,:), export_lines(:,:,:)
    integer:: line_length, max_length, export_length, current_length

    !Extra field line things
    real(num), allocatable:: top_starts(:,:), bottom_starts(:,:)
    integer, allocatable:: top_inds(:,:), bottom_inds(:,:)
    integer:: ntop, nbottom


!*******************************************************************************
END MODULE shared_data
!*******************************************************************************
