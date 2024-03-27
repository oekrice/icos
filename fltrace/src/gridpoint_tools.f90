MODULE gridpoint_tools
    USE shared_data

    IMPLICIT NONE

    contains

    function coeffs(p, t0, t1, t2)

    !Given the coordinates of the triangle vertices t0, t1, t2, finds the linear combination of them used to make the vector p
    !Solving a linear system, using Gauss-Jordan elimination (woop 1st year algebra)
    real(num), dimension(:):: p(0:2), pt(0:2), t0(0:2), t1(0:2), t2(0:2), coeffs(0:2), swap1(0:2)
    real(num), dimension(:,:):: mat(0:2,0:2)
    real(num):: fact1, fact2, fact3, swap2

    pt = p(:)

    fact3 = 0.0
    mat(:,0) = t0(:)
    mat(:,1) = t1(:)
    mat(:,2) = t2(:)

    if (abs(mat(0,0)) < 1e-10 .and. abs(mat(2,0)) > 1e-10) then
        swap1 = mat(0,:)
        mat(0,:) = mat(2,:)
        mat(2,:) = swap1
        swap2 = pt(0)
        pt(0) = pt(2)
        pt(2) = swap2
    end if
    if (abs(mat(0,0)) < 1e-10 .and. abs(mat(1,0)) > 1e-10) then
        swap1 = mat(0,:)
        mat(0,:) = mat(1,:)
        mat(1,:) = swap1
        swap2 = pt(0)
        pt(0) = pt(1)
        pt(1) = swap2
    end if

    fact1 = mat(1,0)/mat(0,0)
    fact2 = mat(2,0)/mat(0,0)

    pt(1) = pt(1) - fact1*pt(0)
    pt(2) = pt(2) - fact2*pt(0)

    mat(1,:) = mat(1,:) - fact1*mat(0,:)
    mat(2,:) = mat(2,:) - fact2*mat(0,:)

    if (abs(mat(1,1)) < 1e-10) then
        swap1 = mat(2,:)
        mat(2,:) = mat(1,:)
        mat(1,:) = swap1
        swap2 = pt(2)
        pt(2) = pt(1)
        pt(1) = swap2
    end if
    fact3 = mat(2,1)/mat(1,1)

    pt(2) = pt(2) - fact3*pt(1)
    mat(2,:) = mat(2,:) - fact3*mat(1,:)

    coeffs(2) = pt(2)/mat(2,2)
    coeffs(1) = (pt(1) - mat(1,2)*pt(2)/mat(2,2))/mat(1,1)
    coeffs(0) = (pt(0) - mat(0,1)*coeffs(1) - mat(0,2)*coeffs(2))/mat(0,0)

    end function coeffs

    function solve3d(A_mat, rhs)
    !Solves 3x3 system Ax = b using Gauss_jordan elimination
    real(num), dimension(:):: rhs(0:2), solve3d(0:2), swap1(0:2)
    real(num), dimension(:,:):: lhs(0:2,0:2), A_mat(0:2,0:2)
    real(num):: fact1, fact2, fact3, swap2

    fact3 = 0.0
    lhs = A_mat(:,:)

    if (abs(lhs(0,0)) < 1e-10 .and. abs(lhs(2,0)) > 1e-10) then
        swap1 = lhs(0,:)
        lhs(0,:) = lhs(2,:)
        lhs(2,:) = swap1
        swap2 = rhs(0)
        rhs(0) = rhs(2)
        rhs(2) = swap2
    end if
    if (abs(lhs(0,0)) < 1e-10 .and. abs(lhs(1,0)) > 1e-10) then
        swap1 = lhs(0,:)
        lhs(0,:) = lhs(1,:)
        lhs(1,:) = swap1
        swap2 = rhs(0)
        rhs(0) = rhs(1)
        rhs(1) = swap2
    end if

    fact1 = lhs(1,0)/lhs(0,0)
    fact2 = lhs(2,0)/lhs(0,0)

    rhs(1) = rhs(1) - fact1*rhs(0)
    rhs(2) = rhs(2) - fact2*rhs(0)

    lhs(1,:) = lhs(1,:) - fact1*lhs(0,:)
    lhs(2,:) = lhs(2,:) - fact2*lhs(0,:)

    if (abs(lhs(1,1)) < 1e-10) then
        swap1 = lhs(2,:)
        lhs(2,:) = lhs(1,:)
        lhs(1,:) = swap1
        swap2 = rhs(2)
        rhs(2) = rhs(1)
        rhs(1) = swap2
    end if
    fact3 = lhs(2,1)/lhs(1,1)

    rhs(2) = rhs(2) - fact3*rhs(1)
    lhs(2,:) = lhs(2,:) - fact3*lhs(1,:)

    solve3d(2) = rhs(2)/lhs(2,2)
    solve3d(1) = (rhs(1) - lhs(1,2)*rhs(2)/lhs(2,2))/lhs(1,1)
    solve3d(0) = (rhs(0) - lhs(0,1)*solve3d(1) - lhs(0,2)*solve3d(2))/lhs(0,0)

    end function solve3d

    function cross1(a,b)
    !Calculates the cross product of the vectors a and b
    real(num), dimension(:):: a(0:2), b(0:2), cross1(0:2)
    cross1(0) = a(1)*b(2) - a(2)*b(1)
    cross1(1) = a(2)*b(0) - a(0)*b(2)
    cross1(2) = a(0)*b(1) - a(1)*b(0)
    end function cross1

    function cross2(a,b)
    !Calculates the cross product of an array of vectors a and b (on gridpoints)
    real(num), dimension(:,:,:,:):: a(0:2,0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max)
    real(num), dimension(:,:,:,:):: b(0:2,0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max)
    real(num), dimension(:,:,:,:):: cross2(0:2,0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max)
    integer:: seg, k ,ri
    do seg = seg_min, seg_max
    do k = 0, ntri+nghosttri-1
    do ri = 0, nr+2
    cross2(0,ri,k,seg) = a(1,ri,k,seg)*b(2,ri,k,seg) - a(2,ri,k,seg)*b(1,ri,k,seg)
    cross2(1,ri,k,seg) = a(2,ri,k,seg)*b(0,ri,k,seg) - a(0,ri,k,seg)*b(2,ri,k,seg)
    cross2(2,ri,k,seg) = a(0,ri,k,seg)*b(1,ri,k,seg) - a(1,ri,k,seg)*b(0,ri,k,seg)
    end do
    end do
    end do
    end function cross2

    function cross3(a,b)
    !Slight modification of the above - used for the upwinding in the outflow term
    real(num), dimension(:,:,:,:):: a(0:2,0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max)
    real(num), dimension(:,:,:,:):: b(0:2,0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max)
    real(num), dimension(:,:,:,:):: cross3(0:2,0:nr+2,0:ntri+nghosttri-1,seg_min:seg_max)
    integer:: seg, k ,ri
    cross3 = 0.0_num
    do seg = seg_min, seg_max
    do k = 0, ntri+nghosttri-1
    do ri = 2, nr+1
    cross3(0,ri,k,seg) = a(1,ri,k,seg)*b(2,ri-1,k,seg) - a(2,ri,k,seg)*b(1,ri-1,k,seg)
    cross3(1,ri,k,seg) = a(2,ri,k,seg)*b(0,ri-1,k,seg) - a(0,ri,k,seg)*b(2,ri-1,k,seg)
    cross3(2,ri,k,seg) = a(0,ri,k,seg)*b(1,ri-1,k,seg) - a(1,ri,k,seg)*b(0,ri-1,k,seg)
    end do
    end do
    end do
    end function cross3



    function vout_fn(r)
    ! The outflow velocity function
    real(num):: r_crit, vout_fn, r

    r_crit = 10.0
    vout_fn = voutfact*((rmax**2*exp(-2.0*r_crit/r))/(r**2*exp(-2.0*r_crit/rmax)))

    end function vout_fn

    function shear_fn(lat)
    ! The shearing velocity as a function of latitude (radians from the equator). In degrees/day, multiplied by the latitude.
    ! So is actually more like R_0 per day...
    real(num):: lat, shear_fn

    shear_fn = 0.25*(shearfact*(0.18-2.3*cos(lat)**2 - 1.62*cos(lat)**4)*cos(lat))

    end function shear_fn

    subroutine printmat2(array,nx,ny)  !prints a matrix with dimensions nx, ny the right way around
        implicit none
        real(num):: array(nx,ny)
        integer:: nx, ny, i, j

        if (proc_num == 0) then
            do i = 1, nx, 1
                write(*,10)(array(i,j), j=1,nx)
            end do
            print*, ''
            10 format(100f10.4)  !Maximum matrix size 100
        end if
        return
    end subroutine


END MODULE gridpoint_tools
