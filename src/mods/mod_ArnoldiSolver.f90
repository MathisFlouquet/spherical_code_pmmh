module mod_ArnoldiSolver
 
    use iso_c_binding
 
    use mod_GlobalVars
    use mod_PrecompSH
    use mod_read
    use mod_TimeStep
    use mod_Output
    use mod_tools
    use mod_ActNewtonSolver
 
    implicit none
 
contains
 
 
 
subroutine output_vector(E_vec_ptr,F_vec_ptr,T_vec_ptr)
 
    double precision, dimension(2 * KK2 * shtns%nlm), target :: E_vec_ptr 
    double precision, dimension(2 * KK4 * shtns%nlm), target :: F_vec_ptr
    double precision, dimension(2 * KK2 * shtns%nlm), target :: T_vec_ptr 
 
    double complex, dimension(:,:), pointer :: E_vec
    double complex, dimension(:,:), pointer :: F_vec
    double complex, dimension(:,:), pointer :: T_vec
 
    call c_f_pointer(c_loc(E_vec_ptr), E_vec, [KK2, shtns%nlm])
    call c_f_pointer(c_loc(F_vec_ptr), F_vec, [KK4, shtns%nlm])
    call c_f_pointer(c_loc(T_vec_ptr), T_vec, [KK2, shtns%nlm])
 
 
    call comp_U_from_EF(E_vec, F_vec, Ur, Ut, Up)
    call ToReal(T_vec, T_real, KK2)
 
    call Output_files(Ur, Up, T_real, step)
 
end subroutine output_vector
 
! ===========================================================================
!  write_mode_restart : write one eigenvector as a restart file.
!  The packed real vector (E, F, T) is reinterpreted as complex spectral
!  fields and written with the exact same record structure as
!  writeRestart, so it can be reread directly by readRestart.
!  A trailing (t0, step_0) = (0, 0) record is written so the file is read
!  back as a fresh timestepping restart.
!  Files are named Restart_mode_<idx>.b in the run directory.
! ===========================================================================
 
subroutine write_mode_restart(E_vec_ptr, F_vec_ptr, T_vec_ptr, idx)
 
    double precision, dimension(2 * KK2 * shtns%nlm), target :: E_vec_ptr
    double precision, dimension(2 * KK4 * shtns%nlm), target :: F_vec_ptr
    double precision, dimension(2 * KK2 * shtns%nlm), target :: T_vec_ptr
    integer, intent(in) :: idx
 
    double complex, dimension(:,:), pointer :: E_vec
    double complex, dimension(:,:), pointer :: F_vec
    double complex, dimension(:,:), pointer :: T_vec
 
    character(len=64) :: fname
 
    call c_f_pointer(c_loc(E_vec_ptr), E_vec, [KK2, shtns%nlm])
    call c_f_pointer(c_loc(F_vec_ptr), F_vec, [KK4, shtns%nlm])
    call c_f_pointer(c_loc(T_vec_ptr), T_vec, [KK2, shtns%nlm])
 
    write(fname, '(A,I0,A)') trim(directory)//"/Restart_mode_", idx, ".b"
    open(unit=71, file=trim(fname), form="unformatted")
    write(71) E_vec
    write(71) F_vec
    write(71) T_vec
    write(71) 0.0d0    ! t0     -> reread as a timestepping restart
    write(71) 0        ! step_0
    close(71)
 
end subroutine write_mode_restart
 
! ===========================================================================
!  apply_propagator : apply the temporal integration NB_loc times.
!  Used in ARPACK with the reverse communication loop (ido = 1 or -1)
! ===========================================================================
 
 
subroutine apply_propagator(x_in, x_out, n_prop, LinNonLinTimeStep_ptr)
 
    integer,          intent(in)  :: n_prop
    double precision, intent(in)  :: x_in(n_prop)
    double precision, intent(out) :: x_out(n_prop)
    procedure(),      pointer     :: LinNonLinTimeStep_ptr
 
    double precision :: EFT(n_prop)
    integer          :: step, NB_loc
 
    NB_loc = 1000
 
    EFT = x_in
 
    do step = 1, NB_loc
        call act_lin_non_lin(LinNonLinTimeStep_ptr,          &
                             EFT(1),                         &
                             EFT(2*KK2*shtns%nlm + 1),      &
                             EFT(2*(KK2+KK4)*shtns%nlm + 1))
    end do
 
    x_out = EFT
 
end subroutine apply_propagator
 
 
! ===========================================================================
!  arnoldi_convective_solver : ARPACK IRAM version (dnaupd / dneupd)
!
!  Real non-symmetric standard eigenvalue problem : P*x = mu*x
!  where P ~ exp(A * NB_loc * delta_t) is the propagator.
!
!  Physical eigenvalues of A are recovered via:
!    lambda = ln(mu) / (NB_loc * delta_t)
! The eigenvectors mu are ordered by magnitude
! ===========================================================================
 
subroutine arnoldi_convective_solver()
 
    ! ----- ARPACK parameters -----
    integer            :: n, nev, ncv, ldv, lworkl, maxiter
    integer            :: ido, info_arpack
    integer            :: iparam(11), ipntr(14)
    character(len=1)   :: bmat
    character(len=2)   :: which
    double precision   :: tol
    double precision   :: sigmar, sigmai
    ! n : 2*(2*KK2+KK4)*shtns%nlm size of the space
    ! nev : number of eigenvalues calculated ( eigenvalues come by pair )
    ! ncv : size of the krylov space, most of the time ncv>2*nev
    ! ldv : 
    ! lworkl : dimension of the work array 
    ! maxiter : maximum number of restart
    ! ido : -1, 1 or 99, use to communicate with ARPACK
    ! info_arpack : ...
    ! iparam : informations given by and for arpack
    ! ipntr : informations given by and for arpack
    ! bmat : 'I' , standard problem ('G' leads to generalized)
    ! which : to indicate the criteria used to choose the eigenvalues : here 'LM' so largest magnitude
    ! tol : tolerance by which the eigenvalues are calculated : norm(λ_i_λ_{i+1})<=tol*norm(λ_i)
    ! sigmar : to store real part
    ! sigmai : to store imaginary part
 
    double precision, allocatable :: resid(:)       ! residual / initial vector         (n)
    double precision, allocatable :: v(:,:)         ! Krylov basis                      (n x ncv)
    double precision, allocatable :: workd(:)       ! reverse communication work array  (3*n)
    double precision, allocatable :: workl(:)       ! internal work array               (lworkl)
    double precision, allocatable :: workev(:)      ! work array for dneupd             (3*ncv)
    logical,          allocatable :: select_v(:)    ! selection array for dneupd        (ncv)
 
    ! ----- eigenvalues and eigenvectors -----
    double precision, allocatable :: dr(:)          ! real parts of eigenvalues         (nev+1)
    double precision, allocatable :: di(:)          ! imaginary parts of eigenvalues    (nev+1)
    double precision, allocatable :: z(:,:)         ! Ritz vectors                      (n x nev+1)
    double precision, allocatable :: mode_R(:,:)    ! real parts of physical modes      (n x nev+1)
    double precision, allocatable :: mode_I(:,:)    ! imaginary parts of physical modes (n x nev+1)
 
    ! ----- a posteriori residual (cf. dnsimp.f) -----
    double precision, allocatable :: ax(:)          ! temporary for P*z                 (n)
    double precision :: norm_Pv1, norm_Pv2          ! ||P*z(:,j)|| and ||P*z(:,j+1)||
    double precision :: norm_r1, norm_r2            ! ||residual_real|| and ||residual_imag||
    double precision :: spectral_resid              ! ||P*v - mu*v|| / ||P*v||
 
    ! ----- sigma and omega (Hermitian inner product estimates) -----
    double precision :: sigma_mode                  ! = 1/(2*T_prop) * ln((v,v)/(u,u))
    double precision :: omega_mode                  ! = 1/T_prop * acos(Re(u,v)/sqrt((u,u)*(v,v)))
    double precision :: norm_u, norm_v              ! Hermitian norms of u and v
    double precision :: uv_real                     ! Re(u,v) = sum of real dot products
    double precision :: cos_arg                     ! argument of acos, clamped to [-1,1]
 
    ! ----- physical eigenvalue recovery -----
    integer          :: NB_loc
    double precision :: T_prop                      ! = NB_loc * delta_t
    double precision :: mu_abs, mu_arg
    double precision :: lambda_r, lambda_i
 
    ! ----- BLAS / LAPACK -----
    double precision, external :: dnrm2, dlapy2, ddot
    external :: daxpy
 
    ! ----- misc -----
    procedure(), pointer :: LinNonLinTimeStep_ptr
    integer          :: i, j, nconv
    integer          :: matvec_count, restart_count
    logical          :: rvec, first
    character(len=32):: filename
    integer          :: iunit
 

    ! ==================================================================
    !  Open output text file
    ! ==================================================================
    iunit = 42
    open(unit=iunit, file='arnoldi_output.txt', status='replace', action='write')

    ! ==================================================================
    !  ARPACK parameters
    ! ==================================================================
    n       = 2*(2*KK2+KK4)*shtns%nlm
    nev     = 2
    ncv     = 70
    maxiter =  4
    tol     = 1.0d-11
    bmat    = 'I'
    which   = 'LM'
    ldv     = n
    lworkl  = 3*ncv**2 + 6*ncv
 
    ! Propagation parameters (must match apply_propagator)
    NB_loc = 1000
    T_prop = dble(NB_loc) * delta_t
 
    allocate(resid(n))
    allocate(v(ldv, ncv))
    allocate(workd(3*n))
    allocate(workl(lworkl))
    allocate(workev(3*ncv))
    allocate(select_v(ncv))
    allocate(dr(nev+1))
    allocate(di(nev+1))
    allocate(z(n, nev+1))
    allocate(mode_R(n, nev+1))
    allocate(mode_I(n, nev+1))
    allocate(ax(n))
 
    ! ==================================================================
    !  ARPACK iparam control 
    ! ==================================================================
    iparam    = 0
    iparam(1) = 1        ! ishift=1 : exact shifts
    iparam(3) = maxiter  ! maximum number of implicit restarts
    iparam(7) = 1        ! mode=1 : standard problem A*x = lambda*x
 
    ! ==================================================================
    !  Precomputation of matrices
    ! ==================================================================
    call precompBuildXY()
 
    if (solver == "arnoldi_convective_explicit") then
 
        select case (time_step)
            case ("cn")
                print*, "Precomputing the Xe and Ye matrices..."
                write(iunit,*) "Precomputing the Xe and Ye matrices..."
                call precompXeYe()
                print*, "Precomputing the Xf and Yf matrices..."
                write(iunit,*) "Precomputing the Xf and Yf matrices..."
                call precompXfYf()
                print*, "Precomputing the XT and YT matrices..."
                write(iunit,*) "Precomputing the XT and YT matrices..."
                call precompXTYT()
            case ("fbe")
                print*, "Precomputing the Xe and Ye matrices..."
                write(iunit,*) "Precomputing the Xe and Ye matrices..."
                call precompXeYe_BDF()
                print*, "Precomputing the Xf and Yf matrices..."
                write(iunit,*) "Precomputing the Xf and Yf matrices..."
                call precompXfYf_BDF()
                print*, "Precomputing the XT and YT matrices..."
                write(iunit,*) "Precomputing the XT and YT matrices..."
                call precompXTYT_BDF()
        end select
 
        LinNonLinTimeStep_ptr => compute_lin_time_step_convective_explicit_CN_FBE
 
    else if (solver == "arnoldi_convective_implicit") then
 
        select case (time_step)
            case ("cn")
                print*, "Precomputing the Xef and Yef matrices..."
                write(iunit,*) "Precomputing the Xef and Yef matrices..."
                call PrecompimplicitXY()
                print*, "Precomputing the XT and YT matrices..."
                write(iunit,*) "Precomputing the XT and YT matrices..."
                call precompXTYT()
                LinNonLinTimeStep_ptr => compute_lin_time_step_convective_implicit_CN
            case ("fbe")
                print*, "Precomputing the Xef, Ye and Yf matrices..."
                write(iunit,*) "Precomputing the Xef, Ye and Yf matrices..."
                call PrecompimplicitXY_BDF()
                print*, "Precomputing the XT and YT matrices..."
                write(iunit,*) "Precomputing the XT and YT matrices..."
                call precompXTYT_BDF()
                LinNonLinTimeStep_ptr => compute_lin_time_step_convective_implicit_FBE
        end select
 
    end if
 
    print*
    write(iunit,*)
    print*, "--------------------- Initial values ---------------------------"
    write(iunit,*) "--------------------- Initial values ---------------------------"
    print*
    write(iunit,*)
    print*, "Volume of the shell = ", Vol
    write(iunit,*) "Volume of the shell = ", Vol
    print*, "Reading restart files"
    write(iunit,*) "Reading restart files"
    call readDim()
    call readRestart()
 
    E_base = E
    F_base = F
    T_base = T


    ! Pour état conductif, fixer E F T C_base à 0!
    E=0
    F=0
    T=0
    C_base=0
    E_base=E
    F_base=F
    T_base=T
    
    call comp_radial_derivatives(E_base, F_base, T_base, Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St)
    call comp_U_cU_gT(Ur_base, Ut_base, Up_base, cUr_base, cUt_base, cUp_base, gTr_base, gTt_base, gTp_base)
 


    call c_f_pointer(c_loc(E_base), E_base_ptr, [2 * KK2 * shtns%nlm])
    call c_f_pointer(c_loc(F_base), F_base_ptr, [2 * KK4 * shtns%nlm])
    call c_f_pointer(c_loc(T_base), T_base_ptr, [2 * KK2 * shtns%nlm])

    step=0
    call output_vector(E_base_ptr,F_base_ptr,T_base_ptr)




    ! ==================================================================
    !  ARPACK reverse communication loop 
    ! ==================================================================
    info_arpack  = 0
    ido          = 0
    matvec_count = 0
    restart_count= 0
 
    do
        call dnaupd(ido, bmat, n, which, nev, tol,  &
                    resid, ncv, v, ldv,              &
                    iparam, ipntr, workd, workl,     &
                    lworkl, info_arpack)
 
        if (ido == -1 .or. ido == 1) then ! Arpack needs a new vector to be created
 
            call apply_propagator(workd(ipntr(1)), workd(ipntr(2)), n, LinNonLinTimeStep_ptr)
            matvec_count = matvec_count + 1 ! Number of vector created
 
            if (matvec_count == ncv) then   ! If we have built ncv vector, then the first Krylov space is created
                print*, "--- Initial Krylov space built ---"
                write(iunit,*) "--- Initial Krylov space built ---"
                print*, "    matvecs : ", matvec_count
                write(iunit,*) "    matvecs : ", matvec_count
                flush(6)
                flush(iunit)
            else if (matvec_count > ncv) then
                if (mod(matvec_count - ncv, ncv - nev) == 0) then ! To monitor the number of restart
                    restart_count = restart_count + 1
                    print*, "--- Restart ", restart_count, " / ", maxiter,  &
                            " | matvecs : ", matvec_count,                  &
                            " | converged : ", iparam(5), " / ", nev
                    write(iunit,*) "--- Restart ", restart_count, " / ", maxiter,  &
                            " | matvecs : ", matvec_count,                          &
                            " | converged : ", iparam(5), " / ", nev
                    flush(6)
                    flush(iunit)
                end if
            end if
            
            print*, "  wavespeed :", C_base
            write(iunit,*) "  wavespeed :", C_base
        else
            exit
        end if
 
    end do
 
    ! ==================================================================
    !  Check dnaupd return code
    ! ==================================================================
    if (info_arpack < 0) then
        print*, "ARPACK dnaupd error: info =", info_arpack
        write(iunit,*) "ARPACK dnaupd error: info =", info_arpack
        close(iunit)
        stop
    else if (info_arpack == 1) then
        print*, "Warning: dnaupd reached maxiter without full convergence."
        write(iunit,*) "Warning: dnaupd reached maxiter without full convergence."
        print*, "Converged Ritz values:", iparam(5), "/", nev
        write(iunit,*) "Converged Ritz values:", iparam(5), "/", nev
    else
        print*, "dnaupd converged. Ritz values computed:", iparam(5)
        write(iunit,*) "dnaupd converged. Ritz values computed:", iparam(5)
    end if
    print*, "Total matvecs  :", matvec_count
    write(iunit,*) "Total matvecs  :", matvec_count
    print*, "Total restarts :", restart_count
    write(iunit,*) "Total restarts :", restart_count
    flush(6)
    flush(iunit)
 
 
    rvec   = .true.
    sigmar = 0.0d0
    sigmai = 0.0d0
 
    call dneupd(rvec, 'A', select_v, dr, di, z, ldv,  &
                sigmar, sigmai, workev,                &
                bmat, n, which, nev, tol,              &
                resid, ncv, v, ldv,                    &
                iparam, ipntr, workd, workl,           &
                lworkl, info_arpack)
 
    if (info_arpack /= 0) then
        print*, "ARPACK dneupd error: info =", info_arpack
        write(iunit,*) "ARPACK dneupd error: info =", info_arpack
        close(iunit)
        stop
    end if
 
    nconv = iparam(5)
 
    ! ==================================================================
    !  A posteriori residual + physical eigenvalues
    !
    !  Residual computed as: ||P*v - mu*v|| / ||P*v||
    !
    !  For complex conjugate pairs (mu, conj(mu)), the residual is
    !  identical. We compute it once using 2 propagator calls for the
    !  pair (P*z_j and P*z_{j+1}) and skip the conjugate.
    !  -> nconv/2 pairs * 2 calls = nconv calls total (not 2*nconv)
    !
    !  Physical eigenvalue:
    !    lambda = ln(mu) / T_prop   with T_prop = NB_loc * delta_t
    !
    !  Hermitian estimates (u = eigenvector, v = P*u) :
    !    sigma = 1/(2*T_prop) * ln((v,v)/(u,u))
    !    omega = 1/T_prop     * acos(Re(u,v) / sqrt((u,u)*(v,v)))
    !
    !  For real vectors (u,u) = ||u||^2, (v,v) = ||v||^2,
    !  Re(u,v) = dot(u,v).
    !  For complex pairs u = z_R + i*z_I, v = v_R + i*v_I :
    !    (u,u) = ||z_R||^2 + ||z_I||^2
    !    (v,v) = ||v_R||^2 + ||v_I||^2
    !    Re(u,v) = dot(z_R,v_R) + dot(z_I,v_I)
    ! ==================================================================
    print*
    write(iunit,*)
    print*, "--- ARPACK IRAM completed ---"
    write(iunit,*) "--- ARPACK IRAM completed ---"
    print*, "NB_loc =", NB_loc, "  delta_t =", delta_t, "  T_prop =", T_prop
    write(iunit,*) "NB_loc =", NB_loc, "  delta_t =", delta_t, "  T_prop =", T_prop
    print*
    write(iunit,*)
    print*, "  (j) | Re(mu)      | Im(mu)      | |mu|        |", &
            " ||r||/||Pv|| | Re(lambda)  | Im(lambda)  "
    write(iunit,*) "  (j) | Re(mu)      | Im(mu)      | |mu|        |", &
            " ||r||/||Pv|| | Re(lambda)  | Im(lambda)  "
 
    first = .true. ! eigenvalues comes by pair, first =.true. 
                   ! if we are dealing with the first of the two element of a pair
    do j = 1, nconv
 
        if (di(j) == 0.0d0) then
            ! -----------------------------------------------------------
            !  Real eigenvalue : 1 propagator call
            !  u = z(:,j)  (real vector)
            !  v = P*u stored in ax
            !
            !  (u,u) = ||u||^2 = norm_u^2
            !  (v,v) = ||v||^2 = norm_Pv1^2
            !  Re(u,v) = dot(u,v) = uv_real
            !
            !  sigma and omega must be computed before ax is modified
            !  by daxpy (residual calculation).
            ! -----------------------------------------------------------
            call apply_propagator(z(1,j), ax, n, LinNonLinTimeStep_ptr)
            norm_Pv1 = dnrm2(n, ax, 1)
 

 
            ! --- spectral residual ---
            call daxpy(n, -dr(j), z(1,j), 1, ax, 1)
            spectral_resid = dnrm2(n, ax, 1) / norm_Pv1
 
            mu_abs   = abs(dr(j))
            lambda_r = log(mu_abs) / T_prop
            lambda_i = 0.0d0
 
            print '(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES12.5,A,ES12.5)', &
                "  mode ", j,                                                        &
                " | ", dr(j),                                                        &
                " | ", di(j),                                                        &
                " | ", mu_abs,                                                       &
                " | ", spectral_resid,                                               &
                " | ", lambda_r,                                                     &
                " | ", lambda_i                                                     
            write(iunit,'(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES12.5,A,ES12.5)') &
                "  mode ", j,                                                        &
                " | ", dr(j),                                                        &
                " | ", di(j),                                                        &
                " | ", mu_abs,                                                       &
                " | ", spectral_resid,                                               &
                " | ", lambda_r,                                                     &
                " | ", lambda_i                                                     
 
 
        else if (first) then
            ! -----------------------------------------------------------
            !  Complex conjugate pair : 2 propagator calls for the pair
            !
            !  u = z(:,j) + i*z(:,j+1)   (z_R + i*z_I)
            !  v = P*z(:,j) + i*P*z(:,j+1)   (v_R + i*v_I)
            !
            !  Residual (existing logic, unchanged):
            !    real part = P*z_R - dr(j)*z_R + di(j)*z_I
            !    imag part = P*z_I - di(j)*z_R - dr(j)*z_I
            !
            !  Hermitian inner products (computed before each daxpy):
            !    After first  apply_propagator : ax = v_R = P*z_R
            !      uv_real  = dot(z_R, v_R)
            !    After second apply_propagator : ax = v_I = P*z_I
            !      uv_real += dot(z_I, v_I)
            !
            !  Then:
            !    (u,u) = ||z_R||^2 + ||z_I||^2
            !    (v,v) = ||v_R||^2 + ||v_I||^2
            !    Re(u,v) = uv_real
            !
            !  Conjugate pair has the same sigma, same |omega|, opposite sign
            !  for Im(lambda) and omega -> printed directly for j+1.
            ! -----------------------------------------------------------
 
            ! Real part of residual
            call apply_propagator(z(1,j), ax, n, LinNonLinTimeStep_ptr)
            norm_Pv1 = dnrm2(n, ax, 1)
 
            ! dot products for Re(u,v) while ax = v_R is intact
            uv_real =  ddot(n, z(1,j), 1, ax, 1)       ! dot(z_R, v_R)
 
            call daxpy(n, -dr(j), z(1,j),   1, ax, 1)
            call daxpy(n,  di(j), z(1,j+1), 1, ax, 1)
            norm_r1 = dnrm2(n, ax, 1)
 
            ! Imag part of residual
            call apply_propagator(z(1,j+1), ax, n, LinNonLinTimeStep_ptr)
            norm_Pv2 = dnrm2(n, ax, 1)
 
            ! dot products for Re(u,v) while ax = v_I is intact
            uv_real = uv_real + ddot(n, z(1,j+1), 1, ax, 1)   ! += dot(z_I, v_I)
 
            call daxpy(n, -di(j), z(1,j),   1, ax, 1)
            call daxpy(n, -dr(j), z(1,j+1), 1, ax, 1)
            norm_r2 = dnrm2(n, ax, 1)
 
            ! Combined complex residual
            spectral_resid = dlapy2(norm_r1, norm_r2) / dlapy2(norm_Pv1, norm_Pv2)
 

 
            ! Physical eigenvalue
            mu_abs   = dlapy2(dr(j), di(j))
            mu_arg   = atan2(di(j), dr(j))
            lambda_r = log(mu_abs) / T_prop
            lambda_i = mu_arg / T_prop
 
            ! Print mode j
            print '(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES12.5,A,ES12.5)', &
                "  mode ", j,                                                        &
                " | ", dr(j),                                                        &
                " | ", di(j),                                                        &
                " | ", mu_abs,                                                       &
                " | ", spectral_resid,                                               &
                " | ", lambda_r,                                                     &
                " | ", lambda_i                                                     
            write(iunit,'(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES12.5,A,ES12.5)') &
                "  mode ", j,                                                        &
                " | ", dr(j),                                                        &
                " | ", di(j),                                                        &
                " | ", mu_abs,                                                       &
                " | ", spectral_resid,                                               &
                " | ", lambda_r,                                                     &
                " | ", lambda_i                                                     
 
            ! Print mode j+1 (conjugate : same sigma, opposite omega)
            print '(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)', &
                "  mode ", j+1,                                                      &
                " | ", dr(j+1),                                                      &
                " | ", di(j+1),                                                      &
                " | ", mu_abs,                                                       &
                " | ", spectral_resid,                                               &
                " | ", lambda_r,                                                     &
                " | ", -lambda_i                                                    
            write(iunit,'(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)') &
                "  mode ", j+1,                                                      &
                " | ", dr(j+1),                                                      &
                " | ", di(j+1),                                                      &
                " | ", mu_abs,                                                       &
                " | ", spectral_resid,                                               &
                " | ", lambda_r,                                                     &
                " | ", -lambda_i                                                    
 
            first = .false.
 
        else
            first = .true.
 
        end if
 
    end do
    flush(6)
    flush(iunit)
 
    ! ==================================================================
    !  Reconstruct real and imaginary parts of the physical modes
    ! ==================================================================
    mode_R = 0.0d0
    mode_I = 0.0d0
 
    i = 1
    do while (i <= nconv)
        if (di(i) == 0.0d0) then
            mode_R(:, i) = z(:, i)
            i = i + 1
        else
            ! Complex conjugate pair : x = z(:,i) +/- i*z(:,i+1)
            mode_R(:, i)   =  z(:, i)
            mode_R(:, i+1) =  z(:, i)
            mode_I(:, i)   =  z(:, i+1)
            mode_I(:, i+1) = -z(:, i+1)
            i = i + 2
        end if
    end do
 
    ! ==================================================================
    !  Print summary
    ! ==================================================================
    print*
    write(iunit,*)
    print*, " ====== ARPACK summary ======"
    write(iunit,*) " ====== ARPACK summary ======"
    print*, " Size of the matrix is               ", n
    write(iunit,*) " Size of the matrix is               ", n
    print*, " Number of Ritz values requested      ", nev
    write(iunit,*) " Number of Ritz values requested      ", nev
    print*, " Number of Arnoldi vectors (NCV)      ", ncv
    write(iunit,*) " Number of Arnoldi vectors (NCV)      ", ncv
    print*, " Portion of the spectrum: which =     ", which
    write(iunit,*) " Portion of the spectrum: which =     ", which
    print*, " Number of converged Ritz values      ", nconv
    write(iunit,*) " Number of converged Ritz values      ", nconv
    print*, " Number of OP*x                       ", iparam(9)
    write(iunit,*) " Number of OP*x                       ", iparam(9)
    print*, " Number of implicit restarts           ", restart_count
    write(iunit,*) " Number of implicit restarts           ", restart_count
    print*, " Convergence tolerance                 ", tol
    write(iunit,*) " Convergence tolerance                 ", tol
    print*
    write(iunit,*)
    flush(6)
    flush(iunit)
 

    ! Write Dim_mode.b once (current grid): shared by all Restart_mode_*.b files.
    ! Dedicated name (NOT Dim.b) so the input Dim.b, which travels with the input
    ! Restart.b, is never overwritten. Same record structure as writeDim.
    open(unit=72, file=trim(directory)//"/Dim_mode.b", form="unformatted")
    write(72) KK
    write(72) LL
    write(72) MM
    write(72) shtns%nlm
    write(72) mres
    close(72)
 
    j = 1
    do while (j <= nconv)
                ! On fait des calculs avec act_lin_non_lin sur des vecteurs de taille n
        ! c'est à dire de taille 2*(2KK2+KK4)*shtns%nlm, ce sont des vecteurs 1D réels représentants
        ! des vecteurs 2D réels/1D complexes.
        ! Or on obtient alors des vecteurs propres de act_lin_non_lin qui sont complexes conjugués.
        ! Un vecteur propre a donc une partie réelle de taille n, et une partie imaginaire de taille n
        ! Mais ce n'est pas grave, act_lin_non_lin c'est un certain opérateur L :
        ! L=exp(NB*deltaT*A) et A et l'opérateur qui nous intéresse.
        ! Or alors, les vecteurs propres de L et A sont les mêmes, et les valeurs propres de L 
        ! sont simplements les valeurs propres de A auxquelles on applique exp(Nb*deltaT*...)
        ! On va donc considérer que la partie réelle et la partie imaginaire des vecteurs propres sont des
        ! solutions à part entières
        if (di(j) == 0.0d0) then
            ! Real eigenvalue: only a real eigenvector, no imaginary part to output
            step = j
            call output_vector(mode_R(1:2*KK2*shtns%nlm, j),                              &
                               mode_R(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_R(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j))
            call write_mode_restart(mode_R(1:2*KK2*shtns%nlm, j),                              &
                               mode_R(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_R(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j), j)
            j = j + 1
        else
            ! Complex conjugate pair: output real part (mode j) and imaginary part (mode j+1)
            step = j
            call output_vector(mode_R(1:2*KK2*shtns%nlm, j),                              &
                               mode_R(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_R(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j))
            call write_mode_restart(mode_R(1:2*KK2*shtns%nlm, j),                              &
                               mode_R(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_R(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j), j)

            step = j + 1
            call output_vector(mode_I(1:2*KK2*shtns%nlm, j),                              &
                               mode_I(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_I(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j))
            call write_mode_restart(mode_I(1:2*KK2*shtns%nlm, j),                              &
                               mode_I(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_I(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j), j + 1)
            j = j + 2
        end if
    end do
 
    close(iunit)

end subroutine arnoldi_convective_solver
 
end module mod_ArnoldiSolver