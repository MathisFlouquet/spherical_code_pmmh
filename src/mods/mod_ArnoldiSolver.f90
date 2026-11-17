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

    ! Number of linearised time steps in one application of the propagator:
    ! P = exp(A * T_prop) with T_prop = NB_loc * delta_t. Single definition,
    ! used both by apply_propagator and to recover lambda = ln(mu) / T_prop.
    ! T_prop_target (option -T_prop) and nev (option -nev) are in mod_GlobalVars.
    integer :: NB_loc = 0   ! = nint(T_prop_target / delta_t), set at solver start

contains
 
 
 
subroutine output_vector(E_vec_ptr,F_vec_ptr,T_vec_ptr, ekin_m)
 
    double precision, dimension(2 * KK2 * shtns%nlm), target :: E_vec_ptr 
    double precision, dimension(2 * KK4 * shtns%nlm), target :: F_vec_ptr
    double precision, dimension(2 * KK2 * shtns%nlm), target :: T_vec_ptr 
    ! Optional: kinetic energy of the field per azimuthal wavenumber m (index m/mres)
    double precision, dimension(0:MM), intent(out), optional :: ekin_m
 
    double complex, dimension(:,:), pointer :: E_vec
    double complex, dimension(:,:), pointer :: F_vec
    double complex, dimension(:,:), pointer :: T_vec
    double complex, dimension(:,:), allocatable :: T_save
 
    call c_f_pointer(c_loc(E_vec_ptr), E_vec, [KK2, shtns%nlm])
    call c_f_pointer(c_loc(F_vec_ptr), F_vec, [KK4, shtns%nlm])
    call c_f_pointer(c_loc(T_vec_ptr), T_vec, [KK2, shtns%nlm])
 
 
    call comp_U_from_EF(E_vec, F_vec, Ur, Ut, Up)
    call ToReal(T_vec, T_real, KK2)
 
    ! Output_files computes T_k_spectral_<step>.dat from the GLOBAL spectral T,
    ! not from its T_real argument: point T to the field being written, then
    ! restore it (otherwise the modes' T_k_spectral files show the base state).
    T_save = T
    T = T_vec
    call Output_files(Ur, Up, T_real, step)
    T = T_save
 
    if (present(ekin_m)) call comp_spectral_KE(ekin_m)   ! from Ur, Ut, Up above
 
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
!  one_step_lambda : eigenvalue seen by ONE linearised time step, for the
!  complex Ritz vector u = zR + i zI :
!     alpha = <u, (Phi - I) u> / <u, u>,   lambda_1 = ln(1 + alpha) / delta_t
!  ln(mu)/T_prop only gives Im(lambda) modulo 2*pi/T_prop (aliasing as soon
!  as |Im lambda| > pi/T_prop). lambda_1 is not aliased (|Im lambda| dt << pi)
!  and is used to pick the right branch.
! ===========================================================================
subroutine one_step_lambda(zR, zI, n_prop, LinNonLinTimeStep_ptr, l1r, l1i)

    integer,          intent(in)  :: n_prop
    double precision, intent(in)  :: zR(n_prop), zI(n_prop)
    procedure(),      pointer     :: LinNonLinTimeStep_ptr
    double precision, intent(out) :: l1r, l1i

    double precision, allocatable :: PR(:), PI(:)
    double precision :: uu, ar, ai, mr, mi

    allocate(PR(n_prop), PI(n_prop))
    PR = zR
    call act_lin_non_lin(LinNonLinTimeStep_ptr, PR(1), PR(2*KK2*shtns%nlm + 1), &
                         PR(2*(KK2+KK4)*shtns%nlm + 1))
    PI = zI
    call act_lin_non_lin(LinNonLinTimeStep_ptr, PI(1), PI(2*KK2*shtns%nlm + 1), &
                         PI(2*(KK2+KK4)*shtns%nlm + 1))
    PR = PR - zR          ! (Phi - I) zR
    PI = PI - zI          ! (Phi - I) zI

    uu = dot_product(zR, zR) + dot_product(zI, zI)
    ar = (dot_product(zR, PR) + dot_product(zI, PI)) / uu
    ai = (dot_product(zR, PI) - dot_product(zI, PR)) / uu
    mr = 1.d0 + ar
    mi = ai
    l1r = log(sqrt(mr*mr + mi*mi)) / delta_t
    l1i = atan2(mi, mr) / delta_t

    deallocate(PR, PI)

end subroutine one_step_lambda

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
    integer          :: step

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
    integer            :: n, ncv, ldv, lworkl, maxiter   ! nev: global (option -nev)
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
    double precision :: T_prop                      ! = NB_loc * delta_t
    double precision :: mu_abs, mu_arg
    double precision :: lambda_r, lambda_i
    double precision :: lam1_r, lam1_i              ! one-step (non aliased) estimate
    integer          :: k_alias                     ! 2*pi/T_prop branch of Im(lambda)

    ! ----- diagnostics -----
    integer(8)       :: clk_start, clk_end, clk_rate
    double precision :: wall_time                   ! wall-clock time of the dnaupd loop (s)
    double precision, allocatable :: rz_r(:), rz_i(:), rz_b(:), rz_a(:)  ! last Ritz values (ncv)
    logical,          allocatable :: rz_done(:)
    integer          :: n_ritz_out
    integer(8)       :: clk_solver                  ! start of the solver (total wall time)

    ! ----- results kept for the summary file (index j = Ritz value) -----
    double precision, allocatable :: ev_lr(:), ev_li(:), ev_mu(:), ev_res(:)   ! lambda, |mu|, residual
    double precision, allocatable :: ev_l1r(:), ev_l1i(:)  ! one-step estimate (complex pairs)
    integer,          allocatable :: ev_kind(:)    ! 0 real, 1 first of a pair, 2 conjugate
    integer,          allocatable :: ev_kal(:)     ! alias branch k (complex pairs)
    double precision, allocatable :: ev_ke(:,:)    ! kinetic energy per m (0:MM, j)
    double precision, allocatable :: ke_tmp(:)     ! (0:MM)
    character(len=100) :: run_dir                  ! -directory of the run
    character(len=32)  :: status_msg
    logical          :: conductive_base            ! base state identically zero
 
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
    !  Outputs: progress goes to the console (run.log); the summary
    !  arnoldi_output.txt is written at the end (write_summary), the
    !  eigenvectors in the sub-folder Eigenvectors/ of the run directory.
    ! ==================================================================
    iunit = 42
    run_dir = directory
    call system_clock(clk_solver)

    ! ==================================================================
    !  ARPACK parameters
    ! ==================================================================
    n       = 2*(2*KK2+KK4)*shtns%nlm
    ! nev (option -nev, default 2 = one complex pair) is the number of
    ! eigenvalues sought. ncv and maxiter are derived from it.
    !
    ! Krylov size ncv: ARPACK needs ncv >= nev + 2 and advises ncv >= 2*nev.
    ! Cost of the search ~ max(ncv + 1, matvecs actually needed) propagator
    ! calls: a large ncv is wasted on an easy problem, a small one may need
    ! many restarts (each costs ncv - nev calls) when eigenvalues are close.
    ! Tests (Ra = 60, m = 4, conductive state, nev = 2, docs/PARAMETRES.md):
    ! ncv = 6 to 10 was optimal; 3*nev + 4 (= 10 for nev = 2) keeps a margin
    ! for harder spectra.
    !
    ! maxiter = maximum number of implicit restarts, only reached when the
    ! search fails; generous on purpose.
    !
    ! tol is relative to |mu|; the error on lambda is ~ tol / T_prop.
    ncv     = 3*nev + 4
    maxiter = 150*nev
    tol     = 1.0d-9
    print*, "Arnoldi parameters: nev =", nev, " ncv =", ncv, " maxiter =", maxiter, &
            " tol =", tol, " T_prop target =", T_prop_target
    bmat    = 'I'
    which   = 'LM'
    ldv     = n
    lworkl  = 3*ncv**2 + 6*ncv
 
    ! Propagation time: NB_loc steps of delta_t, T_prop ~ T_prop_target
    NB_loc = max(1, nint(T_prop_target / delta_t))
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
    allocate(ev_lr(nev+1), ev_li(nev+1), ev_mu(nev+1), ev_res(nev+1))
    allocate(ev_l1r(nev+1), ev_l1i(nev+1), ev_kind(nev+1), ev_kal(nev+1))
    allocate(ev_ke(0:MM, nev+1), ke_tmp(0:MM))
    ev_lr = 0.d0; ev_li = 0.d0; ev_mu = 0.d0; ev_res = 0.d0
    ev_l1r = 0.d0; ev_l1i = 0.d0; ev_kind = 0; ev_kal = 0; ev_ke = 0.d0
    nconv = 0
    matvec_count = 0
    restart_count = 0
    wall_time = 0.d0
 
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
                call precompXeYe()
                print*, "Precomputing the Xf and Yf matrices..."
                call precompXfYf()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT()
            case ("fbe")
                print*, "Precomputing the Xe and Ye matrices..."
                call precompXeYe_BDF()
                print*, "Precomputing the Xf and Yf matrices..."
                call precompXfYf_BDF()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT_BDF()
        end select
 
        LinNonLinTimeStep_ptr => compute_lin_time_step_convective_explicit_CN_FBE
 
    else if (solver == "arnoldi_convective_implicit") then
 
        select case (time_step)
            case ("cn")
                print*, "Precomputing the Xef and Yef matrices..."
                call PrecompimplicitXY()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT()
                LinNonLinTimeStep_ptr => compute_lin_time_step_convective_implicit_CN
            case ("fbe")
                print*, "Precomputing the Xef, Ye and Yf matrices..."
                call PrecompimplicitXY_BDF()
                print*, "Precomputing the XT and YT matrices..."
                call precompXTYT_BDF()
                LinNonLinTimeStep_ptr => compute_lin_time_step_convective_implicit_FBE
        end select
 
    end if
 
    print*
    print*, "--------------------- Initial values ---------------------------"
    print*
    print*, "Volume of the shell = ", Vol
    print*, "Reading restart files"
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
    
    ! Reported in the summary: catches a forgotten conductive-state hardcode
    conductive_base = (maxval(abs(E_base)) == 0.d0) .and. (maxval(abs(F_base)) == 0.d0) &
                      .and. (maxval(abs(T_base)) == 0.d0)

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
    call system_clock(clk_start, clk_rate)
 
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
                print*, "    matvecs : ", matvec_count
                flush(6)
            else if (matvec_count > ncv) then
                if (mod(matvec_count - ncv, ncv - nev) == 0) then ! To monitor the number of restart
                    restart_count = restart_count + 1
                    print*, "--- Restart ", restart_count, " / ", maxiter,  &
                            " | matvecs : ", matvec_count,                  &
                            " | converged : ", iparam(5), " / ", nev
                    flush(6)
                end if
            end if
            
            print*, "  wavespeed :", C_base
        else
            exit
        end if
 
    end do
 
    ! ==================================================================
    !  Check dnaupd return code
    ! ==================================================================
    call system_clock(clk_end)
    wall_time = dble(clk_end - clk_start) / dble(clk_rate)
    if (info_arpack < 0) then
        print*, "ARPACK dnaupd error: info =", info_arpack
        write(status_msg,'(A,I0)') "dnaupd error, info = ", info_arpack
        call write_summary(status_msg)
        stop
    else if (info_arpack == 1) then
        print*, "Warning: dnaupd reached maxiter without full convergence."
        print*, "Converged Ritz values:", iparam(5), "/", nev
    else
        print*, "dnaupd converged. Ritz values computed:", iparam(5)
    end if
    print*, "Total matvecs  :", matvec_count
    print*, "Total restarts :", restart_count
    print*, "Wall time of the Arnoldi search (s) :", wall_time
    flush(6)

    ! ------------------------------------------------------------------
    !  Last Ritz values of the Hessenberg matrix (converged or not), from
    !  workl(ipntr(6:8)), sorted by decreasing |mu|. Useful when the run
    !  did not converge: shows what the search was heading to and how far
    !  (Ritz estimate / |mu|). Im(lambda) here is the principal branch
    !  arg(mu)/T_prop, i.e. known only modulo 2*pi/T_prop.
    ! ------------------------------------------------------------------
    allocate(rz_r(ncv), rz_i(ncv), rz_b(ncv), rz_a(ncv), rz_done(ncv))
    rz_r = workl(ipntr(6):ipntr(6)+ncv-1)
    rz_i = workl(ipntr(7):ipntr(7)+ncv-1)
    rz_b = workl(ipntr(8):ipntr(8)+ncv-1)
    rz_a = sqrt(rz_r**2 + rz_i**2)
    rz_done = .false.
    n_ritz_out = min(ncv, 6)
    print*, "Last Ritz values (|mu| descending): |mu|, Re(lambda), Im(lambda) mod 2pi/T, estimate/|mu|"
    do i = 1, n_ritz_out
        j = maxloc(rz_a, 1, mask = .not. rz_done)
        rz_done(j) = .true.
        print '(I4,4ES14.5)', i, rz_a(j), log(rz_a(j)) / T_prop, &
              atan2(rz_i(j), rz_r(j)) / T_prop, rz_b(j) / rz_a(j)
    end do
    flush(6)
 
 
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
        write(status_msg,'(A,I0)') "dneupd error, info = ", info_arpack
        call write_summary(status_msg)
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
    print*, "--- ARPACK IRAM completed ---"
    print*, "NB_loc =", NB_loc, "  delta_t =", delta_t, "  T_prop =", T_prop
    print*
    print*, "  (j) | Re(mu)      | Im(mu)      | |mu|        |", &
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
            ev_kind(j) = 0
            ev_lr(j) = lambda_r;  ev_li(j) = lambda_i
            ev_mu(j) = mu_abs;    ev_res(j) = spectral_resid
 
 
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

            ! Remove the 2*pi/T_prop aliasing of Im(lambda) using one time step
            call one_step_lambda(z(1,j), z(1,j+1), n, LinNonLinTimeStep_ptr, lam1_r, lam1_i)
            k_alias  = nint((lam1_i * T_prop - mu_arg) / (2.d0 * acos(-1.d0)))
            lambda_i = (mu_arg + 2.d0 * acos(-1.d0) * dble(k_alias)) / T_prop
            print '(A,ES12.5,A,ES12.5,A,I0)', "  one-step lambda = ", lam1_r, " + i ", lam1_i, &
                "   -> alias branch k = ", k_alias
 
            ! Print mode j
            print '(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES12.5,A,ES12.5)', &
                "  mode ", j,                                                        &
                " | ", dr(j),                                                        &
                " | ", di(j),                                                        &
                " | ", mu_abs,                                                       &
                " | ", spectral_resid,                                               &
                " | ", lambda_r,                                                     &
                " | ", lambda_i                                                     
            ev_kind(j) = 1
            ev_lr(j) = lambda_r;  ev_li(j) = lambda_i
            ev_mu(j) = mu_abs;    ev_res(j) = spectral_resid
            ev_l1r(j) = lam1_r;   ev_l1i(j) = lam1_i;   ev_kal(j) = k_alias
 
            ! Print mode j+1 (conjugate : same sigma, opposite omega)
            print '(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES12.5)', &
                "  mode ", j+1,                                                      &
                " | ", dr(j+1),                                                      &
                " | ", di(j+1),                                                      &
                " | ", mu_abs,                                                       &
                " | ", spectral_resid,                                               &
                " | ", lambda_r,                                                     &
                " | ", -lambda_i                                                    
            if (j + 1 <= nconv) then
                ev_kind(j+1) = 2
                ev_lr(j+1) = lambda_r;  ev_li(j+1) = -lambda_i
                ev_mu(j+1) = mu_abs;    ev_res(j+1) = spectral_resid
            end if
 
            first = .false.
 
        else
            first = .true.
 
        end if
 
    end do
    flush(6)
 
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
    print*, " ====== ARPACK summary ======"
    print*, " Size of the matrix is               ", n
    print*, " Number of Ritz values requested      ", nev
    print*, " Number of Arnoldi vectors (NCV)      ", ncv
    print*, " Portion of the spectrum: which =     ", which
    print*, " Number of converged Ritz values      ", nconv
    print*, " Number of OP*x                       ", iparam(9)
    print*, " Number of implicit restarts           ", restart_count
    print*, " Convergence tolerance                 ", tol
    print*
    flush(6)
 

    ! ==================================================================
    !  Eigenvectors -> sub-folder Eigenvectors/ of the run directory:
    !  Restart_mode_<j>.b + Dim_mode.b (readable as restart) and the
    !  visualisation files Data_*_<j>.dat, Ekin_spectral_<j>.dat,
    !  T_k_spectral_<j>.dat, r/theta/phi.dat. The writing routines use the
    !  global "directory", which is redirected here and restored after.
    !  For a complex pair, j = real part and j+1 = imaginary part.
    ! ==================================================================
    if (len_trim(run_dir) + len("/Eigenvectors") <= len(directory)) then
        directory = trim(run_dir)//"/Eigenvectors"
        call execute_command_line('mkdir -p "'//trim(directory)//'"')
    else
        print*, "Warning: directory name too long, eigenvectors written in ", trim(run_dir)
    end if
    call output_coordinates()

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
        ! The eigenvectors of P = exp(T_prop A) are those of A. A complex
        ! eigenvector u = z_R + i z_I is stored as two real fields (real part
        ! j, imaginary part j+1); the physical perturbation is
        ! Re(u exp(lambda t)) = exp(sigma t) (z_R cos(omega t) - z_I sin(omega t)).
        if (di(j) == 0.0d0) then
            ! Real eigenvalue: only a real eigenvector, no imaginary part to output
            step = j
            call output_vector(mode_R(1:2*KK2*shtns%nlm, j),                              &
                               mode_R(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_R(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j), &
                               ev_ke(:, j))
            call write_mode_restart(mode_R(1:2*KK2*shtns%nlm, j),                              &
                               mode_R(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_R(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j), j)
            j = j + 1
        else
            ! Complex conjugate pair: output real part (mode j) and imaginary part (mode j+1)
            step = j
            call output_vector(mode_R(1:2*KK2*shtns%nlm, j),                              &
                               mode_R(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_R(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j), &
                               ev_ke(:, j))
            call write_mode_restart(mode_R(1:2*KK2*shtns%nlm, j),                              &
                               mode_R(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_R(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j), j)

            step = j + 1
            call output_vector(mode_I(1:2*KK2*shtns%nlm, j),                              &
                               mode_I(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_I(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j), &
                               ke_tmp)
            call write_mode_restart(mode_I(1:2*KK2*shtns%nlm, j),                              &
                               mode_I(2*KK2*shtns%nlm + 1 : 2*(KK2+KK4)*shtns%nlm, j),    &
                               mode_I(2*(KK2+KK4)*shtns%nlm + 1 : 2*(2*KK2+KK4)*shtns%nlm, j), j + 1)
            ! Energy of the complex mode = real part + imaginary part
            ev_ke(:, j) = ev_ke(:, j) + ke_tmp
            if (j + 1 <= nev + 1) ev_ke(:, j+1) = ev_ke(:, j)
            j = j + 2
        end if
    end do

    directory = run_dir

    if (info_arpack == 0 .and. nconv >= nev) then
        status_msg = "converged"
    else
        status_msg = "NOT fully converged"
    end if
    call write_summary(status_msg)
    print*, "Summary written in ", trim(run_dir)//"/arnoldi_output.txt"

contains

    ! ------------------------------------------------------------------
    !  write_summary : readable report of the run in arnoldi_output.txt
    !  (parameters, search statistics, eigenvalues, azimuthal content of
    !  each eigenvector, files written, search diagnostics).
    ! ------------------------------------------------------------------
    subroutine write_summary(status)
        character(len=*), intent(in) :: status
        character(len=2048) :: cmd
        character(len=8)    :: d8
        character(len=10)   :: t10
        integer(8)          :: clk_now
        integer             :: jj, ii, m_dom
        double precision    :: ke_tot, total_time
        character(len=100)  :: files

        call system_clock(clk_now)
        total_time = dble(clk_now - clk_solver) / dble(clk_rate)
        call date_and_time(d8, t10)
        call get_command(cmd)

        open(unit=iunit, file=trim(run_dir)//'/arnoldi_output.txt', status='replace', action='write')
        write(iunit,'(A)') "======================================================================"
        write(iunit,'(A)') " Arnoldi (ARPACK) - eigenvalues of the linearised operator"
        write(iunit,'(A)') "======================================================================"
        write(iunit,'(A)') " Date    : "//d8(1:4)//"-"//d8(5:6)//"-"//d8(7:8)//" "// &
                           t10(1:2)//":"//t10(3:4)//":"//t10(5:6)
        write(iunit,'(A)') " Command : "//trim(cmd)
        write(iunit,'(A)')
        write(iunit,'(A)') " --- Physics"
        write(iunit,'(A,ES12.5,A,ES12.5,A,ES12.5)') "  Ra = ", Ra, "   Ek = ", Ek, "   Pr = ", Pr
        if (conductive_base) then
            write(iunit,'(A,ES12.5)') "  Base state : zero (conductive state), C_base = ", C_base
        else
            write(iunit,'(A,ES12.5)') "  Base state : "//trim(restart_filename)//", C_base = ", C_base
        end if
        write(iunit,'(A)') " --- Resolution"
        write(iunit,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') "  KK = ", KK, "   LL = ", LL, "   MM = ", MM, &
            "   mres = ", mres, "   (problem size n = ", n, ")"
        write(iunit,'(A)') " --- Method"
        write(iunit,'(A)') "  solver = "//trim(solver)//"   time_step = "//trim(time_step)
        write(iunit,'(A,ES10.3,A,ES10.3,A,I0,A)') "  delta_t = ", delta_t, "   T_prop = ", T_prop, &
            "   (NB_loc = ", NB_loc, " steps per propagator call)"
        write(iunit,'(A,I0,A,I0,A,I0,A,ES9.2,A)') "  nev = ", nev, "   ncv = ", ncv, &
            "   maxiter = ", maxiter, "   tol = ", tol, "   which = "//which
        write(iunit,'(A)')
        write(iunit,'(A)') " --- Search"
        write(iunit,'(A,I0,A,I0,A)') "  Status            : "//trim(status)//" (", nconv, " / ", nev, &
            " converged Ritz values)"
        write(iunit,'(A,I0,A,I0)') "  Propagator calls  : ", matvec_count, "   restarts : ", restart_count
        write(iunit,'(A,F10.1,A,F10.1,A)') "  Wall time         : ", wall_time, " s (search), ", &
            total_time, " s (total)"
        write(iunit,'(A)')
        write(iunit,'(A)') " --- Eigenvalues   lambda = sigma + i omega = ln(mu) / T_prop"
        write(iunit,'(A)') "  Perturbation ~ Re( u exp(lambda t) ); Im(lambda) branch fixed with one time step."
        write(iunit,'(A)') "  m : azimuthal wavenumber carrying most of the kinetic energy of the eigenvector."
        write(iunit,'(A)')
        write(iunit,'(A)') "   j |     sigma       omega     |    |mu|     residual |  m  KE frac | files in Eigenvectors/"
        do jj = 1, nconv
            ke_tot = sum(ev_ke(:, jj))
            if (ke_tot > 0.d0) then
                m_dom = maxloc(ev_ke(:, jj), 1) - 1
            else
                m_dom = 0
            end if
            select case (ev_kind(jj))
                case (0)
                    write(files,'(A,I0,A)') "Restart_mode_", jj, ".b (real eigenvector)"
                case (1)
                    write(files,'(A,I0,A,I0,A)') "Restart_mode_", jj, ".b = real part, _", jj+1, ".b = imag. part"
                case default
                    write(files,'(A,I0)') "complex conjugate of ", jj-1
            end select
            if (ke_tot > 0.d0) then
                write(iunit,'(I4,A,2ES12.4,A,ES11.4,ES10.2,A,I3,F7.3,A)') jj, " | ", ev_lr(jj), ev_li(jj), &
                    " | ", ev_mu(jj), ev_res(jj), " | ", m_dom*mres, ev_ke(m_dom, jj) / ke_tot, &
                    "  | "//trim(files)
            else
                write(iunit,'(I4,A,2ES12.4,A,ES11.4,ES10.2,A)') jj, " | ", ev_lr(jj), ev_li(jj), &
                    " | ", ev_mu(jj), ev_res(jj), " |     -      | "//trim(files)
            end if
        end do
        if (nconv == 0) write(iunit,'(A)') "  (no converged eigenvalue)"
        write(iunit,'(A)')
        write(iunit,'(A)') " --- Files"
        write(iunit,'(A)') "  Eigenvectors/Restart_mode_<j>.b + Dim_mode.b : eigenvector j, readable as a restart"
        write(iunit,'(A)') "  Eigenvectors/Data_{eq_plane,mer_plane,mid_gap}_<j>.dat, r/theta/phi.dat : visualisation"
        write(iunit,'(A)') "  Eigenvectors/Ekin_spectral_<j>.dat (KE per m), T_k_spectral_<j>.dat (T per Chebyshev k)"
        write(iunit,'(A)') "  Data_*_0.dat, Ekin_spectral_0.dat, ... in the run directory : base state"
        write(iunit,'(A)')
        write(iunit,'(A)') " --- Diagnostics"
        do jj = 1, nconv
            if (ev_kind(jj) == 1) write(iunit,'(A,I0,A,ES12.5,A,ES12.5,A,I0)') "  pair ", jj, &
                ": one-step lambda = ", ev_l1r(jj), " + i ", ev_l1i(jj), "   alias branch k = ", ev_kal(jj)
        end do
        write(iunit,'(A)') "  Last Ritz values of the search (|mu| descending; Im(lambda) modulo 2 pi / T_prop):"
        write(iunit,'(A)') "     i |    |mu|      sigma       omega mod   estimate/|mu|"
        if (allocated(rz_a)) then
            rz_done = .false.
            do ii = 1, min(ncv, 6)
                jj = maxloc(rz_a, 1, mask = .not. rz_done)
                rz_done(jj) = .true.
                write(iunit,'(I6,A,4ES12.4)') ii, " | ", rz_a(jj), log(rz_a(jj)) / T_prop, &
                    atan2(rz_i(jj), rz_r(jj)) / T_prop, rz_b(jj) / rz_a(jj)
            end do
        end if
        close(iunit)
    end subroutine write_summary

end subroutine arnoldi_convective_solver
 
end module mod_ArnoldiSolver