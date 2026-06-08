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

    integer, parameter :: NEV    = 20
    integer, parameter :: NCV    = 60
    integer, parameter :: NB_LOC = 100

contains

! ===========================================================================
!  randn : génère un scalaire gaussien (Box-Muller)
! ===========================================================================
function randn() result(z)
    double precision :: z
    double precision :: u1, u2
    call random_number(u1)
    call random_number(u2)
    if (u1 < 1.0d-12) u1 = 1.0d-12
    z = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * acos(-1.0d0) * u2)
end function randn

! ===========================================================================
!  arnoldi_convective_solver
! ===========================================================================
subroutine arnoldi_convective_solver()

    ! --- taille de l'espace ---
    integer :: sz

    ! --- Base de Krylov et matrice de Hessenberg ---
    double precision, allocatable :: V(:,:)      ! (sz x NCV+1)
    double precision, allocatable :: H(:,:)      ! (NCV+1 x NCV)
    double precision, allocatable :: H_sq(:,:)   ! (NCV x NCV)  copie pour dgeev
    double precision, allocatable :: w(:)        ! (sz)

    ! --- dgeev ---
    double precision, allocatable :: WR(:)       ! (NCV) parties réelles des vp
    double precision, allocatable :: WI(:)       ! (NCV) parties imag des vp
    double precision, allocatable :: EVL(:,:)    ! (NCV x NCV) vp gauches (non utilisés)
    double precision, allocatable :: EVR(:,:)    ! (NCV x NCV) vp droits
    double precision, allocatable :: WORK(:)
    integer :: info_dgeev

    ! --- Sélection des NEV plus grands modules ---
    double precision, allocatable :: mu_abs(:)   ! (NCV)
    integer,          allocatable :: idx(:)      ! (NCV) indices triés
    integer  :: tmp_idx
    double precision :: tmp_abs

    ! --- Vecteurs de Ritz ---
    double precision, allocatable :: ritz_R(:)   ! (sz)
    double precision, allocatable :: ritz_I(:)   ! (sz)

    ! --- Résidu spectral ---
    double precision, allocatable :: Pz_R(:)     ! (sz)
    double precision, allocatable :: Pz_I(:)     ! (sz)
    double precision, allocatable :: res_R(:)    ! (sz)
    double precision, allocatable :: res_I(:)    ! (sz)
    double precision :: norm_Pz, norm_res, spectral_resid

    ! --- Valeurs propres physiques ---
    double precision :: T_prop, mu_r, mu_i, mu_mod, mu_arg
    double precision :: lambda_r, lambda_i

    ! --- Pointeur sur le solveur linéaire ---
    procedure(), pointer :: LinNonLinTimeStep_ptr

    ! --- Divers ---
    integer :: i, j, k, s, rank

    ! ==================================================================
    !  Taille de l'espace et allocations
    ! ==================================================================
    sz     = 2*(2*KK2+KK4)*shtns%nlm
    T_prop = dble(NB_LOC) * delta_t

    allocate(V(sz, NCV+1))
    allocate(H(NCV+1, NCV))
    allocate(H_sq(NCV, NCV))
    allocate(w(sz))
    allocate(WR(NCV), WI(NCV))
    allocate(EVL(NCV, NCV), EVR(NCV, NCV))
    allocate(WORK(4*NCV))
    allocate(mu_abs(NCV), idx(NCV))
    allocate(ritz_R(sz), ritz_I(sz))
    allocate(Pz_R(sz), Pz_I(sz))
    allocate(res_R(sz), res_I(sz))

    ! ==================================================================
    !  Préparations : matrices, lecture, état de base
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

    E_base = 0;  F_base = 0;  T_base = 0
    E = 0;       F = 0;       T = 0
    call comp_radial_derivatives(E_base, F_base, T_base, Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St)
    call comp_U_cU_gT(Ur_base, Ut_base, Up_base, cUr_base, cUt_base, cUp_base, &
                      gTr_base, gTt_base, gTp_base)

    ! ==================================================================
    !  1. Vecteur initial aléatoire normalisé -> V(:,1)
    ! ==================================================================
    V = 0.0d0
    H = 0.0d0

    do i = 1, sz
        V(i,1) = randn()
    end do
    V(:,1) = V(:,1) / sqrt(dot_product(V(:,1), V(:,1)))

    ! ==================================================================
    !  2. Itérations d'Arnoldi (Gram-Schmidt modifié)
    !
    !  À chaque pas j :
    !    w = P * V(:,j)   avec P = NB_LOC applications de act_lin_non_lin
    !    pour i=1..j : H(i,j) = <V(:,i), w>   et   w -= H(i,j)*V(:,i)
    !    H(j+1,j) = ||w||
    !    V(:,j+1) = w / H(j+1,j)
    ! ==================================================================
    print*
    print*, "--- Arnoldi iterations ---"

    do j = 1, NCV

        ! --- Propagation : NB_LOC pas de time step ---
        w = V(:,j)
        do s = 1, NB_LOC
            call act_lin_non_lin(LinNonLinTimeStep_ptr,        &
                                 w(1),                         &
                                 w(2*KK2*shtns%nlm + 1),       &
                                 w(2*(KK2+KK4)*shtns%nlm + 1))
        end do

        ! --- Gram-Schmidt modifié ---
        do i = 1, j
            H(i,j) = dot_product(V(:,i), w)
            w = w - H(i,j) * V(:,i)
        end do

        H(j+1,j) = sqrt(dot_product(w, w))

        if (j < NCV) then
            V(:,j+1) = w / H(j+1,j)
        end if

        print '(A,I3,A,I3,A,ES10.3)', &
            "  Step ", j, " /", NCV, "   beta = ", H(j+1,j)
        flush(6)

    end do

    ! ==================================================================
    !  3. dgeev sur la matrice de Hessenberg carrée H(1:NCV, 1:NCV)
    ! ==================================================================
    H_sq = H(1:NCV, 1:NCV)

    call dgeev('N', 'V', NCV, H_sq, NCV, WR, WI, EVL, NCV, EVR, NCV, &
               WORK, 4*NCV, info_dgeev)

    if (info_dgeev /= 0) then
        print*, "dgeev failed, info =", info_dgeev
        stop
    end if

    ! ==================================================================
    !  4. Tri des valeurs propres par |mu| décroissant
    ! ==================================================================
    do i = 1, NCV
        mu_abs(i) = sqrt(WR(i)**2 + WI(i)**2)
        idx(i)    = i
    end do

    do i = 2, NCV
        tmp_idx = idx(i)
        tmp_abs = mu_abs(tmp_idx)
        k = i - 1
        do while (k >= 1 .and. mu_abs(idx(k)) < tmp_abs)
            idx(k+1) = idx(k)
            k = k - 1
        end do
        idx(k+1) = tmp_idx
    end do

    ! ==================================================================
    !  5. Résidu spectral + valeurs propres physiques
    !
    !  Pour une valeur propre réelle mu :
    !    z     = V(:,1:NCV) * EVR(:,j)
    !    résidu = ||P*z - mu*z|| / ||P*z||
    !
    !  Pour une paire complexe (mu, conj(mu)), indices (j, j+1) :
    !    z = V*EVR(:,j) + i*V*EVR(:,j+1)
    !    résidu = sqrt(||rés_R||²+||rés_I||²) / sqrt(||Pz_R||²+||Pz_I||²)
    !
    !  lambda = ln(mu) / T_prop   avec T_prop = NB_LOC * delta_t
    ! ==================================================================
    print*
    print*, "--- Arnoldi maison : résultats ---"
    print*, "NB_LOC =", NB_LOC, "   delta_t =", delta_t, "   T_prop =", T_prop
    print*
    print '(A)', "  rank | Re(mu)      | Im(mu)      | |mu|        |" // &
                 " résidu      | Re(lambda)  | Im(lambda)"

    do rank = 1, NEV

        j      = idx(rank)
        mu_r   = WR(j)
        mu_i   = WI(j)
        mu_mod = mu_abs(idx(rank))

        ritz_R = matmul(V(:,1:NCV), EVR(:,j))

        if (mu_i == 0.0d0) then
            ! ----------------------------------------------------------
            !  Valeur propre réelle
            ! ----------------------------------------------------------
            Pz_R = ritz_R
            do s = 1, NB_LOC
                call act_lin_non_lin(LinNonLinTimeStep_ptr,        &
                                     Pz_R(1),                      &
                                     Pz_R(2*KK2*shtns%nlm + 1),    &
                                     Pz_R(2*(KK2+KK4)*shtns%nlm + 1))
            end do

            norm_Pz        = sqrt(dot_product(Pz_R, Pz_R))
            res_R          = Pz_R - mu_r * ritz_R
            norm_res       = sqrt(dot_product(res_R, res_R))
            spectral_resid = norm_res / norm_Pz

            lambda_r = log(mu_mod) / T_prop
            lambda_i = 0.0d0

        else
            ! ----------------------------------------------------------
            !  Paire complexe conjuguée : dgeev place conj(mu) en j+1
            !  avec WI(j) > 0 et EVR(:,j+1) = partie imag du vecteur
            ! ----------------------------------------------------------
            ritz_I = matmul(V(:,1:NCV), EVR(:,j+1))

            ! Propagation de la partie réelle
            Pz_R = ritz_R
            do s = 1, NB_LOC
                call act_lin_non_lin(LinNonLinTimeStep_ptr,        &
                                     Pz_R(1),                      &
                                     Pz_R(2*KK2*shtns%nlm + 1),    &
                                     Pz_R(2*(KK2+KK4)*shtns%nlm + 1))
            end do

            ! Propagation de la partie imaginaire
            Pz_I = ritz_I
            do s = 1, NB_LOC
                call act_lin_non_lin(LinNonLinTimeStep_ptr,        &
                                     Pz_I(1),                      &
                                     Pz_I(2*KK2*shtns%nlm + 1),    &
                                     Pz_I(2*(KK2+KK4)*shtns%nlm + 1))
            end do

            norm_Pz = sqrt(dot_product(Pz_R,Pz_R) + dot_product(Pz_I,Pz_I))

            ! P*z - mu*z :
            !   partie réelle = Pz_R - mu_r*ritz_R + mu_i*ritz_I
            !   partie imag   = Pz_I - mu_i*ritz_R - mu_r*ritz_I
            res_R = Pz_R - mu_r*ritz_R + mu_i*ritz_I
            res_I = Pz_I - mu_i*ritz_R - mu_r*ritz_I

            norm_res       = sqrt(dot_product(res_R,res_R) + dot_product(res_I,res_I))
            spectral_resid = norm_res / norm_Pz

            mu_arg   = atan2(mu_i, mu_r)
            lambda_r = log(mu_mod) / T_prop
            lambda_i = mu_arg     / T_prop

        end if

        print '(A,I3,A,ES12.5,A,ES12.5,A,ES12.5,A,ES10.3,A,ES12.5,A,ES12.5)', &
            "  ", rank,             &
            " | ", mu_r,            &
            " | ", mu_i,            &
            " | ", mu_mod,          &
            " | ", spectral_resid,  &
            " | ", lambda_r,        &
            " | ", lambda_i
        flush(6)

    end do

    ! ==================================================================
    !  Libération mémoire
    ! ==================================================================
    deallocate(V, H, H_sq, w)
    deallocate(WR, WI, EVL, EVR, WORK)
    deallocate(mu_abs, idx)
    deallocate(ritz_R, ritz_I, Pz_R, Pz_I, res_R, res_I)

end subroutine arnoldi_convective_solver

end module mod_ArnoldiSolver