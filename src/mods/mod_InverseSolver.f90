! =============================================================================
!  mod_InverseSolver : méthode de la puissance inverse pour
!        lambda u = (L + N_U) u
!
!  Système résolu à chaque itération externe (par GMRES) :
!        A u_{n+1} = B u_n
!  avec
!        A = [ (I - dt L)^-1 (I + dt N_U) ] - I      -> subA_inv
!        B = (I - dt L)^-1 dt                        -> TermeDroite
!
!  A  = intégration temporelle linéarisée moins l'identité.
!  B  = copie de l'intégration temporelle linéarisée SANS les termes explicites
!       (on enlève l'appel à comp_LinNonLin).
!
!  Les lignes commençant par  !>>>  sont mes propositions, DÉSACTIVÉES.
!  Voir en particulier la PROPOSITION 3 (facteur delta_t) : sans elle, la
!  valeur propre renvoyée vaut delta_t * lambda et le shift effectif vaut
!  sigma_r / delta_t.
! =============================================================================

module mod_InverseSolver

    use iso_c_binding

    use mod_GlobalVars
    use mod_PrecompSH
    use mod_precompXY
    use mod_Matrices
    use mod_read
    use mod_TimeStep
    use mod_Output
    use mod_tools
    use mod_ActNewtonSolver
    use mod_IterativeSolvers
    use mod_ArnoldiSolver ! pour output_vector et write_mode_restart

    implicit none

contains

! =============================================================================
!  subA_inv :  A x = (Phi - I) x
!  Phi = une intégration temporelle linéarisée complète.
!  Pas de condition de phase (contrairement au subA du Newton) : c_per = 0.
! =============================================================================
subroutine subA_inv(time_stepping_sub, EFT_in, EFT_out)

    implicit none

    procedure(), pointer, intent(in) :: time_stepping_sub

    double precision, dimension(2*shtns%nlm*(2*KK2+KK4)), intent(inout) :: EFT_in
    double precision, dimension(2*shtns%nlm*(2*KK2+KK4)), intent(out)   :: EFT_out

    c_per = 0.

    ! On sauvegarde l'entrée
    EFT_out = EFT_in

    ! Une intégration temporelle linéarisée : EFT_in <- Phi * EFT_in (en place)
    call act_lin_non_lin(time_stepping_sub, EFT_in(1), EFT_in(2*shtns%nlm*KK2 + 1), &
                        & EFT_in(2*shtns%nlm*(KK2+KK4) + 1))

    ! (Phi - I) x = Phi x - x
    EFT_out = EFT_in - EFT_out

end subroutine subA_inv

! =============================================================================
!  TermeDroite : copie de compute_lin_time_step_convective_implicit_CN (ou _FBE)
!  avec le SEUL changement suivant : on enlève l'appel à comp_LinNonLin.
! =============================================================================
subroutine TermeDroite(time_stepping_sub, EFT_in, EFT_out)

    implicit none

    procedure(), pointer, intent(in) :: time_stepping_sub

    double precision, dimension(2*shtns%nlm*(2*KK2+KK4)), intent(inout), target :: EFT_in
    double precision, dimension(2*shtns%nlm*(2*KK2+KK4)), intent(out)   :: EFT_out

    integer :: lm, k, m, top, stride

    ! E_per / F_per / T_per : pointeurs globaux utilises par le pas de temps
    call c_f_pointer(c_loc(EFT_in(1)),                         E_per, [KK2, shtns%nlm])
    call c_f_pointer(c_loc(EFT_in(2*shtns%nlm*KK2 + 1)),       F_per, [KK4, shtns%nlm])
    call c_f_pointer(c_loc(EFT_in(2*shtns%nlm*(KK2+KK4) + 1)), T_per, [KK2, shtns%nlm])

    ! -------------------------------------------------------------------------
    !  B u = (I - dt L)^-1 * dt * u
    !
    !  Copie du pas de temps linearise dans laquelle :
    !    - on n'appelle PAS comp_LinNonLin,
    !    - on remplace DE, DF, DT par E_per, F_per, T_per,
    !    - on SUPPRIME le terme de masse (Banded_Mult avec Yef en cn,
    !      Ye_mat/Yf_mat en fbe, wT pour T).
    !
    !  Pourquoi passer par Chb2/Chb4 et pas simplement delta_t*E_per :
    !  Chb2(:KK,k) = ChbPoly(k-1,xK)*delta_t  et  Chb2(KK+1,KK+1) = 1.
    !  Autrement dit Chb2 EST "dt fois l'identite" exprime dans la
    !  representation mixte du code : il transforme un champ en coefficients
    !  de Chebyshev vers le second membre du solve implicite, en incluant le
    !  facteur dt et en traitant les lignes de conditions aux limites.
    !  Mettre delta_t*E_per a la place placerait le vecteur dans le mauvais
    !  espace : le second membre aurait alors une composante hors de l'image
    !  de l'operateur, et GMRES stagnerait.
    ! -------------------------------------------------------------------------

    DEr = (0., 0.)
    DFr = (0., 0.)

    select case (time_step)

    case ("cn")

        ! DE <- E_per, DF <- F_per  (dt et CL portes par Chb2 / Chb4)
        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(k, lm)
        do lm = 2, shtns%nlm
            do k = 1, KK2
                DEr(k, lm)%re = Chb2(k, :KK2) .dot. real (E_per(:KK2, lm))
                DEr(k, lm)%im = Chb2(k, :KK2) .dot. aimag(E_per(:KK2, lm))
            end do
        end do
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(k, lm)
        do lm = 2, shtns%nlm
            do k = 1, KK4
                DFr(k, lm)%re = Chb4(k, :KK4) .dot. real (F_per(:KK4, lm))
                DFr(k, lm)%im = Chb4(k, :KK4) .dot. aimag(F_per(:KK4, lm))
            end do
        end do
        !$OMP END PARALLEL DO

        call small2big(DEr, DFr, DA)

        ! Pas de Banded_Mult(Yef, ...) : on veut dt*I*u, pas Yef*u
        stride = 1
        Ap = 0.
        do m = 0, MM*mres, mres
            A1  = 0.
            top = 2*(KK2 + KK4)*(LL - (max(m,1) - 1))
            A1 = Banded_LU_solve(top, 2*(KK2+KK4)-1, Xef(:, stride), 6*(KK2+KK4)-2, DA(stride,1), PIVOT(stride))
            Ap(stride : stride+top-1, 1) = A1(:top, 1)
            stride = stride + top
        end do

        call big2small(Ap, E_per, F_per)

    case ("fbe")

        ! Idem : on garde le terme Chb2*DE avec DE <- E_per, on supprime
        ! le terme de masse Ye_mat*E_per / Yf_mat*F_per.
        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(k, lm)
        do lm = 2, shtns%nlm
            do k = 1, KK2
                DEr(k, lm)%re = Chb2(k, :KK2) .dot. real (E_per(:KK2, lm))
                DEr(k, lm)%im = Chb2(k, :KK2) .dot. aimag(E_per(:KK2, lm))
            end do
        end do
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(k, lm)
        do lm = 2, shtns%nlm
            do k = 1, KK4
                DFr(k, lm)%re = Chb4(k, :KK4) .dot. real (F_per(:KK4, lm))
                DFr(k, lm)%im = Chb4(k, :KK4) .dot. aimag(F_per(:KK4, lm))
            end do
        end do
        !$OMP END PARALLEL DO

        call small2big(DEr, DFr, DA)

        stride = 1
        Ap = 0.
        do m = 0, MM*mres, mres
            A1  = 0.
            top = 2*(KK2 + KK4)*(LL - (max(m,1) - 1))
            A1 = Banded_LU_solve(top, 2*(KK2+KK4)-1, Xef(:, stride), 6*(KK2+KK4)-2, DA(stride,1), PIVOT(stride))
            Ap(stride : stride+top-1, 1) = A1(:top, 1)
            stride = stride + top
        end do

        call big2small(Ap, E_per, F_per)

    end select

    !--- Temperature : DT <- T_per, et on supprime le terme wT (masse).
    !    XT_inv = transpose(inv(XTh) . Chb2) contient deja le facteur dt.
    !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(k, lm)
    do lm = 1, shtns%nlm
        do k = 1, KK2
            wDT(k, lm)%re = XT_inv(:, k, ell(lm)) .dot. real (T_per(:, lm))
            wDT(k, lm)%im = XT_inv(:, k, ell(lm)) .dot. aimag(T_per(:, lm))
        end do
    end do
    !$OMP END PARALLEL DO

    T_per = wDT

    ! E_per / F_per / T_per pointent dans EFT_in, qui contient le resultat
    EFT_out = EFT_in

end subroutine TermeDroite

! =============================================================================
!  subA_inverse : opérateur complet passé à GMRES, avec décalage :
!        A x = (Phi - I) x  -  sigma_r * TermeDroite(x)
!  Si sigma_r = 0, se réduit exactement à subA_inv.
! =============================================================================
subroutine subA_inverse(time_stepping_sub, EFT_in, EFT_out)

    implicit none

    procedure(), pointer, intent(in) :: time_stepping_sub

    double precision, dimension(2*shtns%nlm*(2*KK2+KK4)), intent(inout) :: EFT_in
    double precision, dimension(2*shtns%nlm*(2*KK2+KK4)), intent(out)   :: EFT_out

    double precision, dimension(:), allocatable, save :: shift_vec, in_copy
    integer :: nn

    nn = 2*shtns%nlm*(2*KK2 + KK4)

    if (sigma_r /= 0.d0) then
        if (.not. allocated(shift_vec)) allocate(shift_vec(nn))
        if (.not. allocated(in_copy))   allocate(in_copy(nn))
        in_copy = EFT_in                       ! TermeDroite modifie son entrée
        call TermeDroite(time_stepping_sub, in_copy, shift_vec)
    end if

    call subA_inv(time_stepping_sub, EFT_in, EFT_out)

    if (sigma_r /= 0.d0) then
        EFT_out = EFT_out - sigma_r * shift_vec
    end if

end subroutine subA_inverse

! =============================================================================
!  subA_cplx : opérateur à DÉCALAGE COMPLEXE, écrit en arithmétique réelle.
!
!  On veut appliquer  (A - sigma_c B)  avec  sigma_c = sigma_r + i*sigma_i
!  à un vecteur complexe  x = x_R + i*x_I.  En développant :
!
!     (A - sigma_c B)(x_R + i x_I) =
!         [ A x_R - sigma_r B x_R + sigma_i B x_I ]      <- partie réelle
!       + [ A x_I - sigma_r B x_I - sigma_i B x_R ] i    <- partie imaginaire
!
!  Le vecteur passé à GMRES est la concaténation X = [x_R ; x_I], de taille 2n.
!  Chaque application coûte 2 subA_inv (A) + 2 TermeDroite (B).
!
!  Pourquoi le décalage complexe : la puissance inverse converge vers la valeur
!  propre la plus proche de sigma DANS LE PLAN COMPLEXE. Une paire de Hopf
!  a+ib ne peut jamais être approchée à moins de |b| par un sigma réel, alors
!  que le mode de phase (lambda = 0, présent dès que l'état de base est une
!  onde) est SUR l'axe réel : il gagne toujours. En sortant de l'axe réel on
!  vise directement la paire. Bonus : sigma_c étant proche de lambda et loin de
!  sa conjuguée, la dégénérescence lambda/lambda_barre disparaît et l'itération
!  converge vers un vecteur propre unique (plus besoin de sous-espace 2D).
! =============================================================================
subroutine subA_cplx(time_stepping_sub, X_in, X_out)

    implicit none

    procedure(), pointer, intent(in) :: time_stepping_sub

    double precision, dimension(4*shtns%nlm*(2*KK2+KK4)), intent(inout) :: X_in
    double precision, dimension(4*shtns%nlm*(2*KK2+KK4)), intent(out)   :: X_out

    double precision, dimension(:), allocatable, save :: cp, AxR, AxI, BxR, BxI
    integer :: nn

    nn = 2*shtns%nlm*(2*KK2 + KK4)

    if (.not. allocated(cp))  allocate(cp(nn), AxR(nn), AxI(nn), BxR(nn), BxI(nn))

    ! --- partie réelle de l'entrée
    cp = X_in(1:nn)
    call subA_inv   (time_stepping_sub, cp, AxR)   ! A x_R
    cp = X_in(1:nn)
    call TermeDroite(time_stepping_sub, cp, BxR)   ! B x_R

    ! --- partie imaginaire de l'entrée
    cp = X_in(nn+1:2*nn)
    call subA_inv   (time_stepping_sub, cp, AxI)   ! A x_I
    cp = X_in(nn+1:2*nn)
    call TermeDroite(time_stepping_sub, cp, BxI)   ! B x_I

    X_out(1:nn)      = AxR - sigma_r*BxR + sigma_i*BxI
    X_out(nn+1:2*nn) = AxI - sigma_r*BxI - sigma_i*BxR

end subroutine subA_cplx

! =============================================================================
!  inverse_convective_solver_cplx : puissance inverse à décalage complexe.
!  Utilisée dès que sigma_i /= 0. Réutilise subA_inv, TermeDroite et gmresm.
! =============================================================================
subroutine inverse_convective_solver_cplx(LinNonLinTimeStep_ptr, n, iunit)

    implicit none

    procedure(), pointer, intent(in) :: LinNonLinTimeStep_ptr
    integer, intent(in) :: n, iunit

    integer :: it, gmres_iters, gmres_total, i, seed_size
    integer, dimension(:), allocatable :: seed
    double precision, allocatable :: U(:), B(:), W(:), cp(:)
    double precision :: rr, ri, uu, uwr, uwi, den, ww, err, dlam
    double precision :: lam_r, lam_i, lam_r_old, lam_i_old
    procedure(), pointer :: op_ptr
    logical :: converged

    allocate(U(2*n), B(2*n), W(2*n), cp(n))

    ! ---------------- vecteur initial : aléatoire, lissé par B, normalisé
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = [(1234 + 17*i, i = 1, seed_size)]
    call random_seed(put=seed)
    call random_number(U)
    U = U - 0.5d0

    cp = U(1:n)
    call TermeDroite(LinNonLinTimeStep_ptr, cp, B(1:n))
    cp = U(n+1:2*n)
    call TermeDroite(LinNonLinTimeStep_ptr, cp, B(n+1:2*n))
    call dot(n, B(1),   B(1),   rr)
    call dot(n, B(n+1), B(n+1), ri)
    U = B / sqrt(rr + ri)

    print*
    write(iunit,*)
    print*, "------- Beginning inverse power iteration (COMPLEX shift) -------"
    write(iunit,*) "------- Beginning inverse power iteration (COMPLEX shift) -------"
    print*, " shift sigma = ", sigma_r, " + i * ", sigma_i
    write(iunit,*) " shift sigma = ", sigma_r, " + i * ", sigma_i
    print*
    write(iunit,*)

    op_ptr => subA_cplx
    converged   = .false.
    lam_r = 0.d0 ; lam_i = 0.d0
    gmres_total = 0

    do it = 1, max_inv

        ! --- second membre  B u  (partie réelle puis imaginaire)
        cp = U(1:n)
        call TermeDroite(LinNonLinTimeStep_ptr, cp, B(1:n))
        cp = U(n+1:2*n)
        call TermeDroite(LinNonLinTimeStep_ptr, cp, B(n+1:2*n))

        ! --- résolution  (A - sigma_c B) W = B u   par GMRES, taille 2n
        call gmresm(op_ptr, LinNonLinTimeStep_ptr, restart_gmres, 2*n, max_gmres, &
                    & tol_gmres, B, W, gmres_iters)
        gmres_total = gmres_total + gmres_iters

        ! --- quotient de Rayleigh complexe :  lambda = sigma_c + <u,u>/<u,w>
        !     <u,u> = |u_R|^2 + |u_I|^2                              (réel)
        !     <u,w> = (u_R.w_R + u_I.w_I) + i (u_R.w_I - u_I.w_R)
        lam_r_old = lam_r ; lam_i_old = lam_i

        call dot(n, U(1),   U(1),   rr)
        call dot(n, U(n+1), U(n+1), ri)
        uu = rr + ri

        call dot(n, U(1),   W(1),   rr)
        call dot(n, U(n+1), W(n+1), ri)
        uwr = rr + ri
        call dot(n, U(1),   W(n+1), rr)
        call dot(n, U(n+1), W(1),   ri)
        uwi = rr - ri

        den = uwr*uwr + uwi*uwi
        lam_r = sigma_r + uu*uwr/den
        lam_i = sigma_i   - uu*uwi/den

        ! --- normalisation et résidu vectoriel, à une PHASE COMPLEXE près :
        !     min_theta || w - e^{i theta} u ||  ->  err^2 = 2 - 2|<u,w>|
        call dot(n, W(1),   W(1),   rr)
        call dot(n, W(n+1), W(n+1), ri)
        ww = rr + ri
        W  = W / sqrt(ww)

        call dot(n, U(1),   W(1),   rr)
        call dot(n, U(n+1), W(n+1), ri)
        uwr = rr + ri
        call dot(n, U(1),   W(n+1), rr)
        call dot(n, U(n+1), W(1),   ri)
        uwi = rr - ri
        err = sqrt(max(0.d0, 2.d0 - 2.d0*sqrt(uwr*uwr + uwi*uwi)/uu))

        dlam = sqrt((lam_r - lam_r_old)**2 + (lam_i - lam_i_old)**2)

        print*, "Inverse iteration ", it, " | lambda = ", lam_r, " + i ", lam_i, &
                & " | d(lambda) = ", dlam, " | vec. residual = ", err, " | gmres its = ", gmres_iters
        write(iunit,*) "Inverse iteration ", it, " | lambda = ", lam_r, " + i ", lam_i, &
                & " | d(lambda) = ", dlam, " | vec. residual = ", err, " | gmres its = ", gmres_iters
        flush(6)
        flush(iunit)

        U = W

        if (err < inv_eps) then
            converged = .true.
            exit
        end if

    end do

    print*
    write(iunit,*)
    if (converged) then
        print*, "----------------------- Converged! ------------------------------"
        write(iunit,*) "----------------------- Converged! ------------------------------"
    else
        print*, "--------------- Maximum number of iterations reached ------------"
        write(iunit,*) "--------------- Maximum number of iterations reached ------------"
    end if
    print*, "Eigenvalue closest to (", sigma_r, ",", sigma_i, ") : lambda = ", lam_r, " + i ", lam_i
    write(iunit,*) "Eigenvalue closest to (", sigma_r, ",", sigma_i, ") : lambda = ", lam_r, " + i ", lam_i
    print*, "Total GMRES iterations: ", gmres_total
    write(iunit,*) "Total GMRES iterations: ", gmres_total

    ! Sortie : partie réelle puis partie imaginaire du mode propre
    step = 1
    call output_vector(U(1), U(2*shtns%nlm*KK2 + 1), U(2*shtns%nlm*(KK2+KK4) + 1))
    call write_mode_restart(U(1), U(2*shtns%nlm*KK2 + 1), U(2*shtns%nlm*(KK2+KK4) + 1), 1)
    step = 2
    call output_vector(U(n+1), U(n + 2*shtns%nlm*KK2 + 1), U(n + 2*shtns%nlm*(KK2+KK4) + 1))
    call write_mode_restart(U(n+1), U(n + 2*shtns%nlm*KK2 + 1), U(n + 2*shtns%nlm*(KK2+KK4) + 1), 2)

    deallocate(U, B, W, cp, seed)

end subroutine inverse_convective_solver_cplx

! =============================================================================
!  inverse_convective_solver : programme principal du solveur.
! =============================================================================
subroutine inverse_convective_solver()

    implicit none

    integer :: n, it, gmres_iters, gmres_total
    double precision :: uu, uw, ww, lambda, lambda_old, err, err_p, err_m
    double precision, allocatable :: u(:), b(:), w(:)
    procedure(), pointer :: LinNonLinTimeStep_ptr
    procedure(), pointer :: subA_inv_ptr
    logical :: converged
    integer :: iunit, i, seed_size
    integer, dimension(:), allocatable :: seed

    iunit = 43
    open(unit=iunit, file='inverse_output.txt', status='replace', action='write')

    n = 2*(2*KK2 + KK4)*shtns%nlm

    ! ==================================================================
    !  Précalcul des matrices
    ! ==================================================================
    call precompBuildXY()

    if (solver == "inverse_convective_explicit") then

        select case (time_step)
            case ("cn")
                print*, "Precomputing the Xe, Ye, Xf, Yf, XT and YT matrices..."
                write(iunit,*) "Precomputing the Xe, Ye, Xf, Yf, XT and YT matrices..."
                call precompXeYe()
                call precompXfYf()
                call precompXTYT()
            case ("fbe")
                print*, "Precomputing the Xe, Ye, Xf, Yf, XT and YT matrices..."
                write(iunit,*) "Precomputing the Xe, Ye, Xf, Yf, XT and YT matrices..."
                call precompXeYe_BDF()
                call precompXfYf_BDF()
                call precompXTYT_BDF()
        end select

        LinNonLinTimeStep_ptr => compute_lin_time_step_convective_explicit_CN_FBE

    else if (solver == "inverse_convective_implicit") then

        select case (time_step)
            case ("cn")
                print*, "Precomputing the Xef, Yef, XT and YT matrices..."
                write(iunit,*) "Precomputing the Xef, Yef, XT and YT matrices..."
                call PrecompimplicitXY()
                call precompXTYT()
                LinNonLinTimeStep_ptr => compute_lin_time_step_convective_implicit_CN
            case ("fbe")
                print*, "Precomputing the Xef, Ye, Yf, XT and YT matrices..."
                write(iunit,*) "Precomputing the Xef, Ye, Yf, XT and YT matrices..."
                call PrecompimplicitXY_BDF()
                call precompXTYT_BDF()
                LinNonLinTimeStep_ptr => compute_lin_time_step_convective_implicit_FBE
        end select

    end if

    ! ==================================================================
    !  État de base
    ! ==================================================================
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

    call comp_radial_derivatives(E_base, F_base, T_base, Qu, Su, Tu, Qcu, Scu, Tcu, Tr, St)
    call comp_U_cU_gT(Ur_base, Ut_base, Up_base, cUr_base, cUt_base, cUp_base, gTr_base, gTt_base, gTp_base)

    call c_f_pointer(c_loc(E_base), E_base_ptr, [2*KK2*shtns%nlm])
    call c_f_pointer(c_loc(F_base), F_base_ptr, [2*KK4*shtns%nlm])
    call c_f_pointer(c_loc(T_base), T_base_ptr, [2*KK2*shtns%nlm])

    step = 0
    call output_vector(E_base_ptr, F_base_ptr, T_base_ptr)

    print*, "  Base state wavespeed C_base = ", C_base
    write(iunit,*) "  Base state wavespeed C_base = ", C_base

    ! ==================================================================
    !  Décalage complexe (sigma_i /= 0) : on bascule sur le solveur dédié,
    !  seule façon d'atteindre une paire de Hopf (voir subA_cplx).
    ! ==================================================================
    if (sigma_i /= 0.d0) then
        call inverse_convective_solver_cplx(LinNonLinTimeStep_ptr, n, iunit)
        close(iunit)
        return
    end if

    ! ==================================================================
    !  Vecteur initial : aléatoire, lissé par un TermeDroite, normalisé
    ! ==================================================================
    allocate(u(n), b(n), w(n))

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = [(1234 + 17*i, i = 1, seed_size)]
    call random_seed(put=seed)
    call random_number(u)
    u = u - 0.5d0

    w = u
    call TermeDroite(LinNonLinTimeStep_ptr, w, b)
    call dot(n, b, b, ww)
    u = b / sqrt(ww)

    ! ==================================================================
    !  Itération de la puissance inverse
    ! ==================================================================
    print*
    write(iunit,*)
    print*, "--------------- Beginning inverse power iteration ---------------"
    write(iunit,*) "--------------- Beginning inverse power iteration ---------------"
    print*, " shift sigma_r = ", sigma_r
    write(iunit,*) " shift sigma_r = ", sigma_r
    print*
    write(iunit,*)

    subA_inv_ptr => subA_inverse
    converged   = .false.
    lambda      = 0.d0
    gmres_total = 0

    do it = 1, max_inv

        ! Second membre  b = B u_n = TermeDroite(u)   (copie : modifie son entrée)
        w = u
        call TermeDroite(LinNonLinTimeStep_ptr, w, b)

        ! Résolution  A w = b  par GMRES
        call gmresm(subA_inv_ptr, LinNonLinTimeStep_ptr, restart_gmres, n, max_gmres, &
                    & tol_gmres, b, w, gmres_iters)
        gmres_total = gmres_total + gmres_iters

        ! Quotient de Rayleigh : w ~ u / (lambda - sigma)
        ! NB : si la PROPOSITION 3 (facteur delta_t) est désactivée, la valeur
        !      obtenue vaut delta_t * lambda_physique.
        lambda_old = lambda
        call dot(n, u, u, uu)
        call dot(n, u, w, uw)
        lambda = sigma_r + uu/uw

        ! Normalisation et résidu vectoriel (au signe près)
        call dot(n, w, w, ww)
        w = w / sqrt(ww)
        call dot(n, w - u, w - u, err_p)
        call dot(n, w + u, w + u, err_m)
        err = sqrt(min(err_p, err_m))

        print*, "Inverse iteration ", it, " | lambda = ", lambda, &
                & " | d(lambda) = ", abs(lambda - lambda_old), &
                & " | vec. residual = ", err, " | gmres its = ", gmres_iters
        write(iunit,*) "Inverse iteration ", it, " | lambda = ", lambda, &
                & " | d(lambda) = ", abs(lambda - lambda_old), &
                & " | vec. residual = ", err, " | gmres its = ", gmres_iters
        flush(6)
        flush(iunit)

        u = w

        if (err < inv_eps) then
            converged = .true.
            exit
        end if

    end do

    ! ==================================================================
    !  Résultats
    ! ==================================================================
    print*
    write(iunit,*)
    if (converged) then
        print*, "----------------------- Converged! ------------------------------"
        write(iunit,*) "----------------------- Converged! ------------------------------"
    else
        print*, "--------------- Maximum number of iterations reached ------------"
        write(iunit,*) "--------------- Maximum number of iterations reached ------------"
        print*, "Si lambda oscille pendant que le residu vectoriel stagne, la cible"
        print*, "est probablement une paire complexe (Hopf) : la puissance inverse"
        print*, "reelle ne peut pas converger vers un vecteur unique."
    end if
    print*, "Eigenvalue closest to ", sigma_r, " : lambda = ", lambda
    write(iunit,*) "Eigenvalue closest to ", sigma_r, " : lambda = ", lambda
    print*, "Total GMRES iterations: ", gmres_total
    write(iunit,*) "Total GMRES iterations: ", gmres_total

    step = 1
    call output_vector(u(1), u(2*shtns%nlm*KK2 + 1), u(2*shtns%nlm*(KK2+KK4) + 1))
    call write_mode_restart(u(1), u(2*shtns%nlm*KK2 + 1), u(2*shtns%nlm*(KK2+KK4) + 1), 1)

    deallocate(u, b, w, seed)
    close(iunit)

end subroutine inverse_convective_solver

end module mod_InverseSolver
