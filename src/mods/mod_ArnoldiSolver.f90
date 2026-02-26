module mod_ArnoldiSolver

    use iso_c_binding

    use mod_GlobalVars
    use mod_PrecompSH
    use mod_read
    use mod_TimeStep
    use mod_Output
    use mod_tools

    implicit none

contains

subroutine arnoldi_convective_solver()

    procedure(), pointer :: LinNonLinTimeStep_ptr
    logical :: condition
    double precision :: tolerance

    ! We begin by precomputing the matrices for the resolution of the linear system
    call precompBuildXY() ! Compute building blocks for the matrices

    ! Check for solver type
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

        ! Pointers for subroutines
        LinNonLinTimeStep_ptr => compute_lin_time_step_convective_explicit_CN_FBE

    else if (solver == "arnoldi_convective_implicit") then

        select case (time_step)
        
            case ("cn")

            print*, "Precomputing the Xef and Yef matrices..."
            call PrecompimplicitXY()
            print*, "Precomputing the XT and YT matrices..."
            call precompXTYT()
            ! Pointers for subroutines
            LinNonLinTimeStep_ptr => compute_lin_time_step_convective_implicit_CN

            case ("fbe")

            print*, "Precomputing the Xef, Ye and Yf matrices..."
            call PrecompimplicitXY_BDF()
            print*, "Precomputing the XT and YT matrices..."
            call precompXTYT_BDF()
            ! Pointers for subroutines
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

    ! Save base state
    E_base = E
    F_base = F
    T_base = T

    ! Compute base field explicit terms
    ! ... see calls to comp_radial_derivatives, comp_U_cU_gT and ToReal in comp_RHS_with_rot

    ! Initialize condition
    condition = .true.

    ! Initialize first vector randomly
    ! ....

    do while (condition)

        ! 1# Time integration
        do step = 1, NTS

            ! call to act_lin_non_lin for linearized timestep

            ! add if condition to save vectors

        end do

        ! 2# Call to Arnoldi iteration
        ! See gmeresm subroutine in mod_IterativeSolvers

        ! 3# Produce Hij

        ! 4# Compute eigenvalues and eigenvectors of Hij

        ! 5# Check AV phi with lambda V phi and update condition

    end do

    ! Finalize simulation

end subroutine arnoldi_convective_solver

end module mod_ArnoldiSolver