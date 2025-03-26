program ps_decoherence

! Evolving primordial fluctuations undergoing decoherence effects.
! The code solves single mode evolution and produces deformed power
! spectrum for an arbitrary array of state accidents. Please cite
! arxiv:XXXX.XXXXX if you use it in your research.
!!!!!!!!!!!!!!!!!!
! Instructions:
! To compile (intel fortran is finally free!):
! ifort -r8 -O3 Power_spectrum.f90 -o name_of_exec.out
! To run:
! ./name_of_exec.out
!!!!!!!!!!!!!!!!!!!!!

    implicit none
    ! This is exactly what you think it is...
    real, parameter :: twopi = 6.28318530717958647692528676655900577Q0
    ! Locating the accident in length scale, frac_e is where the accident
    ! and frac_o is where it ends, h_tran is the width in k-space
    real, parameter :: frac_e = 0.7, frac_o = 0.8, h_tran = 20.0
    ! size of the state vector
    integer, parameter :: nvar = 7


    ! Physical (kphys) and comoving (kcom) wavevectors
    real :: kphys, kcom
    ! timestep and counter to inject different modes
    real :: dt, N_ref
    
    ! different state vectors
    real :: state(nvar), state_init(nvar), state_back(nvar)
    ! Dummy human-readable variables
    real :: phi, phi_dot, Hbl, ln_a, L_k, L_k_prime, theta_k_prime

    ! Control variables
    logical :: change, break, elipse, back
   ! File and iteration counter
    integer :: j, iter
   ! Setting up the height of the accident: lO - lE controls the duration
    real :: lE, lO

    ! Setting up the accidents
    ! Case 0: No accidents
    ! Case 1: One accident
    ! Case 2: Hysteresis
    ! Case 1 is set by default
    integer, parameter :: mod_pref = 1
    ! exp(N_star)*k_phys sets the position of the accident in k-space
    ! delta_star sets the width
    real :: N_star, delta_star
    ! Number of modes to produce the power spectrum
    integer, parameter :: N_mod = 100

    ! Integer labels for the Wigner ellipse evolution file and the power spectrum
    integer :: unit_output, unit_summary
    ! Name of the output file
    character(len=100) :: filename  ! Nombre del archivo de salida

    ! Set as .true. to produce Wigner ellipse evolution, and as .false. to skip it
    elipse = .false. !false

    
    ! Evolving backround from N = -5 and saving in state_back
    back = .true.

    ! Initial background field conditions in human-readable form passed to
    ! state vector
    phi = 25.0
    phi_dot = 0.0
    ln_a = -5.0
    state_back(1) = phi
    state_back(2) = phi_dot
    state_back(3) = hubb(state_back)
    state_back(4) = ln_a

! Evolution loop to converge to the attractor (ending after 5 e-folds of evolution)
    do while (back)
        call gl8(state_back, dt, kcom)
        if (state_back(4) >= 0.0) then
            back = .false.
        end if
    end do


    ! Open file for power spectrum
    unit_summary = 20
    open(unit=unit_summary, file="PS_EO78_TEST_1.dat", status="unknown", action="write", form="formatted")

    ! Loop to produce the power spectrum from N_mod modes
    do iter = 1, N_mod
        N_ref = 0.5 * iter
        N_star = 0.5 * N_mod
        delta_star = 0.5 * N_mod / 10

        ! Calling background variables
        state_init = state_back ! Campo escalar φ pre-evolucioanada
        kphys= 1000.0 * state_back(3)
        ! Setting up the timestep
        dt = twopi * (1.0/ (1000.0 * state_init(3))) / 60.0
        ! Setting up the location of the accident
        lE = exp((frac_e - 1.0) * log(kphys/twopi) - frac_e * log(state_back(3)))
        lO = exp((frac_o - 1.0) * log(kphys/twopi) - frac_o * log(state_back(3)))

        ! Turning on the perturbations
        change = .true.
        break = .true.

        ! Open file for dynamics of the Wigner ellipse
        if (elipse) then
            write(filename, '("Elipse_Open_N_ref_", I1, ".dat")') iter
            open(unit=unit_output, file=filename, status="unknown", action="write", form="formatted")
        end if

        ! 🔹 Evolving the background up to N_ref e-folds...
        state = state_init
        do while (change)
            call gl8(state, dt, kcom)
            if (state(4) >= N_ref) then
                ! save state vector
                state_init = state
                ! Compute conformal k and set-up initial perturbations in
                ! Minkowski vacuum
                kcom = kphys * exp(state(4))  ! Cálculo del modo conforme
                L_k = 1.0 / sqrt(2.0 * sqrt(kcom*kcom - zpp_over_z(state_init)))
                L_k_prime = 0.0
                theta_k_prime = sqrt(kcom**2 - zpp_over_z(state_init))
                ! Mukhanov-Sasaki state vector
                state_init(5) = L_k
                state_init(6) = L_k_prime
                state_init(7) = theta_k_prime
                ! End evolution
                change = .false.
            end if
        end do

        ! 🔹 Main evolution loop up to the end of inflation (ε >= 1.0)
        state = state_init
        j = 0
        do while (break)
            ! Reducing time step to evolve squeezed modes
            call gl8(state, dt / 1.0, kcom)
            j = j + 1
            if (epsilon(state) >= 1.0) then
                ! Exit when inflation ends
                break = .false.
            end if

            ! Store ellipse evolution every 1000 steps
            if (elipse .and. mod(j, 1000) == 0) then
                write(unit_output, '(13(ES15.7E3,1X))') state(1), state(2), state(3), state(4), state(5), state(6), state(7), &
                                                        epsilon(state), hubb(state), zpp_over_z(state), kcom, kphys, z_func(state)
            end if
        end do

        ! Save and write state vector
        write(unit_summary, '(13(ES15.7E3,1X))') state(1), state(2), state(3), &
        state(4), state(5), state(6), state(7), epsilon(state), hubb(state), &
        zpp_over_z(state), kcom, kphys, z_func(state)

        ! Close file and write
        if (elipse) then
            close(unit=unit_output)
        end if
        write(*,*) "Simulation ends for N_ref = ", N_ref

        ! Reset background
        state = state_back
        state_init = state_back
    end do

    ! Close power spectrum file
    close(unit=unit_summary)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the main code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

    ! 🔹 Equations of motion
    subroutine evalf(y, dydx, kcom)
        real, intent(in) :: y(nvar)
        real, intent(out) :: dydx(nvar)
        real, intent(in) ::kcom

        ! dphi/dt = phidot
        dydx(1) = y(2)
        ! dphidot/dt = -3H*phidot - V´
        dydx(2) = -3.0 * y(3) * y(2) - Vprime(y)
        ! dH/dt = -phidot**2/2
        dydx(3) = -0.5 * y(2)*y(2)
        ! d ln_a/dt = H
        dydx(4) = y(3)
        ! d L_k/dt = L_k´/a
        dydx(5) = y(6) * exp(-y(4))
        ! d L_k'/dt = L_k''/a = -(k**2 - z''/z - theta_k_prime**2)*L_k/a
        dydx(6) = -(kcom*kcom - zpp_over_z(y) - y(7)*y(7)) * y(5) * exp(-y(4))
        ! d theta_k_prime/dt = -2 L'*theta_k_prime/L/a + source/(2*a*theta_prime**2*L_k**2)
        dydx(7) = -2.0 * (y(7) * y(6) / y(5)) * exp(-y(4)) + &
                  source_open(mod_pref, y, kcom) * exp(-y(4)) / (2.0 * y(7) * y(5)*y(5))  ! Adición de la fuente
    end subroutine evalf

    ! 🔹 8th-order Gauss-Legendre integrator
subroutine gl8(y, dt, kcom)
    integer, parameter :: s = 4, n = 7
    real y(n), g(n,s), dt; integer i, k
    real, intent(in) ::kcom
    
    ! Butcher tableau for 8th order Gauss-Legendre method
    real, parameter :: a(s,s) = reshape((/ &
             0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
             0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
             0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
            -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
             0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
             0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
             0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
             0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
    real, parameter ::   b(s) = (/ &
             0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
             0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)
    
    ! iterate trial steps
    g = 0.0; do k = 1,16
            g = matmul(g,a)
            do i = 1,s
                    call evalf(y + g(:,i)*dt, g(:,i), kcom)
            end do
    end do
    
    ! update the solution
    y = y + matmul(g,b)*dt

end subroutine gl8

    ! 🔹 Potential, derivatives and auxiliary functions
    ! V(phi)
    function Vphi(y)
        real, intent(in) :: y(nvar)
        real Vphi, phi
        phi = y(1)
        Vphi = 1.0d-14 * phi*phi*phi*phi / 4.0
    end function Vphi

    ! V'(phi)
    function Vprime(y)
        real, intent(in) :: y(nvar)
        real Vprime, phi
        phi = y(1)
        Vprime = 1.0d-14 * phi*phi*phi
    end function Vprime

    ! V''(phi)
    function Vprimeprime(y)
        real, intent(in) :: y(nvar)
        real Vprimeprime, phi
        phi = y(1)
        Vprimeprime = 3.0*1.0d-14 * phi*phi
    end function Vprimeprime

    ! 🔹 Hubble parameter
    function hubb(y)
        real, intent(in) :: y(nvar)
        real hubb, phi, phi_dot
        phi = y(1); phi_dot = y(2)
        hubb = sqrt(phi_dot*phi_dot / 6.0 + Vphi(y) / 3.0)
    end function hubb

    ! 🔹 epsilon
    function epsilon(y)
        real, intent(in) :: y(nvar)
        real epsilon, phi_dot, Hubble
        Hubble = y(3); phi_dot = y(2)
        epsilon = phi_dot*phi_dot / (2.0 * Hubble*Hubble)
    end function epsilon

    ! 🔹 Computing z''/z
    function zpp_over_z(y)
        real y(nvar), epsilon1, epsilon2, epsilon3, zpp_over_z, hubbp, hubbpp, hubbppp
        real Vpp, ddot_phi, tdot_phi, eps2dot
        real z_func, phi_dot, Hubble, e_fold
        e_fold = y(4); Hubble = y(3); phi_dot = y(2)
            Vpp = Vprimeprime(y)
            ddot_phi = -3.d0*Hubble*phi_dot-Vprime(y)
            epsilon1 = epsilon(y)
            hubbp = -0.5*phi_dot*phi_dot
            tdot_phi = -3.0*hubbp*phi_dot-3.0*Hubble*ddot_phi-Vpp*phi_dot
            hubbpp = -phi_dot*ddot_phi
            hubbppp = -ddot_phi**2-phi_dot*tdot_phi
            epsilon2 = hubbpp/(Hubble*hubbp)-2.0*hubbp/Hubble**2
            eps2dot = hubbppp/(Hubble*hubbp)-(hubbpp*(hubbp**2+hubbpp*Hubble))/(Hubble*hubbp)**2-2.0*hubbpp/Hubble**2+4.0*hubbp**2/Hubble**3
            epsilon3 = eps2dot/(Hubble*epsilon2)

            zpp_over_z = exp(2.0*e_fold)*Hubble**2*(2.d0-epsilon1+1.5*epsilon2+0.25*epsilon2**2-0.5*epsilon2*epsilon1+0.5*epsilon2*epsilon3)
    end function zpp_over_z

    ! 🔹 Computing z
    function z_func(y)
        real, intent(in) :: y(nvar)
        real z_func, phi_dot, Hubble, ln_a
        ln_a = y(4); Hubble = y(3); phi_dot = y(2)
        z_func = phi_dot * exp(ln_a) / Hubble
    end function z_func

    ! 🔹 Smooth Heaviside unit square for the accident (in k)
    pure function SmoothHeaviside(x) result(H)
        real, intent(in) :: x
        real H
        H = 0.5 + 0.5 * tanh(h_tran * x)
    end function SmoothHeaviside


    ! 🔹 Not-as-smooth Heaviside unit square for the accident (in t)
    pure function SmoothHeaviside_M(x) result(H)
        real, intent(in) :: x
        real H
        H = 0.5 + 0.5 * tanh(h_tran * x / 5.0)
    end function SmoothHeaviside_M



    ! 🔹 Source term in equations of motion (cases)
    function source_open(mod_pref, y, kcom) result(source)
        integer, intent(in) :: mod_pref
        real, intent(in) :: y(nvar)
        real, intent(in) :: kcom
        real e_fold, source
        e_fold = y(4)

        select case (mod_pref)
            case (0)  ! No accident
                source = 0.0

            case (1)  ! One accident
                source = (SmoothHeaviside(1.0 - lE / (e_fold - log(kcom))) - &
                          SmoothHeaviside(1.0 - lO / (e_fold - log(kcom)))) * &
                         (SmoothHeaviside_M(1.0 - (N_star-delta_star)/(log(kcom) - log(kphys))) - &
                          SmoothHeaviside_M(1.0 - (N_star+delta_star)/(log(kcom) - log(kphys))))
            case (2)  ! Hysteresis
                source = (SmoothHeaviside(1.0 - lE / (e_fold - log(kcom))) - &
                         2.0 * SmoothHeaviside(1.0 - lO / (e_fold - log(kcom))) + &
                         SmoothHeaviside(1.0 - (2.0 * lO - lE) / (e_fold - log(kcom))) ) * &
                        (SmoothHeaviside_M(1.0 - (N_star-delta_star)/(log(kcom) - log(kphys))) - &
                         SmoothHeaviside_M(1.0 - (N_star+delta_star)/(log(kcom) - log(kphys))))

            case default
                print *, "¡Error! mod_pref must be 0, 1, 2, 3."
                stop
        end select
    end function source_open

end program ps_decoherence
