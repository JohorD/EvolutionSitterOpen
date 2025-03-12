program simanhosc
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)  ! Quadruple precision
    real(qp), parameter :: twopi = 6.28318530717958647692528676655900577_qp  ! 2π
    integer, parameter :: nvar = 7  ! Number of equations in the system

    ! Main variables
    real(qp) :: kphys, kcom, dt  ! Physical and conformal modes
    real(qp) :: y(nvar), yinit(nvar)  ! System variables
    integer :: unit_output  ! Output file unit
    logical :: change, break  ! Control variables
    integer :: j  ! Iteration counter

    ! Initial setup
    yinit(1) = 12.5_qp  ! Initial scalar field value
    yinit(2) = 0.0_qp  ! Initial derivative of φ
    yinit(3) = hubb(yinit)  ! Initial Hubble parameter H
    yinit(4) = 0.0_qp  ! ln(a), conformal time parameter
    kphys = 1.0q3 * yinit(3)  ! Initial physical mode
    
    dt = twopi * (1.0_qp / (1.0q3 * yinit(3))) / 60.0_qp  ! Initial time step based on `y(3)`
    
    change = .true.  ! Enable background evolution
    break = .true.  ! Enable perturbation evolution

    ! Open output file
    unit_output = 10
    open(unit=unit_output, file="Evolution_Cosmic_Perturbations_Large.dat", status="unknown", action="write", form="formatted")

    ! 🔹 Background evolution until ln(a) >= 5.0
    y = yinit
    do while (change)
        call gl10(y, dt)  ! Integration with 10th-order Gauss-Legendre
        if (y(4) >= 5.0_qp) then
            yinit = y  ! Save adjusted initial state
            kcom = kphys * exp(y(4))  ! Calculate conformal mode

            ! Initialize Sasaki-Mukhanov scalar perturbations
            yinit(5) = 1.0_qp / sqrt(2.0_qp * sqrt(kcom**2 - zpp_over_z(yinit)))  
            yinit(6) = 0.0_qp  
            yinit(7) = sqrt(kcom**2 - zpp_over_z(yinit))  
            change = .false.  ! End initial adjustment
        end if
    end do

    ! 🔹 Main perturbation evolution until ε >= 1.00
    y = yinit
    j = 0
    do while (break)
        call gl10(y, dt / 20.0_qp)  ! Integration with refined step
        j = j + 1
        if (epsilon(y) >= 1.0_qp) then 
            break = .false.  ! Exit condition
        end if

        ! Save results every 1000 steps
        if (mod(j, 1000) == 0) then
            write(unit_output, '(12(ES45.34E4,1X))') y(1), y(2), y(3), y(4), y(5), y(6), y(7), &
                                            epsilon(y), hubb(y), zpp_over_z(y), kcom, kphys
        end if
    end do

    ! Close output file
    close(unit=unit_output)
    write(*,*) "Simulation finished."

contains

    ! 🔹 System equations
    subroutine evalf(y, dydx)
        real(qp), intent(in) :: y(nvar)
        real(qp), intent(out) :: dydx(nvar)

        ! Equations of motion
        dydx(1) = y(2)
        dydx(2) = -3.0_qp * y(3) * y(2) - Vprime(y)
        dydx(3) = -0.5_qp * y(2)**2
        dydx(4) = y(3)
        dydx(5) = y(6) * exp(-y(4))
        dydx(6) = -(kcom**2 - zpp_over_z(y) - y(7)**2) * y(5) * exp(-y(4)) ! Equivalent to -(kcom**2 - zpp_over_z(y) - 0.25_qp * y(5)**(-4)) * y(5) * exp(-y(4))
        dydx(7) = -2.0_qp * (y(7) * y(6) / y(5)) * exp(-y(4))
    end subroutine evalf

    ! 🔹 10th-order Gauss-Legendre integrator
    subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 7
        real(qp), intent(inout) :: y(n)
        real(qp), intent(in) :: dt
        real(qp) :: g(n, s)
        integer :: i, k

        ! Butcher tableau for 10th-order Gauss-Legendre method
        real(qp), parameter :: a(s,s) = reshape([ &
        0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
        1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
        1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
        1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
        1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
        1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
        1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
        4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
        2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
        1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
        1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
        2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
        0.5923172126404727187856601017997934066Q-1], [s,s])
        real(qp), parameter ::   b(s) = [ &
        1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
        2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
        1.1846344252809454375713202035995868132Q-1]
        ! Trial iteration
        g = 0.0_qp
        do k = 1, 16
            g = matmul(g, a)
            do i = 1, s
                call evalf(y + g(:, i) * dt, g(:, i))
            end do
        end do

        ! Update solution
        y = y + matmul(g, b) * dt
    end subroutine gl10

    ! 🔹 Potential-related functions
    function Vphi(y)
        real(qp), intent(in) :: y(nvar)
        real(qp) :: Vphi
        Vphi = 1.0q-14 * y(1)**4 / 4.0_qp
    end function Vphi

    function Vprime(y)
        real(qp), intent(in) :: y(nvar)
        real(qp) :: Vprime
        Vprime = 1.0q-14 * y(1)**3
    end function Vprime
    
    function Vprimeprime(y)
        real(qp), intent(in) :: y(nvar)
        real(qp) :: Vprimeprime
        Vprimeprime = 3.0_qp * 1.0q-14 * y(1)**2
    end function Vprimeprime

    ! 🔹 Hubble function
    function hubb(y)
        real(qp), intent(in) :: y(nvar)
        real(qp) :: hubb
        hubb = sqrt(y(2)**2 / 6.0_qp + Vphi(y) / 3.0_qp)
    end function hubb

    ! 🔹 Slow-roll parameter ε
    function epsilon(y)
        real(qp), intent(in) :: y(nvar)
        real(qp) :: epsilon
        epsilon = y(2)**2 / (2.0_qp * y(3)**2)
    end function epsilon

    ! 🔹 Calculation of z''/z
    function zpp_over_z(y)
        real(qp), intent(in) :: y(nvar)
        real(qp) :: zpp_over_z
        zpp_over_z = exp(2.0_qp * y(4)) * (2.0_qp * y(3)**2 - Vprimeprime(y) - &
                     2.0_qp * y(2) * Vprime(y) / y(3) - 3.5_qp * y(2)**2 * y(3) + &
                     0.5_qp * exp(y(4)) * y(2)**4 / y(3)**2)
    end function zpp_over_z

end program simanhosc
