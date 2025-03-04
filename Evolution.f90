program simanhosc
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: dt = 0.00001  
    real(dp), parameter :: k_gamma_over_k = 10.0_dp
    real(dp), parameter :: x_ref = 1.0_dp
    real(dp), parameter :: p = 6.1_dp
    real(dp), parameter :: lE_H = 0.1_dp
    real(dp), parameter :: lO_H = 0.2_dp
    real(dp), parameter :: NE_H = -3.0_dp
    real(dp), parameter :: NO_H = -2.0_dp
    real(dp), parameter :: h_tran = 20.0_dp  ! Parámetro de transición para Heaviside suave
    integer, parameter :: mod_pref = 1  ! 0 = Close, 1 = Open Mod 1, 2 = Open Mod 2, 3 = Open Mod 3
    integer :: l, l_max, unit_output
    real(dp) :: y(3), Ne, Ne_final

    ! Configuración del rango de evolución
    Ne = -10.0_dp   ! Ne inicial
    Ne_final = 6.0_dp  ! Ne final deseado

    ! Calcular el número de iteraciones necesarias
    l_max = int((Ne_final - Ne) / dt)

    ! Condiciones iniciales
    y(1) = 1.0_dp
    y(2) = 0.0_dp
    y(3) = 1.0_dp

    ! Abrir archivo para guardar los datos
    unit_output = 10
    open(unit=unit_output, file="Evolution_Sitter_Mod_3.dat", status="unknown", action="write", form="formatted")

    ! Bucle de evolución
    do l = 0, l_max
        if (mod(l, 1000) == 0) then
            write(unit_output, '(4ES24.16)') Ne, y(1), y(2), y(3)
        end if
        
        call gl10(y, dt, Ne)
        Ne = Ne + dt   ! Actualizar Ne

        ! Verificar que los valores de y no sean demasiado grandes
        if (maxval(abs(y)) > 1.0d200) then
            print *, "¡Advertencia! Valores demasiado grandes en Ne =", Ne
            exit
        end if
    end do

    ! Cerrar archivo
    close(unit=unit_output)

    write(*,*) "Simulación finalizada. Evolucionó desde Ne =", -10.0_dp, " hasta Ne =", Ne

contains

    ! Función suavizada de Heaviside con tanh
    pure function SmoothHeaviside(x) result(H)
        real(dp), intent(in) :: x
        real(dp) :: H
        H = 0.5_dp + 0.5_dp * tanh(h_tran * x)
    end function SmoothHeaviside

    ! Función que devuelve el término adicional en dy3/dNe según mod_pref
    function source_open(mod_pref, Ne) result(source)
        integer, intent(in) :: mod_pref
        real(dp), intent(in) :: Ne
        real(dp) :: source

        select case (mod_pref)
            case (0)  ! Sitter_Close (sin término adicional)
                source = 0.0_dp

            case (1)  ! Sitter_Open_Mod_1_Smo
                source = 2.0_dp * (k_gamma_over_k**2) * (x_ref**(p-3.0_dp)) * &
                         SmoothHeaviside(1.0_dp - lE_H * exp(-Ne)) * &
                         exp(Ne*(p - 4.0_dp))

            case (2)  ! Sitter_Open_Mod_2_Smo
                source = 2.0_dp * (k_gamma_over_k**2) * (x_ref**(p-3.0_dp)) * &
                         (SmoothHeaviside(1.0_dp - lE_H * exp(-Ne)) - &
                          SmoothHeaviside(1.0_dp - lO_H * exp(-Ne))) * &
                         exp(Ne*(p - 4.0_dp))

            case (3)  ! Sitter_Open_Mod_3_Smo
                source = 2.0_dp * (k_gamma_over_k**2) * (x_ref**(p-3.0_dp)) * &
                         (SmoothHeaviside(1.0_dp - lE_H * exp(-Ne)) - &
                          SmoothHeaviside(1.0_dp - lO_H * exp(-Ne))) * &
                         (SmoothHeaviside(exp(-NE_H) - exp(-Ne)) - &
                          SmoothHeaviside(exp(-NO_H) - exp(-Ne))) * &
                         exp(Ne*(p - 4.0_dp))

            case default
                print *, "¡Error! mod_pref debe ser 0, 1, 2 o 3."
                stop
        end select
    end function source_open

    ! Ecuaciones del sistema
    subroutine evalf(y, dydx, Ne)
        real(dp), intent(in) :: y(3), Ne
        real(dp), intent(out) :: dydx(3)

        dydx(1) = 2.0_dp * exp(-Ne) * y(2)
        dydx(2) = exp(-Ne) * y(3) - (exp(-Ne) - 2.0_dp * exp(Ne)) * y(1)
        dydx(3) = -2.0_dp * (exp(-Ne) - 2.0_dp * exp(Ne)) * y(2) + source_open(mod_pref, Ne)
    end subroutine evalf

    ! Integrador de Gauss-Legendre de orden 10
    subroutine gl10(y, dt, Ne)
        integer, parameter :: s = 5, n = 3
        real(dp), intent(inout) :: y(n), Ne
        real(dp), intent(in) :: dt
        real(dp) :: g(n, s)
        integer :: i, k

        ! Butcher tableau para 10º orden Gauss-Legendre
        real(dp), parameter :: a(s,s) = reshape([ &
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

        real(dp), parameter :: b(s) = [ &
        1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
        2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
        1.1846344252809454375713202035995868132Q-1]

        g = 0.0_dp
        do k = 1, 16
            g = matmul(g, a)
            do i = 1, s
                call evalf(y + g(:, i) * dt, g(:, i), Ne)
            end do
        end do

        y = y + matmul(g, b) * dt
    end subroutine gl10

end program simanhosc
