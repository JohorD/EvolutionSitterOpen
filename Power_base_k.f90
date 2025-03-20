program simanhosc
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)  ! Definir precisión aumentada
    real(dp), parameter :: twopi = 2.0_dp * acos(-1.0_dp)  ! Definir a partir de acos
    integer, parameter :: nvar = 7  ! Número de ecuaciones en el sistema
    integer, parameter :: num_points = 8  ! Número de puntos para ln(a)
    real(dp), parameter :: N_ref_ini = 1.0_dp  !5.0_dp  ! Valor inicial de N_ref
    real(dp), parameter :: N_ref_fin = 5.0_dp  !40.0_dp  ! Valor final de N_ref

    ! Variables principales
    real(dp) :: kphys, kcom, dt, N_ref  ! Modo físico, modo conforme y parámetro N_ref
    real(dp) :: y(nvar), yinit(nvar)  ! Variables del sistema
    integer :: unit_output, unit_summary  ! Unidades para archivos de salida
    logical :: change, break, elipse  ! Variables de control
    integer :: j, iter, cont_N  ! Contador de iteraciones y de archivos
    real(dp) :: kcom_values(num_points), ln_a_values(num_points)
    character(len=100) :: filename  ! Nombre del archivo de salida

    ! Configuración para el entorno
    real(dp), parameter :: a_ref = 1.0_dp  ! 1.0_dp
    real(dp), parameter :: p = 2.1_dp  ! 2.1 6.1 --2.5 6.1
    real(dp), parameter :: h_tran = 20.0_dp  ! Parámetro de transición para Heaviside suave
    integer, parameter :: mod_pref = 5 ! 0 = Close, 1 = Open Mod 1, 2 = Open Mod 2, 3 = Open Mod 3, ...

    real(dp) :: k_gamma, lE, lO, NE, NO

    ! Inicializar la variable elipse
    elipse = .false.  ! Cambiar a .false. para solo crear el archivo de resumen

    ! Calcular kphys a partir de yinit(3)
    yinit(1) = 18.0_dp  ! Campo escalar inicial
    yinit(2) = 0.0_dp  ! Derivada de φ inicial
    yinit(3) = hubb(yinit)  ! Cálculo inicial de Hubble H
    yinit(4) = -5.0_dp   !0.0_dp   ! Establecer ln(a) inicial
    kphys = 1.0d3 * yinit(3)  ! Modo físico inicial

    ! Generar valores de kcom uniformemente espaciados
    kcom_values = [(kphys * exp(N_ref_ini) + (cont_N - 1) * (kphys * exp(N_ref_fin) - &
                  kphys * exp(N_ref_ini)) / (num_points - 1), cont_N = 1, num_points)]

    ! Calcular los valores correspondientes de ln(a)
    ln_a_values = log(kcom_values / kphys)

    ! Abrir archivo de resumen para guardar los últimos valores
    unit_summary = 20
    open(unit=unit_summary, file="Summary_Evolution_Cosmic_Perturbations_N_uni_5.dat", &
    status="unknown", action="write", form="formatted")

    ! Bucle para iterar sobre los valores de N_ref
    do iter = 1, num_points
        N_ref = ln_a_values(iter)

        ! Configuración inicial
        yinit(1) = 18.0_dp  ! Campo escalar inicial 30
        yinit(2) = 0.0_dp  ! Derivada de φ inicial
        yinit(3) = hubb(yinit)  ! Cálculo inicial de Hubble H
        yinit(4) = -5.0_dp   !0.0_dp  ! Establecer ln(a) inicial en 0 para usar el atractor
        kphys = 1.0d3 * yinit(3)  ! Modo físico inicial

        dt = twopi * (1.0_dp / (1.0d3 * yinit(3))) / 60.0_dp  ! 🔹 Paso de tiempo inicial basado en y(3)

        change = .true.  ! Activar evolución del background
        break = .true.  ! Activar evolución de perturbaciones

        ! Generar nombre de archivo dinámico
        if (elipse) then
            write(filename, '("Evolution_Cosmic_Perturbations_18_open_2_N_ref_", I1, ".dat")') iter
            ! Abrir archivo para guardar datos si elipse es verdadero
            open(unit=unit_output, file=filename, status="unknown", action="write", form="formatted")
        end if

        ! 🔹 Evolución del background hasta ln(a) >= N_ref
        y = yinit
        do while (change)
            call gl10(y, dt, kcom)  ! Integración con Gauss-Legendre de orden 10
            if (y(4) >= N_ref) then
                yinit = y  ! Guardar estado inicial ajustado
                kcom = (kphys/twopi) * exp(y(4))  ! Cálculo del modo conforme

                ! Inicialización de perturbaciones escalares de Sasaki-Mukhanov
                yinit(5) = 1.0_dp / sqrt(2.0_dp * sqrt(kcom**2 - zpp_over_z(yinit)))
                yinit(6) = 0.0_dp
                yinit(7) = sqrt(kcom**2 - zpp_over_z(yinit))

                ! Inicialización de la fuente
                k_gamma = 10.0_dp * kcom
                lE = 0.7_dp * log(1 / yinit(3))
                lO = 0.8_dp * log(1 / yinit(3))
                NE = 8.0_dp
                NO = 9.0_dp

                change = .false.  ! Terminar ajuste inicial
            end if
        end do

        ! 🔹 Evolución principal de perturbaciones hasta ε >= 1.0
        y = yinit
        j = 0
        do while (break)
            call gl10(y, dt / 20.0_dp, kcom)  ! Integración con paso refinado 20
            j = j + 1
            if (epsilon(y) >= 1.0_dp) then
                break = .false.  ! Condición de salida
            end if

            ! Guardar resultados cada 1000 pasos si elipse es verdadero
            if (elipse .and. mod(j, 1000) == 0) then
                write(unit_output, '(13(ES15.7E3,1X))') y(1), y(2), y(3), y(4), y(5), y(6), y(7), &
                                                        epsilon(y), hubb(y), zpp_over_z(y), kcom, kphys, z_func(y)
            end if
        end do

        ! Guardar los últimos valores en el archivo de resumen
        write(unit_summary, '(13(ES15.7E3,1X))') y(1), y(2), y(3), y(4), y(5), y(6), y(7), &
                                                   epsilon(y), hubb(y), zpp_over_z(y), kcom, kphys, z_func(y)

        ! Cerrar archivo de salida si elipse es verdadero y limpiar variables
        if (elipse) then
            close(unit=unit_output)
        end if
        write(*,*) "Simulación finalizada para N_ref = ", N_ref

        ! Reinicializar variables para la siguiente iteración
        y = 0.0_dp
        yinit = 0.0_dp
    end do

    ! Cerrar archivo de resumen
    close(unit=unit_summary)

contains

    ! 🔹 Ecuaciones del sistema
    subroutine evalf(y, dydx, kcom)
        real(dp), intent(in) :: y(nvar)
        real(dp), intent(out) :: dydx(nvar)
        real(dp), intent(in) :: kcom

        ! Ecuaciones de movimiento del sistema
        dydx(1) = y(2)
        dydx(2) = -3.0_dp * y(3) * y(2) - Vprime(y)
        dydx(3) = -0.5_dp * y(2)**2
        dydx(4) = y(3)
        dydx(5) = y(6) * exp(-y(4))
        dydx(6) = -(kcom**2 - zpp_over_z(y) - y(7)**2) * y(5) * exp(-y(4))
        dydx(7) = -2.0_dp * (y(7) * y(6) / y(5)) * exp(-y(4)) + &
                  source_open(mod_pref, y, kcom) * exp(-y(4)) / (2.0_dp * y(7) * y(5)**2)  ! Adición de la fuente
    end subroutine evalf

    ! 🔹 Integrador Gauss-Legendre de orden 10
    subroutine gl10(y, dt, kcom)
        integer, parameter :: s = 5, n = 7
        real(dp), intent(inout) :: y(n)
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: kcom
        real(dp) :: g(n, s)
        integer :: i, k

        ! Butcher tableau for 10th order Gauss-Legendre method
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
        ! Iteración de prueba
        g = 0.0_dp
        do k = 1, 16
            g = matmul(g, a)
            do i = 1, s
                call evalf(y + g(:, i) * dt, g(:, i), kcom)
            end do
        end do

        ! Actualizar solución
        y = y + matmul(g, b) * dt
    end subroutine gl10

    ! 🔹 Funciones relacionadas con el potencial
    function Vphi(y)
        real(dp), intent(in) :: y(nvar)
        real(dp) :: Vphi
        Vphi = 1.0d-14 * y(1)**4 / 4.0_dp
    end function Vphi

    function Vprime(y)
        real(dp), intent(in) :: y(nvar)
        real(dp) :: Vprime
        Vprime = 1.0d-14 * y(1)**3
    end function Vprime

    function Vprimeprime(y)
        real(dp), intent(in) :: y(nvar)
        real(dp) :: Vprimeprime
        Vprimeprime = 3.0_dp*1.0d-14 * y(1)**2
    end function Vprimeprime

    ! 🔹 Función de Hubble
    function hubb(y)
        real(dp), intent(in) :: y(nvar)
        real(dp) :: hubb
        hubb = sqrt(y(2)**2 / 6.0_dp + Vphi(y) / 3.0_dp)
    end function hubb

    ! 🔹 Parámetro de slow-roll ε
    function epsilon(y)
        real(dp), intent(in) :: y(nvar)
        real(dp) :: epsilon
        epsilon = y(2)**2 / (2.0_dp * y(3)**2)
    end function epsilon

    ! 🔹 Cálculo de z''/z
    function zpp_over_z(y)
        real(dp), intent(in) :: y(nvar)
        real(dp) :: zpp_over_z
        zpp_over_z = exp(2.0_dp * y(4)) * (2.0_dp * y(3)**2 - Vprimeprime(y) - &
                     2.0_dp * y(2) * Vprime(y) / y(3) - 3.5_dp * y(2)**2 + &
                     0.5_dp * y(2)**4 / y(3)**2)
    end function zpp_over_z

    ! 🔹 Cálculo de z
    function z_func(y)
        real(dp), intent(in) :: y(nvar)
        real(dp) :: z_func
        z_func = y(2) * exp(y(4)) / y(3)
    end function z_func

    ! 🔹 Función suavizada de Heaviside con tanh
    pure function SmoothHeaviside(x) result(H)
        real(dp), intent(in) :: x
        real(dp) :: H
        H = 0.5_dp + 0.5_dp * tanh(h_tran * x)
    end function SmoothHeaviside

    ! 🔹 Función que devuelve el término adicional en dy3/dNe según mod_pref
    function source_open(mod_pref, y, kcom) result(source)
        integer, intent(in) :: mod_pref
        real(dp), intent(in) :: y(nvar)
        real(dp), intent(in) :: kcom
        real(dp) :: source

        select case (mod_pref)
            case (0)  ! Sitter_Close (sin término adicional)
                source = 0.0_dp

            case (1)  ! Open_Mod_1_Smo
                source = (k_gamma**2) * exp((p - 3.0_dp) * y(4)) * (a_ref**(3.0_dp - p)) * &
                         SmoothHeaviside(1.0_dp - lE / (y(4) - log(kcom)))

            case (2)  ! Open_Mod_2_Smo
                source = (k_gamma**2) * exp((p - 3.0_dp) * y(4)) * (a_ref**(3.0_dp - p)) * &
                         (SmoothHeaviside(1.0_dp - lE / (y(4) - log(kcom))) - &
                         SmoothHeaviside(1.0_dp - lO / (y(4) - log(kcom))))

            case (3)  ! Open_Mod_3_Smo
                source = (k_gamma**2) * exp((p - 3.0_dp) * y(4)) * (a_ref**(3.0_dp - p)) * &
                         SmoothHeaviside(1.0_dp - NE / y(4))

            case (4)  ! Open_Mod_4_Smo
                source = (k_gamma**2) * exp((p - 3.0_dp) * y(4)) * (a_ref**(3.0_dp - p)) * &
                         (SmoothHeaviside(1.0_dp - NE / y(4)) - &
                         SmoothHeaviside(1.0_dp - NO / y(4)))

            case (5)  ! Open_Mod_5_Smo
                source = (k_gamma**2) * exp((p - 3.0_dp) * y(4)) * (a_ref**(3.0_dp - p)) * &
                         (SmoothHeaviside(1.0_dp - lE / (y(4) - log(kcom))) - &
                         SmoothHeaviside(1.0_dp - lO / (y(4) - log(kcom)))) * &
                         (SmoothHeaviside(1.0_dp - NE / y(4)) - SmoothHeaviside(1.0_dp - NO / y(4)))

            case default
                print *, "¡Error! mod_pref debe ser 0, 1, 2, 3, 4 o 5."
                stop
        end select
    end function source_open

end program simanhosc
