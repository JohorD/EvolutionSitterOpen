program simanhosc
    implicit none
    real, parameter :: twopi = 2.0 * acos(-1.0)  ! Definir a partir de acos
    integer, parameter :: nvar = 7  ! Número de ecuaciones en el sistema

    ! Variables principales
    real :: kphys, kcom, dt, N_ref  ! Modo físico, modo conforme y parámetro N_ref
    real :: y(nvar), yinit(nvar)  ! Variables del sistema
    integer :: unit_output, unit_summary  ! Unidades para archivos de salida
    logical :: change, break, elipse  ! Variables de control
    integer :: j, iter  ! Contador de iteraciones y de archivos

    ! Configuración para el entorno
    real, parameter :: a_ref = 1.0  ! 1.0_dp
    real, parameter :: p = 2.1  ! 2.1 6.1 --2.5 6.1
    real, parameter :: h_tran = 20.0  ! Parámetro de transición para Heaviside suave
    integer, parameter :: mod_pref = 7  ! 0 = Close, 1 = Open Mod 1, 2 = Open Mod 2, 3 = Open Mod 3, ...

    real :: k_gamma, lE, lO, NE, NO
    character(len=100) :: filename  ! Nombre del archivo de salida

    ! Inicializar la variable elipse
    elipse = .false.  ! Cambiar a .false. para solo crear el archivo de resumen

    ! Abrir archivo de resumen para guardar los últimos valores
    unit_summary = 20
    open(unit=unit_summary, file="Summary_Evolution_Cosmic_Perturbations_N_7_PP_EO78H.dat", status="unknown", action="write", form="formatted")

    ! Bucle para iterar sobre los valores de N_ref
    do iter = 1, 120   !100
        N_ref = 0.6 * iter   !0.6 * iter

        ! Configuración inicial
        yinit(1) = 25.0  ! Campo escalar inicial 30
        yinit(2) = 0.0  ! Derivada de φ inicial
        yinit(3) = hubb(yinit)  ! Cálculo inicial de Hubble H
        yinit(4) = 0.0  ! Establecer ln(a) inicial a N_ref
        kphys = 1000.0 * yinit(3)  ! Modo físico inicial

        dt = twopi * (1.0/ (1000.0 * yinit(3))) / 60.0  ! 🔹 Paso de tiempo inicial basado en y(3)

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
            call gl8(y, dt, kcom)  ! Integración con Gauss-Legendre de orden 10
            if (y(4) >= N_ref) then
                yinit = y  ! Guardar estado inicial ajustado
                kcom = kphys * exp(y(4))  ! Cálculo del modo conforme

                ! Inicialización de perturbaciones escalares de Sasaki-Mukhanov
                yinit(5) = 1.0
                yinit(6) = 0.0
                yinit(7) = 1.0

                ! Inicialización de la fuente
                k_gamma = 10.0 * kcom
                lE = - 0.7 * log(yinit(3))
                lO = - 0.8 * log(yinit(3))
                NE = 25.0
                NO = 45.0

                change = .false.  ! Terminar ajuste inicial
            end if
        end do

        ! 🔹 Evolución principal de perturbaciones hasta ε >= 1.0
        y = yinit
        j = 0
        do while (break)
            call gl8(y, dt / 20.0, kcom)  ! Integración con paso refinado 20
            j = j + 1
            if (epsilon(y) >= 1.0) then
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
        y = 0.0
        yinit = 0.0
    end do

    ! Cerrar archivo de resumen
    close(unit=unit_summary)

contains

    ! 🔹 Ecuaciones del sistema
    subroutine evalf(y, dydx, kcom)
        real, intent(in) :: y(nvar)
        real, intent(out) :: dydx(nvar)
        real, intent(in) ::kcom

        ! Ecuaciones de movimiento del sistema
        dydx(1) = y(2)
        dydx(2) = -3.0 * y(3) * y(2) - Vprime(y)
        dydx(3) = -0.5 * y(2)*y(2)
        dydx(4) = y(3)
        dydx(5) = 2.0*kcom*(y(6) + source_open(mod_pref, y)) * exp(-y(4))
        dydx(6) = (kcom*y(7)-y(5)*(kcom - zpp_over_z(y)/kcom))*  exp(-y(4))
        dydx(7) = -2.0 * (kcom - zpp_over_z(y)/kcom)*exp(-y(4))*y(6)  ! Adición de la fuente
    end subroutine evalf

    ! 🔹 Integrador Gauss-Legendre de orden 10
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

    ! 🔹 Funciones relacionadas con el potencial
    function Vphi(y)
        real, intent(in) :: y(nvar)
        real Vphi
        Vphi = 1.0d-14 * y(1)*y(1)*y(1)*y(1) / 4.0
    end function Vphi

    function Vprime(y)
        real, intent(in) :: y(nvar)
        real Vprime
        Vprime = 1.0d-14 * y(1)*y(1)*y(1)
    end function Vprime

    function Vprimeprime(y)
        real, intent(in) :: y(nvar)
        real Vprimeprime
        Vprimeprime = 3.0*1.0d-14 * y(1)*y(1)
    end function Vprimeprime

    ! 🔹 Función de Hubble
    function hubb(y)
        real, intent(in) :: y(nvar)
        real hubb
        hubb = sqrt(y(2)*y(2) / 6.0 + Vphi(y) / 3.0)
    end function hubb

    ! 🔹 Parámetro de slow-roll ε
    function epsilon(y)
        real, intent(in) :: y(nvar)
        real epsilon
        epsilon = y(2)*y(2) / (2.0 * y(3)*y(3))
    end function epsilon

    ! 🔹 Cálculo de z''/z
    function zpp_over_z(y)
real y(nvar), epsilon1, epsilon2, epsilon3, zpp_over_z, hubbp, hubbpp, hubbppp
real Vpp, ddot_phi, tdot_phi, eps2dot
    Vpp = Vprimeprime(y)
    ddot_phi = -3.d0*y(3)*y(2)-Vprime(y)
    epsilon1 = epsilon(y)
    hubbp = -0.5*y(2)**2
    tdot_phi = -3.0*hubbp*y(2)-3.0*y(3)*ddot_phi-Vpp*y(2)
    hubbpp = -y(2)*ddot_phi
    hubbppp = -ddot_phi**2-y(2)*tdot_phi
    epsilon2 = hubbpp/(y(3)*hubbp)-2.0*hubbp/y(3)**2
    eps2dot = hubbppp/(y(3)*hubbp)-(hubbpp*(hubbp**2+hubbpp*y(3)))/(y(3)*hubbp)**2-2.0*hubbpp/y(3)**2+4.0*hubbp**2/y(3)**3
    epsilon3 = eps2dot/(y(3)*epsilon2)

    zpp_over_z = exp(2.0*y(4))*y(3)**2*(2.d0-epsilon1+1.5*epsilon2+0.25*epsilon2**2-0.5*epsilon2*epsilon1+0.5*epsilon2*epsilon3)
    end function zpp_over_z

    ! 🔹 Cálculo de z
    function z_func(y)
        real, intent(in) :: y(nvar)
        real z_func
        z_func = y(2) * exp(y(4)) / y(3)
    end function z_func

    ! 🔹 Función suavizada de Heaviside con tanh
    pure function SmoothHeaviside(x) result(H)
        real, intent(in) :: x
        real H
        H = 0.5 + 0.5 * tanh(h_tran * x)
    end function SmoothHeaviside

! 🔹 Cálculo de Ventana, Fourier de Heaviside
function ventana(y, lE, kcom)
    real, intent(in) :: y(nvar)
    real, intent(in) :: lE
    real, intent(in) ::kcom
    real ventana
    ventana = (kcom**(-3)) * ( sin(kcom * lE * exp(-y(4))) - kcom * lE * exp(-y(4)) * &
              cos(kcom * lE * exp(-y(4))))
end function ventana

    ! 🔹 Función que devuelve el término adicional en dy3/dNe según mod_pref
    function source_open(mod_pref, y) result(source)
        integer, intent(in) :: mod_pref
        real, intent(in) :: y(nvar)
        real source

        select case (mod_pref)
            case (0)  ! Sitter_Close (sin término adicional)
                source = 0.0

            case (1)  ! Open_Mod_1_Smo
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
                         SmoothHeaviside(1.0 - lE / (y(4) - log(kcom)))

            case (2)  ! Open_Mod_2_Smo
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
                         (SmoothHeaviside(1.0 - lE / (y(4) - log(kcom))) - &
                         SmoothHeaviside(1.0 - lO / (y(4) - log(kcom))))

            case (3)  ! Open_Mod_3_Smo
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
                         SmoothHeaviside(1.0 - NE / y(4))

            case (4)  ! Open_Mod_4_Smo
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
                         (SmoothHeaviside(1.0 - NE / y(4)) - &
                         SmoothHeaviside(1.0 - NO / y(4)))

            case (5)  ! Open_Mod_5_Smo
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
                         (SmoothHeaviside(1.0 - lE / (y(4) - log(kcom))) - &
                         SmoothHeaviside(1.0 - lO / (y(4) - log(kcom)))) * &
                         (SmoothHeaviside(1.0 - NE / y(4)) - SmoothHeaviside(1.0 - NO / y(4)))

           case (6)  ! Open_Mod_6_Smo
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
                         (ventana(y, lE, kcom))
             
           case (7)  ! Open_Mod_7_Smo
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
                         (ventana(y, lE, kcom) - ventana(y, lO, kcom))
           case(8) ! Hysteresis
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
         (ventana(y, lE, kcom) - 2.0*ventana(y, lO, kcom) +  ventana(y, 2.0*lO-lE, kcom))

            case default
                print *, "¡Error! mod_pref debe ser 0, 1, 2, 3, 4 o 5."
                stop
        end select
    end function source_open

end program simanhosc
