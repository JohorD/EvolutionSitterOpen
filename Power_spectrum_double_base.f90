program ps_decoherence
    implicit none
    real, parameter :: twopi = 6.28318530717958647692528676655900577, frac_e = 0.7, frac_o = 0.8!twopi = 2.0 * acos(-1.0)  ! Definir a partir de acos
    integer, parameter :: nvar = 7  ! Número de ecuaciones en el sistema

    ! Variables principales
    real :: kphys_ori, kphys_ac, kcom, dt, N_ref, Nprima  ! Modo físico, modo conforme y parámetro N_ref
    real :: y(nvar), yinit(nvar), yback(nvar)  ! Variables del sistema
    integer :: unit_output, unit_summary  ! Unidades para archivos de salida
    logical :: change, break, elipse, back  ! Variables de control
    integer :: j, iter  ! Contador de iteraciones y de archivos

    ! Configuración para el entorno
    real, parameter :: a_ref = 1.0  ! 1.0_dp
    real, parameter :: p = 2.1 ! 2.1 6.1 --2.5 6.1
    integer, parameter :: mod_pref = 2  ! 0 = Close, 1 = Open Mod 1, 2 = Open Mod 2, 3 = Open Mod 3, ...

    real :: k_gamma, lE, lO
    character(len=100) :: filename  ! Nombre del archivo de salida

    ! Inicializar la variable elipse
    elipse = .false.  ! Cambiar a .false. para solo crear el archivo de resumen
        
    ! Evolucionar el background desde N = -5 y guardar en u vector yback
    
    back = .true. ! Comenzamos la preevlución

    ! Configuración inicial
    yback(1) = 25.0  ! Campo escalar inicial 30
    yback(2) = 0.0  ! Derivada de φ inicial
    yback(3) = hubb(yback)  ! Cálculo inicial de Hubble H
    yback(4) = -5.0  ! Establecer ln(a) inicial a N_ref
    
    Nprima = 5.0
    


    do while (back)
        call gl8(yback, dt, kcom)  ! Integración con Gauss-Legendre de orden 10
        if (yback(4) >= 0.0) then
            back = .false.  ! Terminar ajuste inicial
        end if
    end do

    ! Terminamos la pre-evolución en N = 0 y la llamamos en lo que sigue

!!!! Npi 8--20---40

    
    ! Abrir archivo de resumen para guardar los últimos valores
    unit_summary = 20
    open(unit=unit_summary, file="PS_EO78_Nmid_TEST_.dat", status="unknown", action="write", form="formatted")

    ! Bucle para iterar sobre los valores de N_ref
    do iter = 1, 100!0   !100
        N_ref = 0.1 * iter   !0.6 * iter

        ! Configuración inicial
        yinit(1) = yback(1) ! Campo escalar φ pre-evolucioanada
        yinit(2) = yback(2) ! Derivada de φ pre-evolucioanada
        yinit(3) = yback(3) ! Valor de H pre-evolucioanada
        yinit(4) = yback(4) ! ln(a) pre-evolucioanada  N = 0
        kphys_ori = 1000.0 * yback(3)

        ! 🔹 Modo físico basado en yback(3)
        ! 🔹 Paso de tiempo inicial basado en yback(3)

        if (Nprima > N_ref) then
            kphys_ac = kphys_ori
            dt = twopi * (1.0/ (1000.0 * yback(3))) / 60.0
        else
            kphys_ac = kphys_ori * exp(Nprima)
            dt = twopi * (1.0/ (1000.0 * yback(3))) / 180.0
            !!!!
            lE = (twopi/kphys_ac) * exp(Nprima*(frac_e-1.0))
            lO = (twopi/kphys_ac) * exp(Nprima*(frac_o-1.0))
            !!!!
        end if



        change = .true.  ! Activar evolución del background
        break = .true.  ! Activar evolución de perturbaciones

        ! Generar nombre de archivo dinámico
        if (elipse) then
            write(filename, '("Elipse_Open_N_ref_", I1, ".dat")') iter
            ! Abrir archivo para guardar datos si elipse es verdadero
            open(unit=unit_output, file=filename, status="unknown", action="write", form="formatted")
        end if

        ! 🔹 Evolución del background hasta ln(a) >= N_ref, hasta llegar a los kcom de interés

        y = yinit
        do while (change)
            call gl8(y, dt, kcom)  ! Integración con Gauss-Legendre de orden 8
            if (y(4) >= N_ref) then
                
                yinit = y  ! Guardar estado inicial ajustado
                
                kcom = kphys_ac * exp(y(4))  ! Cálculo del modo conforme

                ! Inicialización de perturbaciones escalares de Sasaki-Mukhanov
                yinit(5) = 1.0 / sqrt(2.0 * sqrt(kcom*kcom - zpp_over_z(yinit)))
                yinit(6) = 0.0
                yinit(7) = sqrt(kcom**2 - zpp_over_z(yinit))

                ! Inicialización de la intensidad de la fuente
                k_gamma = 10.0 * kcom

                change = .false.  ! Terminar ajuste inicial
            end if
        end do

        ! 🔹 Evolución principal de perturbaciones hasta ε >= 1.0
        y = yinit
        j = 0
        do while (break)
            call gl8(y, dt, kcom)  ! Integración con paso refinado 20
            j = j + 1
            if (epsilon(y) >= 1.0) then
                break = .false.  ! Condición de salida
            end if

            ! Guardar resultados cada 1000 pasos si elipse es verdadero
            if (elipse .and. mod(j, 1000) == 0) then
                write(unit_output, '(13(ES15.7E3,1X))') y(1), y(2), y(3), y(4), y(5), y(6), y(7), &
                                                        epsilon(y), hubb(y), zpp_over_z(y), kcom, kphys_ac, z_func(y)
            end if
        end do

        ! Guardar los últimos valores en el archivo de resumen
        write(unit_summary, '(13(ES15.7E3,1X))') y(1), y(2), y(3), y(4), y(5), y(6), y(7), &
                                                   epsilon(y), hubb(y), zpp_over_z(y), kcom, kphys_ac, z_func(y)

        ! Cerrar archivo de salida si elipse es verdadero y limpiar variables
        if (elipse) then
            close(unit=unit_output)
        end if
        write(*,*) "Simulación finalizada para N_ref = ", N_ref

        ! Reinicializar variables para la siguiente iteración
        y = yback
        yinit = yback
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
        dydx(5) = y(6) * exp(-y(4))
        dydx(6) = -(kcom*kcom - zpp_over_z(y) - y(7)*y(7)) * y(5) * exp(-y(4))
        dydx(7) = -2.0 * (y(7) * y(6) / y(5)) * exp(-y(4)) + &
                  source_open(mod_pref, y) * exp(-y(4)) / (2.0 * y(7) * y(5)*y(5))  ! Adición de la fuente
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
                         (ventana(y, lE, kcom))
             
           case (2)  ! Open_Mod_2_Smo
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
                         (ventana(y, lE, kcom) - ventana(y, lO, kcom))
           case(3) ! Hysteresis
                source = (k_gamma*k_gamma) * exp((p - 3.0) * y(4)) * (a_ref**(3.0 - p)) * &
         (ventana(y, lE, kcom) - 2.0*ventana(y, lO, kcom) +  ventana(y, 2.0*lO-lE, kcom))

            case default
                print *, "¡Error! mod_pref debe ser 0, 1, 2, 3, 4 o 5."
                stop
        end select
    end function source_open

end program ps_decoherence
