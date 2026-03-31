module parameters
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: GAMMA = 1.4_dp
    integer, parameter :: NX = 500
    integer, parameter :: NG = 3 
    integer, parameter :: NTOT = NX + 2 * NG
    real(dp), parameter :: DX = 1.0_dp / real(NX, dp)
    real(dp), parameter :: CFL = 0.4_dp
    real(dp), parameter :: T_END = 0.17_dp
end module parameters

module exact_solver
    use parameters
    implicit none

contains
    ! Function to find the root for the Exact Riemann Solver
    function riemann_f(P, p1, rho1, p4, rho4, u1, u4) result(res)
        real(dp), intent(in) :: P, p1, rho1, p4, rho4, u1, u4
        real(dp) :: res, c1, c4, PRL, CRL, MACHLEFT, term
        
        c1 = sqrt(GAMMA * p1 / rho1)
        c4 = sqrt(GAMMA * p4 / rho4)
        PRL = p4 / p1
        CRL = c4 / c1
        MACHLEFT = (u1 - u4) / c1
        
        term = 1.0_dp + MACHLEFT*(GAMMA-1.0_dp)/2.0_dp - &
               (GAMMA-1.0_dp)*CRL*(P-1.0_dp) / &
               sqrt(2.0_dp*GAMMA*(GAMMA-1.0_dp + (GAMMA+1.0_dp)*P))
               
        if (term < 0.0_dp) then
            res = -1.0e10_dp
        else
            res = (term**(2.0_dp*GAMMA/(GAMMA-1.0_dp))) / P - PRL
        end if
    end function riemann_f

    ! Bisection method to solve for P34
    function solve_p34(p1, rho1, p4, rho4, u1, u4) result(P34)
        real(dp), intent(in) :: p1, rho1, p4, rho4, u1, u4
        real(dp) :: P34, a, b, c, fa, fc
        integer :: iter
        
        a = 1.0e-6_dp
        b = 1.0e4_dp
        fa = riemann_f(a, p1, rho1, p4, rho4, u1, u4)
        
        do iter = 1, 100
            c = 0.5_dp * (a + b)
            fc = riemann_f(c, p1, rho1, p4, rho4, u1, u4)
            
            if (abs(fc) < 1.0e-8_dp .or. 0.5_dp*(b-a) < 1.0e-8_dp) exit
            
            if (sign(1.0_dp, fc) == sign(1.0_dp, fa)) then
                a = c
                fa = fc
            else
                b = c
            end if
        end do
        P34 = c
    end function solve_p34

    subroutine get_exact_solution(x, t, rho1, u1, p1, rho4, u4, p4, rho_ex, u_ex, p_ex)
        real(dp), dimension(NX), intent(in) :: x
        real(dp), intent(in) :: t, rho1, u1, p1, rho4, u4, p4
        real(dp), dimension(NX), intent(out) :: rho_ex, u_ex, p_ex
        real(dp) :: p34, p3, rho3, rho2, u2, c1, c2, c4, alpha
        real(dp) :: x0, spos, conpos, pos1, pos2
        integer :: i
        
        alpha = (GAMMA + 1.0_dp) / (GAMMA - 1.0_dp)
        c1 = sqrt(GAMMA * p1 / rho1)
        c4 = sqrt(GAMMA * p4 / rho4)
        
        p34 = solve_p34(p1, rho1, p4, rho4, u1, u4)
        
        p3 = p34 * p4
        rho3 = rho4 * (1.0_dp + alpha * p34) / (alpha + p34)
        rho2 = rho1 * (p34 * p4 / p1)**(1.0_dp / GAMMA)
        u2 = u1 - u4 + (2.0_dp/(GAMMA-1.0_dp))*c1*(1.0_dp - (p34*p4/p1)**((GAMMA-1.0_dp)/(2.0_dp*GAMMA)))
        c2 = sqrt(GAMMA * p3 / rho2)
        
        x0 = 0.5_dp
        spos   = x0 + t * c4 * sqrt((GAMMA-1.0_dp)/(2.0_dp*GAMMA) + (GAMMA+1.0_dp)/(2.0_dp*GAMMA)*p34) + t*u4
        conpos = x0 + u2 * t + t * u4
        pos1   = x0 + (u1 - c1) * t
        pos2   = x0 + (u2 + u4 - c2) * t

        do i = 1, NX
            if (x(i) <= pos1) then
                p_ex(i) = p1; rho_ex(i) = rho1; u_ex(i) = u1
            else if (x(i) <= pos2) then
                p_ex(i) = p1 * (1.0_dp + (pos1 - x(i)) / (c1 * alpha * t))**(2.0_dp * GAMMA / (GAMMA - 1.0_dp))
                rho_ex(i) = rho1 * (1.0_dp + (pos1 - x(i)) / (c1 * alpha * t))**(2.0_dp / (GAMMA - 1.0_dp))
                u_ex(i) = u1 + (2.0_dp / (GAMMA + 1.0_dp)) * (x(i) - pos1) / t
            else if (x(i) <= conpos) then
                p_ex(i) = p3; rho_ex(i) = rho2; u_ex(i) = u2 + u4
            else if (x(i) <= spos) then
                p_ex(i) = p3; rho_ex(i) = rho3; u_ex(i) = u2 + u4
            else
                p_ex(i) = p4; rho_ex(i) = rho4; u_ex(i) = u4
            end if
        end do
    end subroutine get_exact_solution

end module exact_solver

module riemann_schemes
    use parameters
    implicit none

contains

    ! --- Minmod Limiter ---
    elemental function minmod(a, b) result(res)
        real(dp), intent(in) :: a, b
        real(dp) :: res
        if (a * b > 0.0_dp) then
            res = sign(1.0_dp, a) * min(abs(a), abs(b))
        else
            res = 0.0_dp
        end if
    end function minmod

    ! --- WENO5 Reconstruction (Fixed) ---
    ! Calculates the reconstructed state at the forward interface of cell v_0
    elemental function weno5(v_m2, v_m1, v_0, v_p1, v_p2) result(q)
        real(dp), intent(in) :: v_m2, v_m1, v_0, v_p1, v_p2
        real(dp) :: q
        real(dp) :: b0, b1, b2, a0, a1, a2, w0, w1, w2
        real(dp), parameter :: eps = 1.0e-6_dp

        b0 = 13.0_dp/12.0_dp * (v_m2 - 2.0_dp*v_m1 + v_0)**2 + 0.25_dp * (v_m2 - 4.0_dp*v_m1 + 3.0_dp*v_0)**2
        b1 = 13.0_dp/12.0_dp * (v_m1 - 2.0_dp*v_0 + v_p1)**2 + 0.25_dp * (v_m1 - v_p1)**2
        b2 = 13.0_dp/12.0_dp * (v_0 - 2.0_dp*v_p1 + v_p2)**2 + 0.25_dp * (3.0_dp*v_0 - 4.0_dp*v_p1 + v_p2)**2
        
        a0 = 0.1_dp / (eps + b0)**2
        a1 = 0.6_dp / (eps + b1)**2
        a2 = 0.3_dp / (eps + b2)**2
        
        w0 = a0 / (a0 + a1 + a2)
        w1 = a1 / (a0 + a1 + a2)
        w2 = a2 / (a0 + a1 + a2)
        
        q = w0 * (1.0_dp/3.0_dp*v_m2 - 7.0_dp/6.0_dp*v_m1 + 11.0_dp/6.0_dp*v_0) &
          + w1 * (-1.0_dp/6.0_dp*v_m1 + 5.0_dp/6.0_dp*v_0 + 1.0_dp/3.0_dp*v_p1) &
          + w2 * (1.0_dp/3.0_dp*v_0 + 5.0_dp/6.0_dp*v_p1 - 1.0_dp/6.0_dp*v_p2)
    end function weno5

    ! --- Physics & Flux ---
    subroutine get_euler_vars(q, f, u, c)
        real(dp), dimension(3), intent(in) :: q
        real(dp), dimension(3), intent(out) :: f
        real(dp), intent(out) :: u, c
        real(dp) :: rho, p, E

        rho = max(q(1), 1.0e-7_dp)
        u   = q(2) / rho
        E   = q(3)
        p   = max((GAMMA - 1.0_dp) * (E - 0.5_dp * rho * u**2), 1.0e-7_dp)
        c   = sqrt(GAMMA * p / rho)

        f(1) = rho * u
        f(2) = rho * u**2 + p
        f(3) = u * (E + p)
    end subroutine get_euler_vars

    ! --- Compute RHS ---
    subroutine compute_rhs(q, rhs, scheme)
        real(dp), dimension(3, NTOT), intent(in) :: q
        real(dp), dimension(3, NTOT), intent(out) :: rhs
        integer, intent(in) :: scheme
        
        real(dp), dimension(3, NTOT) :: F_int
        real(dp), dimension(3) :: qL, qR, fL, fR
        real(dp) :: uL, cL, uR, cR, alpha, slope_L, slope_R
        real(dp) :: dq_im1, dq_i, dq_ip1
        integer :: i, var

        rhs = 0.0_dp
        F_int = 0.0_dp

        do i = NG, NTOT - NG
            do var = 1, 3
                select case(scheme)
                    case(1) ! First Order
                        qL(var) = q(var, i)
                        qR(var) = q(var, i+1)
                    case(2) ! Minmod
                        dq_im1 = q(var, i) - q(var, i-1)
                        dq_i   = q(var, i+1) - q(var, i)
                        dq_ip1 = q(var, i+2) - q(var, i+1)
                        
                        slope_L = minmod(dq_im1, dq_i)
                        qL(var) = q(var, i) + 0.5_dp * slope_L
                        
                        slope_R = minmod(dq_i, dq_ip1)
                        qR(var) = q(var, i+1) - 0.5_dp * slope_R
                    case(3) ! WENO5 (Fixed Indices)
                        qL(var) = weno5(q(var, i-2), q(var, i-1), q(var, i), q(var, i+1), q(var, i+2))
                        qR(var) = weno5(q(var, i+3), q(var, i+2), q(var, i+1), q(var, i), q(var, i-1))
                end select
            end do

            ! Local Lax-Friedrichs Flux
            call get_euler_vars(qL, fL, uL, cL)
            call get_euler_vars(qR, fR, uR, cR)

            alpha = max(abs(uL) + cL, abs(uR) + cR)
            
            do var = 1, 3
                F_int(var, i) = 0.5_dp * (fL(var) + fR(var)) - 0.5_dp * alpha * (qR(var) - qL(var))
            end do
        end do

        ! Finite Difference
        do i = NG+1, NTOT - NG
            do var = 1, 3
                rhs(var, i) = -(F_int(var, i) - F_int(var, i-1)) / DX
            end do
        end do
    end subroutine compute_rhs

    ! --- Boundary Conditions ---
    subroutine apply_bcs(q)
        real(dp), dimension(3, NTOT), intent(inout) :: q
        integer :: i
        do i = 1, NG
            q(:, i) = q(:, NG+1)
            q(:, NTOT-i+1) = q(:, NTOT-NG)
        end do
    end subroutine apply_bcs

end module riemann_schemes
program euler1d
    use parameters
    use exact_solver
    use riemann_schemes
    implicit none

    real(dp), dimension(3, NTOT) :: Q, Q1, Q2, RHS
    real(dp), dimension(NX) :: x_phys, rho_ex, u_ex, p_ex, e_ex
    real(dp), dimension(NX, 3) :: res_rho, res_u, res_p, res_e
    real(dp), dimension(3) :: exec_times
    
    real(dp) :: start_time, end_time, t, dt, max_wave_speed
    real(dp) :: rhoL, uL, pL, rhoR, uR, pR
    real(dp) :: rho, u, p, E, c
    real(dp) :: err_rho, err_u, err_p, err_e
    real(dp) :: cost_rho, cost_u, cost_p, cost_e
    
    integer :: i, scheme
    character(len=25), dimension(3) :: scheme_names = (/"1st Order (LF)           ", &
                                                        "2nd Order (Minmod)       ", &
                                                        "WENO5                    "/)

    ! Sod's Shock Tube Initial Conditions
    rhoL = 1.0_dp;   uL = 0.0_dp; pL = 1.0_dp
    rhoR = 0.125_dp; uR = 0.0_dp; pR = 0.1_dp

    ! 1. Calculate Exact Solution
    do i = 1, NX
        x_phys(i) = (real(i, dp) - 0.5_dp) * DX
    end do
    call get_exact_solution(x_phys, T_END, rhoL, uL, pL, rhoR, uR, pR, rho_ex, u_ex, p_ex)
    e_ex = p_ex / ((GAMMA - 1.0_dp) * rho_ex)

    ! 2. Run Numerical Schemes
    do scheme = 1, 3
        ! Initialization
        do i = 1, NTOT
            if ((real(i - NG, dp) - 0.5_dp) * DX < 0.5_dp) then
                Q(1, i) = rhoL; Q(2, i) = rhoL * uL; Q(3, i) = pL / (GAMMA - 1.0_dp) + 0.5_dp * rhoL * uL**2
            else
                Q(1, i) = rhoR; Q(2, i) = rhoR * uR; Q(3, i) = pR / (GAMMA - 1.0_dp) + 0.5_dp * rhoR * uR**2
            end if
        end do
        call apply_bcs(Q)
        t = 0.0_dp

        ! Start Timer
        call cpu_time(start_time)

        do while (t < T_END)
            max_wave_speed = 0.0_dp
            do i = NG+1, NTOT-NG
                rho = max(Q(1, i), 1.0e-7_dp)
                u   = Q(2, i) / rho
                p   = max((GAMMA - 1.0_dp) * (Q(3, i) - 0.5_dp * rho * u**2), 1.0e-7_dp)
                c   = sqrt(GAMMA * p / rho)
                max_wave_speed = max(max_wave_speed, abs(u) + c)
            end do
            dt = CFL * DX / max_wave_speed
            if (t + dt > T_END) dt = T_END - t

            call compute_rhs(Q, RHS, scheme); Q1 = Q + dt * RHS; call apply_bcs(Q1)
            call compute_rhs(Q1, RHS, scheme); Q2 = 0.75_dp * Q + 0.25_dp * (Q1 + dt * RHS); call apply_bcs(Q2)
            call compute_rhs(Q2, RHS, scheme); Q = 1.0_dp/3.0_dp * Q + 2.0_dp/3.0_dp * (Q2 + dt * RHS); call apply_bcs(Q)
            t = t + dt
        end do

        ! Stop Timer
        call cpu_time(end_time)
        exec_times(scheme) = end_time - start_time

        ! Extract final physics arrays
        do i = 1, NX
            res_rho(i, scheme) = Q(1, i+NG)
            res_u(i, scheme)   = Q(2, i+NG) / res_rho(i, scheme)
            E                  = Q(3, i+NG)
            res_p(i, scheme)   = (GAMMA - 1.0_dp) * (E - 0.5_dp * res_rho(i, scheme) * res_u(i, scheme)**2)
            res_e(i, scheme)   = res_p(i, scheme) / ((GAMMA - 1.0_dp) * res_rho(i, scheme))
        end do
    end do

    ! 3. Print L1 Errors Output Table
    print '(/,A)', "================================================================="
    print '(A,I0,A,F5.2,A)', "  L1 Errors  (Sod's Shock Tube, N=", NX, ", t=", T_END, ")"
    print '(A)', "================================================================="
    print '(A)', "  Scheme                       Density   Velocity   Pressure  Int. Energy"
    print '(A)', "-----------------------------------------------------------------"
    do scheme = 1, 3
        err_rho = sum(abs(res_rho(:, scheme) - rho_ex)) * DX
        err_u   = sum(abs(res_u(:, scheme) - u_ex)) * DX
        err_p   = sum(abs(res_p(:, scheme) - p_ex)) * DX
        err_e   = sum(abs(res_e(:, scheme) - e_ex)) * DX
        print '(2X,A,4F11.6)', scheme_names(scheme), err_rho, err_u, err_p, err_e
    end do
    print '(A)', "================================================================="

    ! 4. Print Efficiency / Cost Output Table
    print '(/,A)', "============================================================================="
    print '(A,I0,A,F5.2,A)', "  Cost (Time / L1 Error)  (Sod's Shock Tube, N=", NX, ", t=", T_END, ")"
    print '(A)', "============================================================================="
    print '(A)', "  Scheme                       Density     Velocity     Pressure       Energy"
    print '(A)', "-----------------------------------------------------------------------------"
    do scheme = 1, 3
        err_rho = max(sum(abs(res_rho(:, scheme) - rho_ex)) * DX, 1.0e-12_dp)
        err_u   = max(sum(abs(res_u(:, scheme) - u_ex)) * DX, 1.0e-12_dp)
        err_p   = max(sum(abs(res_p(:, scheme) - p_ex)) * DX, 1.0e-12_dp)
        err_e   = max(sum(abs(res_e(:, scheme) - e_ex)) * DX, 1.0e-12_dp)
        
        cost_rho = exec_times(scheme) / err_rho
        cost_u   = exec_times(scheme) / err_u
        cost_p   = exec_times(scheme) / err_p
        cost_e   = exec_times(scheme) / err_e
        
        print '(2X,A,4F13.4)', scheme_names(scheme), cost_rho, cost_u, cost_p, cost_e
    end do
    print '(A)', "============================================================================="

end program euler1d