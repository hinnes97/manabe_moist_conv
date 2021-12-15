! Hamish Innes 12/2021, University of Oxford

module adjust_mod
  
  use phys, only : T_TP => H2O_TriplePointT, P_TP => H2O_TriplePointP, L_sub => H2O_L_sublimation, &
       L_vap => H2O_L_vaporization_TriplePoint, CP_v => H2O_cp, CP_d => N2_cp, &
       mu_d => N2_MolecularWeight, mu_v => H2O_MolecularWeight, Rstar

  implicit none

  ! FOR FSOLVE PARAMS:
  integer, parameter :: rk = kind(1.0D+00)
  real(kind=rk) :: T1, T2, p1, p2, dp1, dp2, q1, q2, pfactor
  
contains

  subroutine adjust(p, dp, T, q, mask)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Performs moist adiabatic adjustment, conserving column-integrated moist
    ! enthalpy: \int(cp*T + Lq)*dp in the dilute limit, as in Manabe 1965. The
    ! adjustment is performed pairwise in layers, until convergence is reached.
    ! For two layers (labelled 1 and 2, with post-adjustment state having
    ! dashes), the conservation of moist enthalpy is:
    !
    ! cp*(T1*dp1+T2*dp2) + L*(q1(T1)*dp1+q2(T2)*dp2)
    !                    = cp*(T1'*dp1+T2'*dp2) + L*(q1'(T1')*dp1+q2'(T2')*dp2)
    !
    ! We can relate T1' to T2' by T1' = T2'*(p1/p2)**(dlnT/dlnp), and then since
    ! qi' is related to Ti' by the Clausius Clapyeron relation, we can solve for
    ! T2' using non-linear equation solver fsolve (More 1980).
    !
    ! Algorithm heavily based on Ray Pierrehumbert's dry convection code

    !==========================================================================
    ! Input variables
    !==========================================================================    
    real, intent(in)   , dimension(:) :: p,dp ! Pressure, pressure thickness

    !==========================================================================
    ! Output variables
    !==========================================================================
    logical, intent(out),dimension(:) :: mask ! True where adiabatic

    !==========================================================================
    ! Mixed input/output
    !==========================================================================
    real, intent(inout), dimension(:) :: T,q ! Temperature and specific humid.

    !==========================================================================
    ! Local variables
    !==========================================================================
    ! Tune these parameters as necessasry
    real, parameter          :: delta = 0.000001 ! Small number speeds up convergence
    real(kind=rk), parameter :: tol = 0.00001 ! Tolerance for non-linear solver
    integer, parameter       :: N_iter = 50   ! Number of up-down iterations
    
    real :: qsat1, qsat2, pfact, grad
    real(kind=rk) :: output(1), f_output(1)

    integer :: n,k
    integer :: info
    integer :: npz 

    !==========================================================================
    ! Main body
    !==========================================================================
    
    npz = size(p)
    info = 0
    do n=1,N_iter

       ! Downwards pass
       do k=1,npz-1
          call sat(p(k), T(k), qsat1)
          call sat(p(k+1), T(k+1), qsat2)

          if ( (q(k) .gt. qsat1*(1.-delta)) .and. (q(k+1) .gt. qsat2*(1.-delta)) ) then

             ! Equivalent to doing large-scale condensation without latent heating, remove
             ! if performed elsewhere

             q(k)   = qsat1
             q(k+1) = qsat2

             call gradient(p(k+1), T(k+1), grad)

             pfact = exp(grad*log(p(k)/p(k+1)))

             ! Test for instability
             if (T(k) .lt. T(k+1)*pfact*(1. + delta) ) then
                ! INSERT ROOT FINDER HERE FOR T(k+1)
                ! HACK: set module variables so that fsolve function can have parameters

                T1 = real(T(k), rk)
                T2 = real(T(k+1), rk)
                p1 = real(p(k), rk)
                p2 = real(p(k+1), rk)
                dp1 = real(dp(k), rk)
                dp2 = real(dp(k+1), rk)
                q1 = real(q(k), rk)
                q2 = real(q(k+1), rk)
                pfactor = real(pfact, rk)

                !Initial guess
                output(1) = 250.0D+00
                call find_my_root(1, output, f_output)

                ! Use fsolve to find root of non-linear equation
                call fsolve(find_my_root, 1, output, f_output, tol, info)

                if (info .ne. 1) write(*,*) 'ERROR IN FSOLVE, CODE: ', info, 'level: ', k
                
                T(k+1) = real(output(1), kind(1.0))
                T(k)   = T(k+1)*pfact

                call sat(p(k), T(k), q(k))
                call sat(p(k+1),T(k+1), q(k+1))

                mask(k) = .true.
                mask(k+1) = .true.
             else
                mask(k) = .false.
                mask(k+1) = .false.
             endif

          endif
       enddo

       ! Upwards pass 
       do k=npz-1,1,-1
          call sat(p(k), T(k), qsat1)
          call sat(p(k+1), T(k+1), qsat2)

          if ( (q(k) .gt. qsat1*(1.-delta)) .and. (q(k+1) .gt. qsat2*(1.-delta)) ) then

             ! Equivalent to doing large-scale condensation without latent heating, remove
             ! if performed elsewhere

             q(k)   = qsat1
             q(k+1) = qsat2

             call gradient(p(k+1), T(k+1), grad)

             pfact =exp(grad*log(p(k)/p(k+1)))

             if (T(k) .lt. T(k+1)*pfact*(1. + delta) ) then
                ! INSERT ROOT FINDER HERE FOR T(k+1)
                ! HACK: set module variables so that fsolve function can have parameters
                T1 = real(T(k), rk)
                T2 = real(T(k+1), rk)
                p1 = real(p(k), rk)
                p2 = real(p(k+1), rk)
                dp1 = real(dp(k), rk)
                dp2 = real(dp(k+1), rk)
                q1 = real(q(k), rk)
                q2 = real(q(k+1), rk)
                pfactor = real(pfact, rk)

                !Initial guess
                output(1) = 250.0D+00
                call find_my_root(1, output, f_output)

                call fsolve(find_my_root, 1, output, f_output, tol, info)

                if (info .ne. 1) write(*,*) 'ERROR IN FSOLVE, CODE: ', info
                
                T(k+1) = real(output(1), kind(1.0))
                T(k) = T(k+1)*pfact

                call sat(p(k), T(k), q(k))
                call sat(p(k+1),T(k+1), q(k+1))

                mask(k) = .true.
                mask(k+1) = .true.
             else
                mask(k) = .false.
                mask(k+1) = .false.
             endif
          endif
       enddo
    enddo
    
    
  end subroutine adjust

  subroutine gradient(p,T, dlnTdlnp)
    !==========================================================================
    ! Description
    !==========================================================================
    ! Calculates the moist adiabatic gradient d(log(T))/d(log(p)) as in
    ! Ding and Pierrehumbert 2016
    
    !==========================================================================
    ! Input variables
    !==========================================================================
    real, intent(in) :: p, T !Pressure and temperature

    !==========================================================================
    ! Output variables
    !==========================================================================
    real, intent(out) :: dlnTdlnp ! Moist adiabatic gradient

    !==========================================================================
    ! Local variables
    !==========================================================================
    real :: eps = mu_v/mu_d
    real :: L, psat, qsat, rsat, num, denom, temp

    !==========================================================================
    ! Main body
    !==========================================================================
    
    if (T > T_TP) then
       L = L_vap
    else
       L = L_sub
    endif

    call sat(p, T, qsat, rsat, psat)

    num   = 1 + (L*mu_d/Rstar/T)*rsat
    denom = 1 + ((cp_v/cp_d) + ((L*mu_v/Rstar/T) - 1)*(L/cp_d/T) )*rsat

    temp = Rstar/mu_d/cp_d * num/denom
    
    dlnTdlnp = 1. / ((psat/p)*L*mu_v/Rstar/T + (p - psat)/p/temp)

  end subroutine gradient
  
  subroutine sat(p,T, qsat, rsat, psat)
    !==========================================================================
    ! Description
    !==========================================================================

    ! Calculates the saturation vapour pressure, mass mixing ratio and mass
    ! concentration of water vapour according to simple Clausius Clapeyron
    ! relation. Change this function to some parametrisation if more accuracy
    ! required.

    !==========================================================================
    ! Input variables
    !==========================================================================
    
    real, intent(in)  :: p, T ! Pressure, temperature

    !==========================================================================
    ! Output variables
    !==========================================================================
    
    real, intent(out) :: qsat ! Mass concentration
    real, intent(out), optional :: rsat, psat ! Mass mixing ratio, sat pressure

    !==========================================================================
    ! Local variables
    !==========================================================================
    real    :: L, eps, psat_temp, rsat_temp

    !==========================================================================
    ! Main body
    !==========================================================================
    
    eps = mu_v/mu_d

    if (T .gt. T_TP) then
       L = L_vap
    else
       L = L_sub
    endif


    psat_temp = P_TP * exp(-L/Rstar*mu_v * (1./T - 1./T_TP))
    rsat_temp = psat_temp / (p-psat_temp) * eps
    qsat = rsat_temp / (1 + rsat_temp)

    if(present(psat)) psat = psat_temp
    if(present(rsat)) rsat = rsat_temp

  end subroutine sat


  subroutine find_my_root(n, T, fvec)

    !==========================================================================
    ! Description
    !==========================================================================
    ! Function passed to fsolve, which will find its roots in T. Represents the
    ! conservation of moist enthalpy in the dilute limit

    !==========================================================================
    ! Input variables
    !==========================================================================
    integer, intent(in) :: n  ! Number of non-linear equations to be solved (1)
    real(kind=rk),    intent(in) :: T(n) ! Temperature to be solved for

    !==========================================================================
    ! Output variables
    !==========================================================================
    real(kind=rk),    intent(out) :: fvec(n) ! Value of function at T

    !==========================================================================
    ! Local variables
    !==========================================================================
    real(kind=rk) :: q1_new, q2_new, L
    real :: q1_temp, q2_temp

    !==========================================================================
    ! Main body
    !==========================================================================
    
    call sat(real(p1,4), real(T(1)*pfactor,4), q1_temp)
    call sat(real(p2,4), real(T(1),4), q2_temp)

    q1_new = real(q1_temp, rk)
    q2_new = real(q2_temp, rk)
    
    if (T2 .gt. real(T_TP, rk)) then
       L = real(L_vap, rk)
    else
       L = real(L_sub, rk)
    endif

    ! Moist enthalpy equation divided by cp*T, to be roughly order unity for solver
    fvec(1) = 1 - (T1*dp1 + T2*dp2 + L/cp_d*(q1 - q1_new)*dp1 + L/cp_d*(q2 - q2_new)*dp2)&
         / (dp2 + dp1*pfactor) / T(1)
    
  end subroutine find_my_root
end module adjust_mod
