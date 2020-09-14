subroutine regcoil_evaluate_coil_surface()

  ! This subroutine takes the arrays rmnc_coil, zmns_coil, etc, and evaluates the position vector r_coil
  ! and its first and second derivatives with respect to theta and zeta.

  use regcoil_Gaussian_quadrature
  use regcoil_variables
  use stel_kinds
  
  implicit none
  
  integer :: imn, m, n, iflag, j, js
  real(dp) :: rmnc, rmns, zmnc, zmns
  real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dcosangledtheta, dsinangledzeta, dcosangledzeta
  real(dp) :: d2sinangledtheta2, d2sinangledthetadzeta, d2sinangledzeta2, d2cosangledtheta2, d2cosangledthetadzeta, d2cosangledzeta2
  real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta, d2sinangle2dzeta2, d2cosangle2dzeta2
  real(dp), dimension(:,:), allocatable :: major_R_squared
  integer :: itheta, izeta
  integer :: tic, toc, countrate, clock
  real(dp), dimension(:,:), allocatable :: sin_m_theta, cos_m_theta, sin_n_zeta, cos_n_zeta
  real(dp), dimension(:,:), allocatable :: fundamental_form_E, fundamental_form_F, fundamental_form_G
  real(dp), dimension(:,:), allocatable :: fundamental_form_L, fundamental_form_M, fundamental_form_P
  real(dp), dimension(:,:), allocatable :: temp_matrix
  real(dp) :: a, temp, d1, d2
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  call system_clock(tic,countrate)

  allocate(sin_m_theta(mnmax_coil,ntheta_coil))
  allocate(cos_m_theta(mnmax_coil,ntheta_coil))
  allocate(sin_n_zeta( mnmax_coil,nzetal_coil))
  allocate(cos_n_zeta( mnmax_coil,nzetal_coil))
  do itheta = 1,ntheta_coil
     sin_m_theta(:,itheta) = sin(xm_coil * theta_coil(itheta))
     cos_m_theta(:,itheta) = cos(xm_coil * theta_coil(itheta))
  end do
  do izeta = 1,nzetal_coil
     sin_n_zeta(:,  izeta) = sin(xn_coil * zetal_coil(izeta))
     cos_n_zeta(:,  izeta) = cos(xn_coil * zetal_coil(izeta))
  end do

  r_coil = 0
  drdtheta_coil = 0
  drdzeta_coil = 0
  d2rdtheta2_coil = 0
  d2rdthetadzeta_coil = 0
  d2rdzeta2_coil = 0
  
  do izeta = 1,nzetal_coil
     angle2 = zetal_coil(izeta)
     sinangle2 = sin(angle2)
     cosangle2 = cos(angle2)
     dsinangle2dzeta = cosangle2
     dcosangle2dzeta = -sinangle2
     d2sinangle2dzeta2 = -sinangle2
     d2cosangle2dzeta2 = -cosangle2
     do imn = 1, mnmax_coil
        m = xm_coil(imn)
        n = xn_coil(imn)
        rmnc = rmnc_coil(imn)
        rmns = rmns_coil(imn)
        zmnc = zmnc_coil(imn)
        zmns = zmns_coil(imn)
        do itheta = 1,ntheta_coil
           !angle = m*theta_coil(itheta) - n*zetal_coil(izeta)
           !sinangle = sin(angle)
           !cosangle = cos(angle)
           ! Trig angle sum formulae for angle = m*theta_coil(itheta) - n*zetal_coil(izeta):
           sinangle = sin_m_theta(imn,itheta) * cos_n_zeta(imn,izeta) - cos_m_theta(imn,itheta) * sin_n_zeta(imn,izeta)
           cosangle = cos_m_theta(imn,itheta) * cos_n_zeta(imn,izeta) + sin_m_theta(imn,itheta) * sin_n_zeta(imn,izeta)
           !if (abs(sinangle - sin(angle)) > 1d-10) stop "Error sin"
           !if (abs(cosangle - cos(angle)) > 1d-10) stop "Error cos"
           dsinangledtheta = cosangle*m
           dcosangledtheta = -sinangle*m
           dsinangledzeta = -cosangle*n
           dcosangledzeta = sinangle*n
           
           r_coil(1,itheta,izeta) = r_coil(1,itheta,izeta) + rmnc * cosangle * cosangle2 + rmns * sinangle * cosangle2
           r_coil(2,itheta,izeta) = r_coil(2,itheta,izeta) + rmnc * cosangle * sinangle2 + rmns * sinangle * sinangle2
           r_coil(3,itheta,izeta) = r_coil(3,itheta,izeta) + zmns * sinangle             + zmnc * cosangle
           
           drdtheta_coil(1,itheta,izeta) = drdtheta_coil(1,itheta,izeta) + rmnc * dcosangledtheta * cosangle2 + rmns * dsinangledtheta * cosangle2
           drdtheta_coil(2,itheta,izeta) = drdtheta_coil(2,itheta,izeta) + rmnc * dcosangledtheta * sinangle2 + rmns * dsinangledtheta * sinangle2
           drdtheta_coil(3,itheta,izeta) = drdtheta_coil(3,itheta,izeta) + zmns * dsinangledtheta + zmnc * dcosangledtheta
           
           drdzeta_coil(1,itheta,izeta) = drdzeta_coil(1,itheta,izeta) + rmnc * (dcosangledzeta * cosangle2 + cosangle * dcosangle2dzeta) &
                + rmns * (dsinangledzeta * cosangle2 + sinangle * dcosangle2dzeta)
           drdzeta_coil(2,itheta,izeta) = drdzeta_coil(2,itheta,izeta) + rmnc * (dcosangledzeta * sinangle2 + cosangle * dsinangle2dzeta) &
                + rmns * (dsinangledzeta * sinangle2 + sinangle * dsinangle2dzeta)
           drdzeta_coil(3,itheta,izeta) = drdzeta_coil(3,itheta,izeta) + zmns * dsinangledzeta + zmnc * dcosangledzeta

           ! 2nd derivatives are only constructed on 1 field period.
           if (izeta > nzeta_coil) cycle

           d2sinangledtheta2  = -m*m*sinangle
           d2sinangledthetadzeta = m*n*sinangle
           d2sinangledzeta2  = -n*n*sinangle
           d2cosangledtheta2  = -m*m*cosangle
           d2cosangledthetadzeta = m*n*cosangle
           d2cosangledzeta2  = -n*n*cosangle

           d2rdtheta2_coil(1,itheta,izeta) = d2rdtheta2_coil(1,itheta,izeta) + rmnc * d2cosangledtheta2 * cosangle2 + rmns * d2sinangledtheta2 * cosangle2
           d2rdtheta2_coil(2,itheta,izeta) = d2rdtheta2_coil(2,itheta,izeta) + rmnc * d2cosangledtheta2 * sinangle2 + rmns * d2sinangledtheta2 * sinangle2
           d2rdtheta2_coil(3,itheta,izeta) = d2rdtheta2_coil(3,itheta,izeta) + zmns * d2sinangledtheta2 + zmnc * d2cosangledtheta2
           
           d2rdthetadzeta_coil(1,itheta,izeta) = d2rdthetadzeta_coil(1,itheta,izeta) + rmnc * (d2cosangledthetadzeta * cosangle2 + dcosangledtheta * dcosangle2dzeta) &
                + rmns * (d2sinangledthetadzeta * cosangle2 + dsinangledtheta * dcosangle2dzeta)
           d2rdthetadzeta_coil(2,itheta,izeta) = d2rdthetadzeta_coil(2,itheta,izeta) + rmnc * (d2cosangledthetadzeta * sinangle2 + dcosangledtheta * dsinangle2dzeta) &
                + rmns * (d2sinangledthetadzeta * sinangle2 + dsinangledtheta * dsinangle2dzeta)
           d2rdthetadzeta_coil(3,itheta,izeta) = d2rdthetadzeta_coil(3,itheta,izeta) + zmns * d2sinangledthetadzeta + zmnc * d2cosangledthetadzeta
           
           d2rdzeta2_coil(1,itheta,izeta) = d2rdzeta2_coil(1,itheta,izeta) + rmnc * (d2cosangledzeta2 * cosangle2 + dcosangledzeta * dcosangle2dzeta &
                + dcosangledzeta * dcosangle2dzeta + cosangle * d2cosangle2dzeta2) &
                + rmns * (d2sinangledzeta2 * cosangle2 + dsinangledzeta * dcosangle2dzeta &
                + dsinangledzeta * dcosangle2dzeta + sinangle * d2cosangle2dzeta2)
           d2rdzeta2_coil(2,itheta,izeta) = d2rdzeta2_coil(2,itheta,izeta) + rmnc * (d2cosangledzeta2 * sinangle2 + dcosangledzeta * dsinangle2dzeta &
                + dcosangledzeta * dsinangle2dzeta + cosangle * d2sinangle2dzeta2) &
                + rmns * (d2sinangledzeta2 * sinangle2 + dsinangledzeta * dsinangle2dzeta &
                + dsinangledzeta * dsinangle2dzeta + sinangle * d2sinangle2dzeta2)
           d2rdzeta2_coil(3,itheta,izeta) = d2rdzeta2_coil(3,itheta,izeta) + zmns * d2sinangledzeta2 + zmnc * d2cosangledzeta2
        end do
     end do
  end do

  deallocate(sin_m_theta, cos_m_theta, sin_n_zeta, cos_n_zeta)

  call system_clock(toc)
  if (verbose) print *,"  Evaluating coil surface & derivatives:",real(toc-tic)/countrate," sec."

!!$  print *,"mnmax_coil:",mnmax_coil
!!$  print *,"xm_coil:",xm_coil
!!$  print *,"xn_coil:",xn_coil
!!$  print *,"rmnc_coil:",rmnc_coil
!!$  print *,"zmns_coil:",zmns_coil
!!$  print *,"r_coil(1,:,:):"
!!$  do j = 1,ntheta_coil
!!$     print *,r_coil(1,j,:)
!!$  end do
!!$  print *,"r_coil(2,:,:):"
!!$  do j = 1,ntheta_coil
!!$     print *,r_coil(2,j,:)
!!$  end do
!!$  print *,"r_coil(3,:,:):"
!!$  do j = 1,ntheta_coil
!!$     print *,r_coil(3,j,:)
!!$  end do

!!$  print *,"d2xdtheta2:"
!!$  do j = 1,ntheta_coil
!!$     print *,d2rdtheta2_coil(1,j,:)
!!$  end do
!!$  print *,"d2xdthetadzeta:"
!!$  do j = 1,ntheta_coil
!!$     print *,d2rdthetadzeta_coil(1,j,:)
!!$  end do
!!$  print *,"d2xdzeta2:"
!!$  do j = 1,ntheta_coil
!!$     print *,d2rdzeta2_coil(1,j,:)
!!$  end do
!!$
!!$  print *,"d2ydtheta2:"
!!$  do j = 1,ntheta_coil
!!$     print *,d2rdtheta2_coil(2,j,:)
!!$  end do
!!$  print *,"d2ydthetadzeta:"
!!$  do j = 1,ntheta_coil
!!$     print *,d2rdthetadzeta_coil(2,j,:)
!!$  end do
!!$  print *,"d2ydzeta2:"
!!$  do j = 1,ntheta_coil
!!$     print *,d2rdzeta2_coil(2,j,:)
!!$  end do
!!$
!!$  print *,"d2zdtheta2:"
!!$  do j = 1,ntheta_coil
!!$     print *,d2rdtheta2_coil(3,j,:)
!!$  end do
!!$  print *,"d2zdthetadzeta:"
!!$  do j = 1,ntheta_coil
!!$     print *,d2rdthetadzeta_coil(3,j,:)
!!$  end do
!!$  print *,"d2zdzeta2:"
!!$  do j = 1,ntheta_coil
!!$     print *,d2rdzeta2_coil(3,j,:)
!!$  end do

  ! Evaluate cross product:
  normal_coil(1,:,:) = drdzeta_coil(2,:,:) * drdtheta_coil(3,:,:) - drdtheta_coil(2,:,:) * drdzeta_coil(3,:,:)
  normal_coil(2,:,:) = drdzeta_coil(3,:,:) * drdtheta_coil(1,:,:) - drdtheta_coil(3,:,:) * drdzeta_coil(1,:,:)
  normal_coil(3,:,:) = drdzeta_coil(1,:,:) * drdtheta_coil(2,:,:) - drdtheta_coil(1,:,:) * drdzeta_coil(2,:,:)
  
  if (allocated(norm_normal_coil)) deallocate(norm_normal_coil)
  allocate(norm_normal_coil(ntheta_coil, nzeta_coil),stat=iflag)
  if (iflag .ne. 0) stop 'Allocation error! regcoil_init_coil_surface 11'
  norm_normal_coil = sqrt(normal_coil(1,:,1:nzeta_coil)**2 + normal_coil(2,:,1:nzeta_coil)**2 &
       +  normal_coil(3,:,1:nzeta_coil)**2)
  
  area_coil = nfp * dtheta_coil * dzeta_coil * sum(norm_normal_coil)

  ! Evaluate coefficients of the first and second fundamental forms:
  allocate(fundamental_form_E(ntheta_coil, nzeta_coil))
  allocate(fundamental_form_F(ntheta_coil, nzeta_coil))
  allocate(fundamental_form_G(ntheta_coil, nzeta_coil))
  allocate(fundamental_form_L(ntheta_coil, nzeta_coil))
  allocate(fundamental_form_M(ntheta_coil, nzeta_coil))
  allocate(fundamental_form_P(ntheta_coil, nzeta_coil))
  fundamental_form_E = drdtheta_coil(1,:,1:nzeta_coil) * drdtheta_coil(1,:,1:nzeta_coil) + drdtheta_coil(2,:,1:nzeta_coil) * drdtheta_coil(2,:,1:nzeta_coil) + drdtheta_coil(3,:,1:nzeta_coil) * drdtheta_coil(3,:,1:nzeta_coil)
  fundamental_form_F = drdtheta_coil(1,:,1:nzeta_coil) *  drdzeta_coil(1,:,1:nzeta_coil) + drdtheta_coil(2,:,1:nzeta_coil) *  drdzeta_coil(2,:,1:nzeta_coil) + drdtheta_coil(3,:,1:nzeta_coil) *  drdzeta_coil(3,:,1:nzeta_coil)
  fundamental_form_G =  drdzeta_coil(1,:,1:nzeta_coil) *  drdzeta_coil(1,:,1:nzeta_coil) +  drdzeta_coil(2,:,1:nzeta_coil) *  drdzeta_coil(2,:,1:nzeta_coil) +  drdzeta_coil(3,:,1:nzeta_coil) *  drdzeta_coil(3,:,1:nzeta_coil)
  fundamental_form_L = (normal_coil(1,:,1:nzeta_coil) *     d2rdtheta2_coil(1,:,:) + normal_coil(2,:,1:nzeta_coil) *     d2rdtheta2_coil(2,:,:) + normal_coil(3,:,1:nzeta_coil) *     d2rdtheta2_coil(3,:,:)) / norm_normal_coil
  fundamental_form_M = (normal_coil(1,:,1:nzeta_coil) * d2rdthetadzeta_coil(1,:,:) + normal_coil(2,:,1:nzeta_coil) * d2rdthetadzeta_coil(2,:,:) + normal_coil(3,:,1:nzeta_coil) * d2rdthetadzeta_coil(3,:,:)) / norm_normal_coil
  fundamental_form_P = (normal_coil(1,:,1:nzeta_coil) *      d2rdzeta2_coil(1,:,:) + normal_coil(2,:,1:nzeta_coil) *      d2rdzeta2_coil(2,:,:) + normal_coil(3,:,1:nzeta_coil) *      d2rdzeta2_coil(3,:,:)) / norm_normal_coil

  mean_curvature_coil = (fundamental_form_L * fundamental_form_G + fundamental_form_P * fundamental_form_E - 2 * fundamental_form_M * fundamental_form_F) &
       / (2 * (fundamental_form_E * fundamental_form_G - fundamental_form_F * fundamental_form_F))

  ! Initialize s arrays
  allocate(s_integration(ns_integration))
  allocate(s_weights(ns_integration))
  allocate(s_magnetization(ns_magnetization))
  allocate(temp_matrix(ns_integration,ns_integration))
  select case (trim(s_integration_option))
  case (s_integration_option_Chebyshev)
     call regcoil_Chebyshev_grid(ns_integration,   0.0_dp, 1.0_dp, s_integration, s_weights, temp_matrix)
  case (s_integration_option_Gaussian)
     call get_legendre_grids(0.0_dp, 1.0_dp, s_integration, s_weights)
  case default
     stop "Unrecognized s_integration_option"
  end select
  deallocate(temp_matrix)
  allocate(temp_matrix(ns_magnetization,ns_magnetization))
  allocate(s_magnetization_weights(ns_magnetization))
  call regcoil_Chebyshev_grid(ns_magnetization, 0.0_dp, 1.0_dp, s_magnetization, s_magnetization_weights, temp_matrix)
  deallocate(temp_matrix)
  allocate(interpolate_magnetization_to_integration(ns_integration,ns_magnetization))
  call regcoil_Chebyshev_interpolation_matrix(ns_magnetization, ns_integration, s_magnetization, s_integration, interpolate_magnetization_to_integration)
  print *,"s_integration:",s_integration
  print *,"s_weights:",s_weights
  print *,"s_magnetization",s_magnetization
  print *,"s_magnetization_weights",s_magnetization_weights
  print *,"interpolate_magnetization_to_integration:"
  do j = 1,ns_integration
     print *,interpolate_magnetization_to_integration(j,:)
  end do

  ! Initialize d, the magnet thickness
  allocate(d(ntheta_coil, nzeta_coil))
  d = d_initial
  d0 = d_initial
  if (random_d_initial) then
     ! Set random seed based on system time
     ! https://gcc.gnu.org/onlinedocs/gcc-4.6.4/gfortran/RANDOM_005fSEED.html

     CALL RANDOM_SEED(size = n)
     ALLOCATE(seed(n))
          
     CALL SYSTEM_CLOCK(COUNT=clock)
          
     seed = clock + 37 * (/ (j - 1, j = 1, n) /)
     CALL RANDOM_SEED(PUT = seed)
     print *,"Setting random seed:", seed
          
     DEALLOCATE(seed)

     ! Each point of d is random in the interval [0, d_initial]:
!!$     do izeta = 1, nzeta_coil
!!$        do itheta = 1, ntheta_coil
!!$           call random_number(temp)
!!$           d(itheta, izeta) = temp * d_initial
!!$        end do
!!$     end do

     d = 0
     ! Add some low-order Fourier modes with random amplitudes
     do m = 0, 4
        do n = -4, 4
           call random_number(temp)
           do izeta = 1, nzeta_coil
              do itheta = 1, ntheta_coil
                 d(itheta, izeta) = d(itheta, izeta) + (temp - 0.5) * cos(m * theta_coil(itheta) - n * nfp * zeta_coil(izeta))
              end do
           end do
        end do
     end do
     ! Shift and scale d so it lies between [delta*d_initial, d_initial]:
     d = d - minval(d)
     d = (d / maxval(d) * 0.99 + 0.01) * d_initial
  end if


  allocate(Jacobian_ssquared_term(ntheta_coil, nzeta_coil))
  ! We wll need this quantity later to generate the Jacobian of the (s, theta, zeta) coordinates in the magnetization region:
  Jacobian_ssquared_term = (fundamental_form_M * fundamental_form_M - fundamental_form_L * fundamental_form_P) / norm_normal_coil

  allocate(max_d_before_singularity(ntheta_coil, nzeta_coil))
  do izeta = 1, nzeta_coil
     do itheta = 1, ntheta_coil
        ! Use quadratic formula to find when Jacobian=0 for s=1:
        a = (fundamental_form_M(itheta,izeta) * fundamental_form_M(itheta,izeta) - fundamental_form_L(itheta,izeta) * fundamental_form_P(itheta,izeta)) &
             / (norm_normal_coil(itheta,izeta) * norm_normal_coil(itheta,izeta))
        if (abs(a) < 1e-12) then
           ! Handle the a=0 case separately
           d1 = 1 / (2 * sign_normal * mean_curvature_coil(itheta,izeta))
           if (d1 > 0) then
              max_d_before_singularity(itheta,izeta) = d1
           else
              max_d_before_singularity(itheta,izeta) = 1.0d+200
           end if
           !print *,"a=0:",a,"  mean_curv:",mean_curvature_coil(itheta,izeta),"  temp:",temp,"  d1:",d1,"  d2:",d2,"  max_d:",max_d_before_singularity(itheta,izeta)
           cycle
        end if
        temp = sqrt(4 * mean_curvature_coil(itheta,izeta) * mean_curvature_coil(itheta,izeta) + 4 * a)
        d1 = (-2 * sign_normal * mean_curvature_coil(itheta,izeta) + temp) / (2 * a)
        d2 = (-2 * sign_normal * mean_curvature_coil(itheta,izeta) - temp) / (2 * a)
        if (d1 < 0) then
           if (d2 >= 0) then
              max_d_before_singularity(itheta,izeta) = d2
           else
              max_d_before_singularity(itheta,izeta) = 1.0d+200
           end if
        else
           if ((d2 < d1) .and. (d2 >= 0)) then
              max_d_before_singularity(itheta,izeta) = d2
           else
              max_d_before_singularity(itheta,izeta) = d1
           end if
        end if
        !print "(3(a,es10.2))","d1:",d1,"  d2:",d2,"  max_d:",max_d_before_singularity(itheta,izeta)
        !print *,"a:",a,"  mean_curv:",mean_curvature_coil(itheta,izeta),"  temp:",temp,"  d1:",d1,"  d2:",d2,"  max_d:",max_d_before_singularity(itheta,izeta)
     end do
  end do

  print "(a,es12.4)"," Maximum allowed uniform d before singularity:",minval(max_d_before_singularity)
  if (d_initial >= minval(max_d_before_singularity)) stop "Error! d_initial is too large. Jacobian is singular."

  deallocate(fundamental_form_E, fundamental_form_F, fundamental_form_G, fundamental_form_L, fundamental_form_M, fundamental_form_P)
  
  ! Compute coil surface volume using \int (1/2) R^2 dZ dzeta.
  ! These quantities will be evaluated on the half theta grid, which is the natural grid for dZ,
  ! but we will need to interpolate R^2 from the full to half grid.
  allocate(major_R_squared(ntheta_coil,nzetal_coil))
  major_R_squared = r_coil(1,:,:)*r_coil(1,:,:) + r_coil(2,:,:)*r_coil(2,:,:)
  ! First handle the interior of the theta grid:
  volume_coil = sum((major_R_squared(1:ntheta_coil-1,:) + major_R_squared(2:ntheta_coil,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
       * (r_coil(3,2:ntheta_coil,:)-r_coil(3,1:ntheta_coil-1,:))) ! dZ
  ! Add the contribution from the ends of the theta grid:
  volume_coil = volume_coil + sum((major_R_squared(1,:) + major_R_squared(ntheta_coil,:)) * (0.5d+0) & ! R^2, interpolated from full to half grid
       * (r_coil(3,1,:)-r_coil(3,ntheta_coil,:))) ! dZ
  volume_coil = abs(volume_coil * dzeta_coil / 2) ! r includes all nfp periods already, so no factor of nfp needed.
  deallocate(major_R_squared)
  if (verbose) print "(a,es10.3,a,es10.3,a)"," Coil surface area:",area_coil," m^2. Volume:",volume_coil," m^3."
  
end subroutine  regcoil_evaluate_coil_surface
