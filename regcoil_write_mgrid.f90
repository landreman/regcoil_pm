subroutine regcoil_write_mgrid()

  use regcoil_variables
  use stel_constants
  use stel_kinds
  use omp_lib
  use ezcdf
  USE mgrid_mod, ONLY: vn_nextcur, vn_mgmode, vn_ir, vn_jz, vn_kp, vn_nfp, vn_rmin, vn_rmax, vn_zmin, vn_zmax, vn_coilgrp, vn_coilcur, vn_br0, vn_bz0, vn_bp0
  
  implicit none

  integer :: nextcur = 1


  CHARACTER(LEN=*), PARAMETER :: coildim(2) = (/'stringsize          ', 'external_coil_groups'/),  groupdim(1)= (/'external_coils'/), cylcoord(3)= (/'rad','zee','phi'/)

  INTEGER     :: iextc, istat
  INTEGER     :: kp2, jz2, kp_odd, jz_odd
  REAL(dp), ALLOCATABLE, DIMENSION(:, :, :) :: br, bz, bp
  REAL(dp), ALLOCATABLE, DIMENSION(:) :: extcur
  REAL(dp) :: rmin, rmax, zmin, zmax
  REAL(dp) :: fperiod, delr, delz, delp
  LOGICAL     :: lstell_sym=.true.
  CHARACTER(LEN=1) :: mgrid_mode='S'
  CHARACTER(LEN=70) :: mgrid_file, coil_file
  CHARACTER(LEN=60) :: mgrid_ext

  INTEGER :: ngrid, ig
  CHARACTER(LEN=100) :: temp
  CHARACTER(LEN=100), ALLOCATABLE, DIMENSION(:) ::  vn_br, vn_bp, vn_bz
  character(len=*), parameter :: coil_group_name = 'TF_and_PM'
  REAL(dp) :: rad, zee, phi
  INTEGER :: k, j, i

  real(dp) :: constants
  integer :: tic, toc, countrate
  integer :: ir, jz, kp
  integer :: l_coil, izeta_coil, izetal_coil, itheta_coil

  ir = mgrid_ir
  jz = mgrid_jz
  kp = mgrid_kp
  rmin = mgrid_rmin
  rmax = mgrid_rmax
  zmin = mgrid_zmin
  zmax = mgrid_zmax
  ig = 1

  call system_clock(tic,countrate)
  if (verbose) print *,"Beginning to write mgrid file."

  allocate(cos_zetal(nzetal_coil))
  allocate(sin_zetal(nzetal_coil))
  do izetal_coil = 1, nzetal_coil
     cos_zetal(izetal_coil) = cos(zetal_coil(izetal_coil))
     sin_zetal(izetal_coil) = sin(zetal_coil(izetal_coil))
  end do
  ! In contrast to the analogous loop in build_matrice, here "constants" needs dtheta_coil * dzeta_coil
  constants =  dtheta_coil * dzeta_coil * mu0 / (4 * pi)

  allocate(d_times_unit_normal_coil(3,ntheta_coil,nzetal_coil))
  do l_coil = 0, (nfp-1)
     do izeta_coil = 1, nzeta_coil
        izetal_coil = izeta_coil + l_coil*nzeta_coil
        do itheta_coil = 1, ntheta_coil
           d_times_unit_normal_coil(:,itheta_coil,izetal_coil) = sign_normal * d(itheta_coil,izeta_coil) * normal_coil(:,itheta_coil,izetal_coil) / norm_normal_coil(itheta_coil,izeta_coil)
        end do
     end do
  end do

  ! Begin stuff copied from makegrid

  allocate(extcur(1))
  extcur = 1

  print *,"mgrid_ir:",mgrid_ir,"  mgrid_jz:",mgrid_jz,"  mgrid_kp:",mgrid_kp
  allocate(br(mgrid_ir, mgrid_jz, mgrid_kp))
  allocate(bp(mgrid_ir, mgrid_jz, mgrid_kp))
  allocate(bz(mgrid_ir, mgrid_jz, mgrid_kp))
  br = 0
  bp = 0
  bz = 0


  IF (lstell_sym) THEN
     kp2 = kp/2;  jz2 = jz/2
     kp_odd = MOD(kp,2)
     jz_odd = MOD(jz,2)
!        Must be sure zmax = -zmin
     IF (ABS(zmax) > ABS(zmin)) THEN
        zmax = ABS(zmax)
        zmin = -zmax
     ELSE
        zmin = -ABS(zmin)
        zmax = -zmin
     END IF
  ELSE
     kp2 = kp;    jz2 = jz
     kp_odd = 0;  jz_odd = 0
  END IF

  fperiod = (8*ATAN(one))/nfp
  delr = (rmax-rmin)/(ir-1)
  delz = (zmax-zmin)/(jz-1)
  delp = fperiod/kp

  !$omp parallel do default(none) private(i,j,k,rad,phi,zee) shared(br,bp,bz,delr,delp,delz,ir,jz,kp,kp2,jz2,kp_odd,jz_odd,rmin,zmin,lstell_sym,d_times_unit_normal_coil,cos_zetal,sin_zetal,constants)
  DO i=1,ir
     rad = rmin + (i-1)*delr
     !print "(a,i4,a,i5)","  thread",omp_get_thread_num()," handling i=",i
     k = 1                     ! this is always a symmetry plane
     phi = (k-1)*delp
     DO j=1,jz2 + jz_odd
        zee = zmin + (j-1)*delz
        CALL regcoil_mgrid_bfield (rad, phi, zee, br(i,j,k), bp(i,j,k), bz(i,j,k), constants)
        
        IF (lstell_sym) THEN
           br(i,jz+1-j,k) = -br(i,j,k)
           bz(i,jz+1-j,k) =  bz(i,j,k)
           bp(i,jz+1-j,k) =  bp(i,j,k)
        END IF
     END DO
     
     DO k=2,kp2+kp_odd
        phi = (k-1)*delp
        DO j=1,jz
           zee = zmin + (j-1)*delz
           CALL regcoil_mgrid_bfield (rad, phi, zee, br(i,j,k), bp(i,j,k), bz(i,j,k), constants)
           IF (lstell_sym) THEN
              br(i,jz+1-j,kp+2-k) = -br(i,j,k)
              bz(i,jz+1-j,kp+2-k) =  bz(i,j,k)
              bp(i,jz+1-j,kp+2-k) =  bp(i,j,k)
           END IF
        END DO
     END DO

     IF ((kp_odd == 0) .and. lstell_sym) THEN       ! another symmetry plane
        k = kp2 + 1
        phi = (k-1)*delp
        DO j=1,jz2 + jz_odd
           zee = zmin + (j-1)*delz
           CALL regcoil_mgrid_bfield (rad, phi, zee, br(i,j,k), bp(i,j,k), bz(i,j,k), constants)
           br(i,jz+1-j,k) = -br(i,j,k)
           bz(i,jz+1-j,k) =  bz(i,j,k)
           bp(i,jz+1-j,k) =  bp(i,j,k)
        END DO
     END IF

  end do
  !$omp end parallel do

  deallocate(cos_zetal, sin_zetal, d_times_unit_normal_coil)

  ! ---------------------------------------------------------------
  ! Done evaluating B. Now write .nc file.
  ! ---------------------------------------------------------------

  mgrid_file = 'mgrid.nc'

  CALL cdf_open(ngrid,mgrid_file,'w',istat)
  if (istat .ne. 0) STOP 'Error opening mgrid output file'

!
!     DEFINE DATA VARIABLES, DIMENSION NAMES
!
  CALL cdf_define(ngrid, vn_ir, mgrid_ir)
  CALL cdf_define(ngrid, vn_jz, mgrid_jz)
  CALL cdf_define(ngrid, vn_kp, mgrid_kp)
  CALL cdf_define(ngrid, vn_nfp, nfp)
  CALL cdf_define(ngrid, vn_nextcur, nextcur)
  CALL cdf_define(ngrid, vn_rmin, mgrid_rmin)
  CALL cdf_define(ngrid, vn_zmin, mgrid_zmin)
  CALL cdf_define(ngrid, vn_rmax, mgrid_rmax)
  CALL cdf_define(ngrid, vn_zmax, mgrid_zmax)
  IF (nextcur .eq. 1) THEN
     !CALL cdf_define(ngrid, vn_coilgrp,coil_group(1)%s_name) 
     CALL cdf_define(ngrid, vn_coilgrp,coil_group_name) 
  ELSE
     !CALL cdf_define(ngrid, vn_coilgrp,coil_group(1:nextcur)%s_name, dimname=coildim)
     CALL cdf_define(ngrid, vn_coilgrp,coil_group_name, dimname=coildim)
  END IF
  CALL cdf_define(ngrid, vn_mgmode, mgrid_mode)
  CALL cdf_define(ngrid, vn_coilcur, extcur(1:nextcur), dimname=groupdim)
!
!     STORED AS 3D ARRAYS (ACTUALLY 4D, BUT CUT THROUGH IG)
!
  ALLOCATE (vn_br(nextcur), vn_bz(nextcur), vn_bp(nextcur))

  DO ig = 1, nextcur
     write (temp, '(a,i3.3)') "_",ig
     vn_br(ig) = vn_br0 // temp
     vn_bp(ig) = vn_bp0 // temp
     vn_bz(ig) = vn_bz0 // temp
     CALL cdf_define(ngrid, vn_br(ig), br, dimname=cylcoord)
     CALL cdf_define(ngrid, vn_bp(ig), bp, dimname=cylcoord)
     CALL cdf_define(ngrid, vn_bz(ig), bz, dimname=cylcoord)
  END DO


!
!     WRITE OUT DATA
!
  CALL cdf_write(ngrid, vn_ir, mgrid_ir)
  CALL cdf_write(ngrid, vn_jz, mgrid_jz)
  CALL cdf_write(ngrid, vn_kp, mgrid_kp)
  CALL cdf_write(ngrid, vn_nfp, nfp)
  CALL cdf_write(ngrid, vn_nextcur, nextcur)
  CALL cdf_write(ngrid, vn_rmin, mgrid_rmin)
  CALL cdf_write(ngrid, vn_zmin, mgrid_zmin)
  CALL cdf_write(ngrid, vn_rmax, mgrid_rmax)
  CALL cdf_write(ngrid, vn_zmax, mgrid_zmax)
  IF (nextcur .eq. 1) THEN
     !CALL cdf_write(ngrid, vn_coilgrp, coil_group(1)%s_name)
     CALL cdf_write(ngrid, vn_coilgrp, coil_group_name)
  ELSE
     !CALL cdf_write(ngrid, vn_coilgrp, coil_group(1:nextcur)%s_name)
     CALL cdf_write(ngrid, vn_coilgrp, coil_group_name)
  END IF

!
!     SET UP CYLINDRICAL COMPONENTS OF MAGNETIC FIELD ON GRID
!     SUM OVER CURRENT GROUPS IG = 1,NEXTCUR
!     NOTE: USER MUST SUPPLY SUBROUTINE "BFIELD" TO COMPUTE THESE VALUES
!
!  GROUPS: DO ig = 1,nextcur

!     CALL compute_bfield(ig)

  ig = 1  ! We only have 1 group
  CALL cdf_write(ngrid, vn_br(ig), br)
  CALL cdf_write(ngrid, vn_bp(ig), bp)
  CALL cdf_write(ngrid, vn_bz(ig), bz)

!  END DO GROUPS

  CALL cdf_write(ngrid, vn_mgmode, mgrid_mode)
  CALL cdf_write(ngrid, vn_coilcur, extcur(1:nextcur))
  
  CALL cdf_close(ngrid)

  DEALLOCATE (vn_br, vn_bz, vn_bp)

  call system_clock(toc)
  if (verbose) print "(a,f8.3,a)"," Done writing mgrid. Took ",real(toc-tic)/countrate,' sec.'

end subroutine regcoil_write_mgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine regcoil_mgrid_bfield(rr, pp, z, b_r, b_phi, b_z, constants)
    
    use regcoil_variables
    use stel_constants

    implicit none

    real(dp), intent(in) :: rr, pp, z, constants
    real(dp), intent(inout) :: b_r, b_phi, b_z

    integer :: l_coil, izeta_coil, izetal_coil, itheta_coil
    real(dp) :: x, y, dx, dy, dz, dr2inv, dr32inv, factor, cosphi, sinphi
    integer :: js, ks

    b_r = 0
    b_phi = 0
    b_z = 0

    sinphi = sin(pp)
    cosphi = cos(pp)

    x = rr * cosphi
    y = rr * sinphi
    do ks = 1, ns_magnetization
       do izeta_coil = 1, nzeta_coil
          do l_coil = 0, (nfp-1)
             izetal_coil = izeta_coil + l_coil*nzeta_coil
             do itheta_coil = 1, ntheta_coil
                do js = 1, ns_integration
                   dx = x - (r_coil(1,itheta_coil,izetal_coil) + s_integration(js) * d_times_unit_normal_coil(1,itheta_coil,izetal_coil))
                   dy = y - (r_coil(2,itheta_coil,izetal_coil) + s_integration(js) * d_times_unit_normal_coil(2,itheta_coil,izetal_coil))
                   dz = z - (r_coil(3,itheta_coil,izetal_coil) + s_integration(js) * d_times_unit_normal_coil(3,itheta_coil,izetal_coil))
                   
                   dr2inv = 1/(dx*dx + dy*dy + dz*dz)
                   dr32inv = dr2inv*sqrt(dr2inv)
                   factor = dr32inv * Jacobian_coil(itheta_coil, izeta_coil, js) * s_weights(js) &
                        * interpolate_magnetization_to_integration(js, ks) * constants
                   
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ! R component of magnetization
                   ! e_R = e_X * cos(zeta) + e_Y * sin(zeta)
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$                   inductance(index_plasma,index_coil,ks,1) = inductance(index_plasma,index_coil,ks,1) - &
!!$                        (normal_plasma(1,itheta_plasma,izeta_plasma)*cos_zetal(izetal_coil) &
!!$                        +normal_plasma(2,itheta_plasma,izeta_plasma)*sin_zetal(izetal_coil) &
!!$                        - (3*dr2inv) * normal_plasma_dot_dr * &
!!$                        (cos_zetal(izetal_coil)*dx &
!!$                        +sin_zetal(izetal_coil)*dy )) * factor

                   ! For B_R, 
                   ! normal_plasma(1,...) = e_x dot e_R = cosphi.
                   ! normal_plasma(2,...) = e_y dot e_R = sinphi.
                   b_r = b_r - &
                        (cosphi*cos_zetal(izetal_coil) &
                        +sinphi*sin_zetal(izetal_coil) &
                        - (3*dr2inv) * (cosphi * dx + sinphi * dy) * &
                        (cos_zetal(izetal_coil)*dx &
                        +sin_zetal(izetal_coil)*dy )) * factor * magnetization_vector(itheta_coil, izeta_coil, ks, 1, nsaved)

                   ! For B_phi, 
                   ! normal_plasma(1,...) = e_x dot e_phi = -sinphi.
                   ! normal_plasma(2,...) = e_y dot e_phi = cosphi.
                   b_phi = b_phi - &
                        (-sinphi*cos_zetal(izetal_coil) &
                        +cosphi*sin_zetal(izetal_coil) &
                        - (3*dr2inv) * (-sinphi * dx + cosphi * dy) * &
                        (cos_zetal(izetal_coil)*dx &
                        +sin_zetal(izetal_coil)*dy )) * factor * magnetization_vector(itheta_coil, izeta_coil, ks, 1, nsaved)

                   b_z = b_z - &
                        (- (3*dr2inv) * dz * &
                        (cos_zetal(izetal_coil)*dx &
                        +sin_zetal(izetal_coil)*dy )) * factor * magnetization_vector(itheta_coil, izeta_coil, ks, 1, nsaved)
                   
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ! zeta component of magnetization
                   ! e_zeta = e_X * (-sin(zeta)) + e_Y * cos(zeta)
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$                   inductance(index_plasma,index_coil,ks,2) = inductance(index_plasma,index_coil,ks,2) - &
!!$                        (normal_plasma(1,itheta_plasma,izeta_plasma)*(-sin_zetal(izetal_coil)) &
!!$                        +normal_plasma(2,itheta_plasma,izeta_plasma)*  cos_zetal(izetal_coil) &
!!$                        - (3*dr2inv) * normal_plasma_dot_dr * &
!!$                        (-sin_zetal(izetal_coil)*dx &
!!$                        + cos_zetal(izetal_coil)*dy )) * factor

                   ! For B_R, 
                   ! normal_plasma(1,...) = e_x dot e_R = cosphi.
                   ! normal_plasma(2,...) = e_y dot e_R = sinphi.
                   b_r = b_r - &
                        (cosphi*(-sin_zetal(izetal_coil)) &
                        +sinphi*  cos_zetal(izetal_coil) &
                        - (3*dr2inv) * (cosphi * dx + sinphi * dy) * &
                        (-sin_zetal(izetal_coil)*dx &
                        + cos_zetal(izetal_coil)*dy )) * factor * magnetization_vector(itheta_coil, izeta_coil, ks, 2, nsaved)
                   
                   ! For B_phi, 
                   ! normal_plasma(1,...) = e_x dot e_phi = -sinphi.
                   ! normal_plasma(2,...) = e_y dot e_phi = cosphi.
                   b_phi = b_phi - &
                        (-sinphi*(-sin_zetal(izetal_coil)) &
                        +cosphi*  cos_zetal(izetal_coil) &
                        - (3*dr2inv) * (-sinphi * dx + cosphi * dy) * &
                        (-sin_zetal(izetal_coil)*dx &
                        + cos_zetal(izetal_coil)*dy )) * factor * magnetization_vector(itheta_coil, izeta_coil, ks, 2, nsaved)
                   
                   b_z = b_z - &
                        (- (3*dr2inv) * dz * &
                        (-sin_zetal(izetal_coil)*dx &
                        + cos_zetal(izetal_coil)*dy )) * factor * magnetization_vector(itheta_coil, izeta_coil, ks, 2, nsaved)
                   
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   ! Z component of magnetization
                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$                   inductance(index_plasma,index_coil,ks,3) = inductance(index_plasma,index_coil,ks,3) - &
!!$                        (normal_plasma(3,itheta_plasma,izeta_plasma) &
!!$                        - (3*dr2inv) * normal_plasma_dot_dr * &
!!$                        dz) * factor
                   
                   ! For B_R, 
                   ! normal_plasma(1,...) = e_x dot e_R = cosphi.
                   ! normal_plasma(2,...) = e_y dot e_R = sinphi.
                   b_r = b_r - &
                        ( - (3*dr2inv) * (cosphi * dx + sinphi * dy) * &
                        dz) * factor * magnetization_vector(itheta_coil, izeta_coil, ks, 3, nsaved)
                   
                   ! For B_phi, 
                   ! normal_plasma(1,...) = e_x dot e_phi = -sinphi.
                   ! normal_plasma(2,...) = e_y dot e_phi = cosphi.
                   b_phi = b_phi - &
                        (- (3*dr2inv) * (-sinphi * dx + cosphi * dy) * &
                        dz) * factor * magnetization_vector(itheta_coil, izeta_coil, ks, 3, nsaved)
                   
                   b_z = b_z - &
                        (1 &
                        - (3*dr2inv) * dz * &
                        dz) * factor * magnetization_vector(itheta_coil, izeta_coil, ks, 3, nsaved)
                   
                end do
             end do
          end do
       end do
    end do

    ! Perhaps the TF field in this next line should only be included some of the time, not always?
    b_phi = b_phi + mu0 * net_poloidal_current_Amperes / (2 * pi * rr)

!!$    if (abs(b_r) < 1e-14) then
!!$       print "(4(a,es12.5))","  Small b_r=",b_r," for rr=",rr," pp=",pp," z=",z
!!$    end if
!!$    if (abs(b_r) >1000) then
!!$       print "(4(a,es12.5))","  Large b_r=",b_r," for rr=",rr," pp=",pp," z=",z
!!$    end if
!!$    if (abs(b_phi) >1000) then
!!$       print "(4(a,es12.5))","  Large b_phi=",b_phi," for rr=",rr," pp=",pp," z=",z
!!$    end if
!!$    if (abs(b_z) >1000) then
!!$       print "(4(a,es12.5))","  Large b_z=",b_z," for rr=",rr," pp=",pp," z=",z
!!$    end if

  end subroutine regcoil_mgrid_bfield


