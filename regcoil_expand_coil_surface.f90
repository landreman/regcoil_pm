subroutine regcoil_expand_coil_surface(theta, zeta, x,y,z)
  
  use regcoil_variables, only: nfp, xm_coil, xn_coil, mnmax_coil, rmnc_coil, zmns_coil, rmns_coil, zmnc_coil, lasym, d0, dmns, dmnc, &
       mnmax_magnetization, xm_magnetization, xn_magnetization, sign_normal
  use stel_kinds
  use stel_constants

  implicit none

  real(dp), intent(in) :: theta,zeta
  real(dp), intent(out) :: x,y,z

  integer :: imn
  real(dp) :: dxdtheta, dxdzeta, dydtheta, dydzeta, dzdtheta, dzdzeta
  real(dp) :: angle, sinangle, cosangle
  real(dp) :: dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta
  real(dp) :: phi, sinphi, cosphi
  real(dp) :: dsinphidzeta, dcosphidzeta
  real(dp) :: normal_x, normal_y, normal_z, normal_abs
  real(dp) :: separation

  x=0
  y=0
  z=0
  dxdtheta=0
  dydtheta=0
  dzdtheta=0
  dxdzeta=0
  dydzeta=0
  dzdzeta=0

  phi = zeta
  sinphi = sin(phi)
  cosphi = cos(phi)
  dsinphidzeta = cosphi
  dcosphidzeta = -sinphi

  do imn = 1,mnmax_coil
     angle = xm_coil(imn)*theta - xn_coil(imn)*zeta
     cosangle = cos(angle)
     sinangle = sin(angle)
     dsinangledtheta = cosangle*xm_coil(imn)
     dcosangledtheta = -sinangle*xm_coil(imn)
     dsinangledzeta = -cosangle*xn_coil(imn)
     dcosangledzeta = sinangle*xn_coil(imn)

     x = x + rmnc_coil(imn) * cosangle * cosphi
     y = y + rmnc_coil(imn) * cosangle * sinphi
     z = z + zmns_coil(imn) * sinangle
     
     dxdtheta = dxdtheta + rmnc_coil(imn) * dcosangledtheta * cosphi
     dydtheta = dydtheta + rmnc_coil(imn) * dcosangledtheta * sinphi
     dzdtheta = dzdtheta + zmns_coil(imn) * dsinangledtheta
     
     dxdzeta = dxdzeta + rmnc_coil(imn) * (dcosangledzeta * cosphi + cosangle * dcosphidzeta)
     dydzeta = dydzeta + rmnc_coil(imn) * (dcosangledzeta * sinphi + cosangle * dsinphidzeta)
     dzdzeta = dzdzeta + zmns_coil(imn) * dsinangledzeta

     if (lasym) then
        x = x + rmns_coil(imn) * sinangle * cosphi
        y = y + rmns_coil(imn) * sinangle * sinphi
        z = z + zmnc_coil(imn) * cosangle
     
        dxdtheta = dxdtheta + rmns_coil(imn) * dsinangledtheta * cosphi
        dydtheta = dydtheta + rmns_coil(imn) * dsinangledtheta * sinphi
        dzdtheta = dzdtheta + zmnc_coil(imn) * dcosangledtheta
     
        dxdzeta = dxdzeta + rmns_coil(imn) * (dsinangledzeta * cosphi + sinangle * dcosphidzeta)
        dydzeta = dydzeta + rmns_coil(imn) * (dsinangledzeta * sinphi + sinangle * dsinphidzeta)
        dzdzeta = dzdzeta + zmnc_coil(imn) * dcosangledzeta
     end if
     
  end do

  ! Evaluate the (non-unit-magnitude) surface normal vector from N = (dr/dv) cross (dr/du):
  normal_x = dydzeta * dzdtheta - dydtheta * dzdzeta
  normal_y = dzdzeta * dxdtheta - dzdtheta * dxdzeta
  normal_z = dxdzeta * dydtheta - dxdtheta * dydzeta

  ! Now normalize the normal vector so it has unit magnitude:
  normal_abs = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z)
  normal_x = normal_x / normal_abs
  normal_y = normal_y / normal_abs
  normal_z = normal_z / normal_abs

  separation = d0
  do imn = 1,mnmax_magnetization
     angle = xm_magnetization(imn)*theta - xn_magnetization(imn)*zeta
     cosangle = cos(angle)
     separation = separation + dmnc(imn) * cosangle
     if (lasym) then
        sinangle = sin(angle)
        separation = separation + dmns(imn) * sinangle
     end if
  end do
  separation = separation * sign_normal

  ! Move in the normal direction away from the inner magnetization surface:
  x = x + normal_x * separation
  y = y + normal_y * separation
  z = z + normal_z * separation

end subroutine regcoil_expand_coil_surface
