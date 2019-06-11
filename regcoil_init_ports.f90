subroutine regcoil_init_ports()

  use regcoil_variables
  use stel_kinds
  
  implicit none

  integer :: itheta, izeta, j_port, j_sign, sign
  real(dp) :: f

  allocate(ports_weight(ntheta_coil, nzeta_coil))

  ports_weight = 0
  do itheta = 1, ntheta_coil
     do izeta = 1, nzeta_coil
        do j_sign = 1,2 ! Adds a stellarator-symmetric copy of the port
           sign = j_sign*2 - 3 ! So sign is either -1 or +1
           do j_port = 1, nports
              f = 1 - 2*(1 - cos(theta_coil(itheta) - sign*ports_theta0(j_port)))/(ports_theta_width(j_port) ** 2) &
                   - 2*(1 - cos(nfp*(zeta_coil(izeta) - sign*ports_zeta0(j_port))))/(nfp*nfp*ports_zeta_width(j_port) ** 2)
              ports_weight(itheta,izeta) = ports_weight(itheta,izeta) + (tanh(ports_sharpness*f)+1)/2
           end do
        end do
     end do
  end do

  ports_weight = 1 + ports_weight * ports_magnitude

end subroutine  regcoil_init_ports
