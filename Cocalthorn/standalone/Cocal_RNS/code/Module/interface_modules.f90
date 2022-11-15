module interface_modules
  use phys_constant, only : long
  implicit none
  interface 
!______________________________________________________________________
    subroutine poisson_solver(sou, pot)
      real(8),pointer  ::  sou(:,:,:),pot(:,:,:)
    end subroutine poisson_solver
!    subroutine poisson_solver_no_dipole(sou, pot)
!      real(8),pointer  ::  sou(:,:,:),pot(:,:,:)
!    end subroutine poisson_solver_no_dipole
    subroutine poisson_solver_no_lm1(sou, pot)
      real(8),pointer  ::  sou(:,:,:),pot(:,:,:)
    end subroutine poisson_solver_no_lm1
    subroutine grgrad_midpoint(fnc,dfdx,dfdy,dfdz)
      real(8), pointer :: fnc(:,:,:)
      real(8), pointer :: dfdx(:,:,:)
      real(8), pointer :: dfdy(:,:,:)
      real(8), pointer :: dfdz(:,:,:)
    end subroutine grgrad_midpoint
    subroutine grgrad_4th(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8), intent(out) :: dfdx, dfdy, dfdz
      integer :: irg, itg, ipg
    end subroutine grgrad_4th
    subroutine interpo_linear_type0(val,fnc,ir,it,ip)
      real(8) :: val
      real(8), pointer :: fnc(:,:,:)
      integer :: ir, it, ip
    end subroutine interpo_linear_type0
    subroutine interpo_linear_type0_2Dsurf(val,fnc,it,ip)
      real(8) :: val
      real(8), pointer :: fnc(:,:)
      integer :: it, ip
    end subroutine interpo_linear_type0_2Dsurf
    subroutine interpo_lag4th_2Dsurf(val,fnc,tv,pv)
      real(8), intent(out) :: val
      real(8), intent(in)  :: tv, pv
      real(8), pointer :: fnc(:,:)
    end subroutine interpo_lag4th_2Dsurf
    subroutine interpo_gr2fl(a,b)
      real(8), pointer :: a(:,:,:), b(:,:,:)
    end subroutine interpo_gr2fl
    subroutine interpo_fl2gr(a,b)
      real(8), pointer :: a(:,:,:), b(:,:,:)
    end subroutine interpo_fl2gr
    subroutine interpo_radial1p_grav(grv,val,rv,it,ip)
      real(8), pointer :: grv(:,:,:)
      real(8), intent(out) :: val
      real(8), intent(in)  :: rv
      integer, intent(in)  :: it, ip
    end subroutine interpo_radial1p_grav
!______________________________________________________________________
    subroutine sourceterm_HaC(sou)    
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_HaC
    subroutine sourceterm_trG(sou)
      real(8), pointer :: sou(:,:,:)
    end subroutine sourceterm_trG
    subroutine sourceterm_MoC(souvec,sou)
      real(8), pointer :: sou(:,:,:) 
      real(8), pointer :: souvec(:,:,:,:)
    end subroutine sourceterm_MoC
    subroutine compute_shift(potx,poty,potz,gvec,bfnc)      
      real(8), pointer :: bfnc(:,:,:) 
      real(8), pointer :: gvec(:,:,:,:)
      real(8), pointer :: potx(:,:,:), poty(:,:,:), potz(:,:,:)
    end subroutine compute_shift
    subroutine compute_alps2alph(pot,psi)
      real(8), pointer     :: pot(:,:,:), psi(:,:,:)
    end subroutine compute_alps2alph
!______________________________________________________________________
    subroutine hydrostatic_eq(emd)
      real(8), pointer :: emd(:,:,:)
    end subroutine hydrostatic_eq
    subroutine calc_surface(rsnew,emd)
      real(8), pointer :: emd(:,:,:)
      real(8), pointer :: rsnew(:,:)
    end subroutine calc_surface
    subroutine error_matter(pot,pot_bak,error,flag)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
    end subroutine error_matter
!______________________________________________________________________ 
    subroutine update_grfield(pot,grfield,convf)
      real(8), pointer :: pot(:,:,:)
      real(8), pointer :: grfield(:,:,:)
      real(8), intent(in) :: convf
    end subroutine update_grfield
    subroutine update_matter(potf,mtfield,convf)
      real(8), pointer :: potf(:,:,:)
      real(8), pointer :: mtfield(:,:,:)
      real(8), intent(in) :: convf
    end subroutine update_matter
    subroutine update_surface(potrs,rsnew,convf)
      real(8), pointer :: potrs(:,:)
      real(8), pointer :: rsnew(:,:)
      real(8), intent(in) :: convf
    end subroutine update_surface
    subroutine update_parameter(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter
    subroutine update_parameter_axisym(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_axisym
    subroutine update_parameter_triaxial(convf)
      real(8), intent(in) :: convf
    end subroutine update_parameter_triaxial
!______________________________________________________________________
    subroutine source_rest_mass(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_rest_mass
    subroutine source_adm_mass(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_adm_mass
    subroutine source_komar_mass(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_komar_mass
    subroutine source_komar_mass_compact(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_komar_mass_compact
    subroutine source_proper_mass(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_proper_mass
    subroutine source_ang_mom(souf)
      real(8), pointer     :: souf(:,:,:)
    end subroutine source_ang_mom
    subroutine source_ang_mom_asymp(sousf,irg)
      real(8), pointer     :: sousf(:,:)
      integer, intent(in)  :: irg
    end subroutine source_ang_mom_asymp
    subroutine source_mp_minus_madm(soug,souf)
      real(8), pointer     :: soug(:,:,:), souf(:,:,:)
    end subroutine source_mp_minus_madm
    subroutine vol_int_fluid(souf,vol)
      real(8), pointer     :: souf(:,:,:)
      real(8), intent(out) :: vol
    end subroutine vol_int_fluid
    subroutine vol_int_grav(soug,vol)
      real(8), pointer     :: soug(:,:,:)
      real(8), intent(out) :: vol
    end subroutine  vol_int_grav
    subroutine surf_int_grav(sousf,surf,irg)
      real(8), pointer     :: sousf(:,:)
      real(8), intent(out) :: surf
      integer, intent(in)  :: irg
    end subroutine  surf_int_grav
    subroutine radial_int_fluid(sou,radius,it,ip)
      real(8), pointer     :: sou(:)
      real(8), intent(out) :: radius
      integer, intent(in)  :: it, ip
    end subroutine radial_int_fluid
!______________________________________________________________________
    subroutine minv(aa,bb,nn,nnz)
      integer :: nn, nnz
      real(8) :: aa(nnz,nnz),bb(nnz)
    end subroutine minv
!______________________________________________________________________
  end interface
end module interface_modules
