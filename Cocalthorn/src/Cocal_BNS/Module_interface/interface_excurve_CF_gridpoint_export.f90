module interface_excurve_CF_gridpoint_export
  implicit none
  interface 
    subroutine excurve_CF_gridpoint_export(alph,bvxd,bvyd,bvzd,tfkxx,tfkxy,tfkxz,tfkyy,tfkyz,tfkzz)
      real(8), pointer ::  alph(:,:,:), bvxd(:,:,:), bvyd(:,:,:), bvzd(:,:,:), tfkxx(:,:,:),  &
        &    tfkxy(:,:,:), tfkxz(:,:,:), tfkyy(:,:,:), tfkyz(:,:,:), tfkzz(:,:,:)
    end subroutine excurve_CF_gridpoint_export
  end interface
end module interface_excurve_CF_gridpoint_export
