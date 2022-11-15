module def_emfield_cartesian
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  vaca(:,:,:), alvaca(:,:,:)
  real(long), pointer  ::  vaxdca(:,:,:), vaydca(:,:,:), vazdca(:,:,:)
  real(long), pointer  ::  vaxuca(:,:,:), vayuca(:,:,:), vazuca(:,:,:)
  real(long), pointer  ::  jtdca(:,:,:),jxdca(:,:,:),jydca(:,:,:),jzdca(:,:,:)
  real(long), pointer  ::  jtuca(:,:,:),jxuca(:,:,:),jyuca(:,:,:),jzuca(:,:,:)
end module def_emfield_cartesian
