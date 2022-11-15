module def_emfield
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  va(:,:,:), alva(:,:,:)
  real(long), pointer  ::  vaxd(:,:,:), vayd(:,:,:), vazd(:,:,:)
  real(long), pointer  ::  vaxu(:,:,:), vayu(:,:,:), vazu(:,:,:)
  real(long), pointer  ::  jtd(:,:,:), jxd(:,:,:), jyd(:,:,:), jzd(:,:,:)
  real(long), pointer  ::  jtu(:,:,:), jxu(:,:,:), jyu(:,:,:), jzu(:,:,:)
  real(long), pointer  ::  jtuf(:,:,:),jxuf(:,:,:),jyuf(:,:,:),jzuf(:,:,:)
end module def_emfield
