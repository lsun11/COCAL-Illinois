module interface_IO_input_CF_irrot_BNS
  implicit none
  interface 
    subroutine IO_input_CF_irrot_BNS(impt, emd,vep,ome,ber,radi,rs,psi,alph,bvxd,bvyd,bvzd)
      real(8), pointer :: emd(:,:,:), vep(:,:,:), rs(:,:), psi(:,:,:), alph(:,:,:), &
         &                bvxd(:,:,:), bvyd(:,:,:), bvzd(:,:,:)
      real(8) ::  ome,ber,radi
    end subroutine IO_input_CF_irrot_BNS 
  end interface
end module interface_IO_input_CF_irrot_BNS
