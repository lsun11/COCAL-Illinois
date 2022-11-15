sed -i "s/!  call grid_r_bhex/  call grid_r_bhex/g" \
       ../code/Subroutine_mpatch/coordinate_patch_kit_grav_mpt.f90
sed -i "s/!  call grid_r_bhex/  call grid_r_bhex/g" \
       ../code/Analysis/Subroutine/coordinate_patch_kit_grav_noGreen_mpt.f90
sed -i "s/!  call grid_r_bhex/  call grid_r_bhex/g" \
       ../code/Subroutine/coordinate_patch_kit_bhex.f90
sed -i "s/!  call grid_r_bhex/  call grid_r_bhex/g" \
       ../code/Analysis/Subroutine/coordinate_patch_kit_grav_noGreen.f90
sed -i "s/pBH/eBH/g" \
       ../code/Subroutine_mpatch/coordinate_patch_kit_grav_mpt.f90
sed -i "s/pBH/eBH/g" \
       ../code/Analysis/Subroutine/coordinate_patch_kit_grav_noGreen_mpt.f90
sed -i "s/pBH/eBH/g" \
       ../code/Subroutine/coordinate_patch_kit_bhex.f90
sed -i "s/pBH/eBH/g" \
       ../code/Analysis/Subroutine/coordinate_patch_kit_grav_noGreen.f90
