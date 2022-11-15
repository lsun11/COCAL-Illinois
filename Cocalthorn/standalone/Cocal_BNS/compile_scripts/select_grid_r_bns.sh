sed -i "s/  call grid_r/  call grid_r_bns/g" \
       ../code/Subroutine_mpatch/coordinate_patch_kit_grav_mpt.f90
sed -i "s/  call grid_r/  call grid_r_bns/g" \
       ../code/Analysis/Subroutine/coordinate_patch_kit_grav_noGreen_mpt.f90
sed -i "s/  call grid_r/  call grid_r_bns/g" \
       ../code/Subroutine/coordinate_patch_kit_bhex.f90
sed -i "s/  call grid_r/  call grid_r_bns/g" \
       ../code/Analysis/Subroutine/coordinate_patch_kit_grav_noGreen.f90

