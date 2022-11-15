cp ../code/Subroutine/calc_omega_drot_Newton.f90_type0 \
   ../code/Subroutine/calc_omega_drot_Newton.f90
cp ../code/Subroutine/calc_omega_drot_bisection.f90_type0 \
   ../code/Subroutine/calc_omega_drot_bisection.f90
cp ../code/Subroutine/update_parameter_axisym_peos_drot.f90_type0 \
   ../code/Subroutine/update_parameter_axisym_peos_drot.f90

sed -i "s/!rotlaw_type0//g" ../code/Subroutine/calc_omega_drot.f90
sed -i "s/!rotlaw_type012//g" ../code/Subroutine/calc_omega_drot_secant.f90
sed -i "s/!rotlaw_type01//g" ../code/Subroutine/calc_omega_drot_secant.f90
sed -i "s/!rotlaw_type0//g" ../code/Subroutine/calc_omega_drot_secant.f90
#sed -i "s/!rotlaw_type12OJ//g" ../code/Subroutine/iteration_peos.f90
