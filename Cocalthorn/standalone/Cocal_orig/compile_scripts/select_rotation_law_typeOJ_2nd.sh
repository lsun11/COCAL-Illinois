cp ../code/Subroutine/calc_omega_drot_bisection.f90_typeOJ \
   ../code/Subroutine/calc_omega_drot_bisection.f90
cp ../code/Subroutine/update_parameter_axisym_peos_drot.f90_typeOJ_2nd \
   ../code/Subroutine/update_parameter_axisym_peos_drot.f90

sed -i "s/!rotlaw_OJ2nd//g" ../code/Subroutine/calc_omega_drot.f90
sed -i "s/!rotlaw_typeOJ//g" ../code/Subroutine/calc_omega_drot_secant.f90
sed -i "s/!rotlaw_typeOJ//g" ../code/Subroutine/calc_omega_drot_bisection.f90
sed -i "s/!rotlaw_typeOJ//g" ../code/Subroutine/iteration_peos.f90
