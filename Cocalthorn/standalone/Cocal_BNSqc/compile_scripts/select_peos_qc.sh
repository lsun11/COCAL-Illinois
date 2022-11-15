cp ../code/EOS/Subroutine_mpatch/peos_qc_initialize_mpt.f90  ../code/EOS/Subroutine_mpatch/peos_initialize_mpt.f90


cp ../code/EOS/Subroutine/peos_qc_initialize.f90  ../code/EOS/Subroutine/peos_initialize.f90

cp ../code/EOS/Subroutine/peos_qc_lookup.f90      ../code/EOS/Subroutine/peos_lookup.f90

cp ../code/EOS/Subroutine/peos_qc_h2qprho.f90     ../code/EOS/Subroutine/peos_h2qprho.f90

cp ../code/EOS/Subroutine/peos_qc_q2hprho.f90     ../code/EOS/Subroutine/peos_q2hprho.f90
sed -i "s/!s1//g"                                 ../code/EOS/Subroutine/peos_q2hprho.f90
#sed -i "s/!sg//g"                                 ../code/EOS/Subroutine/peos_q2hprho.f90

