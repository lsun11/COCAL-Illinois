# This may sometimes cause a non convergence of a solution
sed -i "s/!cmout//g" ../code/Subroutine/correct_matter_source_midpoint.f90
sed -i -e '19,23d'   ../code/Subroutine/correct_matter_source_midpoint.f90
