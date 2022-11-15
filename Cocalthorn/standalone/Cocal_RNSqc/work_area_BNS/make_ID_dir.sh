rm -rf ID_dir
mkdir ID_dir
cp *.dat *.ini *.las ./ID_dir/
cd ID_dir
sed -i "s/d/e/g" peos_parameter*.dat
sed -i "s/E/e/g" peos_parameter*.dat
