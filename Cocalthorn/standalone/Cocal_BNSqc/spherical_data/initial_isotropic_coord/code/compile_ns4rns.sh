#! /usr/bin/tcsh -f
#pgf77 -Bstatic -fastsse -tp athlonxp -r8 -o ahotest GR_IRDR_ver4.f
#g77 -fno-automatic -O -o ahotest GR_IRDR_ver4.f
#pgf77 -Bstatic -fastsse -tp k8-64 -o ahotest GR_KEH_NS_ver1.f
/opt/intel/fce/9.1.040/bin/ifort -fast -o exe_1Dini GR_KEH_NS_ver1_4rns.f

mv exe_1Dini ../work_area_1D/.
