/*
 * Reads binary neutron star initial data, exporting Cocal structures
 * onto standard C arrays  on a Cartesian grid.
*/


// C headers                                                                                      
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <time.h>
#include <ctype.h>
#include <new>

using namespace std;


// Isolated Rotating star
extern "C" void coc2pri_rs_(char *dir_path, double Xb[], double Yb[], double Zb[],
			    int *Nx, int *Ny, int *Nz,
			    double *gtt, double *gtx, double *gty,double *gtz, 
			    double *gxx, double *gxy,double *gxz,
			    double *gyy, double *gyz,double *gzz,
			    double *Axx,double *Axy, double *Axz,
			    double *Ayy,double *Ayz, double *Azz,
			    double *Pb,double *rhob,
			    double *vbx ,  double *vby, double *vbz,size_t *path_len);


int main(int argc, char** argv){

  if (argv[1]==NULL || argv[10]==NULL) {
    printf("Sorry.  I was expecting 10 command line arguments:\n");
    printf("Usage: ./read_bin_ns <xmin> <ymin> <zmin> <dx> <dy> <dz> <nx> <ny> <nz> <processor number>\n");
    exit(1);
  }


  //  Reads Cartesian grid parameters                                                                                                                                                                
        
  int Nx, Ny, Nz ;
  double x_min, x_max, y_min, y_max, z_min, z_max ;


  x_min = strtod(argv[1],NULL) ;
  y_min = strtod(argv[2],NULL) ;
  z_min = strtod(argv[3],NULL) ;

  double dx = strtod(argv[4],NULL) ;
  double dy = strtod(argv[5],NULL) ;
  double dz = strtod(argv[6],NULL) ;
  Nx = (int)strtod(argv[7],NULL) ;
  Ny = (int)strtod(argv[8],NULL) ;
  Nz = (int)strtod(argv[9],NULL) ;
  int procnumber = (int)strtod(argv[10],NULL) ;
  int ntotal = Nx*Ny*Nz;


  // define size of the grid

  x_max = x_min+(Nx-1)*dx;
  y_max = y_min+(Ny-1)*dy;
  z_max = z_min+(Nz-1)*dz;


  // screen output
  cout << "---------------- " << endl ;
  cout << "Cartesian grid : " << endl ;
  cout << "---------------- " << endl ;
  cout << "   Number of points in the x direction : " << Nx << endl ;
  cout << "   Number of points in the y direction : " << Ny << endl ;
  cout << "   Number of points in the z direction : " << Nz << endl ;
  cout << "   x_min, x_max : " << x_min << " , " << x_max << endl ;
  cout << "   y_min, y_max : " << y_min << " , " << y_max << endl ;
  cout << "   z_min, z_max : " << z_min << " , " << z_max << endl ;



  // define 1D arrays
  printf("Allocating ID arrays\n");

  //Define size of the grids
  double Xb[Nx], Yb[Ny], Zb[Nz];

  //Define full metric
  double *gtt, *gtx,*gty,*gtz;
  double *gxx,*gxy,*gxz;
  double *gyy,*gyz,*gzz;

  //Define traceless curvature extrinsic
  double *Axx,*Axy,*Axz;
  double *Ayy,*Ayz,*Azz;

  //Define fluid variables
  double *rhob, *Pb;
  double *vbx,*vby,*vbz;  

  // Allocate memory
  gtt = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gtt
  gtx = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gtx
  gty = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gty
  gtz = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gtz
 
  gxx = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gxx
  gxy = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gxy
  gxz = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gxz
  gyy = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gyy
  gyz = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gyz
  gzz = new (nothrow) double [Nx*Ny*Nz];  // metric coefficient gzz


  Axx = new (nothrow) double [Nx*Ny*Nz];  // traceless extrinsic curv coefficient Axx
  Axy = new (nothrow) double [Nx*Ny*Nz];  // traceless extrinsic curv coefficient Axy
  Axz = new (nothrow) double [Nx*Ny*Nz];  // traceless extrinsic curv coefficient Axz
  Ayy = new (nothrow) double [Nx*Ny*Nz];  // traceless extrinsic curv coefficient Ayy
  Ayz = new (nothrow) double [Nx*Ny*Nz];  // traceless extrinsic curv coefficient Ayz
  Azz = new (nothrow) double [Nx*Ny*Nz];  // traceless extrinsic curv coefficient Azz

  rhob = new (nothrow) double [Nx*Ny*Nz];  // baryon density
  Pb   = new (nothrow) double [Nx*Ny*Nz];  // pressure

  vbx = new (nothrow) double [Nx*Ny*Nz];   // three velocity v^x measured by Euler observers
  vby = new (nothrow) double [Nx*Ny*Nz];   // three velocity v^y measured by Euler observers
  vbz = new (nothrow) double [Nx*Ny*Nz];   // three velocity v^z measured by Euler observers



  if (gtt == NULL){ cout << "Stopping: memory could not be allocated for a1" << endl;}

  // Construction of the Cartesian grid

  for (int k=0; k<Nz; k++) {Zb[k] =  z_min + dz*k;}
  for (int j=0; j<Ny; j++) {Yb[j] =  y_min + dy*j;}
  for (int i=0; i<Nx; i++) {Xb[i] =  x_min + dx*i;}


  //  We have defined 1D C arrays for the metric and extrinsic curvature
  //  Thus, we need to define a relation for how to translate from indices i,j,k
  //  to the 1D c arrays. The standard way to do this 

  // get the path for the ID
  printf("Getting the path for the ID\n");
  ifstream fparam("read_bin_path.par");

  char path[400],comment[400],id_type[3];
  char *dir_path, *file_path;

  fparam.getline(comment, 400);
  fparam.getline(comment, 400);
  fparam.getline(path, 400);
  fparam.close();

  printf("path is %s\n",path);

  dir_path=path;

  size_t path_len=strlen(dir_path);
  
  file_path=strcat(path,"/rnspar_mpt1.dat");
  size_t path_len_file=strlen(file_path);


  // // Read ID type
  // int ierr=1;
  // read_id_type_(file_path,id_type,&path_len_file,&ierr);
  // id_type[2]='\0'; // NULL terminate the string
  
  // if(ierr==0){
  //   printf("COCAL_ID:: success reading rnspar_mpt1.dat.\n");
  // }
  // else{
  //   printf("COCAL_ID:: Can't open the rnspar_mpt1.dat.");
  //   exit(1);
  // }


  // Interpolate ID
  printf("COCAL_ID:: Reading ROTATING STAR ID\n");
  coc2pri_rs_(dir_path, Xb, Yb, Zb,
	      &Nx, &Ny, &Nz,
	      gtt, gtx, gty,gtz,
	      gxx, gxy,gxz,
	      gyy, gyz,gzz,
	      Axx, Axy, Axz,
	      Ayy, Ayz, Azz,
	      Pb,rhob,
	      vbx, vby, vbz, &path_len);

  
  
  /*Use Processor Number to name the output file

    Note that procnumber below is an integer that uniquely specifies the chunk
    of grid we want to output.  For an N MPI-processes run, with M levels, there
    are N*M initial data files and unique procnumber's.  We have found that
    this works best when setting up initial data on grids.  All the initial
    data binary file output routine should need is: 
    xmin,ymin,zmin,dx,dy,dz,nx,ny,nz,procnumber*/
  

  ofstream outfilepsi;
  char fname_psi[100];

  printf("********************************\n");
  printf("Dumping ID data file = %i\n",procnumber);
  printf("********************************\n");


  sprintf(fname_psi,"CTS_bin-proc%d.d",procnumber);
  outfilepsi.open(fname_psi,ios::out | ios::ate);


  // Grid parameters
  outfilepsi.write((char *) &(Nx),sizeof(int));
  outfilepsi.write((char *) &(Ny),sizeof(int));
  outfilepsi.write((char *) &(Nz),sizeof(int));
  outfilepsi.write((char *) &(ntotal),sizeof(int));


  // 3D Arrays 

  int intdum=0;
  for (int k=0; k<Nz; k++) {
    for (int j=0; j<Ny; j++) {
      for (int i=0; i<Nx; i++) {

	double gtt_val = gtt[intdum];
	double gtx_val = gtx[intdum];
	double gty_val = gty[intdum];
	double gtz_val = gtz[intdum];

	double gxx_val = gxx[intdum];
	double gxy_val = gxy[intdum];
	double gxz_val = gxz[intdum];
	double gyy_val = gyy[intdum];
	double gyz_val = gyz[intdum];
	double gzz_val = gzz[intdum];

	double Axx_val = Axx[intdum];
        double Axy_val = Axy[intdum];
        double Axz_val = Axz[intdum];
        double Ayy_val = Ayy[intdum];
        double Ayz_val = Ayz[intdum];
        double Azz_val = Azz[intdum];

	double rho0_val = rhob[intdum];

        double ux_val   = vbx[intdum];
        double uy_val   = vby[intdum];
        double uz_val   = vbz[intdum];

	//	printf("%f %f %f %f %i\n",Xb[intdum],Yb[intdum],Zb[intdum],rhob[intdum],procnumber);

        outfilepsi.write((char *) &gtt_val,sizeof(double));
        outfilepsi.write((char *) &gtx_val,sizeof(double));
        outfilepsi.write((char *) &gty_val,sizeof(double));
        outfilepsi.write((char *) &gtz_val,sizeof(double));

        outfilepsi.write((char *) &gxx_val,sizeof(double));
        outfilepsi.write((char *) &gxy_val,sizeof(double));
        outfilepsi.write((char *) &gxz_val,sizeof(double));
        outfilepsi.write((char *) &gyy_val,sizeof(double));
        outfilepsi.write((char *) &gyz_val,sizeof(double));
        outfilepsi.write((char *) &gzz_val,sizeof(double));

        outfilepsi.write((char *) &Axx_val,sizeof(double));
        outfilepsi.write((char *) &Axy_val,sizeof(double));
        outfilepsi.write((char *) &Axz_val,sizeof(double));
        outfilepsi.write((char *) &Ayy_val,sizeof(double));
        outfilepsi.write((char *) &Ayz_val,sizeof(double));
        outfilepsi.write((char *) &Azz_val,sizeof(double));

        outfilepsi.write((char *) &rho0_val,sizeof(double));

        outfilepsi.write((char *) &ux_val,sizeof(double));
        outfilepsi.write((char *) &uy_val,sizeof(double));
        outfilepsi.write((char *) &uz_val,sizeof(double));

        intdum++;
      }
    }
  }


  // EoS parameters

  double gamma = 1.0;                // polytropic index n
//  double kappa = 22654787.871; // polytropic constant
  double kappa = 123.54631384041791; // polytropic constant

  outfilepsi.write((char *) &(gamma),sizeof(double));
  outfilepsi.write((char *) &(kappa),sizeof(double));
  outfilepsi.close();


  printf("********************************\n");
  printf("Done with file =%i\n",procnumber);
  printf("********************************\n");


  //release memory
  delete[]  gtt,gtx,gty,gtz;
  delete[]  gxx,gxy,gxz; 
  delete[]  gyy,gyz,gzz; 

  delete[] Axx,Axy,Axz;
  delete[] Ayy,Ayz,Azz; 

  delete[] rhob,Pb;  
  delete[] vbx,vby,vbz;  

  return 0;
}

