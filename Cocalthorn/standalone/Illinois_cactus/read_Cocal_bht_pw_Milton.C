/*
 * Reads COCAL black hole-torus initial data.
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

int NR_END=1;
//int EXTRAP_POLY_NUM_POINTS=8;
//double POLY_EXTRAP_dr=0.1;
char FREE_ARG; 
//keep the below definitions
//int extra_points=8;

using namespace std;


void polint_statpunc(double xa[], double ya[], int n, double x, double *y, double *dy);
double *vector_statpunc(long nl, long nh);
//void nrerror_statpunc(char error_text[]);
//void free_vector_statpunc(double *v, long nl, long nh);

// function declaration 


// Interporlater  Black hole-torus
extern "C" void coc2pri_bht_wl_(char *dir_path, double Xb[], double Yb[], double Zb[],
		                int *Nx, int *Ny, int *Nz,
			        double *gtt, double *gtx, double *gty,double *gtz, 
			        double *gxx, double *gxy,double *gxz,
			        double *gyy, double *gyz,double *gzz,
			        double *Kxx,double *Kxy, double *Kxz,
			        double *Kyy,double *Kyz, double *Kzz,
			        double *Pb,double *rhob,
			        double *vbx ,  double *vby, double *vbz,size_t *path_len);


 // error meseage
int main(int argc, char** argv){

  if (argv[1]==NULL || argv[10]==NULL) {
    printf("Sorry.  I was expecting 10 command line arguments:\n");
    printf("Usage: ./read_bin_ns <xmin> <ymin> <zmin> <dx> <dy> <dz> <nx> <ny> <nz> <processor number>\n");
    exit(1);
  }


  //  Reads Cartesian grid parameters                                                                                                                                                                
        
  int Nx, Ny, Nz;
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


  // screen output
  cout << "---------------- " << endl ;
  cout << "Cartesian grid : " << endl ;
  cout << "---------------- " << endl ;
  cout << "   Number of points in the x direction : " << Nx     << endl ;
  cout << "   Number of points in the y direction : " << Ny     << endl ;
  cout << "   Number of points in the z direction : " << Nz     << endl ;
  cout << "   x_min                               : " << x_min  << endl ;
  cout << "   y_min                               : " << y_min  << endl ;
  cout << "   z_min                               : " << z_min  << endl ;
  cout << "   x: grid width                       : " << dx     << endl ;
  cout << "   y: grid width                       : " << dy     << endl ;
  cout << "   z: grid width                       : " << dz     << endl ;
  cout << "   procnumer                           : " << procnumber << endl ;


  // get the AH info (i.e., excision radius)
  // For a general case I will read in the 
  //  position of the AH centoid

  ifstream AHdat("read_bin_AH.par");
  printf("Getting the AH info\n");
  
  double xx_bh;
  double yy_bh;
  double zz_bh;
  double r_AH;
  int    extra_points;
  double extra_points_dr;

  AHdat >> xx_bh;
  AHdat >> yy_bh;
  AHdat >> zz_bh;
  AHdat >> r_AH;
  AHdat >> extra_points;
  AHdat >> extra_points_dr;

  cout << "xx_bh = " << xx_bh << endl;
  cout << "yy_bh = " << yy_bh << endl;
  cout << "zz_bh = " << zz_bh << endl;
  cout << "r_AH =  " << r_AH << endl;
  cout << "extra_points = " << extra_points << endl;
  cout << "extra_points_dr = " << extra_points_dr << endl;

  AHdat.close();


  // Define the total number of points. Notice that the junk filled data
  // requires points inside the horizon 

  int ntotal = Nx*Ny*Nz; // Total number of points

  int tot_points_inside_AH=0;
  for (int k=0; k<Nz; k++){
    double zz = z_min + dz*k;
    for (int j=0; j<Ny; j++){
      double yy = y_min + dy*j;
      for (int i=0; i<Nx; i++) {
	double xx = x_min + dx*i;

	// add grid points inside the horizon for the junk filled data  
	if(sqrt((xx - xx_bh)*(xx - xx_bh) + (yy - yy_bh)*(yy - yy_bh) + 
		(zz - zz_bh)*(zz - zz_bh)) < r_AH){
	  ntotal = ntotal + extra_points;
	  tot_points_inside_AH +=1;
	}
      }
    }
  }


  // ***********************************
  // ***  Build the cartesian grids  ***
  // ***********************************

  // define 1D cartesian arrays
  printf("Allocating Cartesian coordinates arrays\n");

  //Define size of the grids
  double* const Xb = new double[ntotal]; // Cartesian grid in the x direction
  double* const Yb = new double[ntotal]; // Cartesian grid in the y direction
  double* const Zb = new double[ntotal]; // Cartesian grid in the z direction


  double *lx = Xb;
  double *ly = Yb;
  double *lz = Zb;

    for (int k=0; k<Nz; k++){
      double zz = z_min + dz*k;
      for (int j=0; j<Ny; j++){
	double yy = y_min + dy*j;
	for (int i=0; i<Nx; i++){
	  *lx = x_min + dx*i;
	  *ly = yy ;
	  *lz = zz ;

	  lx++ ;
	  ly++ ;
	  lz++ ;
	}
      }
    }


    // grid inside the horizon

    for (int k=0; k<Nz; k++){
      for (int j=0; j<Ny; j++){
	for (int i=0; i<Nx; i++) {
	  double zz = z_min + dz * k ;
	  double yy = y_min + dy * j ;
	  double xx = x_min + dx * i ;
	  if(sqrt((xx - xx_bh)*(xx - xx_bh) + (yy - yy_bh)*(yy - yy_bh) + 
		  (zz - zz_bh)*(zz - zz_bh)) < r_AH) {
	    
	    //center the coordinates in the BH centroid
	    xx -= xx_bh;
	    yy -= yy_bh;
	    zz -= zz_bh;

	    double theta = atan2(yy,xx);
	    double phi   = atan2(sqrt(xx*xx + yy*yy),zz);

	    for(int ri=1;ri<=extra_points;ri++) {
	      double rr = r_AH + extra_points_dr*((double)ri-1.0+1e-8);

	      //	      cout << "i="<< i << "lx=" << lx  << endl ;
	      //Recenter coordinates to origin:
	      *lx = rr*sin(phi)*cos(theta) + xx_bh;
	      *ly = rr*sin(phi)*sin(theta) + yy_bh;
	      *lz = rr*cos(phi)            + zz_bh;

	      lx++;
	      ly++;
	      lz++;
	    }
	  }
	}
      }
    }

    // done with Cartesian coordiantes
  printf("Done with Cartesian Coordinates, now ...\n");


  // define 1D arrays
  printf("allocating ID arrays for the dynamical fields\n");

  //Define full metric
  double *gtt, *gtx,*gty,*gtz;
  double *gxx,*gxy,*gxz;
  double *gyy,*gyz,*gzz;

  //Define curvature extrinsic
  double *Kxx,*Kxy,*Kxz;
  double *Kyy,*Kyz,*Kzz;

  //Define fluid variables
  double *rhob, *Pb;
  double *vbx,*vby,*vbz;  

  //  We have defined 1D C arrays for the metric and extrinsic curvature
  //  Thus, we need to define a relation for how to translate from indices i,j,k
  //  to the 1D c arrays. 

  // Allocate memory
  gtt = new (nothrow) double [ntotal];  // metric coefficient gtt
  gtx = new (nothrow) double [ntotal];  // metric coefficient gtx
  gty = new (nothrow) double [ntotal];  // metric coefficient gty
  gtz = new (nothrow) double [ntotal];  // metric coefficient gtz
 
  gxx = new (nothrow) double [ntotal];  // metric coefficient gxx
  gxy = new (nothrow) double [ntotal];  // metric coefficient gxy
  gxz = new (nothrow) double [ntotal];  // metric coefficient gxz
  gyy = new (nothrow) double [ntotal];  // metric coefficient gyy
  gyz = new (nothrow) double [ntotal];  // metric coefficient gyz
  gzz = new (nothrow) double [ntotal];  // metric coefficient gzz


  Kxx = new (nothrow) double [ntotal];  // extrinsic curv coefficient Kxx
  Kxy = new (nothrow) double [ntotal];  // extrinsic curv coefficient Kxy
  Kxz = new (nothrow) double [ntotal];  // extrinsic curv coefficient Kxz
  Kyy = new (nothrow) double [ntotal];  // extrinsic curv coefficient Kyy
  Kyz = new (nothrow) double [ntotal];  // extrinsic curv coefficient Kyz
  Kzz = new (nothrow) double [ntotal];  // extrinsic curv coefficient Kzz

  rhob = new (nothrow) double [ntotal];  // baryon density
  Pb   = new (nothrow) double [ntotal];  // pressure

  vbx = new (nothrow) double [ntotal];   // three velocity v^x measured by Euler observers
  vby = new (nothrow) double [ntotal];   // three velocity v^y measured by Euler observers
  vbz = new (nothrow) double [ntotal];   // three velocity v^z measured by Euler observers

  if (gtt == NULL){ cout << "Stopping: memory could not be allocated for a1" << endl;}


  // get the path for the ID
  printf("Getting the path for the ID\n");
  ifstream fparam("read_bin_path.par");

  char path[400],comment[400];
  char *dir_path, *file_path, *EoS_path;

  fparam.getline(comment, 400);
  fparam.getline(comment, 400);
  fparam.getline(path,    400);
  fparam.close();

  printf("path is %s\n",path);


  dir_path=path;
  size_t path_len=strlen(dir_path);



  // *********************************
  // **** Find the EoS paramters  ****
  // *********************************
  
  // Define EoS parameters

  double gamma1,gamma2,gamma3;   // polytropic index n_i
  double kappa1,kappa2,kappa3;   // polytropic constant K_i
  double rho0, Press0, rho1, rho2, rho3;

  // path for EoS parameters 
  file_path=strcat(path,"/peos_parameter.dat");  
  ifstream pdat(file_path);
  if(!pdat.is_open())  
  {
    cout << " Path is :" << file_path << endl;
    cout << " **** Couldn't open the file peos_parameter.dat ****" << endl;
    exit(0);
  }

  char buf[100];
  int    nphase, ii, iphase;
  double rho_0, pre_0, facrho, facpre, det;
  double rhoini_cgs, fac2;
  double *rhocgs, *abi, *rhoi, *abccgs, *abc, *qi, *hi, *tedi ;
  double g_coc = 6.67428e-08 ;
  double c_coc = 2.99792458e+10;
  double msol_coc = 1.98892e+33;

  // Get number of phases, initial central rho0
  pdat >> nphase >> rhoini_cgs ;
  pdat.get(buf,100);

  // Get reference restmass density, pressure
  pdat >> rho_0 >> pre_0 ;
  pdat.get(buf,100);

  rhocgs = new (nothrow) double [nphase+1];  // Interface densities
  abi    = new (nothrow) double [nphase+1];  // Gammas
  rhoi   = new (nothrow) double [nphase+1];  // Interface densities at G=c=Msun=1
  abccgs = new (nothrow) double [nphase+1];  // Kappas in cgs
  abc    = new (nothrow) double [nphase+1];  // Kappas in G=c=Msun=1
  qi     = new (nothrow) double [nphase+1];  // Kappas in G=c=Msun=1
  hi     = new (nothrow) double [nphase+1];  // Kappas in G=c=Msun=1
  tedi   = new (nothrow) double [nphase+1];  // Kappas in G=c=Msun=1

  // Get all pieces: interface restmass density, gamma
  for(int j=nphase; j>=0; j--)
  { 
    pdat >> rhocgs[j] >> abi[j] ;
    pdat.get(buf,100);
  }


  facrho = pow(g_coc/pow(c_coc,2), 3)*pow(msol_coc,2) ;
  facpre = pow(g_coc,3)*pow(msol_coc,2)/pow(c_coc,8)  ;

  for(int j=0; j<=nphase; j++)   rhoi[j] = facrho*rhocgs[j]  ;

  iphase = 1 ;
  for(int ii=1; ii<=nphase; ii++)
  {
    det = (rho_0-rhocgs[ii])*(rho_0-rhocgs[ii-1]) ;
    if (det <= 0.0e0)
    {
      iphase = ii ;
      break ;
    }
  }

  abc[iphase] = pre_0/pow(rho_0,abi[iphase]) ;
  abc[iphase] = facpre/pow(facrho,abi[iphase])*abc[iphase] ;
  abccgs[iphase] = pre_0/(pow(rho_0,abi[iphase])) ;

  if (iphase > 0)
  {
    for(int j=iphase-1; j>=0; j--)
    {
      abc[   j] = pow(rhoi[j],   abi[j+1]-abi[j]) * abc[j+1] ;
      abccgs[j] = pow(rhocgs[j], abi[j+1]-abi[j]) * abccgs[j+1] ;
    }
  }
  if (iphase < nphase)
  {
    for(int j=iphase+1; j<=nphase; j++)
    {
      abc[   j] = pow(rhoi[  j-1], abi[j-1]-abi[j]) * abc[j-1] ;
      abccgs[j] = pow(rhocgs[j-1], abi[j-1]-abi[j]) * abccgs[j-1] ;
    }
  }
 
  for(int j=0; j<=nphase; j++)
  {
    qi[j] = abc[j]*pow(rhoi[j], abi[j]-1.0);
  }

  hi[0] = 1.0;    //  Surface enthalpy.
  for(int j=1; j<=nphase; j++)
  {
    fac2 = abi[j]/(abi[j] - 1.0)  ;
    hi[j] = hi[j-1] + fac2*(qi[j] - qi[j-1]);
  }
  
  for(int j=0; j<=nphase; j++)
  {
    tedi[j] = rhoi[j]*(hi[j] - qi[j]) ;   // total energy density = rho*h - P 
  }

  cout << "=============================================================================================" << endl;
  cout << " Number of phases :" << nphase << endl;
  cout << setw(18) << "Kappa: [cgs]" << setw(18) << "[G=c=Msun=1]" << setw(18) << "Gamma" << setw(18) << "Density: [cgs]" << setw(18) << "[G=c=Msun=1]" << endl;
  for(int j=0; j<=nphase; j++)  
    cout << scientific << setprecision(10) << setw(18) << abccgs[j] << setw(18) << abc[j] << setw(18) <<
            abi[j] << setw(18) << rhocgs[j] << setw(18) << rhoi[j] << endl;
  cout << "=============================================================================================" << endl;
  cout << setw(34) << "Total energy density: [G=c=Msun=1]" << setw(25) << "Enthalpy [G=c=Msun=1]" << endl;
  for(int j=0; j<=nphase; j++)  
    cout << scientific << setprecision(10) << setw(34) << tedi[j] << setw(25) << hi[j] << endl;
  cout << "=============================================================================================" << endl;

  // Finally convert Kappa's into geometrized units where Msun =1
  // The conversion for K's is given by                                                              
  //  k_geo=kcgs*(c^2/G)^(1/n)*(1/c^2)

  // c^2/G in cgs 
  double c2oG = 1.34663531e+28;
  //c^2 in cgs
  double c2   = 8.98755179e+20;
  //Msun in cm             
  double Msun  =1.47708885e+5;



  // ************************
  // ***  Interpolate ID  ***
  // ************************

  // Difine total number of points to take into accoung
  // the extra points

  int NNx =  Nx + tot_points_inside_AH*extra_points;
  int NNy =  Ny + tot_points_inside_AH*extra_points;
  int NNz =  Nz + tot_points_inside_AH*extra_points;

  printf("COCAL_ID:: Reading black hole-torus ID, total # points without AH Nx = %d and with AH, NNx=%d, ntotal= %d\n",Nx, NNx, ntotal);

    coc2pri_bht_wl_(dir_path, Xb, Yb, Zb,
	          &Nx, &Ny, &Nz,
	          gtt, gtx, gty,gtz,
	          gxx, gxy,gxz,
	          gyy, gyz,gzz,
	          Kxx, Kxy, Kxz,
	          Kyy, Kyz, Kzz,
	          Pb,rhob,
	          vbx, vby, vbz, &path_len);

  printf("COCAL_ID:: Done with Cocal, now we need to fill out the excised region");

  // *************************************
  //  ***  Fill the excised BH region  ***
  // *************************************

    int n=0;
    int num_points_inside_ah = 0;
    double *r_pts = (double *)malloc(sizeof(double)*extra_points);
    for (int k=0; k < Nz; k++)
      for(int j=0; j < Ny; j++) 
	for(int i=0; i < Nx; i++) { 
	  double zz1 = z_min + dz * k;
	  double yy1 = y_min + dy * j;
	  double xx1 = x_min + dx * i;

	  if(sqrt((xx1 -xx_bh)*(xx1 - xx_bh) + (yy1 - yy_bh)*(yy1 - yy_bh) + (zz1 - zz_bh)*(zz1 - zz_bh)) < r_AH) {
	    //Center coordinates at BH. 
	    xx1 -= xx_bh;
	    yy1 -= yy_bh;
	    zz1 -= zz_bh;

	    double rr = sqrt(xx1*xx1 + yy1*yy1 + zz1*zz1);

	    
	    for(int ri=1;ri<=extra_points;ri++) {
	      r_pts[ri] = r_AH + extra_points_dr*((double)ri-1.0+1e-8);}

	    double dydummy=0;

	    polint_statpunc(r_pts, &gtt[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gtt[n], &dydummy);
	    polint_statpunc(r_pts, &gtx[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gtx[n], &dydummy);
	    polint_statpunc(r_pts, &gty[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gty[n], &dydummy);
	    polint_statpunc(r_pts, &gtz[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gtz[n], &dydummy);


	    polint_statpunc(r_pts, &gxx[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gxx[n], &dydummy);
	    polint_statpunc(r_pts, &gxy[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gxy[n], &dydummy);
	    polint_statpunc(r_pts, &gxz[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gxz[n], &dydummy);
	    polint_statpunc(r_pts, &gyy[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gyy[n], &dydummy);
	    polint_statpunc(r_pts, &gyz[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gyz[n], &dydummy);
	    polint_statpunc(r_pts, &gzz[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &gzz[n], &dydummy);


	    polint_statpunc(r_pts, &Kxx[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &Kxx[n], &dydummy);
	    polint_statpunc(r_pts, &Kxy[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &Kxy[n], &dydummy);
	    polint_statpunc(r_pts, &Kxz[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &Kxz[n], &dydummy);
	    polint_statpunc(r_pts, &Kyy[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &Kyy[n], &dydummy);
	    polint_statpunc(r_pts, &Kyz[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &Kyz[n], &dydummy);
	    polint_statpunc(r_pts, &Kzz[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &Kzz[n], &dydummy);

	    polint_statpunc(r_pts, &rhob[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &rhob[n], &dydummy);

	    polint_statpunc(r_pts, &vbx[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &vbx[n], &dydummy);
	    polint_statpunc(r_pts, &vby[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &vby[n], &dydummy);
	    polint_statpunc(r_pts, &vbz[Nx*Ny*Nz + num_points_inside_ah*extra_points-1], extra_points, rr, &vbz[n], &dydummy);

	    num_points_inside_ah++;	    

	  }
	  n++;
	}

    //dumping data  

  ofstream outfilepsi;
  char fname_psi[100];

  printf("********************************\n");
  printf("Dumping ID data file = %i\n",procnumber);
  printf("********************************\n");


  sprintf(fname_psi,"CTS_bin-proc%d.d",procnumber);
  outfilepsi.open(fname_psi,ios::out | ios::ate);


  // Grid parameters
  outfilepsi.write((char *) &(ntotal),sizeof(int));
  outfilepsi.write((char *) &(ntotal),sizeof(int));
  outfilepsi.write((char *) &(ntotal),sizeof(int));
  outfilepsi.write((char *) &(ntotal),sizeof(int));


  // 3D Arrays 
  int intdum=0;
  for (int k=0; k<1; k++) {
    for (int j=0; j<Ny; j++) {
      for (int i=0; i<1; i++) {

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

	double Kxx_val = Kxx[intdum];
        double Kxy_val = Kxy[intdum];
        double Kxz_val = Kxz[intdum];
        double Kyy_val = Kyy[intdum];
        double Kyz_val = Kyz[intdum];
        double Kzz_val = Kzz[intdum];

	double rho0_val = rhob[intdum];

        double ux_val   = vbx[intdum];
        double uy_val   = vby[intdum];
        double uz_val   = vbz[intdum];

	cout << Xb[intdum] << " " << gtt_val << " " << gxx_val <<  " " << Kxx_val  << " " << rho0_val <<  " " << ux_val << endl ;


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

        outfilepsi.write((char *) &Kxx_val,sizeof(double));
        outfilepsi.write((char *) &Kxy_val,sizeof(double));
        outfilepsi.write((char *) &Kxz_val,sizeof(double));
        outfilepsi.write((char *) &Kyy_val,sizeof(double));
        outfilepsi.write((char *) &Kyz_val,sizeof(double));
        outfilepsi.write((char *) &Kzz_val,sizeof(double));

        outfilepsi.write((char *) &rho0_val,sizeof(double));

        outfilepsi.write((char *) &ux_val,sizeof(double));
        outfilepsi.write((char *) &uy_val,sizeof(double));
        outfilepsi.write((char *) &uz_val,sizeof(double));

        intdum++;
      }
    }
  }


  // EoS parameters
  outfilepsi.write((char *) &(nphase),sizeof(int));
  for(int j=0; j<=nphase; j++)
  {
    outfilepsi.write((char *) &abc[j] ,sizeof(double));
    outfilepsi.write((char *) &abi[j] ,sizeof(double));
    outfilepsi.write((char *) &rhoi[j],sizeof(double));
    outfilepsi.write((char *) &tedi[j],sizeof(double));
  }

//  outfilepsi.write((char *) &(gamma1),sizeof(double));
//  outfilepsi.write((char *) &(gamma2),sizeof(double));
//  outfilepsi.write((char *) &(gamma3),sizeof(double));
//  outfilepsi.write((char *) &(kappa1),sizeof(double));
//  outfilepsi.write((char *) &(kappa2),sizeof(double));
//  outfilepsi.write((char *) &(kappa3),sizeof(double));
//  outfilepsi.write((char *) &(rho1),  sizeof(double));
//  outfilepsi.write((char *) &(rho2),  sizeof(double));
  outfilepsi.close();


  printf("********************************\n");
  printf("Done with file =%i\n",procnumber);
  printf("********************************\n");


  //release memory
  delete[]  gtt,gtx,gty,gtz;
  delete[]  gxx,gxy,gxz; 
  delete[]  gyy,gyz,gzz; 

  delete[] Kxx,Kxy,Kxz;
  delete[] Kyy,Kyz,Kzz; 

  delete[] rhob,Pb;  
  delete[] vbx,vby,vbz;  

  delete[] rhocgs, abi, rhoi, abccgs, abc, qi, hi, tedi ;

  return 0;
}



void polint_statpunc(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;

  dif=fabs(x-xa[1]);

  c=vector_statpunc(1,n);
  d=vector_statpunc(1,n);
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) {
	fprintf(stderr,"Numerical Recipes run-time error in stationary puncture...\n");
	fprintf(stderr,"Error in routine polint_statpunc");
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
      }

      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  //free_vector_statpunc(d,1,n);
  //free_vector_statpunc(c,1,n);
  delete[] d,c;  
}


double *vector_statpunc(long nl, long nh)
// allocate a double vector with subscript range v[nl..nh] 
{
  double *v;
  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  if (!v) {
	fprintf(stderr,"Numerical Recipes run-time error in stationary puncture...\n");
	fprintf(stderr,"Error in allocation failure in vector_statpunc()");
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
      }
  return v-nl+NR_END;
}


