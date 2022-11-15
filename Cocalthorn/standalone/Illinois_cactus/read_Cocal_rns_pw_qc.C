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

  char path[400],comment[400];
  char *dir_path, *file_path, *EoS_path;

  fparam.getline(comment, 400);
  fparam.getline(comment, 400);
  fparam.getline(path,    400);
  fparam.close();

  printf("path is %s\n",path);


  dir_path=path;
  size_t path_len=strlen(dir_path);


  //***********************************************************************
  //***********************************************************************
  //***********************************************************************

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
  double rhoini_cgs, fac2  ;
  double sgma, cbar,eneqc,preqc,constqc  ;
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

  rhocgs = new (nothrow) double [nphase+2];  // Interface densities
  abi    = new (nothrow) double [nphase+2];  // Gammas
  rhoi   = new (nothrow) double [nphase+2];  // Interface densities at G=c=Msun=1
  abccgs = new (nothrow) double [nphase+2];  // Kappas in cgs
  abc    = new (nothrow) double [nphase+2];  // Kappas in G=c=Msun=1
  qi     = new (nothrow) double [nphase+2];  // Kappas in G=c=Msun=1
  hi     = new (nothrow) double [nphase+2];  // Kappas in G=c=Msun=1
  tedi   = new (nothrow) double [nphase+2];  // Kappas in G=c=Msun=1

  // Get core density and sigma : P=sigma*(epsilon-epsilon0)-P0
  pdat >> rhocgs[nphase+1] >> sgma ;
  pdat.get(buf,100);

  // Get all pieces: interface restmass density, gamma
  for(int j=nphase; j>=0; j--)
  { 
    pdat >> rhocgs[j] >> abi[j] ;
    pdat.get(buf,100);
  }
//  cout << " Number of phases :" << nphase << endl;
//  for(int j=nphase; j>=0; j--)  cout << setw(20) << rhocgs[j] << setw(20) << abi[j] << endl;

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


  cbar           = hi[nphase]/pow(rhoi[nphase],sgma) ; 

  eneqc          = (hi[nphase] - qi[nphase])*rhoi[nphase] ;  // h = ene/rho + q

  preqc          = qi[nphase]*rhoi[nphase] ;

  constqc        = preqc - sgma*eneqc ;

  rhoi[nphase+1] = facrho*rhocgs[nphase+1] ;

  hi[nphase+1]   = cbar*pow(rhoi[nphase+1],sgma) ;

  qi[nphase+1]   = (sgma*cbar*pow(rhoi[nphase+1],sgma) + constqc/rhoi[nphase+1])/(sgma+1.0) ;

  tedi[nphase+1] = rhoi[nphase+1]*(hi[nphase+1] - qi[nphase+1]) ;

  abi[nphase+1]  = sgma ;

  abc[nphase+1]  = cbar ;

  abccgs[nphase+1] = 0.0 ;

  cout << "=============================================================================================" << endl;
  cout << " Number of phases :" << nphase << endl;
  cout << setw(18) << "Kappa: [cgs]" << setw(18) << "[G=c=Msun=1]" << setw(18) << "Gamma" << setw(18) << "Density: [cgs]" << setw(18) << "[G=c=Msun=1]" << endl;
  for(int j=0; j<=nphase; j++)  
    cout << scientific << setprecision(10) << setw(18) << abccgs[j] << setw(18) << abc[j] << setw(18) <<
            abi[j] << setw(18) << rhocgs[j] << setw(18) << rhoi[j] << endl;

  cout   << scientific << setprecision(10) << setw(18) << abccgs[nphase+1] << setw(18) << abc[nphase+1] << setw(18) <<
            abi[nphase+1] << setw(18) << rhocgs[nphase+1] << setw(18) << rhoi[nphase+1] << endl;
  cout << "=============================================================================================" << endl;
  cout << "Quark EOS: P=sigma*(epsilon-epsilon0)-P0     with sigma = " << scientific << setprecision(10) << setw(18) << sgma << endl;
  cout << "=============================================================================================" << endl;
  cout << setw(34) << "Total energy density: [G=c=Msun=1]" << setw(25) << "Enthalpy [G=c=Msun=1]" << endl;

  for(int j=0; j<=nphase+1; j++)  
    cout << scientific << setprecision(10) << setw(34) << tedi[j] << setw(25) << hi[j] << endl;
  cout << "=============================================================================================" << endl;


  // fist and second lines
//  pdat >> temp1 >> temp2;
//  pdat.get(buf,100);

//  pdat >> rho0 >> Press0;
//  pdat.get(buf,100);

  // next lines

//  pdat >> rho3 >> gamma3;
//  pdat.get(buf,100);

//  pdat >> rho2 >> gamma2;
//  pdat.get(buf,100);

//  pdat >> rho1 >> gamma1;
//  pdat.get(buf,100);

  // From Vasilis' awk script:
  /* The gamma's in Antonis peos_parameter.dat file so far are such that                                                                                                                                          
     if rho[1] <= rho < rho[2], gamma=gam[2], k=kap[2]                                                                                                                                       
     if rho[2] <= rho < rho[3], gamma=gam[3], k=kap[3] 
 
     Antonis said that the last gamma (gam[1]) in peos_parameter.dat                                                                                                                                         
     is not used but here I am assuming that rho[1] can also be a transition                                                                                                                       
     density. So, far the actual transition density is actually only rho[2], but                                                                                                                
     I will treat rho[1] here also as a transition density perhaps for                                                                                                                        
     future applications 
     strictly speaking there is no reason why the EOS should be capped by                                                                                                                            
     the densities rho[1] and rho[3], so I will treat the EOS as follows                                                                                                                                        
     if     rho < rho[1], P=kap[1]*rho^gam[1] 
     if     rho[1] <= rho < rho[2], P=kap[2]*rho^gam[2]
     if     rho[2] <= rho         , P=kap[3]*rho^gam[3] */
 
//  if(rho0<rho1){
//    kappa1 = Press0/pow(rho0,gamma1);

    // use continuity to set the other Kappas
//    kappa2 = kappa1*pow(rho1,gamma1-gamma2);
//    kappa3 = kappa2*pow(rho2,gamma2-gamma3);

//    cout << "condition rho0<rho1"  << endl;

//  }

//  if(rho0<rho2 && rho0 >= rho1){
//    kappa2 = Press0/pow(rho0,gamma2);

    // use continuity to set the other Kappas
//    kappa1 = kappa2*pow(rho1,gamma2-gamma1);
//    kappa3 = kappa2*pow(rho2,gamma2-gamma3);

//    cout << "condition rho0<rho2"  << endl;
//  }


//  if(rho0>=rho2){
//    kappa3 = Press0/pow(rho0,gamma3);

    // use continuity to set the other Kappas
//    kappa2 = kappa3*pow(rho2,gamma3-gamma2);
//    kappa1 = kappa2*pow(rho1,gamma2-gamma1);
//  }



  // Finally convert Kappa's into geometrized units where Msun =1
  // The conversion for K's is given by                                                              
  //  k_geo=kcgs*(c^2/G)^(1/n)*(1/c^2)

//  double n1 = 1.0/(gamma1 - 1);
//  double n2 = 1.0/(gamma2 - 1);
//  double n3 = 1.0/(gamma3 - 1);


  // c^2/G in cgs 
  double c2oG = 1.34663531e+28;
  //c^2 in cgs
  double c2   = 8.98755179e+20;
  //Msun in cm             
  double Msun  =1.47708885e+5;


  //  kappa1 = kappa1;//;*pow(c2oG,1.0/n1) / ( c2* pow(Msun,(2.0/n1)));
//  kappa1 = kappa1*pow(c2oG,1.0/n1)/(c2*pow(Msun,2.0/n1));
//  kappa2 = kappa2*pow(c2oG,1.0/n2)/(c2*pow(Msun,2.0/n2));
//  kappa3 = kappa3*pow(c2oG,1.0/n3)/(c2*pow(Msun,2.0/n3));

  // Now convert the transition densities to geometrized units in which Msun=1                   
  // rho_geo = rho_cgs*G/c^2 and then Multiply by Msun^2                                             
                                                                                                    
//  rho1  =  rho1*Msun*Msun/c2oG;
//  rho2  =  rho2*Msun*Msun/c2oG;

//  cout << "rho1   =" << rho1   << endl;
//  cout << "rho2   =" << rho2   << endl;
//  cout << "kappa1 =" << kappa1 << endl; 
//  cout << "kappa2 =" << kappa2 << endl; 
//  cout << "kappa3 =" << kappa3 << endl; 
//  cout << "gamma1 =" << gamma1 << endl; 
//  cout << "gamma2 =" << gamma2 << endl; 
//  cout << "gamma3 =" << gamma3 << endl; 


  // ************************
  // ***  Interpolate ID  ***
  // ************************

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
  outfilepsi.write((char *) &(nphase),sizeof(int));
  for(int j=0; j<=nphase+1; j++)
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

  delete[] Axx,Axy,Axz;
  delete[] Ayy,Ayz,Azz; 

  delete[] rhob,Pb;  
  delete[] vbx,vby,vbz;  

  delete[] rhocgs, abi, rhoi, abccgs, abc, qi, hi, tedi ;

  return 0;
}

