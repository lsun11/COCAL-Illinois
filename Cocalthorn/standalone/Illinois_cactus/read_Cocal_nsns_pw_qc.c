/*
 * Reads binary neutron star initial data, exporting Lorene structures
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


// Select intial data type
extern "C" void read_id_type_(char *file_path, char *id_type, size_t *path_len_file, int *ierr);


// Irrotational binary NSNS
extern "C" void coc2pri_ir_(char *dir_path, double Xb[], double Yb[], double Zb[],
			    int *Nx, int *Ny, int *Nz,
			    double *gtt, double *gtx, double *gty,double *gtz, 
			    double *gxx, double *gxy,double *gxz,
			    double *gyy, double *gyz,double *gzz,
			    double *Axx,double *Axy, double *Axz,
			    double *Ayy,double *Ayz, double *Azz,
			    double *Pb,double *rhob,
			    double *vbx ,  double *vby, double *vbz,size_t *path_len);


// Corotational binary NSNS
extern "C" void coc2pri_co_(char *dir_path, double Xb[], double Yb[], double Zb[],
			    int *Nx, int *Ny, int *Nz,
			    double *gtt, double *gtx, double *gty,double *gtz, 
			    double *gxx, double *gxy,double *gxz,
			    double *gyy, double *gyz,double *gzz,
			    double *Axx,double *Axy, double *Axz,
			    double *Ayy,double *Ayz, double *Azz,
			    double *Pb,double *rhob,
			    double *vbx ,  double *vby, double *vbz,size_t *path_len);
// Spinning binary NSNS
extern "C" void coc2pri_sp_(char *dir_path, double Xb[], double Yb[], double Zb[],
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

  char path[400],pathII[400],comment[400],id_type[3];
  char *dir_path, *file_path;

  fparam.getline(comment, 400);
  fparam.getline(comment, 400);
  fparam.getline(path, 400);
  fparam.close();

  printf("path is %s\n",path);

  dir_path=path;

  size_t path_len=strlen(dir_path);
  
  // read the parameter of the EOS
  //--------------------------

  file_path=strcat(path,"/peos_parameter_mpt1.dat");  
  size_t path_len_file=strlen(file_path);
  ifstream pdat(file_path);
  if(!pdat.is_open())  
  {
    cout << " Path is :" << file_path << endl;
    cout << " **** Couldn't open the file peos_parameter_mpt1.dat ****" << endl;
    exit(0);
  }

  char buf[100];
  int    nphase, ii, iphase;
  double rho_0, pre_0, facrho, facpre, det;
  double rhoini_cgs, fac2;
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


  //--------------------------

  //file_path=strcat(path,"/rnspar_mpt1.dat");
  //size_t path_len_file=strlen(file_path);


  ifstream fparamII("read_bin_path.par");
  fparamII.getline(comment, 400);
  fparamII.getline(comment, 400);
  fparamII.getline(pathII, 400);
  fparamII.close();

  printf("pathII is %s\n",pathII);

  dir_path=pathII;
  path_len=strlen(dir_path);


  file_path=strcat(pathII,"/rnspar_mpt1.dat");  
  path_len_file=strlen(file_path);

  //  ifstream pdat(file_path);
  if(!pdat.is_open())  
  {
    cout << " Path is :" << file_path << endl;
    cout << " **** Couldn't open the file peos_parameter_mpt1.dat ****" << endl;
    exit(0);
  }



  // Read ID type
  int ierr=1;
  read_id_type_(file_path,id_type,&path_len_file,&ierr);
  id_type[2]='\0'; // NULL terminate the string
  
  if(ierr==0){
    printf("COCAL_ID:: success reading rnspar_mpt1.dat.\n");
  }
  else{
    printf("COCAL_ID:: Can't open the rnspar_mpt1.dat.");
    exit(1);
  }


  // Interpolate ID
  if(strcmp(id_type,"CO")==0){
      printf("COCAL_ID:: Reading corotating BNS ID\n");
      
      coc2pri_co_(dir_path, Xb, Yb, Zb,
		  &Nx, &Ny, &Nz,
		  gtt, gtx, gty,gtz,
		  gxx, gxy,gxz,
		  gyy, gyz,gzz,
		  Axx, Axy, Axz,
		  Ayy, Ayz, Azz,
		  Pb,rhob,
		  vbx, vby, vbz, &path_len);}
    
    else if(strcmp(id_type,"IR")==0){
      printf("COCAL_ID:: Reading irrotational BNS ID\n");
      
      coc2pri_ir_(dir_path, Xb, Yb, Zb,
		  &Nx, &Ny, &Nz,
		  gtt, gtx, gty,gtz,
		  gxx, gxy,gxz,
		  gyy, gyz,gzz,
		  Axx, Axy, Axz,
		  Ayy, Ayz, Azz,
		  Pb,rhob,
		  vbx, vby, vbz, &path_len);}
    else if(strcmp(id_type,"SP")==0){
      printf("COCAL_ID:: Reading spinning BNS ID\n");
      
      coc2pri_sp_(dir_path, Xb, Yb, Zb,
		  &Nx, &Ny, &Nz,
		  gtt, gtx, gty,gtz,
		  gxx, gxy,gxz,
		  gyy, gyz,gzz,
		  Axx, Axy, Axz,
		  Ayy, Ayz, Azz,
		  Pb,rhob,
		  vbx, vby, vbz, &path_len);}
  
  
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

