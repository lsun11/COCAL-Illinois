
/*
 *  Methods of class Data_bhns
 *
 *    (see file bhns_data.h for documentation).
 *
 */

/*
 *   Copyright (c) 2008 Keisuke Taniguchi
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

char bhns_data_aux_C[] = "$Header: /cvsroot/Lorene/Devel/template.C,v 1.4 2003/10/19 20:01:10 e_gourgoulhon Exp $" ;

/*
 * $Id: template.C,v 1.4 2003/10/19 20:01:10 e_gourgoulhon Exp $
 * $Log: template.C,v $
 *
 * $Header: /cvsroot/Lorene/Devel/template.C,v 1.4 2003/10/19 20:01:10 e_gourgoulhon Exp $
 *
 */

// C++ headers
//#include <>

// C headers
#include <math.h>
#include "nr.h"
#define NR_END 1
#define FREE_ARG char*

// New Polynomial extrapolation junk definition
#include "poly_extrap_junk_settings.h"

// Lorene headers
#include "bhns_data.h"
#include "scalar.h"
#include "bin_bhns.h"
#include "eos.h"
#include "unites.h"

               //--------------------------------------//
               //     Constructor from LORENE data     //
               //--------------------------------------//

Data_bhns::Data_bhns(int nbpoints, const double* xi, const double* yi,
		     const double* zi, const char* filename,
		     const double *dx,const double *dy,const double *dz,
		     const double *xmin,const double *ymin,const double *zmin,
		     const int *nx,const int *ny,const int *nz)
    : np(nbpoints) {

    using namespace Unites ;

    // Reading of data
    // ---------------

    FILE* fich = fopen(filename, "r") ;
    if (fich == 0x0) {
        cout << "Problem in opening the file " << filename
	     << " ! " << endl ;
	perror(" reason") ;
	abort() ;
    }

    int mer ;
    fread(&mer, sizeof(int), 1, fich) ; // mer

    Mg3d mg_bh(fich) ;
    Map_af mp_bh(mg_bh, fich) ;

    Mg3d mg_ns(fich) ;
    Map_et mp_ns(mg_ns, fich) ;
    Eos* peos = Eos::eos_from_file(fich) ;

    Bin_bhns bhns(mp_bh, mp_ns, *peos, fich) ;

    fclose(fich) ;

    //NEXT BLOCK OF LINES: TESTING ONLY
    /*
    // Unit of length
    const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( peos ) ;
    kappa_poly = p_eos_poly->get_kap() ;
    gamma_poly = p_eos_poly->get_gam() ;
    double kap_ns2 = pow( kappa_poly,  0.5 /(gamma_poly-1) ) ;
    // Polytropic unit of length in terms of r_unit :
    r_poly = kap_ns2 / sqrt(ggrav) ;
    xx_bh = mp_bh.get_ori_x() / r_poly ;
    yy_bh = 0. ;
    printf("HI. %e\t%e\t%e\n",kappa_poly,xx_bh,yy_bh);
    */

    // Construction of the binary system
    // ---------------------------------

    // Updation of the metric
    (bhns.set_bh()).update_metric_bhns(bhns.get_ns(), bhns.get_bh(), 1.) ;
    (bhns.set_ns()).update_metric_bhns(bhns.get_bh(), bhns.get_ns(), 1.) ;
    (bhns.set_bh()).update_met_der_comp_bhns(bhns.get_ns()) ;
    (bhns.set_ns()).update_met_der_comp_bhns(bhns.get_bh()) ;
    (bhns.set_bh()).extr_curv_bhns(bhns.get_omega(), bhns.get_x_rot(),
				   bhns.get_y_rot()) ;
    (bhns.set_ns()).extr_curv_bhns() ;

    // Updation of the hydro quantities
    (bhns.set_ns()).equation_of_state() ;
    (bhns.set_ns()).kinema_bhns((bhns.get_bh()).is_kerrschild(),
				(bhns.get_bh()).get_mass_bh(),
				bhns.get_separ(), bhns.get_omega(),
				bhns.get_x_rot(), bhns.get_y_rot()) ;
    (bhns.set_ns()).fait_d_psi_bhns() ;
    (bhns.set_ns()).hydro_euler_bhns((bhns.get_bh()).is_kerrschild(),
				     (bhns.get_bh()).get_mass_bh(),
				     bhns.get_separ()) ;

    // Initialisation of member data
    // -----------------------------

    // Unit of length
    const Eos_poly* p_eos_poly = dynamic_cast<const Eos_poly*>( peos ) ;

    double t_poly = 1. ;
    double m_poly = 1. ;
    double j_poly = 1. ;

    if (p_eos_poly != 0x0) {

        kappa_poly = p_eos_poly->get_kap() ;
	gamma_poly = p_eos_poly->get_gam() ;
	double kap_ns2 = pow( kappa_poly,  0.5 /(gamma_poly-1) ) ;

	// Polytropic unit of length in terms of r_unit :
	r_poly = kap_ns2 / sqrt(ggrav) ;

	// Polytropic unit of time in terms of t_unit :
	t_poly = r_poly ;

	// Polytropic unit of mass in terms of m_unit :
	m_poly = r_poly / ggrav ;

	// Polytropic unit of angular momentum in terms of j_unit :
	j_poly = r_poly * r_poly / ggrav ;
    }

    mass_adm = bhns.mass_adm_bhns_vol() / m_poly ;

    double separ = bhns.get_separ() ;
    double y_ns = (bhns.get_ns()).get_mp().get_ori_y() ;
    dist_coord = sqrt(separ*separ + y_ns*y_ns) / r_poly ;

    omega_orb = bhns.get_omega() * t_poly ;

    omega_spin = (bhns.get_bh()).get_omega_spin() * t_poly ;

    angu_mom = bhns.angu_mom_bhns()(2) / j_poly ;

    // Parameters for the Killing vector
    Tbl xi_i(3) ;
    xi_i.set_etat_qcq() ;
    xi_i.set(0) = 0. ; // xi_t0 : xi_bar{theta}
    xi_i.set(1) = 1. ; // xi_p0 : xi_bar{phi}
    xi_i.set(2) = 0. ; // xi_l0 : L

    double phi_i = 0. ; // pp_0
    double theta_i = 0.5 * M_PI ;
    int nrk_phi = 10 ;
    int nrk_theta = 10 ;

    spin_angu_mom =
      (bhns.get_bh()).spin_am_bhns(xi_i, phi_i, theta_i, nrk_phi, nrk_theta)
      / j_poly ;

    mass_bh_irr = (bhns.get_bh()).mass_irr_bhns() / m_poly ;

    rad_ah = (bhns.get_bh()).rad_ah() / r_poly ;

    mass_ns_b = (bhns.get_ns()).mass_b() / m_poly ;

    nbar_c = pow( kappa_poly,  1./(gamma_poly-1) ) *
      (bhns.get_ns()).get_nbar().val_grid_point(0,0,0,0) ;

    xx_bh = mp_bh.get_ori_x() / r_poly ;
    yy_bh = 0. ;

    xx_ns = mp_ns.get_ori_x() / r_poly ;
    yy_ns = mp_ns.get_ori_y() / r_poly ;

    rad_ns_x_comp = (bhns.get_ns()).ray_eq_pi() / r_poly ;
    rad_ns_y_tow = (bhns.get_ns()).ray_eq_pis2() / r_poly ;
    rad_ns_x_opp = (bhns.get_ns()).ray_eq() / r_poly ;
    rad_ns_y_opp = (bhns.get_ns()).ray_eq_3pis2() / r_poly ;
    rad_ns_z = (bhns.get_ns()).ray_pole() / r_poly ;

    cout.precision(16) ;
    cout << "Black hole-neutron star binary :" << endl ;
    cout << "--------------------------------" << endl ;
    cout << "  Adiabatic index            : " << gamma_poly << endl ;
    cout << "  Polytropic constant        : " << kappa_poly
	 << " [rho_nuc c^2 / n_nuc^gamma]" << endl ;
    cout << "  Length unit R_poly := kappa^{1/2(gamma-1)} : "
	 << r_poly / km << " km" << endl ;
    cout << "  ADM mass of the system     : " << mass_adm
	 << " R_poly" << endl ;
    cout << "  Coordinate separation d    : " << dist_coord
	 << " R_poly" << endl ;
    cout << "  Orbital angular velocity   : " << omega_orb
	 << " / R_poly" << endl ;
    cout << "  Spin angular velocity      : " << omega_spin
	 << " / R_poly" << endl ;
    cout << "  Total angular momentum     : " << angu_mom
	 << " R_poly^2" << endl ;
    cout << "  Spin angular momentum     : " << spin_angu_mom
	 << " R_poly^2" << endl ;
    cout << "  Irreducible mass of the BH : " << mass_bh_irr
	 << " R_poly" << endl ;
    cout << "  Radius of the AH           : " << rad_ah
	 << " R_poly" << endl ;
    cout << "  Baryon rest mass of the NS : " << mass_ns_b
	 << " R_poly" << endl ;
    cout << "  Cent.bary.dens. of the NS  : " << nbar_c
	 << " R_poly^{-2}" << endl ;
    cout << "  Position of the BH (X)     : " << xx_bh
	 << " R_poly" << endl ;
    cout << "  Position of the BH (Y)     : " << yy_bh
	 << " R_poly" << endl ;
    cout << "  Position of the NS (X)     : " << xx_ns
	 << " R_poly" << endl ;
    cout << "  Position of the NS (Y)     : " << yy_ns
	 << " R_poly" << endl ;
    cout << "  Radius of the NS (x_comp)  : " << rad_ns_x_comp
	 << " R_poly" << endl ;
    cout << "  Radius of the NS (x_opp)   : " << rad_ns_x_opp
	 << " R_poly" << endl ;
    cout << "  Radius of the NS (y_tow)   : " << rad_ns_y_tow
	 << " R_poly" << endl ;
    cout << "  Radius of the NS (y_opp)   : " << rad_ns_y_opp
	 << " R_poly" << endl ;
    cout << "  Radius of the NS (z)       : " << rad_ns_z
	 << " R_poly" << endl ;

    // Creation of the various arrays on the Cartesian grid
    // ----------------------------------------------------

    alloc_memory() ;

    // Initialisation of the Cartesian grid
    // ------------------------------------

    for (int i=0; i<np; i++) {
        xx[i] = xi[i] ;
    }
    for (int i=0; i<np; i++) {
        yy[i] = yi[i] ;
    }
    for (int i=0; i<np; i++) {
        zz[i] = zi[i] ;
    }

    // Computation of the values at the points of the Cartesian grid
    // -------------------------------------------------------------

    // (Lapse*Conformal) function
    // --------------------------
    const Scalar& lapconf_bh = (bhns.get_bh()).get_lapconf_auto() ;
    const Scalar& lapconf_ns = (bhns.get_ns()).get_lapconf_auto() ;

    Valeur vlapconf_bh = lapconf_bh.get_spectral_va() ;
    Valeur vlapconf_ns = lapconf_ns.get_spectral_va() ;
    vlapconf_bh.coef() ;
    vlapconf_ns.coef() ;

    // Shift vector
    // ------------
    const Vector& shift_bh = (bhns.get_bh()).get_shift_auto() ;
    const Vector& shift_ns = (bhns.get_ns()).get_shift_auto() ;

    Valeur vshift_bhx = shift_bh(1).get_spectral_va() ;
    Valeur vshift_bhy = shift_bh(2).get_spectral_va() ;
    Valeur vshift_bhz = shift_bh(3).get_spectral_va() ;

    Valeur vshift_nsx = shift_ns(1).get_spectral_va() ;
    Valeur vshift_nsy = shift_ns(2).get_spectral_va() ;
    Valeur vshift_nsz = shift_ns(3).get_spectral_va() ;

    vshift_bhx.coef() ;
    vshift_bhy.coef() ;
    vshift_bhz.coef() ;
    vshift_nsx.coef() ;
    vshift_nsy.coef() ;
    vshift_nsz.coef() ;

    // Conformal factor
    // ----------------
    const Scalar& confo_bh = (bhns.get_bh()).get_confo_auto() ;
    const Scalar& confo_ns = (bhns.get_ns()).get_confo_auto() ;

    Valeur vconfo_bh = confo_bh.get_spectral_va() ;
    Valeur vconfo_ns = confo_ns.get_spectral_va() ;

    vconfo_bh.coef() ;
    vconfo_ns.coef() ;

    // Extrinsic curvature : K_{ij}
    // ----------------------------

    // BH part
    const Tensor& d_shift_bh = (bhns.get_bh()).get_d_shift_auto() ;
    Tensor dshift_bh = d_shift_bh ;
    dshift_bh.std_spectral_base() ;
    for (int i=1; i<=3; i++) {
      for (int j=1; j<=3; j++) {
	dshift_bh.set(i,j).dec_dzpuis(2) ;
      }
    }

    Valeur vdshift_bh_xx = dshift_bh(1,1).get_spectral_va() ;
    Valeur vdshift_bh_xy = dshift_bh(1,2).get_spectral_va() ;
    Valeur vdshift_bh_xz = dshift_bh(1,3).get_spectral_va() ;
    Valeur vdshift_bh_yx = dshift_bh(2,1).get_spectral_va() ;
    Valeur vdshift_bh_yy = dshift_bh(2,2).get_spectral_va() ;
    Valeur vdshift_bh_yz = dshift_bh(2,3).get_spectral_va() ;
    Valeur vdshift_bh_zx = dshift_bh(3,1).get_spectral_va() ;
    Valeur vdshift_bh_zy = dshift_bh(3,2).get_spectral_va() ;
    Valeur vdshift_bh_zz = dshift_bh(3,3).get_spectral_va() ;

    vdshift_bh_xx.set_base(((bhns.get_bh()).get_d_shift_auto()(1,1)).get_spectral_va().base) ;
    vdshift_bh_xy.set_base(((bhns.get_bh()).get_d_shift_auto()(1,2)).get_spectral_va().base) ;
    vdshift_bh_xz.set_base(((bhns.get_bh()).get_d_shift_auto()(1,3)).get_spectral_va().base) ;
    vdshift_bh_yx.set_base(((bhns.get_bh()).get_d_shift_auto()(2,1)).get_spectral_va().base) ;
    vdshift_bh_yy.set_base(((bhns.get_bh()).get_d_shift_auto()(2,2)).get_spectral_va().base) ;
    vdshift_bh_yz.set_base(((bhns.get_bh()).get_d_shift_auto()(2,3)).get_spectral_va().base) ;
    vdshift_bh_zx.set_base(((bhns.get_bh()).get_d_shift_auto()(3,1)).get_spectral_va().base) ;
    vdshift_bh_zy.set_base(((bhns.get_bh()).get_d_shift_auto()(3,2)).get_spectral_va().base) ;
    vdshift_bh_zz.set_base(((bhns.get_bh()).get_d_shift_auto()(3,3)).get_spectral_va().base) ;

    vdshift_bh_xx.coef() ;
    vdshift_bh_xy.coef() ;
    vdshift_bh_xz.coef() ;
    vdshift_bh_yx.coef() ;
    vdshift_bh_yy.coef() ;
    vdshift_bh_yz.coef() ;
    vdshift_bh_zx.coef() ;
    vdshift_bh_zy.coef() ;
    vdshift_bh_zz.coef() ;

    // NS part
    const Tensor& d_shift_ns = (bhns.get_ns()).get_d_shift_auto() ;
    Tensor dshift_ns = d_shift_ns ;
    dshift_ns.std_spectral_base() ;
    for (int i=1; i<=3; i++) {
      for (int j=1; j<=3; j++) {
	dshift_ns.set(i,j).dec_dzpuis(2) ;
      }
    }

    Valeur vdshift_ns_xx = dshift_ns(1,1).get_spectral_va() ;
    Valeur vdshift_ns_xy = dshift_ns(1,2).get_spectral_va() ;
    Valeur vdshift_ns_xz = dshift_ns(1,3).get_spectral_va() ;
    Valeur vdshift_ns_yx = dshift_ns(2,1).get_spectral_va() ;
    Valeur vdshift_ns_yy = dshift_ns(2,2).get_spectral_va() ;
    Valeur vdshift_ns_yz = dshift_ns(2,3).get_spectral_va() ;
    Valeur vdshift_ns_zx = dshift_ns(3,1).get_spectral_va() ;
    Valeur vdshift_ns_zy = dshift_ns(3,2).get_spectral_va() ;
    Valeur vdshift_ns_zz = dshift_ns(3,3).get_spectral_va() ;

    vdshift_ns_xx.set_base(((bhns.get_ns()).get_d_shift_auto()(1,1)).get_spectral_va().base) ;
    vdshift_ns_xy.set_base(((bhns.get_ns()).get_d_shift_auto()(1,2)).get_spectral_va().base) ;
    vdshift_ns_xz.set_base(((bhns.get_ns()).get_d_shift_auto()(1,3)).get_spectral_va().base) ;
    vdshift_ns_yx.set_base(((bhns.get_ns()).get_d_shift_auto()(2,1)).get_spectral_va().base) ;
    vdshift_ns_yy.set_base(((bhns.get_ns()).get_d_shift_auto()(2,2)).get_spectral_va().base) ;
    vdshift_ns_yz.set_base(((bhns.get_ns()).get_d_shift_auto()(2,3)).get_spectral_va().base) ;
    vdshift_ns_zx.set_base(((bhns.get_ns()).get_d_shift_auto()(3,1)).get_spectral_va().base) ;
    vdshift_ns_zy.set_base(((bhns.get_ns()).get_d_shift_auto()(3,2)).get_spectral_va().base) ;
    vdshift_ns_zz.set_base(((bhns.get_ns()).get_d_shift_auto()(3,3)).get_spectral_va().base) ;

    vdshift_ns_xx.coef() ;
    vdshift_ns_xy.coef() ;
    vdshift_ns_xz.coef() ;
    vdshift_ns_yx.coef() ;
    vdshift_ns_yy.coef() ;
    vdshift_ns_yz.coef() ;
    vdshift_ns_zx.coef() ;
    vdshift_ns_zy.coef() ;
    vdshift_ns_zz.coef() ;

    // Baryon rest-mass density
    // ------------------------
    const Scalar& nbar_ns = (bhns.get_ns()).get_nbar() ;

    Valeur vnbar = nbar_ns.get_spectral_va() ;
    vnbar.coef() ;

    // Spatial part of the fluid 4-velocity (covariant: u_i)
    // -----------------------------------------------------
    // u_i = gamma_euler * u_euler

    const Vector& ueuler = (bhns.get_ns()).get_u_euler() ;
    const Scalar& gameu = (bhns.get_ns()).get_gam_euler() ;
    const Scalar& psi4_tot = (bhns.get_ns()).get_psi4() ;

    Vector u4 = gameu * psi4_tot * ueuler ;
    u4.std_spectral_base() ;
    u4.annule( (bhns.get_ns()).get_nzet(), mg_ns.get_nzone()-1 ) ;

    Valeur vu4x = u4(1).get_spectral_va() ;
    Valeur vu4y = u4(2).get_spectral_va() ;
    Valeur vu4z = u4(3).get_spectral_va() ;
    vu4x.coef() ;
    vu4y.coef() ;
    vu4z.coef() ;

    double norm = pow( kappa_poly, 1./(gamma_poly-1.) ) ;

    for (int i=0; i<np; i++) {

        double x0 = xx[i] * r_poly ;
	double y0 = yy[i] * r_poly ;
	double z0 = zz[i] * r_poly ;

	// Values of (ll_bh, xi_bh, theta_bh, phi_bh) (grid BH)
	// corresponding to (x,y,z):
	// ----------------------------------------------------
	double r_bh, theta_bh, phi_bh ;  // polar coordinates centered on BH
	mp_bh.convert_absolute(x0, y0, z0, r_bh, theta_bh, phi_bh) ;

	int ll_bh ;     // domain index
	double xi_bh ;  // radial coordinate xi in [0,1] or [-1,1]
	mp_bh.val_lx(r_bh, theta_bh, phi_bh, ll_bh, xi_bh) ;

	// Values of (ll_ns, xi_ns, theta_ns, phi_ns) (grid NS)
	// corresponding to (x,y,z):
	// ----------------------------------------------------
	double r_ns, theta_ns, phi_ns ;  // polar coordinates centered on NS
	mp_ns.convert_absolute(x0, y0, z0, r_ns, theta_ns, phi_ns) ;

	int ll_ns ;     // domain index
	double xi_ns ;  // radial coordinate xi in [0,1] or [-1,1]
	mp_ns.val_lx(r_ns, theta_ns, phi_ns, ll_ns, xi_ns) ;


	// Conformal factor
	// ----------------

	psi_conf[i] =
	  vconfo_bh.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	  + vconfo_ns.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns) ;

	// Lapse function
	// --------------

	lapse[i] =
	  (vlapconf_bh.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	   + vlapconf_ns.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns))
	  / psi_conf[i] ;

	// Shift vector
	// ------------

	shift_x[i] =
	  vshift_bhx.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	  + vshift_nsx.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns) ;

	shift_y[i] =
	  vshift_bhy.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	  + vshift_nsy.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns) ;

	shift_z[i] =
	  vshift_bhz.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	  + vshift_nsz.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns) ;

	// Extrinsic curvature : K_{ij}
	// ----------------------------

	kij_xx[i] = r_poly * (pow(psi_conf[i], 4.) / lapse[i] / 3.)
	  *(2.*vdshift_bh_xx.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    +2.*vdshift_ns_xx.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)
	    -vdshift_bh_yy.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    -vdshift_ns_yy.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)
	    -vdshift_bh_zz.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    -vdshift_ns_zz.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)) ;

	kij_xy[i] = r_poly * 0.5 * (pow(psi_conf[i], 4.) / lapse[i])
	  *(vdshift_bh_xy.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    +vdshift_ns_xy.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)
	    +vdshift_bh_yx.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    +vdshift_ns_yx.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)) ;

	kij_xz[i] = r_poly * 0.5 * (pow(psi_conf[i], 4.) / lapse[i])
	  *(vdshift_bh_xz.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    +vdshift_ns_xz.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)
	    +vdshift_bh_zx.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    +vdshift_ns_zx.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)) ;

	kij_yy[i] = r_poly * (pow(psi_conf[i], 4.) / lapse[i] / 3.)
	  *(2.*vdshift_bh_yy.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    +2.*vdshift_ns_yy.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)
	    -vdshift_bh_xx.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    -vdshift_ns_xx.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)
	    -vdshift_bh_zz.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    -vdshift_ns_zz.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)) ;

	kij_yz[i] = r_poly * 0.5 * (pow(psi_conf[i], 4.) / lapse[i])
	  *(vdshift_bh_yz.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    +vdshift_ns_yz.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)
	    +vdshift_bh_zy.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    +vdshift_ns_zy.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)) ;

	kij_zz[i] = r_poly * (pow(psi_conf[i], 4.) / lapse[i] / 3.)
	  *(2.*vdshift_bh_zz.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    +2.*vdshift_ns_zz.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)
	    -vdshift_bh_xx.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    -vdshift_ns_xx.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)
	    -vdshift_bh_yy.c_cf->val_point(ll_bh, xi_bh, theta_bh, phi_bh)
	    -vdshift_ns_yy.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns)) ;

	// Baryon rest-mass density
	// ------------------------

	nbar[i] =
	  norm * vnbar.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns) ;

	// Spatial part of the fluid 4-velocity (covariant: u_i)
	// -----------------------------------------------------

	u4_x[i] = vu4x.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns) ;

	u4_y[i] = vu4y.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns) ;

	u4_z[i] = vu4z.c_cf->val_point(ll_ns, xi_ns, theta_ns, phi_ns) ;

    }ACA // End of loop on the points

    // Finally, fill in excised BH regions:
    //double dr=0.1;
    //    int EXTRAP_POLY_NUM_POINTS=8;
    double r_AH = rad_ah;
    int n=0;
    int num_points_inside_ah = 0;
    double *r_pts = (double *)malloc(sizeof(double)*EXTRAP_POLY_NUM_POINTS);
    for (int k=0; k<*nz; k++)
      for(int j=0; j<*ny; j++) 
	for(int i=0; i<*nx`; i++) { 
	  //double zz1 = *zmin + *dz * k ;
	  //double yy1 = *ymin + *dy * j ;
	  //double xx1 = *xmin + *dx * i ;
          int index = i + j*(*nx) + k*(*nx)*(*ny);
	  double zz1 = zi[index];
	  double yy1 = yi[index];
	  double xx1 = xi[index];

	  //printf("blah %e %e %e %e\n",xx[n],yy[n],zz[n],psi_conf[n]);
	  if(sqrt((xx1-xx_bh)*(xx1-xx_bh) + (yy1-yy_bh)*(yy1-yy_bh) + zz1*zz1) < r_AH) {
	    //Center coordinates at BH.  Note that zz_bh==yy_bh==0 always!
	    xx1 -= xx_bh;
	    double rr=sqrt(xx1*xx1+yy1*yy1+zz1*zz1);

	    //double theta=atan2(yy1,xx1);
	    //double phi=atan2(sqrt(xx1*xx1+yy1*yy1),zz1);

	    for(int ri=1;ri<=EXTRAP_POLY_NUM_POINTS;ri++) {
	      r_pts[ri] = r_AH + POLY_EXTRAP_dr*((double)ri-1.0+1e-8);
	    }
	    
	    //double r1=r_AH + 1e-8*POLY_EXTRAP_dr;
	    //double rr2=r_AH + 1*POLY_EXTRAP_dr;
	    
	    //double r_pts[2];
	    
	    //r_pts[0]=rr1;
	    //r_pts[1]=rr2;
	    
	    double dydummy=0;
	    polint_statpunc(r_pts, &lapse[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &lapse[n], &dydummy);
	    polint_statpunc(r_pts, &shift_x[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &shift_x[n], &dydummy);
	    polint_statpunc(r_pts, &shift_y[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &shift_y[n], &dydummy);
	    polint_statpunc(r_pts, &shift_z[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &shift_z[n], &dydummy);
	    
	    polint_statpunc(r_pts, &psi_conf[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &psi_conf[n], &dydummy);

	    polint_statpunc(r_pts, &kij_xx[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &kij_xx[n], &dydummy);
	    polint_statpunc(r_pts, &kij_xy[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &kij_xy[n], &dydummy);
	    polint_statpunc(r_pts, &kij_xz[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &kij_xz[n], &dydummy);
	    polint_statpunc(r_pts, &kij_yy[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &kij_yy[n], &dydummy);
	    polint_statpunc(r_pts, &kij_yz[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &kij_yz[n], &dydummy);
	    polint_statpunc(r_pts, &kij_zz[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS-1], EXTRAP_POLY_NUM_POINTS, rr, &kij_zz[n], &dydummy);

	    /*
	    lapse[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&lapse[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.3,rr);
	    
	    shift_x[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&shift_x[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.0,rr);
	    double by0=0.25;
	    if(yy<0) by0*=-1;
	    shift_y[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&shift_y[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],by0,rr);
	    shift_z[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&shift_z[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.0,rr);
	    
	    psi_conf[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&psi_conf[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],2.0,rr);
	    
	    kij_xx[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&kij_xx[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],-6.0,rr);
	    kij_xy[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&kij_xy[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.0,rr);
	    kij_xz[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&kij_xz[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.0,rr);
	    kij_yy[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&kij_yy[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],3.0,rr);
	    kij_yz[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&kij_yz[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.0,rr);
	    kij_zz[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&kij_zz[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],3.0,rr);
	    */	    

	    /*  Don't need to overwrite matter variables inside ah, since these are zero anyway.
	      nbar[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&nbar[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.0,rr);
	      
	      u4_x[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&u4_x[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.0,rr);
	      u4_y[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&u4_y[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.0,rr);
	      u4_z[n] = fourth_order_smooth_extrap_fitpoint(r_pts,&u4_z[(*nx)*(*ny)*(*nz)+num_points_inside_ah*EXTRAP_POLY_NUM_POINTS],0.0,rr);
	    */	    

	    num_points_inside_ah++;
	  }
	  /*
	  printf("lap %e %e %e %e\n",xx[n],yy[n],zz[n],lapse[n]);

	  printf("sx %e %e %e %e\n",xx[n],yy[n],zz[n],shift_x[n]);
	  printf("sy %e %e %e %e\n",xx[n],yy[n],zz[n],shift_y[n]);
	  printf("sz %e %e %e %e\n",xx[n],yy[n],zz[n],shift_z[n]);

	  printf("kxx %e %e %e %e\n",xx[n],yy[n],zz[n],kij_xx[n]);
	  printf("kxy %e %e %e %e\n",xx[n],yy[n],zz[n],kij_xy[n]);
	  printf("kxz %e %e %e %e\n",xx[n],yy[n],zz[n],kij_xz[n]);
	  printf("kyy %e %e %e %e\n",xx[n],yy[n],zz[n],kij_yy[n]);
	  printf("kyz %e %e %e %e\n",xx[n],yy[n],zz[n],kij_yz[n]);
	  printf("kzz %e %e %e %e\n",xx[n],yy[n],zz[n],kij_zz[n]);
	  */
	  n++;
	}
}

