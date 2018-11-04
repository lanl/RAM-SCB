//============================================================================
//    Copyright (c) 2016, Los Alamos National Security, LLC
//    All rights reserved.
//============================================================================

#include <omp.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>

/*===============================================================================================*/
/*Initialize all needed structures and error handlers*/
void ramgsl_initialization_c(int a)
   {
/*    Custom error messages will go here when I can figure them out -ME
      void debug_handler(const char * reason, const char * file, int line, int gsl_errno)
         {
            if (gsl_errno == GSL_EDOM) {

            } else if (gsl_errno == GSL_EDIVERGE) {

            } else if (gsl_errno == GSL_EMAXITER) {

            } else {

            }

            return;
         }
      a = 0 turns off all error messages
      a = 1 turns on custom error messages for debugging
      if (a == 0) (
         gsl_set_error_handler_off();
      } else {
         gsl_error_handler_t *new_handler
             = debug_handler(reason, file, line, gsl_errno)
      }
*/
      a = 0;
      gsl_set_error_handler_off();
      return;
   }

/*===============================================================================================*/
/* Interpolate in 1D and return interpolated values and derivatives */
void interpolation_1d_c(int n, int i1, int i2,
                        double *xa, double *fa, 
                        double *xb, double *fb, int *status)
   {
      int i, err;
      const gsl_interp_type *t;

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();

      /* Spline List
         n = 0 => Linear Spline
         n = 1 => Cubic Spline
         n = 2 => Akima Spline
         n = 3 => Steffen Spline
      */
      
      if (n == 0) {
         t = gsl_interp_linear;
      }else if (n == 1) {  
         if (fa[0] == fa[i1-1] ) {
             t = gsl_interp_cspline_periodic;
         }else{
             t = gsl_interp_cspline;
         }
      }else if (n == 2) {
         if (fa[0] == fa[i1-1] ) {
             t = gsl_interp_akima_periodic;
         }else{
             t = gsl_interp_akima;
         }
      }else if (n == 3) {
         t = gsl_interp_steffen;
      }else{
         t = gsl_interp_linear;
      }
      gsl_spline *spline = gsl_spline_alloc (t, i1);
      err = gsl_spline_init (spline, xa, fa, i1);
      if (err) { 
         //printf("interp error: Array not in ascending order ");
         //for ( i=1; i<i1; i++ ) {
         //    printf("%f ",xa[i]);
         //} 
         //printf("\n");
         status[0] = err;
         return;
      }

      if (xb[0] <= xa[0]) { xb[0] = xa[0] + 1e-10; }
      if (xb[i2-1] >= xa[i1-1]) { xb[i2-1] = xa[i1-1] - 1e-10; }
      for ( i = 0; i<i2; i++ ) {
         // We need to add some sort of extrapolation. For now do linear extrapolation -ME
         if (xb[i] <= xa[0]) {
            fb[i] = fa[0] + (xb[i]-xa[0])/(xa[1]-xa[0])*(fa[1]-fa[0]);
            err = GSL_EDOM;
         }else if (xb[i] >= xa[i1-1]) {
            fb[i] = fa[i-1] + (xb[i]-xa[i1-1])/(xa[i1-2]-xa[i1-1])*(fa[i1-2]-fa[i1-1]);
            err = GSL_EDOM;
         }else{
            err = gsl_spline_eval_e(spline, xb[i], acc, &fb[i]);
            }

//         if (err) {
//            if (err == GSL_EDOM) {
//               printf("interpolation_1d_c warning, GSL_EDOM: x is outside range of xa, x=%f, xa[0]=%f, xa[-1]=%f\n",xb[i],xa[0],xa[i1-1]);
//               printf("  Using extrapolated value, y=%f, ya[0]=%f, ya[-1]=%f\n",fb[i],fa[0],fa[i1-1]);
//            }else{
//               printf("interpolation_1d_c failed: Generic failure\n");
//            }
//         }

      }
      status[0] = err;

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);
      return;
   }

/*===============================================================================================*/
/* Interpolate in 1D smoothing spline and return N values and derivatives */
void interpolation_smooth_c(int i1, int i2,
                            double *xa, double *fa,
                            double *xb, double *fb, int *status)
   {
      int err;
      const int n = i1;
      const int ncoeffs = 24;
      const int nbreak = ncoeffs-2;
      int i, j;
      gsl_bspline_workspace *bw;
      gsl_vector *B;
      gsl_vector *c;
      gsl_vector *y;
      gsl_matrix *X, *cov;
      gsl_multifit_linear_workspace *mw;
      double chisq, yerr, Bj;

      bw = gsl_bspline_alloc(4, nbreak);
      X = gsl_matrix_alloc(n, ncoeffs);
      B = gsl_vector_alloc(ncoeffs);
      c = gsl_vector_alloc(ncoeffs);
      y = gsl_vector_alloc(n);
      cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
      mw = gsl_multifit_linear_alloc(n, ncoeffs);

      gsl_bspline_knots_uniform(xa[0],xa[n-1],bw);
      for (i = 0; i<n; i++) {
          gsl_bspline_eval(xa[i],B,bw);
          gsl_vector_set(y,i,fa[i]);
          for (j =0; j<ncoeffs; j++) {
              Bj = gsl_vector_get(B,j);
              gsl_matrix_set(X, i, j, Bj);
          }
      }

      err = gsl_multifit_linear(X, y, c, cov, &chisq, mw);

      for (i = 0; i<i2; i++) {
          xb[i] = xa[0] + (i+1)*(xa[n-1]-xa[0])/(i2+2);
          gsl_bspline_eval(xb[i],B,bw);
          gsl_multifit_linear_est(B, c, cov, &fb[i], &yerr);
      }
      status[0] = err;

      gsl_bspline_free(bw);
      gsl_vector_free(B);
      gsl_matrix_free(X);
      gsl_vector_free(c);
      gsl_vector_free(y);
      gsl_matrix_free(cov);
      gsl_multifit_linear_free(mw);
      return;
   }

/*===============================================================================================*/
/* Interpolate in 2D and return interpolated values and derivatives */
void interpolation_2d_c(int i1, int j1, int i2, int j2, 
                        double *xa, double *ya, double *fa,
                        double *xb, double *yb, double *fb, int *status)
   {
      int i;
      int j;
      double *za = malloc(i1 * j1 * sizeof(double));

      gsl_interp_accel *xacc = gsl_interp_accel_alloc();
      gsl_interp_accel *yacc = gsl_interp_accel_alloc();
      gsl_interp2d *interp
        = gsl_interp2d_alloc (gsl_interp2d_bilinear, i1, j1);
      /*
        For BiCubic spline, use gsl_interp2d_bicubic
        For BiLinear spline, use gsl_interp2d_bilinear
      */
      for ( j = 0; j<j1; j++ ) {
         for ( i = 0; i<i1; i++ ) {
            gsl_interp2d_set(interp, za, i, j, fa[j*i1+i]);
         }
      }
      gsl_interp2d_init (interp, xa, ya, za, i1, j1);

      for ( j = 0; j<j2; j++ ) {
         for ( i = 0; i<i2; i++ ) {
            fb[j*i2+i] = gsl_interp2d_eval_extrap(interp, xa, ya, za, xb[j*i2+i], yb[j*i2+i], xacc, yacc);
         }
      }
      status[0] = 0;

      gsl_interp2d_free(interp);
      gsl_interp_accel_free(xacc);
      gsl_interp_accel_free(yacc);
      free(za);
 
      return;
   }

/*===============================================================================================*/
/* Interpolate in 3D and return interpolated values and derivatives */
//void interpolation_3d_c(int i1, int j1, int k1, int i2, int j2, int k2,
//                        double *xa, double *ya, double *za, double *fa,
//                        double *xb, double *yb, double *zb, double *fb,
//                        double *dx, double *dy, double *dz, int *status)
//   {
//      return;
//   }
//
/*===============================================================================================*/
/* Return 1D derivatives */
void interpolation_derivs_c(int i1, double *xa, double *fa, double *dx, int *status)
   {
      int i, err;
      const gsl_interp_type *t;

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      t = gsl_interp_steffen;
      gsl_spline *spline = gsl_spline_alloc (t, i1);
      err = gsl_spline_init (spline, xa, fa, i1);
      if (err) { printf("deriv error %f %f %f %f %f %f\n",xa[0],xa[1],xa[2],xa[i1-3],xa[i1-2],xa[i1-1]); }

      for ( i = 0; i<i1; i++ ) {
         err = gsl_spline_eval_deriv_e(spline,xa[i],acc,&dx[i]);
         if (dx[i]==0.0) { dx[i] = dx[i] + 1e-31; }
         if (err) {
            if (err == GSL_EDOM) {
               printf("interpolation_derivs_c warning, GSL_EDOM: x is outside range of xa, x=%f, xa[0]=%f, xa[-1]=%f\n",xa[i],xa[0],xa[i1-1]);
            } else {
               printf("interpolation_derivs_c failed: Generic failure\n");
            }
         }
      }
      status[0] = err;

      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);

      return;
   }

/*===============================================================================================*/
/* Interpolation for Integration */
void integral_interpolate(int i1,
                          double *xa, double *fa,
                          double *xb, double *fb, int *status)
   {
      int err;
      const gsl_interp_type *t;

      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      
      t = gsl_interp_linear;
      gsl_spline *spline = gsl_spline_alloc (t, i1);
      gsl_spline_init (spline, xa, fa, i1);

      err = gsl_spline_eval_e(spline, xb[0], acc, &fb[0]);
      if (err) {
         if (err == GSL_EDOM) {
            printf("integral_interpolate warning, GSL_EDOM: x is outside range of xa, x=%f, xa[0]=%f, xa[-1]=%f\n",xb[0],xa[0],xa[i1-1]);
         } else {
            printf("integral_interpolate failed: Generic failure\n");
         }
      }
      status[0] = err;

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);
      
      return;
   }

/*===============================================================================================*/
/*Evaluate I integral between points a and b*/
void integrate_i_c(int i1, double a, double b, double *cVal, double *bf,
                   double mirror, double *result, double *error, int *status)
   {
      struct f_params { int i; double *theta; double *field; double bmirror; };
      struct f_params params = { i1, cVal, bf, mirror};

      double f_i(double thetaLocal, void * p)
         {
            struct f_params * params = (struct f_params *)p;
            int err[1];
            double bLocal;
      
            integral_interpolate(params->i, params->theta, params->field, &thetaLocal, &bLocal, err);
            if (bLocal >= params->bmirror ) {
               return 0.0;
            } else {
               return sqrt(params->bmirror - bLocal);
            }
         }
      unsigned long int nevals = 0;
      gsl_integration_cquad_workspace *w
         = gsl_integration_cquad_workspace_alloc(100);

      gsl_function F;
      F.function = &f_i;
      F.params = &params;

      int err = gsl_integration_cquad(&F, a, b, 1e-3, 1e-3, w, &result[0], &error[0], &nevals);

      gsl_integration_cquad_workspace_free(w);

/*
      gsl_integration_workspace *w
         = gsl_integration_workspace_alloc(10000);

      gsl_function F;
      F.function = &f_i;
      F.params = &params;

      int err = gsl_integration_qag(&F, a, b, 1e-4, 1e-3, 3, 10000, w, &result[0], &error[0]);

      gsl_integration_workspace_free(w);
*/
      if (err) {
         if (err == GSL_EDIVERGE) {
            printf("integrate_i_c warning, GSL_EDIVERGE: Integral is divergent, or too slowly convergent to be integrated numerically\n");
         } else if (err == GSL_EMAXITER) {
            printf("integrate_i_c warning, GSL_EMAXITER: The maximum number of subdivisions was exceeded\n");
         } else if (err == GSL_EROUND) {
            printf("integrate_i_c warning, GSL_EROUND: Cannot reach tolerance because of roundoff error\n");
         } else if (err == GSL_ESING) {
            printf("integrate_i_c warning, GSL_ESING: A non-integrable singularity or other bad integrand behavior was found in the interval\n");
         } else {
            printf("integrate_i_c warning, GSL_EDOM: Generic failure warning\n");
         }
      }
      status[0] = err;

      return;
   }

/*===============================================================================================*/
/*Evaluate h integral between points a and b*/
void integrate_h_c(int i1, double a, double b, double *cVal, double *bf,
                   double mirror, double *result, double *error, int *status)
   {
      struct f_params { int i; double *theta; double *field; double bmirror; };
      struct f_params params = { i1, cVal, bf, mirror};

      double f_h(double thetaLocal, void * p)
         {
            struct f_params * params = (struct f_params *)p;
            int err[1];
            double bLocal;
      
            integral_interpolate(params->i, params->theta, params->field, &thetaLocal, &bLocal, err);
      
            if (bLocal >= params->bmirror ) {
               return 0.0;
            } else {
               return sqrt(1.0/(params->bmirror - bLocal));
            }
         }
      unsigned long int nevals = 0;
      gsl_integration_cquad_workspace *w
         = gsl_integration_cquad_workspace_alloc(100);

      gsl_function F;
      F.function = &f_h;
      F.params = &params;

      int err = gsl_integration_cquad(&F, a, b, 1e-3, 1e-3, w, &result[0], &error[0], &nevals);

      gsl_integration_cquad_workspace_free(w);
/*
      gsl_integration_workspace *w
         = gsl_integration_workspace_alloc(10000);

      gsl_function F;
      F.function = &f_h;
      F.params = &params;

      int err = gsl_integration_qags(&F, a, b, 1e-4, 1e-3, 10000, w, &result[0], &error[0]);

      gsl_integration_workspace_free(w);
*/
      if (err) {
         if (err == GSL_EDIVERGE) {
            printf("integrate_h_c warning, GSL_EDIVERGE: Integral is divergent, or too slowly convergent to be integrated numerically\n");
         } else if (err == GSL_EMAXITER) {
            printf("integrate_h_c warning, GSL_EMAXITER: The maximum number of subdivisions was exceeded\n");
         } else if (err == GSL_EROUND) {
            printf("integrate_h_c warning, GSL_EROUND: Cannot reach tolerance because of roundoff error\n");
         } else if (err == GSL_ESING) {
            printf("integrate_h_c warning, GSL_ESING: A non-integrable singularity or other bad integrand behavior was found in the interval\n");
         } else {
            printf("integrate_h_c warning, GSL_EDOM: Generic failure warning\n");
         }
      }
      status[0] = err;

      return;
   }

/*===============================================================================================*/
/*Evaluate hDens integral between points a and b*/
void integrate_hdens_c(int i1, double a, double b, double *cVal, double *bf,
                       double mirror, double *dist, double *result, double *error, int *status)
   {
      struct f_params { int i; double *theta; double *field; double bmirror; double *distance; };
      struct f_params params = { i1, cVal, bf, mirror, dist};

      double f_hdens(double thetaLocal, void * p)
         {
            struct f_params * params = (struct f_params *)p;
            int err[1];
            double bLocal, radius;

            double rairden(double r)
               {
                  double rad_local = fmin(r,6.4);
                  double rairden = 1e-5 * pow(10,13.326 
                                                  - 3.6908*rad_local 
                                                  + 1.1362*pow(rad_local,2) 
                                                  - 0.16984*pow(rad_local,3) 
                                                  + 0.009552*pow(rad_local,4));
            
                  return rairden;
               }
      
            integral_interpolate(params->i, params->theta, params->field, &thetaLocal, &bLocal, err);
            integral_interpolate(params->i, params->theta, params->distance, &thetaLocal, &radius, err);

            if (bLocal >= params->bmirror) {
               return 0.0;
            } else {
               return rairden(radius)*sqrt(1.0/(params->bmirror - bLocal));
            }
         }
      unsigned long int nevals = 0;
      gsl_integration_cquad_workspace *w
         = gsl_integration_cquad_workspace_alloc(100);

      gsl_function F;
      F.function = &f_hdens;
      F.params = &params;

      int err = gsl_integration_cquad(&F, a, b, 1e-3, 1e-3, w, &result[0], &error[0], &nevals);

      gsl_integration_cquad_workspace_free(w);
/*
      gsl_integration_workspace *w
         = gsl_integration_workspace_alloc(10000);

      gsl_function F;
      F.function = &f_hdens;
      F.params = &params;

      int err = gsl_integration_qags(&F, a, b, 1e-4, 1e-3, 10000, w, &result[0], &error[0]);
      gsl_integration_workspace_free(w);
*/
      if (err) {
         if (err == GSL_EDIVERGE) {
            printf("integrate_hdens_c warning, GSL_EDIVERGE: Integral is divergent, or too slowly convergent to be integrated numerically\n");
         } else if (err == GSL_EMAXITER) {
            printf("integrate_hdens_c warning, GSL_EMAXITER: The maximum number of subdivisions was exceeded\n");
         } else if (err == GSL_EROUND) {
            printf("integrate_hdens_c warning, GSL_EROUND: Cannot reach tolerance because of roundoff error\n");
         } else if (err == GSL_ESING) {
            printf("integrate_hdens_c warning, GSL_ESING: A non-integrable singularity or other bad integrand behavior was found in the interval\n");
         } else {
            printf("integrate_hdens_c warning, GSL_EDOM: Generic failure warning\n");
         }
      }
      status[0] = err;

      return;
   }

/*===============================================================================================*/
void integrator_c(int nT, int nPa, double *mirror, double *cVal, double *bf,
                  double *dist, double *yI, double *yH, double *yD)
    {
        int err, i, L;
        int LH[nPa], RH[nPa];
        double a[nPa], b[nPa];
        double error;

        // Compute mirror points without relying on a strictly increasing B field
        for (L=1; L<nPa-1; L++) {
            a[L] = 0;
            b[L] = 0;
            LH[L] = 0;
            RH[L] = 0;
            if (mirror[L] < bf[(nT+1)/2]) {
                a[L] = cVal[(nT+1)/2];
                b[L] = cVal[(nT+1)/2];
                LH[L] = (nT+1)/2;
                RH[L] = (nT+1)/2;
                continue;
            }
            for (i=1; i<(nT+1)/2+2; i++) {
                if (mirror[L] <= bf[i-1] && mirror[L] >= bf[i]) {
                    a[L] = cVal[i-1];
                    LH[L] = i-1;
                    break;
                }
            }
            for (i=nT-2; i>(nT+1)/2-2; i--) {
                if (mirror[L] >= bf[i-1] && mirror[L] <= bf[i]) {
                    b[L] = cVal[i];
                    RH[L] = i;
                    break;
                }
            }
        }

        // Next integrate over the whole array to get a 'baseline'
        integrate_i_c(nT, cVal[0], cVal[nT-1], cVal, bf, mirror[nPa-1], &yI[nPa-1], &error, &err);
        integrate_h_c(nT, cVal[0], cVal[nT-1], cVal, bf, mirror[nPa-1], &yH[nPa-1], &error, &err);
        integrate_hdens_c(nT, cVal[0], cVal[nT-1], cVal, bf, mirror[nPa-1], dist, &yD[nPa-1], &error, &err);

        // Next loop over the pitch angles (nPa) starting with nPa
        for (L=nPa-2; L>0; L--) {
            // Need to check if it would mirror outside our domain
            if (mirror[L] >= bf[1]    || a[L] == 0) { a[L] = cVal[0]; }
            if (mirror[L] >= bf[nT-1] || b[L] == 0) { b[L] = cVal[nT-1]; }
            if (a[L] <= cVal[0] || b[L] >= cVal[nT-1]) {
                mirror[L] = mirror[nPa-1];
                yI[L] = yI[nPa-1];
                yH[L] = yH[nPa-1];
                yD[L] = yD[nPa-1];
            }else{
                if ((RH[L]-LH[L])<=4) {
                    yI[L] = yI[L+1];
                    yH[L] = yH[L+1];
                    yD[L] = yD[L+1];
                    continue;
                }
                // Finally we can compute the integrals
                integrate_i_c(nT, a[L], b[L], cVal, bf, mirror[L], &yI[L], &error, &err);
                if (err) { yI[L] = -1; }
                integrate_h_c(nT, a[L], b[L], cVal, bf, mirror[L], &yH[L], &error, &err);
                if (err) { yH[L] = -1; }
                integrate_hdens_c(nT, a[L], b[L], cVal, bf, mirror[L], dist, &yD[L], &error, &err);
                if (err) { yD[L] = -1; }

                // If there was trouble getting the integral, use the previous L's integral
                if (yI[L] <= 0) { yI[L] = yI[L+1]; }
                if (yH[L] <= 0) { yH[L] = yH[L+1]; }
                if (yD[L] <= 0) { yD[L] = yD[L+1]; }
            }
        }
        // The first L (90 degree pitch angle) isn't calculated
        yI[0] = 0;
        yH[0] = yH[1];
        yD[0] = yD[1];

        return;
    }
