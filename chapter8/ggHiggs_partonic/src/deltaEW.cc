//
//  Created by Luca Rottoli on 18/09/2015.
//

// #include <vector>
// #include <cstring>
// #include <string>
// #include <fstream>
// #include <iomanip>
// #include <sstream>
// #include <sys/time.h>
// #include <sys/stat.h>
// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_spline.h>


double deltaEW(double mH, double mt)
{
  const int npoints = 60;
  double x[npoints] = {  100, 110, 120, 130, 140, 145, 150, 151, 152, 153, 154,
			 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165,
			 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176,
			 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187,
			 188, 189, 190, 195, 197, 200, 210, 220, 230, 240, 250,
			 260, 270, 280, 290, 300};
  double y1[npoints] = {  4.18358, 4.54598, 4.94167, 5.33852, 5.71055, 5.87131, 5.9808, 5.98933,
			  5.98977, 5.97938, 5.95375, 5.9065, 5.82784, 5.7022, 5.50743, 5.22133,
			  4.86185, 4.52687, 4.2209, 3.87773, 3.526, 3.20017, 2.9111, 2.6577,
			  2.4351, 2.23819, 2.06185, 1.90148, 1.75341, 1.6143, 1.48097, 1.349,
			  1.2158, 1.07441, 0.919019, 0.740333, 0.530749, 0.298401, 0.077024,
			  -0.102261, -0.270062, -0.453871, -0.64347, -0.822932, -0.987444, -1.13869,
			  -1.27542, -1.77099, -1.91007, -2.07378, -2.38711, -2.50679, -2.50779,
			  -2.45153, -2.37917, -2.28018, -2.17313, -2.05227, -1.95334, -1.87019};
  double y2[npoints] = { 4.18772, 4.54989, 4.94532, 5.3419, 5.71357, 5.87408, 5.98325, 5.99169, 
			 5.99204, 5.98155, 5.95579, 5.90839, 5.82955, 5.70369, 5.5086, 5.22206, 
			 4.86193, 4.52607, 4.21914, 3.87513, 3.52271, 3.19631, 2.90677, 2.65299, 
			 2.43007, 2.23286, 2.05622, 1.8956, 1.74731, 1.60796, 1.47449, 1.3423, 
			 1.20892, 1.06731, 0.911812, 0.73288, 0.522871, 0.290366, 0.0684936, 
			 -0.111577, -0.279963, -0.464333, -0.654448, -0.834337, -0.999445, -1.15087, 
			 -1.28768, -1.78469, -1.92451, -2.08908, -2.40298, -2.52461, -2.52539, 
			 -2.47042, -2.39951, -2.29895, -2.19246, -2.0697, -1.97319, -1.8809 };
  //
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  const gsl_interp_type *t1 = gsl_interp_cspline_periodic;
  gsl_spline *spline1 = gsl_spline_alloc (t1, npoints);
  const gsl_interp_type *t2 = gsl_interp_cspline_periodic;
  gsl_spline *spline2 = gsl_spline_alloc (t2, npoints);
  gsl_spline_init (spline1, x, y1, npoints);
  gsl_spline_init (spline2, x, y2, npoints);
  //
  double mt_c = 172.64;
  double mt_p = 174.22;
  //
  double d=0;
  double d_c, d_p;
  //
  if (mH < 100 || mH > 300) {
    cout << "\033[0;31m" << "Warning: Unable to compute EW correction for mH=" << mH <<"GeV. Allowed range is 100GeV<mH<300GeV" << "\033[0m" << endl;
    return __builtin_nan("");
  } else {
    d_c = gsl_spline_eval (spline1, mH, acc);
    d_p = gsl_spline_eval (spline2, mH, acc);
    d = d_c + (mt-mt_c)*(d_p-d_c)/(mt_p-mt_c);
    d = d/100.;
  }
  //
  return d;
}


