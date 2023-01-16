#include "sz_cp_preserve_utils.hpp"
#include "sz_compress_3d.hpp"
#include "sz_compress_cp_preserve_2d.hpp"
#include "sz_def.hpp"
#include "sz_compression_utils.hpp"
#include <complex>
#include "../test/utils.hpp"

// maximal error bound to keep the sign of A*(1 + e_1) + B*(1 + e_2) + C
template<typename T>
static inline double max_eb_to_keep_sign_2d_online(const T A, const T B, const T C=0){
	double fabs_sum = (fabs(A) + fabs(B));
	if(fabs_sum == 0) return 0;
	return fabs(A + B + C) / fabs_sum;
}

// maximal error bound to keep the sign of 
// a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f
// when f > 0
static inline double max_eb_to_keep_sign_2d_online_corners_gt0(const double a, const double b, const double c, const double d, const double e, const double f){
	double eb = 1;
	{
		// [-t, -t]
		// a*t^2 + b*t^2 + c*t^2 - d*t - e*t + f
		double tmp_a = a + b + c;
		double tmp_b = - d - e;
		double tmp_c = f;
		double x = (-tmp_b - sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a);
		if(x > 0) eb = MINF(eb, x);
	}
	{
		// [-t, t]
		double tmp_a = a + b - c;
		double tmp_b = - d + e;
		double tmp_c = f;
		double x = (-tmp_b - sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a);
		if(x > 0) eb = MINF(eb, x);
	}
	{
		// [t, -t]
		double tmp_a = a + b - c;
		double tmp_b = d - e;
		double tmp_c = f;
		double x = (-tmp_b - sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a);
		if(x > 0) eb = MINF(eb, x);
	}
	{
		// [t, t]
		double tmp_a = a + b + c;
		double tmp_b = d + e;
		double tmp_c = f;
		double x = (-tmp_b - sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a);
		if(x > 0) eb = MINF(eb, x);
	}					
	return eb;
}

// maximal error bound to keep the sign of 
// a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f
// when f < 0
static inline 
double max_eb_to_keep_sign_2d_online_corners_lt0(const double a, const double b, const double c, const double d, const double e, const double f, bool verbose=false){
	double eb = 1;
	{
		// [-t, -t]
		// a*t^2 + b*t^2 + c*t^2 - d*t - e*t + f
		double tmp_a = a + b + c;
		double tmp_b = - d - e;
		double tmp_c = f;
		double x = (-tmp_b + sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a);
		if(x > 0) eb = MINF(eb, x);
        if(verbose) printf("[-t, -t]: eb = %.7f, f(x) = %.7f x^2 + %.7f x + %.7f\n", eb, tmp_a, tmp_b, tmp_c);
	}
	{
		// [-t, t]
		double tmp_a = a + b - c;
		double tmp_b = - d + e;
		double tmp_c = f;
		double x = (-tmp_b + sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a);
		if(x > 0) eb = MINF(eb, x);
        if(verbose) printf("[-t, t]: eb = %.7f, f(x) = %.7f x^2 + %.7f x + %.7f\n", eb, tmp_a, tmp_b, tmp_c);
	}
	{
		// [t, -t]
		double tmp_a = a + b - c;
		double tmp_b = d - e;
		double tmp_c = f;
		double x = (-tmp_b + sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a);
		if(x > 0) eb = MINF(eb, x);
        if(verbose) printf("[t, -t]: eb = %.7f, f(x) = %.7f x^2 + %.7f x + %.7f\n", eb, tmp_a, tmp_b, tmp_c);
	}
	{
		// [t, t]
		double tmp_a = a + b + c;
		double tmp_b = d + e;
		double tmp_c = f;
		double x = (-tmp_b + sqrt(tmp_b*tmp_b - 4*tmp_a*tmp_c))/(2*tmp_a);
		if(x > 0) eb = MINF(eb, x);
        if(verbose) printf("[t, t]: eb = %.7f, f(x) = %.7f x^2 + %.7f x + %.7f\n", eb, tmp_a, tmp_b, tmp_c);
	}							
    if(verbose) printf("Returned eb lt0 = %.7f\n", eb);
	return eb;
}

// maximal error bound to keep the sign of 
// a*e1^2 + b*e1 + c
static inline 
double max_eb_to_keep_sign_quadratic(const double a, const double b, const double c){
    double eb = 1;
    if(c > 0){
        double x = (-b - sqrt(b*b - 4*a*c))/(2*a);
        if(x > 0) eb = MINF(eb, x);
    }
    else{
        double x = (-b + sqrt(b*b - 4*a*c))/(2*a);
        if(x > 0) eb = MINF(eb, x);        
    }
    return eb;
}

// maximal error bound to keep the sign of 
// a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f
// if the qudratic form forms a parabolic cylinder (i.e., 4ab - c^2 = 0)
static inline double max_eb_to_keep_sign_2d_online_quadratic_cylinder(const double a, const double b, const double c, const double d, const double e, const double f, bool verbose=false){
    if(verbose){
        printf("%.7f x*x + %.7f y*y + %.7f x*y + %.7f x +%.7f y + %.7f\n", a, b, c, d, e, f);
    }
    double eb = 1;
    // return 0;
    if(f > 0){
        // f > 0
        // Min(a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f) > 0
        // double eb = MINF(eb, max_eb_to_keep_sign_2d_online_corners_gt0(a, b, c, d, e, f));
        eb = MINF(eb, f/(fabs(d) + fabs(e)));
    }
    else{
        // f < 0
        // Max(a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f) < 0
        eb = MINF(eb, max_eb_to_keep_sign_2d_online_corners_lt0(a, b, c, d, e, f, verbose));
    }
    return eb;
}

// maximal error bound to keep the sign of 
// a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f
static inline double max_eb_to_keep_sign_2d_online_general(const double a, const double b, const double c, const double d, const double e, const double f, bool verbose=false){
	// solve cp
	// 2a * e1 + c * e2 + d = 0
	// 2b * e2 + c * e1 + e = 0
	// double e1 = (-2*b*d + c*e) / (4*a*b - c*c);
	// double e2 = (-2*a*e + c*d) / (4*a*b - c*c);
    // Hessian
    // 2a c
    // c  2b
    double det_H = 4*a*b - c*c;
    double eb = max_eb_to_keep_sign_2d_online_quadratic_cylinder(a, b, c, d, e, f);
    if(det_H >= 0) return eb;
    // det_H < 0
    // saddle point: need to check edges
    // if(verbose){
    //     printf("f = %.7f\ne1 = %.7f, e2 = %.7f\nf(e1, e2) = %.7f\n", f, e1, e2, a*e1*e1 + b*e2*e2 + c*e1*e2 + d*e1 + e*e2 + f);
    //     printf("gradient: %.7f, %.7f\n", 2 * a * e1 + c * e2 + d, 2 * b * e2 + c * e1 + e);
    //     printf("%.7f x*x + %.7f y*y + %.7f x*y + %.7f x +%.7f y + %.7f\n", a, b, c, d, e, f);
    // }
	if(f > 0){
      	// f > 0
		// Min(a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f) > 0
        if(b > 0){
            double eb1 = 0, eb2 = 0;
            // x = t
            // check y0 = - (ct + e)/(2b)
            // y0 < -t
            if(2*b - c > 0){
                eb1 = MAX(0, e / (2*b - c));
            }
            else{
                if(e > 0) eb1 = 1;
            }
            // y0 > t
            if(2*b + c > 0){
                eb1 = MAX(eb1, -e / (2*b + c));
            }
            else{
                if(e < 0) eb1 = 1;
            }
            // y0 in [-t, t]
            if(f - (e*e)/(4*b) > 0){
                eb2 = max_eb_to_keep_sign_quadratic(a - (c*c)/(4*b), d - (c*e)/(2*b), f - (e*e)/(4*b));
            }
            eb = MIN(eb, MAX(eb1, eb2));
            // x = -t
            // check y = - (-ct + e)/(2b)
            // y0 < -t
            if(2*b + c > 0){
                eb1 = MAX(0, e / (2*b + c));
            }
            else{
                if(e > 0) eb1 = 1;
            }
            // y0 > t
            if(2*b - c > 0){
                eb1 = MAX(eb1, -e / (2*b - c));
            }
            else{
                if(e < 0) eb1 = 1;
            }
            // y0 in [-t, t]
            if(f - (e*e)/(4*b) > 0){
                eb2 = max_eb_to_keep_sign_quadratic(a - (c*c)/(4*b), - d + (c*e)/(2*b), f - (e*e)/(4*b));
            }
            eb = MIN(eb, MAX(eb1, eb2));
        }
        if(a > 0){
            double eb1 = 0, eb2 = 0;
            // y = t
            // check x0 = - (ct + d)/(2a)
            // x0 < -t
            if(2*a - c > 0){
                eb1 = MAX(0, d / (2*a - c));
            }
            else{
                if(d > 0) eb1 = 1;
            }
            // x0 > t
            if(2*a + c > 0){
                eb1 = MAX(eb1, -d / (2*a + c));
            }
            else{
                if(d < 0) eb1 = 1;
            }
            // x0 in [-t, t]
            if(f - (d*d)/(4*a) > 0){
                eb2 = max_eb_to_keep_sign_quadratic(b - (c*c)/(4*a), e - (c*d)/(2*a), f - (d*d)/(4*a));
            }
            eb = MIN(eb, MAX(eb1, eb2));
            // y = -t
            // check x0 = - (-ct + d)/(2a)
            // x0 < -t
            if(2*a + c > 0){
                eb1 = MAX(0, d / (2*a + c));
            }
            else{
                if(d > 0) eb1 = 1;
            }
            // x0 > t
            if(2*a - c > 0){
                eb1 = MAX(eb1, -d / (2*a + c));
            }
            else{
                if(d < 0) eb1 = 1;
            }
            // x0 in [-t, t]
            if(f - (d*d)/(4*a) > 0){
                eb2 = max_eb_to_keep_sign_quadratic(b - (c*c)/(4*a), -e + (c*d)/(2*a), f - (d*d)/(4*a));
            }
        }
	}
	else{
		// f < 0
		// Max(a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f) < 0
        if(b < 0){
            double eb1 = 0, eb2 = 0;
            // x = -t
            if(- 2*b - c > 0){
                eb1 = MAX(0, e / (- 2*b - c));
            }
            else{
                if(e > 0) eb1 = 1;
            }
            if(- 2*b + c > 0){
                eb1 = MAX(eb1, - e / (- 2*b + c));
            }
            else{
                if(e < 0) eb1 = 1;
            }
            if(f + (e*e)/(-4*b) < 0){
                eb2 = max_eb_to_keep_sign_quadratic(a + (c*c)/(-4*b), d + (c*e)/(-2*b), f + (e*e)/(-4*b));
            }
            eb = MIN(eb, MAX(eb1, eb2));
            // x = t
            if(- 2*b + c > 0){
                eb1 = MAX(0, e / (- 2*b + c));
            }
            else{
                if(e > 0) eb1 = 1;
            }
            if(- 2*b - c > 0){
                eb1 = MAX(eb1, -e / (- 2*b - c));
            }
            else{
                if(e < 0) eb1 = 1;
            }
            if(f + (e*e)/(-4*b) < 0){
                eb2 = max_eb_to_keep_sign_quadratic(a + (c*c)/(-4*b), - d - (c*e)/(-2*b), f + (e*e)/(-4*b));
            }
            eb = MIN(eb, MAX(eb1, eb2));
        }
        if(a < 0){
            double eb1 = 0, eb2 = 0;
            // y = t
            // check x0 = - (ct + d)/(2a)
            // x0 < -t
            if(- 2*a - c > 0){
                eb1 = MAX(0, d / (- 2*a - c));
            }
            else{
                if(d > 0) eb1 = 1;
            }
            // x0 > t
            if(- 2*a + c > 0){
                eb1 = MAX(eb1, -d / (- 2*a + c));
            }
            else{
                if(d < 0) eb1 = 1;
            }
            // x0 in [-t, t]
            if(f + (d*d)/(-4*a) < 0){
                eb2 = max_eb_to_keep_sign_quadratic(b + (c*c)/(-4*a), e + (c*d)/(-2*a), f + (d*d)/(-4*a));
            }
            eb = MIN(eb, MAX(eb1, eb2));
            // y = -t
            // check x0 = - (-ct + d)/(2a)
            // x0 < -t
            if(- 2*a + c > 0){
                eb1 = MAX(0, d / (- 2*a + c));
            }
            else{
                if(d > 0) eb1 = 1;
            }
            // x0 > t
            if(- 2*a - c > 0){
                eb1 = MAX(eb1, -d / (- 2*a + c));
            }
            else{
                if(d < 0) eb1 = 1;
            }
            // x0 in [-t, t]
            if(f + (d*d)/(-4*a) < 0){
                eb2 = max_eb_to_keep_sign_quadratic(b + (c*c)/(-4*a), -e - (c*d)/(-2*a), f + (d*d)/(-4*a));
            }
            eb = MIN(eb, MAX(eb1, eb2));     
        }
		// else return eb;
	}
    return eb;
}

/*
x1 - x2    rotate    x2 - x3
|    |     ------->  |    |
x0 - x3              x1 - x0
bilinear cell x0, x1, x2, x3, derive cp-preserving eb for x3 given x0, x1, x2
rotate to make (u3, v3) appear only once in A, B, C, D
*/
static double 
derive_cp_eb_bilinear_online(const double u0, const double u1, const double u2, const double u3, const double v0, const double v1, const double v2, const double v3, bool verbose=false){
	// return 1;
	// f(u) = A0 xy + B0 x + C0 y + D0
	// f(v) = A1 xy + B1 x + C1 y + D1
	// solve 
	double A0 = u3 - u0 - u2 + u1;
	double B0 = u0 - u1;
	double C0 = u2 - u1;
	double D0 = u1;

	double A1 = v3 - v0 - v2 + v1;
	double B1 = v0 - v1;
	double C1 = v2 - v1;
	double D1 = v1;

	// M[x, 1]^T = lambda [x, 1]^T
	// 	1/(A0C1 - A1C0)	[-B0C1 + B1C0, -C1D0 + C0D1][x]	= y[x]
	// 									[ A1B0 - A0B1,  A1D0 - A0D1][1]		 [1]
	double A0C1_minus_A1C0 = A0*C1 - A1*C0;
	if(A0C1_minus_A1C0 == 0) return 0;
	double M[4] = {- B0*C1 + B1*C0, - C1*D0 + C0*D1, A1*B0 - A0*B1,  A1*D0 - A0*D1};
	for(int i=0; i<4; i++){
		M[i] /= A0C1_minus_A1C0;
	}
	// original determinant
	double dM = M[0] * M[3] - M[1] * M[2];
	// original trace
	double tM = M[0] + M[3];

	if((dM == 0) || (tM == 0)) return 0;
	// A0C1 - A1C0 = - u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3

	// note: not dividing A0C1 - A1C0 in M
	// M[0, 0] = - u1v0 + u2v0 + u0v1 - u2v1 - u0v2 + u1v2
	// M[0, 1] = u2v1 - u1v2
	// M[1, 0] = u2v0 - u3v0 - u2v1 + u3v1 - u0v2 + u1v2 + u0v3 - u1v3
	// M[1, 1] = -u1v0 + u0v1 + u2v1 -u3v1 

	// tr(M) = -2 u1v0 + u2v0 + 2 u0v1 - u3v1 - u0v2 + u1v3
    // det(M) = u1u1v0v0 - u1u2v0v0 - 2 u0u1v0v1 + u0u2v0v1 + u1u3v0v1 + u0u0v1v1 
    //					- u0u3v1v1 + u0u1v0v2 - u1u3v0v2 - u0u0v1v2 + u0u3v1v2 - u1u1v0v3 
    //						+ u1u2v0v3 + u0u1v1v3 - u0u2v1v3
	double eb = 1;
    if(verbose){
        printf("check delta: tM * tM - 4*dM = %.7f, A0C1 - A1C0 = %.7f\n", tM * tM - 4*dM, A0C1_minus_A1C0);
        printf("unormalized delta = %.7f\n", (tM * tM - 4*dM)*A0C1_minus_A1C0*A0C1_minus_A1C0);
    }
	if(tM * tM - 4*dM == 0){
		// 1 double root
        if(verbose){
            printf("1 double root: eb = 0\n");
        }
		return 0;
	}
	else if(tM * tM - 4*dM < 0) {
		// keep sign of delta
		// keep sign of trace_M * trace_M - 4 * determinant_M
		// tr(M)^2 - 4 det(M) = 
		// 	u2u2v0v0 - 2 u2u3v0v1 + u3u3v1v1 - 2 u0u2v0v2 + 4 u1u3v0v2 - 2 u0u3v1v2
		// +	u0u0v2v2 - 2 u1u2v0v3 + 4 u0u2v1v3 - 2 u1u3v1v3 - 2 u0u1v2v3 + u1u1v3v3
		// = a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f = 0
		double a = u3*u3*v1*v1; // (1+e1)u3(1+e1)u3v1v1
		double b = u1*u1*v3*v3; // u1u1(1+e2)v3(1+e2)v3
		double c = - 2*u1*u3*v1*v3; // -2 u1u3(1+e1)v1v3(1+e2)
		double d = - 2*u2*u3*v0*v1 + 2*u3*u3*v1*v1 + 4*u1*u3*v0*v2 - 2*u0*u3*v1*v2 - 2*u1*u3*v1*v3; // - 2 u2u3v0v1 + (1+e1)u3(1+e1)u3v1v1 + 4 u1u3v0v2 - 2 u0u3v1v2 - 2 u1(1+e1)u3v1(1+e2)v3 
		double e = - 2*u1*u2*v0*v3 + 4*u0*u2*v1*v3 - 2*u1*u3*v1*v3 - 2*u0*u1*v2*v3 + 2*u1*u1*v3*v3; // - 2 u1u2v0v3 + 4 u0u2v1v3 - 2 u1u3v1v3 - 2 u0u1v2v3 + u1u1v3v3
		double f = u2*u2*v0*v0 - 2*u2*u3*v0*v1 + u3*u3*v1*v1 - 2*u0*u2*v0*v2 + 4*u1*u3*v0*v2 - 2*u0*u3*v1*v2
		 +	u0*u0*v2*v2 - 2*u1*u2*v0*v3 + 4*u0*u2*v1*v3 - 2*u1*u3*v1*v3 - 2*u0*u1*v2*v3 + u1*u1*v3*v3;
		// keep sign of A0C1 - A1C0 = - u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3
		{
			eb = MINF(eb, max_eb_to_keep_sign_2d_online(-u3*v1 + u3*v2, u1*v3 - u2*v3, -u1*v0 + u2*v0 + u0*v1 - u0*v2));
		}
        auto tmp_eb = max_eb_to_keep_sign_2d_online_quadratic_cylinder(a, b, c, d, e, f, verbose);
		eb = MINF(eb, tmp_eb);
        if(verbose){
            printf("no roots: eb = %.7f, tmp_eb = %.7f\n", eb, tmp_eb);
        }
	}		
	else{
		// tM * tM - 4*dM > 0
		if((dM > 0) && (1 - tM + dM > 0)){
            if(verbose){
                printf("(dM > 0) && (1 - tM + dM > 0)\n");
            }
			// include 2 cases
			if((tM < 0) || (tM >= 2)){
				// 1) have roots but no root in [0, 1] case 1
				// f(0) > 0, f(1) > 0, Tr(M)/2 not in [0, 1]
			  // keep sign of f(0) = dM
			  {
                    double a = u1*u3*v0*v1 - u0*u3*v1*v1 - u1*u3*v0*v2 + u0*u3*v1*v2;
                    double b = - u1*u1*v0*v3 + u1*u2*v0*v3 + u0*u1*v1*v3 - u0*u2*v1*v3;
                    double c = u1*u1*v0*v0 - u1*u2*v0*v0 - 2*u0*u1*v0*v1 + u0*u2*v0*v1 + u0*u0*v1*v1 + u0*u1*v0*v2 - u0*u0*v1*v2; 
                    eb = MINF(eb, max_eb_to_keep_sign_2d_online(a, b, c));		  	
			  }				
			  // keep sign of f(1) = 1 - tM + dM
			  // (- u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3) ** 2 -
			  // (-2 u1v0 + u2v0 + 2 u0v1 - u3v1 - u0v2 + u1v3) * (- u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3)
			  // + dM
			  // = - u1u3v0v2 + u2u3v0v2 + u0u3v1v2 - u3u3v1v2 - u0u3v2v2 + u3u3v2v2 + u1u2v0v3 - u2u2v0v3 - u0u2v1v3
				// 	+ u2u3v1v3 + u0u2v2v3 + u1u3v2v3 - 2 u2u3v2v3 - u1u2v3v3 + u2u2v3v3
			  {
                    double a = u3*u3*v2*v2 - u3*u3*v1*v2;
                    double b = - u1*u2*v3*v3 + u2*u2*v3*v3;
                    double c = u2*u3*v1*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3;
                    double d = - u1*u3*v0*v2 + u2*u3*v0*v2 + u0*u3*v1*v2 - 2*u3*u3*v1*v2 - u0*u3*v2*v2 + 2*u3*u3*v2*v2 + u2*u3*v1*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3;
                    double e = u1*u2*v0*v3 - u2*u2*v0*v3 - u0*u2*v1*v3	+ u2*u3*v1*v3 + u0*u2*v2*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3 - 2*u1*u2*v3*v3 + 2*u2*u2*v3*v3;
                    double f = - u1*u3*v0*v2 + u2*u3*v0*v2 + u0*u3*v1*v2 - u3*u3*v1*v2 - u0*u3*v2*v2 + u3*u3*v2*v2 + u1*u2*v0*v3 - u2*u2*v0*v3 - u0*u2*v1*v3
                    						+ u2*u3*v1*v3 + u0*u2*v2*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3 - u1*u2*v3*v3 + u2*u2*v3*v3;
					eb = MINF(eb, max_eb_to_keep_sign_2d_online_general(a, b, c, d, e, f));
			  }
			  {
					// keep sign of A0C1 - A1C0 = - u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3
					eb = MINF(eb, max_eb_to_keep_sign_2d_online(-u3*v1 + u3*v2, u1*v3 - u2*v3, -u1*v0 + u2*v0 + u0*v1 - u0*v2));
                    // keep sign of tM 
                    double a = - u3*v1;
                    double b = u1*v3;
                    double c = -2*u1*v0 + u2*v0 + 2*u0*v1 - u0*v2;
                    double eb1 = max_eb_to_keep_sign_2d_online(a, b, c);
                    // keep sign of tr(M) - 2(A0C1 - A1C0)
                    a -= 2*(- u3*v1 + u3*v2);
                    b -= 2*(u1*v3 - u2*v3);
                    c -= 2*(- u1*v0 + u2*v0 + u0*v1 - u0*v2);
                    double eb2 = max_eb_to_keep_sign_2d_online(a, b, c);
                    // only need to keep either 
                    eb = MINF(eb, MAX(eb1, eb2));
			  }
			}
			else{
				// 2) two roots in [0, 1]
				// f(0) > 0, f(1) > 0, Tr(M)/2 in [0, 1]
				// keep sign of delta
				{
					double a = u3*u3*v1*v1; // (1+e1)u3(1+e1)u3v1v1
					double b = u1*u1*v3*v3; // u1u1(1+e2)v3(1+e2)v3
					double c = - 2*u1*u3*v1*v3; // -2 u1u3(1+e1)v1v3(1+e2)
					double d = - 2*u2*u3*v0*v1 + 2*u3*u3*v1*v1 + 4*u1*u3*v0*v2 - 2*u0*u3*v1*v2 - 2*u1*u3*v1*v3; // - 2 u2u3v0v1 + (1+e1)u3(1+e1)u3v1v1 + 4 u1u3v0v2 - 2 u0u3v1v2 - 2 u1(1+e1)u3v1(1+e2)v3 
					double e = - 2*u1*u2*v0*v3 + 4*u0*u2*v1*v3 - 2*u1*u3*v1*v3 - 2*u0*u1*v2*v3 + 2*u1*u1*v3*v3; // - 2 u1u2v0v3 + 4 u0u2v1v3 - 2 u1u3v1v3 - 2 u0u1v2v3 + u1u1v3v3
					double f = u2*u2*v0*v0 - 2*u2*u3*v0*v1 + u3*u3*v1*v1 - 2*u0*u2*v0*v2 + 4*u1*u3*v0*v2 - 2*u0*u3*v1*v2
					 +	u0*u0*v2*v2 - 2*u1*u2*v0*v3 + 4*u0*u2*v1*v3 - 2*u1*u3*v1*v3 - 2*u0*u1*v2*v3 + u1*u1*v3*v3;
					// = a*e1^2 + b*e2^2 + c*e1e2 + d*e1 + e*e2 + f > 0
					eb = MINF(eb, max_eb_to_keep_sign_2d_online_quadratic_cylinder(a, b, c, d, e, f));
				}
                // keep sign of f(0) = dM
                {
                    double a = u1*u3*v0*v1 - u0*u3*v1*v1 - u1*u3*v0*v2 + u0*u3*v1*v2;
                    double b = - u1*u1*v0*v3 + u1*u2*v0*v3 + u0*u1*v1*v3 - u0*u2*v1*v3;
                    double c = u1*u1*v0*v0 - u1*u2*v0*v0 - 2*u0*u1*v0*v1 + u0*u2*v0*v1 + u0*u0*v1*v1 + u0*u1*v0*v2 - u0*u0*v1*v2; 
                    eb = MINF(eb, max_eb_to_keep_sign_2d_online(a, b, c));		  	
                }
                // tr(M) = (-2 u1v0 + u2v0 + 2 u0v1 - u3v1 - u0v2 + u1v3) / (- u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3)
                {
					// keep sign of A0C1 - A1C0 = - u1v0 + u2v0 + u0v1 - u3v1 - u0v2 + u3v2 + u1v3 - u2v3
					eb = MINF(eb, max_eb_to_keep_sign_2d_online(-u3*v1 + u3*v2, u1*v3 - u2*v3, -u1*v0 + u2*v0 + u0*v1 - u0*v2));
                    // keep sign of tM 
                    double a = - u3*v1;
                    double b = u1*v3;
                    double c = -2*u1*v0 + u2*v0 + 2*u0*v1 - u0*v2;
                    eb = MINF(eb, max_eb_to_keep_sign_2d_online(a, b, c));
                    // keep sign of tr(M) - 2(A0C1 - A1C0)
                    a -= 2*(- u3*v1 + u3*v2);
                    b -= 2*(u1*v3 - u2*v3);
                    c -= 2*(- u1*v0 + u2*v0 + u0*v1 - u0*v2);
                    eb = MINF(eb, max_eb_to_keep_sign_2d_online(a, b, c));
                }
                // keep sign of f(1) = 1 - tM + dM
                {
                    double a = u3*u3*v2*v2 - u3*u3*v1*v2;
                    double b = - u1*u2*v3*v3 + u2*u2*v3*v3;
                    double c = u2*u3*v1*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3;
                    double d = - u1*u3*v0*v2 + u2*u3*v0*v2 + u0*u3*v1*v2 - 2*u3*u3*v1*v2 - u0*u3*v2*v2 + 2*u3*u3*v2*v2 + u2*u3*v1*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3;
                    double e = u1*u2*v0*v3 - u2*u2*v0*v3 - u0*u2*v1*v3	+ u2*u3*v1*v3 + u0*u2*v2*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3 - 2*u1*u2*v3*v3 + 2*u2*u2*v3*v3;
                    double f = - u1*u3*v0*v2 + u2*u3*v0*v2 + u0*u3*v1*v2 - u3*u3*v1*v2 - u0*u3*v2*v2 + u3*u3*v2*v2 + u1*u2*v0*v3 - u2*u2*v0*v3 - u0*u2*v1*v3
                    						+ u2*u3*v1*v3 + u0*u2*v2*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3 - u1*u2*v3*v3 + u2*u2*v3*v3;
                    eb = MINF(eb, max_eb_to_keep_sign_2d_online_general(a, b, c, d, e, f));
                }
            }
		}
		else if((dM < 0) && (1 - tM + dM < 0)){
            if(verbose){
                printf("(dM < 0) && (1 - tM + dM < 0)\n");
            }
			// have roots but no root in [0, 1] case 2
			// f(0) < 0, f(1) < 0
			// keep sign of f(0) = dM
            {
                double a = u1*u3*v0*v1 - u0*u3*v1*v1 - u1*u3*v0*v2 + u0*u3*v1*v2;
                double b = - u1*u1*v0*v3 + u1*u2*v0*v3 + u0*u1*v1*v3 - u0*u2*v1*v3;
                double c = u1*u1*v0*v0 - u1*u2*v0*v0 - 2*u0*u1*v0*v1 + u0*u2*v0*v1 + u0*u0*v1*v1 + u0*u1*v0*v2 - u0*u0*v1*v2; 
                auto tmp_eb = max_eb_to_keep_sign_2d_online(a, b, c);
                eb = MINF(eb, tmp_eb);		
                if(verbose) printf("tmp_eb = %.7f\n", tmp_eb);  	
            }
			// keep sign of f(1) = 1 - tM + dM
            {
                double a = u3*u3*v2*v2 - u3*u3*v1*v2;
                double b = - u1*u2*v3*v3 + u2*u2*v3*v3;
                double c = u2*u3*v1*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3;
                double d = - u1*u3*v0*v2 + u2*u3*v0*v2 + u0*u3*v1*v2 - 2*u3*u3*v1*v2 - u0*u3*v2*v2 + 2*u3*u3*v2*v2 + u2*u3*v1*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3;
                double e = u1*u2*v0*v3 - u2*u2*v0*v3 - u0*u2*v1*v3	+ u2*u3*v1*v3 + u0*u2*v2*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3 - 2*u1*u2*v3*v3 + 2*u2*u2*v3*v3;
                double f = - u1*u3*v0*v2 + u2*u3*v0*v2 + u0*u3*v1*v2 - u3*u3*v1*v2 - u0*u3*v2*v2 + u3*u3*v2*v2 + u1*u2*v0*v3 - u2*u2*v0*v3 - u0*u2*v1*v3
                						+ u2*u3*v1*v3 + u0*u2*v2*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3 - u1*u2*v3*v3 + u2*u2*v3*v3;
                // f < 0
                auto tmp_eb = max_eb_to_keep_sign_2d_online_general(a, b, c, d, e, f, verbose);
                eb = MINF(eb, tmp_eb);
                if(verbose) printf("tmp_eb = %.7f\n", tmp_eb);      
            }
		}
		else{
            if(verbose){
                printf("dM * (1 - tM + dM) < 0\n");
            }
			// dM * (1 - tM + dM) < 0
			// one root in [0, 1]
			// keep sign of f(0) = dM
            {
                double a = u1*u3*v0*v1 - u0*u3*v1*v1 - u1*u3*v0*v2 + u0*u3*v1*v2;
                double b = - u1*u1*v0*v3 + u1*u2*v0*v3 + u0*u1*v1*v3 - u0*u2*v1*v3;
                double c = u1*u1*v0*v0 - u1*u2*v0*v0 - 2*u0*u1*v0*v1 + u0*u2*v0*v1 + u0*u0*v1*v1 + u0*u1*v0*v2 - u0*u0*v1*v2; 
                eb = MINF(eb, max_eb_to_keep_sign_2d_online(a, b, c));		  	
            }
            // keep sign of f(1) = (1 - tM + dM)
            {
                double a = u3*u3*v2*v2 - u3*u3*v1*v2;
                double b = - u1*u2*v3*v3 + u2*u2*v3*v3;
                double c = u2*u3*v1*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3;
                double d = - u1*u3*v0*v2 + u2*u3*v0*v2 + u0*u3*v1*v2 - 2*u3*u3*v1*v2 - u0*u3*v2*v2 + 2*u3*u3*v2*v2 + u2*u3*v1*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3;
                double e = u1*u2*v0*v3 - u2*u2*v0*v3 - u0*u2*v1*v3	+ u2*u3*v1*v3 + u0*u2*v2*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3 - 2*u1*u2*v3*v3 + 2*u2*u2*v3*v3;
                double f = - u1*u3*v0*v2 + u2*u3*v0*v2 + u0*u3*v1*v2 - u3*u3*v1*v2 - u0*u3*v2*v2 + u3*u3*v2*v2 + u1*u2*v0*v3 - u2*u2*v0*v3 - u0*u2*v1*v3
                						+ u2*u3*v1*v3 + u0*u2*v2*v3 + u1*u3*v2*v3 - 2*u2*u3*v2*v3 - u1*u2*v3*v3 + u2*u2*v3*v3;
                eb = MINF(eb, max_eb_to_keep_sign_2d_online_general(a, b, c, d, e, f));
            }
		}
        if(verbose){
            printf("1 root: eb = %.7f\n", eb);
        }
	}
	// need to check x and type
	return eb;
}

int inverse_bilinear_interpolation(const double A0, const double B0, const double C0, const double D0,
	const double A1, const double B1, const double C1, const double D1, double pos[2][2], double J[2][2][2])
{
  double M0[4] = {-B0, -D0, -B1, -D1}, // stored in row major
    M1[4] = {A0, C0, A1, C1}; // (yM1 - M0)v = 0, v = {x, 1}^T

  double detM1 = A0*C1 - A1*C0; // TODO: check if detM1==0
  double invM1[4] = {C1/detM1, -C0/detM1, -A1/detM1, A0/detM1};  
  // Q = invM1*M0
  double Q[4] = {
    invM1[0]*M0[0] + invM1[1]*M0[2], 
    invM1[0]*M0[1] + invM1[1]*M0[3], 
    invM1[2]*M0[0] + invM1[3]*M0[2], 
    invM1[2]*M0[1] + invM1[3]*M0[3]
  };
  // compute y=eig(Q)
  double trace = Q[0] + Q[3];
  double det = Q[0]*Q[3] - Q[1]*Q[2];

  if(trace*trace/4 - det < 0) return 0;

  double lambda[2] = {
    static_cast<double>(trace/2 + std::sqrt(trace*trace/4 - det)), 
    static_cast<double>(trace/2 - std::sqrt(trace*trace/4 - det))
  }; 

  double x[2] = {
    (lambda[0]-Q[3])/Q[2], 
    (lambda[1]-Q[3])/Q[2]
  }; 
  double y[2] = {
    lambda[0], 
    lambda[1]
  };
  int nroots = 0;
  for (int i=0; i<2; i++) // check the two roots 
    if (x[i]>=0 && x[i]<1 && y[i]>=0 && y[i]<1) {
      pos[nroots][0] = x[i];
      pos[nroots][1] = y[i];
      J[nroots][0][0] = A0 * y[i] + B0;
      J[nroots][0][1] = A0 * x[i] + C0;
      J[nroots][1][0] = A1 * y[i] + B1;
      J[nroots][1][1] = A1 * x[i] + C1;
      nroots ++;
    }
  return nroots;
}

// copied from ftk
template <typename T>
static std::complex<T> complex_sqrt(const std::complex<T> z)
{
  return pow(z, T(1)/T(2));
}
template <typename T>
inline T solve_quadratic(const T P[3], std::complex<T> x[2])
{
  const T delta = P[1]*P[1] - 4*P[2]*P[0];
  if (delta >= 0) {
    x[0] = (-P[1] + sqrt(delta)) / (2 * P[2]);
    x[1] = (-P[1] - sqrt(delta)) / (2 * P[2]);
  } else {
    x[0] = (-P[1] + complex_sqrt<T>(delta)) / (2 * P[2]);
    x[1] = (-P[1] - complex_sqrt<T>(delta)) / (2 * P[2]);
  }
  return delta;
}
template <typename T>
inline T trace2(T A[2][2])
{
  return A[0][0] + A[1][1];
}
template <typename T>
inline T det2(const T A[2][2])
{
  return A[0][0] * A[1][1] - A[1][0] * A[0][1];
}
template <typename T>
void characteristic_polynomial_2x2(const T A[2][2], T P[3])
{
  P[2] = T(1);
  P[1] = -trace2(A);
  P[0] = det2(A);
}
template <typename T>
inline T solve_eigenvalues2x2(const T M[2][2], std::complex<T> eig[2])
{
  T P[3];
  characteristic_polynomial_2x2(M, P);
  return solve_quadratic(P, eig); // returns delta
}
// copied from ftk end

int get_cp_type(double delta, std::complex<double> eig[2]){
	int cp_type = 0;
  if (delta >= 0) { // two real roots
    if (eig[0].real() * eig[1].real() < 0) {
      cp_type = 3;
    } else if (eig[0].real() < 0) {
      cp_type = 1;
    }
    else if (eig[0].real() > 0){
      cp_type = 2;
    }
    else cp_type = 0;
  } else { // two conjugate roots
    if (eig[0].real() < 0) {
      cp_type = 4;
    } else if (eig[0].real() > 0) {
      cp_type = 5;
    } else 
      cp_type = 6;
  }
  return cp_type;
}

/*
x1 - x2    
|    |     
x0 - x3   
No rotation to avoid wrong types
*/
static int 
bilinear_extract_critical_point(const double u0, const double u1, const double u2, const double u3,
			const double v0, const double v1, const double v2, const double v3, double J[2][2][2], bool verbose=false){
	// solve original
  double 	A0 = u0 - u1 - u3 + u2,  // Axy + Bx + Cy + D = 0
    			B0 = u1 - u0, 
    			C0 = u3 - u0, 
    			D0 = u0,
    			A1 = v0 - v1 - v3 + v2, 
    			B1 = v1 - v0, 
    			C1 = v3 - v0, 
    			D1 = v0; 

	double pos[2][2];
	int num_root = inverse_bilinear_interpolation(A0, B0, C0, D0, A1, B1, C1, D1, pos, J);
	return num_root;
}

static bool 
bilinear_verify_critical_point(int nroots, const double J[2][2][2], double J_[2][2][2], bool verbose=false){
	for(int i=0; i<nroots; i++){
		std::complex<double> eig[2];
		double delta = solve_eigenvalues2x2(J[i], eig);
		std::complex<double> eig_[2];
		double delta_ = solve_eigenvalues2x2(J_[i], eig_);
		if(verbose){
			std::cout << get_cp_type(delta, eig) << " " << get_cp_type(delta_, eig_) << std::endl;
		}
		if(get_cp_type(delta, eig) != get_cp_type(delta_, eig_)) return false;
		// if(delta * delta_ < 0) return false;
		// if(eig[0].real() * eig_[0].real() < 0) return false;
		// if(eig[1].real() * eig_[1].real() < 0) return false;	
	}
	return true;
}

template<typename T>
unsigned char *
sz_compress_cp_preserve_2d_bilinear_online_log(const T * U, const T * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb){

	size_t num_elements = r1 * r2;
	size_t sign_map_size = (num_elements - 1)/8 + 1;
	unsigned char * sign_map_compressed = (unsigned char *) malloc(2*sign_map_size);
	unsigned char * sign_map_compressed_pos = sign_map_compressed;
	unsigned char * sign_map = (unsigned char *) malloc(num_elements*sizeof(unsigned char));
	// Note the convert function has address auto increment
	T * log_U = log_transform(U, sign_map, num_elements);
	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map, num_elements, sign_map_compressed_pos);
	T * log_V = log_transform(V, sign_map, num_elements);
	convertIntArray2ByteArray_fast_1b_to_result_sz(sign_map, num_elements, sign_map_compressed_pos);
	free(sign_map);

    double * eb = (double *) malloc(num_elements*sizeof(double));
	T * decompressed_U = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_U, U, num_elements*sizeof(T));
	T * decompressed_V = (T *) malloc(num_elements*sizeof(T));
	memcpy(decompressed_V, V, num_elements*sizeof(T));

	int * eb_quant_index = (int *) malloc(num_elements*sizeof(int));
	int * data_quant_index = (int *) malloc(2*num_elements*sizeof(int));
	int * eb_quant_index_pos = eb_quant_index;
	int * data_quant_index_pos = data_quant_index;
	// next, row by row
	const int base = 2;
	const double log_of_base = log2(base);
	const int capacity = 65536;
	const int intv_radius = (capacity >> 1);
	unpred_vec<T> unpred_data;
	// offsets to get six adjacent square mesh indices
	// T8 rolls back to T0
	/*
	|	T3	T4	T5
	y	T2	X 	T6
	|	T1	T0  T7
		-	x 	-
	*/
	const int offsets[9] = {
		-(int)r2, -(int)r2 - 1, -1, (int)r2 - 1, (int)r2, (int)r2+1, 1, -(int)r2 + 1, -(int)r2
	};
	// relative indices for the four squares
	// we always put the current index as (x3, y3)
	const T x[4][4] = {
		{1, 0, 0, 1}, // T0_T1_T2_X
		{0, 0, 1, 1}, // T2_T3_T4_X
		{0, 1, 1, 0}, // T4_T5_T6_X
		{1, 1, 0, 0}, // T6_T7_T0_X
	};
	const T y[4][4] = {
		{0, 0, 1, 1},
		{0, 1, 1, 0},
		{1, 1, 0, 0},
		{1, 0, 0, 1},
	};
  // compute offset of index for the surrounding cells	
	int index_offset[4][3][2];
	for(int i=0; i<4; i++){
		for(int j=0; j<3; j++){
			index_offset[i][j][0] = x[i][j] - x[i][3];
			index_offset[i][j][1] = y[i][j] - y[i][3];
		}
	}
	T * cur_log_U_pos = log_U;
	T * cur_log_V_pos = log_V;
	T * cur_U_pos = decompressed_U;
	T * cur_V_pos = decompressed_V;

    int index_eb = 0;
    int failed_verification = 0;
	for(int i=0; i<r1; i++){
		// printf("start %d row\n", i);
		for(int j=0; j<r2; j++){
			double required_eb = max_pwr_eb;
			// derive eb given four adjacent cells
			for(int k=0; k<4; k++){
				bool in_mesh = true;
				for(int p=0; p<3; p++){
					// reserved order!
					// note: x and j are on r2, y and i are on r1
					if(!(in_range(i + index_offset[k][p][1], (int)r1) && in_range(j + index_offset[k][p][0], (int)r2))){
						in_mesh = false;
						break;
					}
				}
				if(in_mesh){
                    bool verbose = false;
                    auto derived_eb = derive_cp_eb_bilinear_online(cur_U_pos[offsets[2*k]], cur_U_pos[offsets[2*k+1]], cur_U_pos[offsets[2*k+2]], cur_U_pos[0],
                        cur_V_pos[offsets[2*k]], cur_V_pos[offsets[2*k+1]], cur_V_pos[offsets[2*k+2]], cur_V_pos[0], verbose);
					required_eb = MINF(required_eb, derived_eb);
				}
			}
			eb[index_eb++] = required_eb;
			if((required_eb > 0) && (*cur_U_pos != 0) && (*cur_V_pos != 0)){
				bool unpred_flag = false;
				T decompressed[2];
				double abs_eb = log2(1 + required_eb);
				*eb_quant_index_pos = eb_exponential_quantize(abs_eb, base, log_of_base);
				// *eb_quant_index_pos = eb_linear_quantize(abs_eb, 1e-2);
				if(*eb_quant_index_pos > 0){
					// compress U and V
					for(int k=0; k<2; k++){
						T * cur_data_pos = (k == 0) ? cur_log_U_pos : cur_log_V_pos;
						T cur_data = *cur_data_pos;
						// get adjacent data and perform Lorenzo
						/*
							d2 X
							d0 d1
						*/
						T d0 = (i && j) ? cur_data_pos[-1 - r2] : 0;
						T d1 = (i) ? cur_data_pos[-r2] : 0;
						T d2 = (j) ? cur_data_pos[-1] : 0;
						T pred = d1 + d2 - d0;
						double diff = cur_data - pred;
						double quant_diff = fabs(diff) / abs_eb + 1;
						if(quant_diff < capacity){
							quant_diff = (diff > 0) ? quant_diff : -quant_diff;
							int quant_index = (int)(quant_diff/2) + intv_radius;
							data_quant_index_pos[k] = quant_index;
							decompressed[k] = pred + 2 * (quant_index - intv_radius) * abs_eb; 
							// check original data
							if(fabs(decompressed[k] - cur_data) >= abs_eb){
								unpred_flag = true;
								break;
							}
						}
						else{
							unpred_flag = true;
							break;
						}
					}
					// do verification for x and critical point type
					if(!unpred_flag){
						bool verbose = false;
						// if((fabs(i - 1688.66) < 1) && (fabs(j - 966.522) < 1)) verbose = true;
						double decompressed_u = (cur_U_pos[0] > 0) ? exp2(decompressed[0]) : -exp2(decompressed[0]);
						double decompressed_v = (cur_V_pos[0] > 0) ? exp2(decompressed[1]) : -exp2(decompressed[1]);
						{
							// k = 0
                            // need to deal with k separately due to different relative vertex positions
							{
								int k = 0;
								// if(verbose) std::cout << k << ":\n";
								bool in_mesh = true;
								for(int p=0; p<3; p++){
									// reserved order!
									// note: x and j are on r2, y and i are on r1
									if(!(in_range(i + index_offset[k][p][1], (int)r1) && in_range(j + index_offset[k][p][0], (int)r2))){
										in_mesh = false;
										break;
									}
								}
								if(in_mesh){
									double J[2][2][2];
									double J_[2][2][2];
									int nroots = bilinear_extract_critical_point(cur_U_pos[-1 - r2], cur_U_pos[-r2], cur_U_pos[0], cur_U_pos[-1],
										cur_V_pos[-1 - r2], cur_V_pos[-r2], cur_V_pos[0], cur_V_pos[-1], J);
									int nroots_ = bilinear_extract_critical_point(cur_U_pos[-1 - r2], cur_U_pos[-r2], decompressed_u, cur_U_pos[-1],
										cur_V_pos[-1 - r2], cur_V_pos[-r2], decompressed_v, cur_V_pos[-1], J_);
									if((nroots != nroots_) || (!bilinear_verify_critical_point(nroots, J, J_, verbose))) unpred_flag = true;
									if(unpred_flag){
                                        // printf("i = %d, j = %d, k = %d\n", i, j, k);
                                        // std::cout << "#roots: " << nroots << " " << nroots_ << ", unpred_flag = " << unpred_flag << std::endl;
                                        // verbose = true;
                                        // printf("~~~~~~~~~~\n");
                                        derive_cp_eb_bilinear_online(cur_U_pos[offsets[2*k]], cur_U_pos[offsets[2*k+1]], cur_U_pos[offsets[2*k+2]], decompressed_u,
                                        cur_V_pos[offsets[2*k]], cur_V_pos[offsets[2*k+1]], cur_V_pos[offsets[2*k+2]], decompressed_v, verbose);
                                        // printf("U: %.7f, %.7f, %.7f, %.7f\n", cur_U_pos[offsets[2*k]], cur_U_pos[offsets[2*k+1]], cur_U_pos[offsets[2*k+2]], decompressed_u);
                                        // printf("V: %.7f, %.7f, %.7f, %.7f\n", cur_V_pos[offsets[2*k]], cur_V_pos[offsets[2*k+1]], cur_V_pos[offsets[2*k+2]], decompressed_v);
                                        // printf("x = %.7f, y = %.7f\n", decompressed_u/cur_U_pos[0] - 1, decompressed_v/cur_V_pos[0] - 1);
                                    }
								}
							}
							if(!unpred_flag){
								// k = 1
								int k = 1;
								// if(verbose) std::cout << k << ":\n";
								bool in_mesh = true;
								for(int p=0; p<3; p++){
									// reserved order!
									// note: x and j are on r2, y and i are on r1
									if(!(in_range(i + index_offset[k][p][1], (int)r1) && in_range(j + index_offset[k][p][0], (int)r2))){
										in_mesh = false;
										break;
									}
								}
								if(in_mesh){
									double J[2][2][2];
									double J_[2][2][2];
									int nroots = bilinear_extract_critical_point(cur_U_pos[-1], cur_U_pos[0], cur_U_pos[r2], cur_U_pos[r2 - 1],
										cur_V_pos[-1], cur_V_pos[0], cur_V_pos[r2], cur_V_pos[r2 - 1], J);
									int nroots_ = bilinear_extract_critical_point(cur_U_pos[-1], decompressed_u, cur_U_pos[r2], cur_U_pos[r2 - 1],
										cur_V_pos[-1], decompressed_v, cur_V_pos[r2], cur_V_pos[r2 - 1], J_);
									if((nroots != nroots_) || (!bilinear_verify_critical_point(nroots, J, J_, verbose))) unpred_flag = true;
                                    // if(unpred_flag){
                                    //     printf("i = %d, j = %d, k = %d\n", i, j, k);
                                    //     std::cout << "#roots: " << nroots << " " << nroots_ << ", unpred_flag = " << unpred_flag << std::endl;
                                    // }
								}
							}
							if(!unpred_flag){
								// k = 2
								int k = 2;
								// if(verbose) std::cout << k << ":\n";
								bool in_mesh = true;
								for(int p=0; p<3; p++){
									// reserved order!
									// note: x and j are on r2, y and i are on r1
									if(!(in_range(i + index_offset[k][p][1], (int)r1) && in_range(j + index_offset[k][p][0], (int)r2))){
										in_mesh = false;
										break;
									}
								}
								if(in_mesh){
									double J[2][2][2];
									double J_[2][2][2];
									int nroots = bilinear_extract_critical_point(cur_U_pos[0], cur_U_pos[1], cur_U_pos[r2 + 1], cur_U_pos[r2],
										cur_V_pos[0], cur_V_pos[1], cur_V_pos[r2 + 1], cur_V_pos[r2], J);
									int nroots_ = bilinear_extract_critical_point(decompressed_u, cur_U_pos[1], cur_U_pos[r2 + 1], cur_U_pos[r2],
										decompressed_v, cur_V_pos[1], cur_V_pos[r2 + 1], cur_V_pos[r2], J_);
									if((nroots != nroots_) || (!bilinear_verify_critical_point(nroots, J, J_, verbose))) unpred_flag = true;
                                    // if(unpred_flag){
                                    //     printf("i = %d, j = %d, k = %d\n", i, j, k);
                                    //     std::cout << "#roots: " << nroots << " " << nroots_ << ", unpred_flag = " << unpred_flag << std::endl;
                                    // }
								}
							}
							if(!unpred_flag){
								// k = 3
								int k = 3;
								// if(verbose) std::cout << k << ":\n";
								bool in_mesh = true;
								for(int p=0; p<3; p++){
									// reserved order!
									// note: x and j are on r2, y and i are on r1
									if(!(in_range(i + index_offset[k][p][1], (int)r1) && in_range(j + index_offset[k][p][0], (int)r2))){
										in_mesh = false;
										break;
									}
								}
								if(in_mesh){
									double J[2][2][2];
									double J_[2][2][2];
									int nroots = bilinear_extract_critical_point(cur_U_pos[- r2], cur_U_pos[-r2 + 1], cur_U_pos[1], cur_U_pos[0],
										cur_V_pos[- r2], cur_V_pos[-r2 + 1], cur_V_pos[1], cur_V_pos[0], J);
									int nroots_ = bilinear_extract_critical_point(cur_U_pos[- r2], cur_U_pos[-r2 + 1], cur_U_pos[1], decompressed_u,
										cur_V_pos[- r2], cur_V_pos[-r2 + 1], cur_V_pos[1], decompressed_v, J_);
									if((nroots != nroots_) || (!bilinear_verify_critical_point(nroots, J, J_, verbose))) unpred_flag = true;
                                    // if(unpred_flag){
                                    //     printf("i = %d, j = %d, k = %d\n", i, j, k);
                                    //     std::cout << "#roots: " << nroots << " " << nroots_ << ", unpred_flag = " << unpred_flag << std::endl;
                                    // }
								}
							}
                            // record number of failed verification
                            if(unpred_flag) failed_verification ++;
						}
					}
				}
				else unpred_flag = true;
				if(unpred_flag){
					// recover quant index
					*(eb_quant_index_pos ++) = 0;
					*(data_quant_index_pos ++) = intv_radius;
					*(data_quant_index_pos ++) = intv_radius;
					unpred_data.push_back(*cur_U_pos);
					unpred_data.push_back(*cur_V_pos);
				}
				else{
					eb_quant_index_pos ++;
					data_quant_index_pos += 2;
					// assign decompressed data
					*cur_log_U_pos = decompressed[0];
					*cur_log_V_pos = decompressed[1];
					*cur_U_pos = (*cur_U_pos > 0) ? exp2(*cur_log_U_pos) : -exp2(*cur_log_U_pos);
					*cur_V_pos = (*cur_V_pos > 0) ? exp2(*cur_log_V_pos) : -exp2(*cur_log_V_pos);
				}
			}
			else{
				// record as unpredictable data
				*(eb_quant_index_pos ++) = 0;
				*(data_quant_index_pos ++) = intv_radius;
				*(data_quant_index_pos ++) = intv_radius;
				unpred_data.push_back(*cur_U_pos);
				unpred_data.push_back(*cur_V_pos);
			}
			cur_log_U_pos ++, cur_log_V_pos ++;
			cur_U_pos ++, cur_V_pos ++;
		}
	}
    // printf("#Failed verification = %d\n", failed_verification);
	// writefile("eb_2d.dat", eb, num_elements);
	// free(eb);
	free(log_U);
	free(log_V);
	free(decompressed_U);
	free(decompressed_V);
	// printf("offsets eb_q, data_q, unpred: %ld %ld %ld\n", eb_quant_index_pos - eb_quant_index, data_quant_index_pos - data_quant_index, unpred_data.size());
	unsigned char * compressed = (unsigned char *) malloc(3*num_elements*sizeof(T));
	unsigned char * compressed_pos = compressed;
	write_variable_to_dst(compressed_pos, base);
	write_variable_to_dst(compressed_pos, intv_radius);
	write_array_to_dst(compressed_pos, sign_map_compressed, 2*sign_map_size);
	free(sign_map_compressed);
	size_t unpredictable_count = unpred_data.size();
	// printf("unpredictable_count=%d\n", unpredictable_count);
	write_variable_to_dst(compressed_pos, unpredictable_count);
	write_array_to_dst(compressed_pos, (T *)&unpred_data[0], unpredictable_count);	
	Huffman_encode_tree_and_data(2*1024, eb_quant_index, num_elements, compressed_pos);
	free(eb_quant_index);
	Huffman_encode_tree_and_data(2*capacity, data_quant_index, 2*num_elements, compressed_pos);
	free(data_quant_index);
	compressed_size = compressed_pos - compressed;
	return compressed;	
}

template
unsigned char *
sz_compress_cp_preserve_2d_bilinear_online_log(const float * U, const float * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);

template
unsigned char *
sz_compress_cp_preserve_2d_bilinear_online_log(const double * U, const double * V, size_t r1, size_t r2, size_t& compressed_size, bool transpose, double max_pwr_eb);
