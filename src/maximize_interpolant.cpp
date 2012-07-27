#include "edgeR.h"
#include "fmm_spline.h"

/**************************************
 * Quadratic solver for a*x^2+b*x+c = 0.
 * 
 **************************************/

struct solution {
    double sol1, sol2;
    bool solvable;
};

solution quad_solver (const double& a, const double& b, const double& c) {
    double back=b*b-4*a*c;
    solution cur_sol;
    if (back<0) {
        cur_sol.solvable=false;
        return(cur_sol);
    }
    double front=-b/(2*a);
    back=std::sqrt(back)/(2*a);
    cur_sol.sol1=front-back;
    cur_sol.sol2=front+back;
    cur_sol.solvable=true;
    return(cur_sol);
}

/************************************
 * maximize_interpolant: mark II
 *
 * It fits the spline and grabs the coefficients of each segment.
 * It then calculates the maxima at the segments neighbouring
 * the maximally highest point. This avoids NR optimisation 
 * as well as the need to call R's splinefun's from within C.
 *
 ***********************************/

extern "C" {

SEXP maximize_interpolant(SEXP spline_pts, SEXP likelihoods) {
    PROTECT(spline_pts=AS_NUMERIC(spline_pts));
    PROTECT(likelihoods=AS_NUMERIC(likelihoods));

    // Loading in the spline x-axis values.
    const int num_pts=LENGTH(spline_pts);
    double* spline_ptr=NUMERIC_POINTER(spline_pts);
    if (num_pts<2) { 
        error("Insufficient spline points provided for interpolation.");
    }

    // Setting up coefficients for the spline interpolation.
    double* x=new double [num_pts];
    double* y=new double [num_pts];
    for (int i=0; i<num_pts; ++i) { x[i]=spline_ptr[i]; }
    double* b=new double [num_pts];
    double* c=new double [num_pts];
    double* d=new double [num_pts];
     
    // Counter for the number of divergences i.e. maximum is out of interpolation range.
    long diverged=0, first=0;
    
    // Setting up the output object.
    const long nrows=LENGTH(likelihoods)/((double) num_pts)+0.5;
    SEXP output;
    PROTECT(output=allocVector(REALSXP, nrows));
    double* out_ptr=NUMERIC_POINTER(output);

    // Now running through the likelihood matrix. 
    double** l_ptrs=new double* [num_pts];
    for (int i=0; i<num_pts; ++i) { 
        if (i==0) { l_ptrs[i]=NUMERIC_POINTER(likelihoods); }
        else { l_ptrs[i]=l_ptrs[i-1]+nrows; }
    }
    long count=0;
    while (count < nrows) {       
        // Defining the y-values to be fitted to the spline. 
        double maxed=-1;
        int maxed_at=-1;
        for (int i=0; i<num_pts; ++i) { 
            y[i]=*(l_ptrs[i]++); 
            // Also getting a good initial guess for the MLE.
            if (maxed_at < 0 || y[i] > maxed) { 
                maxed=y[i];
                maxed_at=i;
            }
        }
        double x_max=x[maxed_at];

        // Calculating the spline coefficients. b, c, and d are overwritten and don't need to be reset.
        fmm_spline(num_pts, x, y, b, c, d);

        /* We want to find the local maximum around our initial guess (specifically, we want the x-value
         * at which this maximum occurs i.e. the MLE). We solve the derivatives of the surrounding cubic 
         * segments on either side of the guessed datapoint. First, the left side (though no cubic
         * exists if the initial guess is the first spline point).
         */
        if (maxed_at>0) {
            const double& ld=d[maxed_at-1];
            const double& lc=c[maxed_at-1];
            const double& lb=b[maxed_at-1];
            solution sol_left=quad_solver(3*ld, 2*lc, lb);
            if (sol_left.solvable) {
                /* Using the solution with the maximum (not minimum). If the curve is mostly increasing, the 
                 * maximal point is located at the smaller solution (i.e. sol1 for a>0). If the curve is mostly
                 * decreasing, the maximal is located at the larger solution (i.e., sol1 for a<0).
                 */
                const double& chosen_sol=sol_left.sol1; 

                /* The spline coefficients are designed such that 'x' in 'y + b*x + c*x^2 + d*x^3' is
                 * equal to 'x_t - x_l' where 'x_l' is the left limit of that spline segment and 'x_t'
                 * is where you want to get an interpolated value. This is necessary in 'splinefun' to 
                 * ensure that you get 'y' (i.e. the original data point) when 'x=0'. For our purposes, 
                 * the actual MLE corresponds to 'x_t' and is equal to 'solution + x_0'.
                 */
                if (chosen_sol > 0 && chosen_sol < x[maxed_at]-x[maxed_at-1]) {
                    double temp=((ld*chosen_sol+lc)*chosen_sol+lb)*chosen_sol+y[maxed_at-1];
                    if (temp > maxed) { 
                        maxed=temp; 
                        x_max=chosen_sol+x[maxed_at-1];
                    }
                }
            }       
        } 
        // Repeating for the other side (not possible if the initial guess is the last spline point).
        if (maxed_at<num_pts-1) {
            const double& rd=d[maxed_at];
            const double& rc=c[maxed_at];
            const double& rb=b[maxed_at];
            solution sol_right=quad_solver(3*rd, 2*rc, rb);
            if (sol_right.solvable) {
                const double& chosen_sol=sol_right.sol1; // see arguments above.
                if (chosen_sol > 0 && chosen_sol < x[maxed_at+1]-x[maxed_at]) {
                    double temp=((rd*chosen_sol+rc)*chosen_sol+rb)*chosen_sol+y[maxed_at];
                    if (temp>maxed) { 
                        maxed=temp; 
                        x_max=chosen_sol+x[maxed_at];
                    }
                }
            }       
        }

        // Cleaning up; noting any divergences and incrementing count and pointers.
        if (maxed_at==0 || maxed_at==num_pts-1) { 
            if (diverged==0) { first=count; }
            ++diverged; 
        }
        ++count;
        *(out_ptr++)=x_max; 
    }

//    if (diverged>=1) { 
//        std::stringstream wrn;
//        wrn << diverged << " instances of divergence starting at row " << first << ".";
//        warning(wrn.str().c_str()); 
//    }

    // Cleaning up:
    delete [] x;
    delete [] y;
    delete [] b;
    delete [] c;
    delete [] d;
    delete [] l_ptrs;
    UNPROTECT(3);

    return(output);    
}


}
