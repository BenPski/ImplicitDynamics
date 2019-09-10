/*
This does the integration and boundary condition computation for the implicit scheme without the equation solving
*/

#include "mex.h"

#include "implicitDynamics.h"
#include "lie_theory.h"

#include <math.h>
#include <gsl/gsl_matrix.h>


void printVector(const gsl_vector *v) {
    int size = v->size;

    printf("Vector:\n");

    int i;
    for (i=0;i<size;i++) {
        printf("%f\n",gsl_vector_get(v,i));
    }
    printf("\n");
}

void printMatrix(const gsl_matrix *m) {
    int height = m->size1;
    int width = m->size2;

    printf("Matrix: \n");

    int i,j;
    for (i=0;i<height;i++) {
        for (j=0;j<width-1;j++) {
            printf("%f ", gsl_matrix_get(m,i,j));
        }
        printf("%f\n", gsl_matrix_get(m,i,width-1));
    }
    printf("\n");
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /*
    Interface to the stepdynamics function for the variaitonal lie integrator of the beam

    for now just passing g, xi ,eta, mu, and lambda back and forth
    leaving the properties in the step dynamics function

    matlab is transposed

    */

    int i,j;

    int n; //assuming everything is the right size
    const mwSize *dim;
    dim = mxGetDimensions(prhs[0]);
    n = (int)dim[1];
    // printf("%d\n", n);
    double *g_in;
    g_in = mxGetPr(prhs[0]);
    gsl_matrix *g = gsl_matrix_alloc(12,n);
    for (i=0; i<12; i++) {
        for (j=0; j<n; j++) {
            gsl_matrix_set(g, i, j, g_in[i+12*j]);
        }
    }
    // printMatrix(g);



    double *xi_in;
    double *eta_in;
    xi_in = mxGetPr(prhs[1]);
    eta_in = mxGetPr(prhs[2]);
    gsl_matrix *y = gsl_matrix_alloc(12,n);
    for (i=0; i<6; i++) {
        for (j=0; j<n; j++) {
            gsl_matrix_set(y, i, j, xi_in[i+6*j]); //xi
            gsl_matrix_set(y, i+6, j, eta_in[i+6*j]); //eta
        }
    }
    // printMatrix(y);

    double *q_in;
    q_in = mxGetPr(prhs[3]);
    const mwSize *q_dim;
    q_dim = mxGetDimensions(prhs[3]);
    int N_act = (int)(q_dim[1]*q_dim[0]);
    gsl_vector *q = gsl_vector_alloc(N_act);
    for (i=0; i<N_act; i++) {
        gsl_vector_set(q, i, q_in[i]);
    }

    double *xi0_in;
    xi0_in = mxGetPr(prhs[4]);
    gsl_vector *xi0 = gsl_vector_alloc(6);
    for (i=0; i<6; i++) {
        gsl_vector_set(xi0, i, xi0_in[i]);
    }

    double dt = mxGetPr(prhs[5])[0];




    gsl_matrix *g_next = gsl_matrix_alloc(12,n);
    gsl_matrix *y_next = gsl_matrix_alloc(12,n);
    gsl_vector *condition = gsl_vector_alloc(6);


    integrate(dt, g, y, q, xi0, g_next, y_next);
    tipCondition(g_next, y_next, q, condition);



    //put it all back in matlab
    double *g_out_vec, *xi_out_vec, *eta_out_vec, *condition_out_vec;
    mxArray *g_out, *xi_out, *eta_out, *condition_out;
    g_out = plhs[0] = mxCreateDoubleMatrix(12,n,mxREAL);
    g_out_vec = mxGetPr(g_out);
    xi_out = plhs[1] = mxCreateDoubleMatrix(6,n,mxREAL);
    xi_out_vec = mxGetPr(xi_out);
    eta_out = plhs[2] = mxCreateDoubleMatrix(6,n,mxREAL);
    eta_out_vec = mxGetPr(eta_out);
    condition_out = plhs[3] = mxCreateDoubleMatrix(6,1,mxREAL);
    condition_out_vec = mxGetPr(condition_out);


    for (i=0;i<n;i++) {
        for (j=0; j<12; j++) {
            g_out_vec[i*12+j] = gsl_matrix_get(g_next,j,i);
        }
        for (j=0;j<6;j++) {
            xi_out_vec[i*6+j] = gsl_matrix_get(y_next, j, i);
            eta_out_vec[i*6+j] = gsl_matrix_get(y_next, j+6, i);

        }
    }
    for (i=0;i<6;i++) {
        condition_out_vec[i] = gsl_vector_get(condition, i);
    }

    gsl_vector_free(q);

    gsl_matrix_free(g);
    gsl_matrix_free(y);
    gsl_matrix_free(g_next);
    gsl_matrix_free(y_next);
    gsl_vector_free(xi0);
    gsl_vector_free(condition);


}
