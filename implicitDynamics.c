/*
Defining the important functions for the implicit dynamics equations

The main ones are:
    the system of odes for integrating
    boundary condition
    equation solving routine

*/


#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multiroots.h>

#include "lie_theory.h"
#include "implicitDynamics.h"
//#include "rkmk.h"


void printVectorHelp(const gsl_vector *v) {
    int size = v->size;

    printf("Vector:\n");

    int i;
    for (i=0;i<size;i++) {
        printf("%f\n",gsl_vector_get(v,i));
    }
    printf("\n");
}

void printMatrixHelp(const gsl_matrix *m) {
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

void rodDynamics(double s, double dt, const gsl_matrix *g, const gsl_vector *y, const gsl_vector *y_prev, const gsl_vector *q, gsl_vector *y_der) {
    //compute the space derivatives of xi and eta for being integrated with rkmk

    gsl_matrix *K = gsl_matrix_calloc(6,6);
    gsl_matrix *V = gsl_matrix_calloc(6,6);
    gsl_matrix *M = gsl_matrix_calloc(6,6);

    gsl_vector *xi_ref = gsl_vector_calloc(6);
    gsl_vector *dxi = gsl_vector_alloc(6);
    gsl_vector *grav = gsl_vector_calloc(3);

    gsl_vector *xi_dot = gsl_vector_alloc(6);
    gsl_vector *eta_dot = gsl_vector_alloc(6);
    gsl_vector *W_vis = gsl_vector_alloc(6);
    gsl_vector *W_grav = gsl_vector_alloc(6);
    gsl_vector *W_act = gsl_vector_calloc(6);
    gsl_vector *W = gsl_vector_calloc(6);
    gsl_vector *r = gsl_vector_alloc(3);
    gsl_vector *pa_der = gsl_vector_alloc(3);

    gsl_matrix *r_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *omega_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *pa_der_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *A = gsl_matrix_alloc(6,6);
    gsl_matrix *A1 = gsl_matrix_alloc(3,3);
    gsl_matrix *A2 = gsl_matrix_alloc(3,3);
    gsl_matrix *A3 = gsl_matrix_alloc(3,3);
    gsl_matrix *ad_xi = gsl_matrix_alloc(6,6);
    gsl_matrix *ad_eta = gsl_matrix_alloc(6,6);
    gsl_matrix *P = gsl_matrix_alloc(3,3);
    gsl_permutation *Perm = gsl_permutation_alloc(6);

    gsl_vector *temp1 = gsl_vector_alloc(3);
    gsl_vector *temp2 = gsl_vector_alloc(3);
    gsl_vector *temp3 = gsl_vector_alloc(3);
    gsl_vector *temp4 = gsl_vector_alloc(6);
    gsl_matrix *tempM1 = gsl_matrix_alloc(3,3);
    gsl_matrix *tempM2 = gsl_matrix_alloc(3,3);

    int sig; //not terribly important

    double E = 1e6;
    double G = E/3;
    double mu = 30000;
    double rho = 1000;
    //double L = 10e-2;
    double D = 1e-2;
    double Area = M_PI/4*pow(D,2);
    double I = M_PI/64*pow(D,4);
    double J = 2*I;

    int i,j,k;

    int n = q->size;

    //setup constants
    gsl_matrix_set(K, 0, 0, E*I);
    gsl_matrix_set(K, 1, 1, E*I);
    gsl_matrix_set(K, 2, 2, G*J);
    gsl_matrix_set(K, 3, 3, G*Area);
    gsl_matrix_set(K, 4, 4, G*Area);
    gsl_matrix_set(K, 5, 5, E*Area);

    gsl_matrix_memcpy(A,K);

    gsl_matrix_set(V, 0, 0, mu*3*I);
    gsl_matrix_set(V, 1, 1, mu*3*I);
    gsl_matrix_set(V, 2, 2, mu*J);
    gsl_matrix_set(V, 3, 3, mu*Area);
    gsl_matrix_set(V, 4, 4, mu*Area);
    gsl_matrix_set(V, 5, 5, mu*3*Area);

    gsl_matrix_set(M, 0, 0, rho*I);
    gsl_matrix_set(M, 1, 1, rho*I);
    gsl_matrix_set(M, 2, 2, rho*J);
    gsl_matrix_set(M, 3, 3, rho*Area);
    gsl_matrix_set(M, 4, 4, rho*Area);
    gsl_matrix_set(M, 5, 5, rho*Area);

    gsl_vector_set(xi_ref, 5, 1);
    gsl_vector_set(grav, 2, -9.81);

    //only computing derivatives at the current point


    //first step compute time derivatives
    //j = round(s/ds)
    for (i=0; i<6; i++) {
        gsl_vector_set(xi_dot, i, 2*(gsl_vector_get(y,i) - gsl_vector_get(y_prev, i))/dt);
        gsl_vector_set(eta_dot, i, 2*(gsl_vector_get(y,i+6) - gsl_vector_get(y_prev, i+6))/dt);
    }

    //external loads
    //viscosity
    gsl_blas_dgemv(CblasNoTrans, -1, V, xi_dot, 0, W_vis);

    //gravity
    gsl_matrix_const_view R = gsl_matrix_const_submatrix(g, 0, 0, 3, 3);
    gsl_blas_dgemv(CblasTrans, rho*Area, &R.matrix, grav, 0, temp1);
    for (i=0;i<3;i++) {
        gsl_vector_set(W_grav, i, 0);
        gsl_vector_set(W_grav, i+3, gsl_vector_get(temp1,i));
    }
    // printf("R\n");
    // printMatrixHelp(&R.matrix);

    //actuator
    gsl_vector_const_view nu = gsl_vector_const_subvector(y,3,3);
    gsl_vector_const_view omega = gsl_vector_const_subvector(y,0,3);
    gsl_vector_const_view xi = gsl_vector_const_subvector(y, 0, 6);
    gsl_vector_const_view eta = gsl_vector_const_subvector(y, 6, 6);

    skew(&omega.vector, omega_skew);
    for (i=0; i<n; i++) {
        gsl_vector_set(r, 0, D/4*cos(2*M_PI*i/n));
        gsl_vector_set(r, 1, D/4*sin(2*M_PI*i/n));
        gsl_vector_set(r, 2, 0);

        skew(r, r_skew);

        //nu+skew(omega)*r
        gsl_vector_memcpy(temp1, &nu.vector);
        gsl_blas_dgemv(CblasNoTrans, 1, omega_skew, r, 1.0, temp1);

        //pa_der
        gsl_blas_dgemv(CblasNoTrans, 1, &R.matrix, temp1, 0.0, pa_der);
        // printf("pa_der\n");
        // printVectorHelp(pa_der);

        //P
        skew(pa_der,pa_der_skew);

        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1/pow(gsl_blas_dnrm2(pa_der),3), &R.matrix, pa_der_skew, 0.0, tempM1);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, tempM1, pa_der_skew, 0.0, tempM2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, tempM2, &R.matrix, 0.0, P);

        //A = q(i)*[-A1,A2;-A3,A3], A3 = Pskew(r), A2 = skew(r)P, A1 = skew(r)*A3
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, gsl_vector_get(q,i), P, r_skew, 0, A3);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, gsl_vector_get(q,i), r_skew, P, 0, A2);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, r_skew, A3, 0, A1);
        // printf("A1\n");
        // printMatrixHelp(A1);
        // printf("A2\n");
        // printMatrixHelp(A2);
        // printf("A3\n");
        // printMatrixHelp(A3);
        for (j=0; j<3; j++) {
            for (k=0; k<3; k++) {
                gsl_matrix_set(A, j, k, gsl_matrix_get(A,j,k) - gsl_matrix_get(A1, j, k));
                gsl_matrix_set(A, j, k+3, gsl_matrix_get(A,j,k+3) + gsl_matrix_get(A2, j, k));
                gsl_matrix_set(A, j+3, k, gsl_matrix_get(A,j+3,k) - gsl_matrix_get(A3, j, k));
                gsl_matrix_set(A, j+3, k+3, gsl_matrix_get(A,j+3,k+3 ) + gsl_matrix_get(A3, j, k));
            }
        }

        //W_act
        gsl_blas_dgemv(CblasNoTrans, 1, omega_skew, temp1, 0, temp2);
        gsl_blas_dgemv(CblasNoTrans, gsl_vector_get(q, i), P, temp2, 0, temp3);
        gsl_blas_dgemv(CblasNoTrans, 1, r_skew, temp3, 0, temp2);
        for (j=0;j<3;j++) {
            gsl_vector_set(W_act, j, gsl_vector_get(W_act, j) + gsl_vector_get(temp2,j));
            gsl_vector_set(W_act, j+3, gsl_vector_get(W_act, j+3) + gsl_vector_get(temp3,j));
        }

        //A\stuff
        adjointSE3(&xi.vector, ad_xi);
        adjointSE3(&eta.vector, ad_eta);

        for (j=0;j<6;j++) {
            gsl_vector_set(dxi, j, gsl_vector_get(&xi.vector, j)-gsl_vector_get(xi_ref,j));
            gsl_vector_set(W, j, gsl_vector_get(W_grav, j) + gsl_vector_get(W_vis, j) + gsl_vector_get(W_act, j));
        }

        gsl_blas_dgemv(CblasNoTrans, 1, M, eta_dot, -1, W);
        gsl_blas_dgemv(CblasNoTrans, 1, M, &eta.vector, 0, temp4);
        gsl_blas_dgemv(CblasTrans, -1, ad_eta, temp4, 1, W);
        gsl_blas_dgemv(CblasNoTrans, 1, K, dxi, 0, temp4);
        gsl_blas_dgemv(CblasTrans, 1, ad_xi, temp4, 1, W);


    }
    // printf("A\n");
    // printMatrixHelp(A);
    gsl_linalg_LU_decomp(A, Perm, &sig);
    gsl_linalg_LU_solve(A, Perm, W, temp4);

    gsl_blas_dgemv(CblasNoTrans, -1, ad_xi, &eta.vector, 1, xi_dot);
    for (i=0; i<6; i++) {
        gsl_vector_set(y_der, i, gsl_vector_get(temp4, i));
        gsl_vector_set(y_der, i+6, gsl_vector_get(xi_dot, i));
    }


    gsl_matrix_free(K);
    gsl_matrix_free(V);
    gsl_matrix_free(M);

    gsl_vector_free(xi_ref);
    gsl_vector_free(dxi);
    gsl_vector_free(grav);

    gsl_vector_free(xi_dot);
    gsl_vector_free(eta_dot);
    gsl_vector_free(W_vis);
    gsl_vector_free(W_grav);
    gsl_vector_free(W_act);
    gsl_vector_free(W);
    gsl_vector_free(r);
    gsl_vector_free(pa_der);

    gsl_matrix_free(r_skew);
    gsl_matrix_free(omega_skew);
    gsl_matrix_free(pa_der_skew);
    gsl_matrix_free(A);
    gsl_matrix_free(A1);
    gsl_matrix_free(A2);
    gsl_matrix_free(A3);
    gsl_matrix_free(ad_xi);
    gsl_matrix_free(ad_eta);
    gsl_matrix_free(P);
    gsl_permutation_free(Perm);

    gsl_vector_free(temp1);
    gsl_vector_free(temp2);
    gsl_vector_free(temp3);
    gsl_vector_free(temp4);
    gsl_matrix_free(tempM1);
    gsl_matrix_free(tempM2);

}

//just hardcoding material parameters and such for now
//unfortunately not quite the right form for rkmk and don't really want to try and generalize it, so just manually do the integration
//also slightly easier to pass previous values in
void integrate(double dt, const gsl_matrix *g, const gsl_matrix *y, const gsl_vector *q, const gsl_vector *xi0, gsl_matrix *g_next, gsl_matrix *y_next) {
    // printf("xi0\n");
    // printVectorHelp(xi0);
    int N = g->size2;

    gsl_vector *y_der = gsl_vector_alloc(12);
    gsl_matrix *y_half = gsl_matrix_calloc(12,N);
    gsl_matrix *g_half = gsl_matrix_calloc(12,N);
    gsl_vector *g_in = gsl_vector_calloc(12);
    gsl_vector *y_in = gsl_vector_alloc(12);
    gsl_matrix *g_curr = gsl_matrix_alloc(4,4);
    gsl_matrix *g_temp = gsl_matrix_alloc(4,4);
    gsl_matrix *y_exp = gsl_matrix_alloc(4,4);
    gsl_vector *eta_temp = gsl_vector_alloc(6);
    gsl_vector *y_temp = gsl_vector_alloc(12);

    double L = 10e-2;
    double ds = L/(N-1);

    int i;

    //setup initial conditions
    for (i=0;i<6;i++) {
        // gsl_matrix_set(y_next, i, 0, gsl_vector_get(xi0,i));
        // gsl_matrix_set(y_next, i+6, 0, 0);

        gsl_vector_set(y_in, i, (gsl_vector_get(xi0,i) + gsl_matrix_get(y,i,0))/2); //xi
        gsl_vector_set(y_in, i+6, 0); //eta
    }
    gsl_matrix_set_col(y_half, 0, y_in);

    gsl_vector_set(g_in, 0, 1);
    gsl_vector_set(g_in, 4, 1);
    gsl_vector_set(g_in, 8, 1);
    gsl_matrix_set_col(g_half,0,g_in);

    //gsl_matrix_set_col(g_next, 0, g_in);


    //integrate at the half step
    for (i=1; i<N; i++) {
        //grab states for input
        unflattenConfig(g_in, g_curr);
        gsl_vector_const_view y_prev = gsl_matrix_const_column(y, i-1);
        gsl_vector_const_view y_in = gsl_matrix_const_column(y_half,i-1);
        gsl_vector_const_view g_in = gsl_matrix_const_column(g_half, i-1);

        unflattenConfig(&g_in.vector, g_curr);

        //derivative
        rodDynamics(ds*i, dt, g_curr, &y_in.vector, &y_prev.vector, q, y_der);

        //integration step
        gsl_vector_view y_step = gsl_matrix_column(y_half,i); //where y gets put
        gsl_vector_view g_step = gsl_matrix_column(g_half,i);

        //y_next = y + ds*y_der
        gsl_vector_memcpy(&y_step.vector, &y_in.vector);
        gsl_vector_scale(y_der,ds);
        gsl_vector_add(&y_step.vector,y_der);

        //g_next = g*expm(ds*y_in)
        gsl_vector_memcpy(y_temp, &y_in.vector);
        gsl_vector_scale(y_temp, ds);
        expSE3(y_temp, y_exp);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, g_curr, y_exp, 0, g_temp);
        // printMatrixHelp(g_curr);
        // printMatrixHelp(y_exp);
        // printMatrixHelp(g_temp);
        flattenConfig(g_temp, &g_step.vector);


        // gsl_vector_memcpy(y_in, &y_step.vector); //step to next state
        // gsl_vector_memcpy(g_in, &g_step.vector);
    }
    // printMatrixHelp(g_half);
    // printMatrixHelp(y_half);

    //g_i+1 = g_i*exp(dt*eta_half)
    for (i=0; i<N; i++) {
        gsl_vector_const_view g_i = gsl_matrix_const_column(g, i);
        gsl_vector_const_view eta_half = gsl_matrix_const_subcolumn(y_half, i, 6, 6);

        unflattenConfig(&g_i.vector, g_curr);
        gsl_vector_memcpy(eta_temp, &eta_half.vector);
        gsl_vector_scale(eta_temp, dt);

        expSE3(eta_temp, y_exp);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, g_curr, y_exp, 0, g_temp);
        gsl_vector_view g_step = gsl_matrix_column(g_next, i);
        flattenConfig(g_temp, &g_step.vector);
    }

    //y_(i+1) = 2*y_half-y
    gsl_matrix_scale(y_half,2);
    gsl_matrix_memcpy(y_next, y_half);
    gsl_matrix_sub(y_next,y);

    gsl_vector_free(y_der);
    gsl_matrix_free(y_half);
    gsl_matrix_free(g_half);
    gsl_vector_free(g_in);
    gsl_vector_free(y_in);
    gsl_matrix_free(g_curr);
    gsl_matrix_free(g_temp);
    gsl_matrix_free(y_exp);
    gsl_vector_free(eta_temp);
    gsl_vector_free(y_temp);


}



void tipCondition(const gsl_matrix *g, const gsl_matrix *y, const gsl_vector *q, gsl_vector *condition) {
    //checks that the wrench balance at the tip is met
    //takes the already integrated values

    int N = y->size2;
    int n = q->size;

    gsl_vector *W = gsl_vector_calloc(6);
    gsl_vector *r = gsl_vector_alloc(3);
    gsl_vector *ta = gsl_vector_alloc(3);
    gsl_vector *xi_ref = gsl_vector_calloc(6);
    gsl_vector *dxi = gsl_vector_alloc(6);

    gsl_matrix *r_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *K = gsl_matrix_calloc(6,6);

    gsl_vector *temp1 = gsl_vector_calloc(3);

    double E = 1e6;
    double G = E/3;
    double D = 1e-2;
    double A = M_PI/4*pow(D,2);
    double I = M_PI/64*pow(D,4);
    double J = 2*I;

    gsl_vector_set(xi_ref,5,1);

    gsl_matrix_set(K, 0, 0, E*I);
    gsl_matrix_set(K, 1, 1, E*I);
    gsl_matrix_set(K, 2, 2, G*J);
    gsl_matrix_set(K, 3, 3, G*A);
    gsl_matrix_set(K, 4, 4, G*A);
    gsl_matrix_set(K, 5, 5, E*A);

    int i,k;
    gsl_vector_const_view nu = gsl_matrix_const_subcolumn(y, N-1, 3, 3);
    gsl_vector_const_view omega = gsl_matrix_const_subcolumn(y, N-1, 0, 3);
    // printMatrixHelp(y);
    for (i=0; i<n; i++) {
        gsl_vector_set(r, 0, D/4*cos(2*M_PI*i/n));
        gsl_vector_set(r, 1, D/4*sin(2*M_PI*i/n));
        gsl_vector_set(r, 2, 0);
        // printVectorHelp(r);

        skew(r, r_skew);

        gsl_vector_memcpy(ta, &nu.vector);

        gsl_blas_dgemv(CblasNoTrans, -1, r_skew, &omega.vector, 1, ta);
        gsl_blas_dgemv(CblasNoTrans, 1, r_skew, ta, 0, temp1);
        // printVectorHelp(ta);
        for (k=0; k<3;k++) {
            gsl_vector_set(W, k, gsl_vector_get(W,k) + gsl_vector_get(q,i)*gsl_vector_get(temp1,k)/gsl_blas_dnrm2(ta));
            gsl_vector_set(W, k+3, gsl_vector_get(W,k+3) + gsl_vector_get(q,i)*gsl_vector_get(ta,k)/gsl_blas_dnrm2(ta));
        }


    }
    // printVectorHelp(W);
    gsl_vector_const_view xi = gsl_matrix_const_subcolumn(y, N-1, 0, 6);
    gsl_vector_memcpy(dxi, &xi.vector);
    gsl_vector_sub(dxi,xi_ref);
    gsl_vector_memcpy(condition, W);
    gsl_blas_dgemv(CblasNoTrans, 1, K, dxi, -1, condition);


    gsl_vector_free(W);
    gsl_vector_free(r);
    gsl_vector_free(ta);
    gsl_vector_free(xi_ref);
    gsl_vector_free(dxi);

    gsl_matrix_free(r_skew);
    gsl_matrix_free(K);

    gsl_vector_free(temp1);

}

// struct dynParams {
//
//     double dt;
//     const gsl_matrix *g;
//     gsl_matrix *g_next;
//     const gsl_matrix *y;
//     gsl_matrix *y_next;
//
//     const gsl_vector *q;
//
//
// };

int dynamics(const gsl_vector *xi0, void *params, gsl_vector *condition) {
    //takes xi0 and the previous steps values of g and y
    //integrates them
    //and computes the boundary condition

    //first pull out stuff from params and initialize matrices
    double dt = ((struct dynParams *) params)->dt;
    const gsl_matrix *g = ((struct dynParams *) params)->g;
    const gsl_matrix *y = ((struct dynParams *) params)->y;
    const gsl_vector *q = ((struct dynParams *) params)->q;

    gsl_matrix *g_next = ((struct dynParams *) params)->g_next;
    gsl_matrix *y_next = ((struct dynParams *) params)->y_next;


    //integration
    integrate(dt, g, y, q, xi0, g_next, y_next);

    //boundary
    tipCondition(g_next, y_next, q, condition);

    return GSL_SUCCESS;

}

void solver(double dt, const gsl_matrix *g, const gsl_matrix *y, const gsl_vector *q, gsl_matrix *g_next, gsl_matrix *y_next) {
    //solves for the next state

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    struct dynParams p = {dt, g, g_next, y, y_next, q};

    int status;
    size_t i, iter = 0;
    const size_t n = 6; //always 6 things to solve for

    gsl_vector *xi0 = gsl_vector_alloc(6);

    gsl_multiroot_function f = {&dynamics, n, &p};

    //initialize guess to last steps value
    for (i=0; i<6; i++) {
        gsl_vector_set(xi0, i, gsl_matrix_get(y, i, 0));
    }

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T, 6);
    gsl_multiroot_fsolver_set(s, &f, xi0);

    do
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        // printf("status: %s\n", gsl_strerror (status));
        if (status)
            break;
        // printf("Iteration: %d\n",iter);
        // printVectorHelp(s->x);
        // printVectorHelp(s->f);

        status = gsl_multiroot_test_residual(s->f, 1e-10);
        // printf("status: %s\n", gsl_strerror (status));
    }
    while (status == GSL_CONTINUE && iter<1000);

    gsl_vector_memcpy(xi0, s->x);
    integrate(dt, g, y, q, xi0, g_next, y_next);


    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(xi0);

}


//
// int main(void) {
//     //just a test
//
//     double dt = 0.1;
//     int N = 10;
//     gsl_matrix *g = gsl_matrix_calloc(12,N);
//     gsl_matrix *y = gsl_matrix_calloc(12,N);
//     gsl_matrix *g_next = gsl_matrix_alloc(12,N);
//     gsl_matrix *y_next = gsl_matrix_alloc(12,N);
//     gsl_vector *q = gsl_vector_alloc(3);
//
//     gsl_vector *condition = gsl_vector_alloc(6);
//     gsl_vector *xi0 = gsl_vector_calloc(6);
//     gsl_vector_set(xi0,5,1);
//
//     int i;
//     for (i=0; i<N; i++) {
//         gsl_matrix_set(y, 5, i, 1);
//         gsl_matrix_set(g, 0, i, 1);
//         gsl_matrix_set(g, 4, i, 1);
//         gsl_matrix_set(g, 8, i, 1);
//         gsl_matrix_set(g, 11, i, 10e-2/(N-1)*i);
//     }
//
//     gsl_vector_set(q, 0, -0.5);
//     gsl_vector_set(q, 1, -0);
//     gsl_vector_set(q, 2, -0);
//
//     printMatrixHelp(g);
//     printMatrixHelp(y);
//
//     solver(dt, g, y, q, g_next, y_next);
//     // struct dynParams p = {dt, g, g_next, y, y_next, q};
//     // dynamics(xi0, &p, condition);
//     // integrate(dt, g, y, q, xi0, g_next, y_next);
//
//     printMatrixHelp(g_next);
//     printMatrixHelp(y_next);
//
//     gsl_matrix_free(g);
//     gsl_matrix_free(y);
//     gsl_matrix_free(g_next);
//     gsl_matrix_free(y_next);
//     gsl_vector_free(q);
//
//     gsl_vector_free(condition);
//     gsl_vector_free(xi0);
// }
