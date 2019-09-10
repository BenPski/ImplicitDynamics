/*
Some useful functions for lie theory

Main focus is the exponenetial map and cayley transform
defined for SO(3) and SE(3)

Work with matrices and aim to get in operator form

prefer the vector form of the algebra
*/

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


/*
Other convenience functions
*/
void flattenConfig(const gsl_matrix *in, gsl_vector *out) {
    int i, j;
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            gsl_vector_set(out, i+j*3, gsl_matrix_get(in,i,j));
        }
        gsl_vector_set(out, 9+i, gsl_matrix_get(in,i,3));
    }
}

void unflattenConfig(const gsl_vector *in, gsl_matrix *out) {
    int i, j;
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            gsl_matrix_set(out, i, j, gsl_vector_get(in,i+j*3));
        }
        gsl_matrix_set(out, i, 3, gsl_vector_get(in,9+i));
    }
    gsl_matrix_set(out, 3,0,0);
    gsl_matrix_set(out, 3,1,0);
    gsl_matrix_set(out, 3,2,0);
    gsl_matrix_set(out, 3,3,1);
}



/*
First the algebra isomorphisms
*/

void skew(const gsl_vector *in, gsl_matrix *out) {
    gsl_matrix_set(out, 0, 0, 0);
    gsl_matrix_set(out, 1, 1, 0);
    gsl_matrix_set(out, 2, 2, 0);

    gsl_matrix_set(out, 0, 1, -gsl_vector_get(in, 2));
    gsl_matrix_set(out, 0, 2, gsl_vector_get(in, 1));
    gsl_matrix_set(out, 1, 0, gsl_vector_get(in, 2));
    gsl_matrix_set(out, 1, 2, -gsl_vector_get(in, 0));
    gsl_matrix_set(out, 2, 0, -gsl_vector_get(in, 1));
    gsl_matrix_set(out, 2, 1, gsl_vector_get(in, 0));
}

void unskew(const gsl_matrix *in, gsl_vector *out) {
    gsl_vector_set(out, 0, gsl_matrix_get(in, 2, 1));
    gsl_vector_set(out, 1, gsl_matrix_get(in, 0, 2));
    gsl_vector_set(out, 2, gsl_matrix_get(in, 1, 0));
}

void se(const gsl_vector *in, gsl_matrix *out) {
    //convenient, though a bit of a hack
    skew(in, out);
    gsl_matrix_set(out, 0, 3, gsl_vector_get(in, 3));
    gsl_matrix_set(out, 1, 3, gsl_vector_get(in, 4));
    gsl_matrix_set(out, 2, 3, gsl_vector_get(in, 5));
    gsl_matrix_set(out, 3, 0, 0);
    gsl_matrix_set(out, 3, 1, 0);
    gsl_matrix_set(out, 3, 2, 0);
    gsl_matrix_set(out, 3, 3, 1);
}

void unse(const gsl_matrix *in, gsl_vector *out) {
    unskew(in, out);
    gsl_vector_set(out, 3, gsl_matrix_get(in, 0, 3));
    gsl_vector_set(out, 4, gsl_matrix_get(in, 1, 3));
    gsl_vector_set(out, 5, gsl_matrix_get(in, 2, 3));
}


/*
SO(3) exponential
*/

void expSO3(const gsl_vector *in, gsl_matrix *out) {
    /*
    unrolling for slight efficiency and its easy anyways
    [      ((cos(t) - 1)*(b^2 + c^2))/t^2 + 1, - (c*sin(t))/t - (a*b*(cos(t) - 1))/t^2,   (b*sin(t))/t - (a*c*(cos(t) - 1))/t^2]
    [   (c*sin(t))/t - (a*b*(cos(t) - 1))/t^2,      ((cos(t) - 1)*(a^2 + c^2))/t^2 + 1, - (a*sin(t))/t - (b*c*(cos(t) - 1))/t^2]
    [ - (b*sin(t))/t - (a*c*(cos(t) - 1))/t^2,   (a*sin(t))/t - (b*c*(cos(t) - 1))/t^2,      ((cos(t) - 1)*(a^2 + b^2))/t^2 + 1]

    */

    double a = gsl_vector_get(in,0);
    double b = gsl_vector_get(in,1);
    double c = gsl_vector_get(in,2);
    double t = sqrt(a*a+b*b+c*c);
    if (t == 0) {
        gsl_matrix_set(out, 0, 0, 1);
        gsl_matrix_set(out, 1, 1, 1);
        gsl_matrix_set(out, 2, 2, 1);

        gsl_matrix_set(out, 0, 1, 0);
        gsl_matrix_set(out, 0, 2, 0);
        gsl_matrix_set(out, 1, 0, 0);
        gsl_matrix_set(out, 1, 2, 0);
        gsl_matrix_set(out, 2, 0, 0);
        gsl_matrix_set(out, 2, 1, 0);
    } else {
        gsl_matrix_set(out, 0, 0, ((cos(t) - 1)*(pow(b,2) + pow(c,2)))/pow(t,2) + 1);
        gsl_matrix_set(out, 0, 1, - (c*sin(t))/t - (a*b*(cos(t) - 1))/pow(t,2));
        gsl_matrix_set(out, 0, 2, (b*sin(t))/t - (a*c*(cos(t) - 1))/pow(t,2));
        gsl_matrix_set(out, 1, 0, (c*sin(t))/t - (a*b*(cos(t) - 1))/pow(t,2));
        gsl_matrix_set(out, 1, 1, ((cos(t) - 1)*(pow(a,2) + pow(c,2)))/pow(t,2) + 1);
        gsl_matrix_set(out, 1, 2, - (a*sin(t))/t - (b*c*(cos(t) - 1))/pow(t,2));
        gsl_matrix_set(out, 2, 0, - (b*sin(t))/t - (a*c*(cos(t) - 1))/pow(t,2));
        gsl_matrix_set(out, 2, 1, (a*sin(t))/t - (b*c*(cos(t) - 1))/pow(t,2));
        gsl_matrix_set(out, 2, 2, ((cos(t) - 1)*(pow(a,2) + pow(b,2)))/pow(t,2) + 1);
    }
}

void tanSO3(const gsl_vector *in, gsl_matrix *out) {
    /*
    [              ((sin(t)/t - 1)*(b^2 + c^2))/t^2 + 1, - (c*(cos(t) - 1))/t^2 - (a*b*(sin(t)/t - 1))/t^2,   (b*(cos(t) - 1))/t^2 - (a*c*(sin(t)/t - 1))/t^2]
    [   (c*(cos(t) - 1))/t^2 - (a*b*(sin(t)/t - 1))/t^2,              ((sin(t)/t - 1)*(a^2 + c^2))/t^2 + 1, - (a*(cos(t) - 1))/t^2 - (b*c*(sin(t)/t - 1))/t^2]
    [ - (b*(cos(t) - 1))/t^2 - (a*c*(sin(t)/t - 1))/t^2,   (a*(cos(t) - 1))/t^2 - (b*c*(sin(t)/t - 1))/t^2,              ((sin(t)/t - 1)*(a^2 + b^2))/t^2 + 1]
    */
    double a = gsl_vector_get(in,0);
    double b = gsl_vector_get(in,1);
    double c = gsl_vector_get(in,2);
    double t = sqrt(a*a+b*b+c*c);
    if (t == 0) {
        gsl_matrix_set(out, 0, 0, 1);
        gsl_matrix_set(out, 1, 1, 1);
        gsl_matrix_set(out, 2, 2, 1);

        gsl_matrix_set(out, 0, 1, 0);
        gsl_matrix_set(out, 0, 2, 0);
        gsl_matrix_set(out, 1, 0, 0);
        gsl_matrix_set(out, 1, 2, 0);
        gsl_matrix_set(out, 2, 0, 0);
        gsl_matrix_set(out, 2, 1, 0);
    } else {
        gsl_matrix_set(out, 0, 0, ((sin(t)/t - 1)*(pow(b,2) + pow(c,2)))/pow(t,2) + 1);
        gsl_matrix_set(out, 0, 1, - (c*(cos(t) - 1))/pow(t,2) - (a*b*(sin(t)/t - 1))/pow(t,2));
        gsl_matrix_set(out, 0, 2, (b*(cos(t) - 1))/pow(t,2) - (a*c*(sin(t)/t - 1))/pow(t,2));
        gsl_matrix_set(out, 1, 0, (c*(cos(t) - 1))/pow(t,2) - (a*b*(sin(t)/t - 1))/pow(t,2));
        gsl_matrix_set(out, 1, 1, ((sin(t)/t - 1)*(pow(a,2) + pow(c,2)))/pow(t,2) + 1);
        gsl_matrix_set(out, 1, 2, - (a*(cos(t) - 1))/pow(t,2) - (b*c*(sin(t)/t - 1))/pow(t,2));
        gsl_matrix_set(out, 2, 0, - (b*(cos(t) - 1))/pow(t,2) - (a*c*(sin(t)/t - 1))/pow(t,2));
        gsl_matrix_set(out, 2, 1, (a*(cos(t) - 1))/pow(t,2) - (b*c*(sin(t)/t - 1))/pow(t,2));
        gsl_matrix_set(out, 2, 2, ((sin(t)/t - 1)*(pow(a,2) + pow(b,2)))/pow(t,2) + 1);
    }
}

void itanSO3(const gsl_vector *in, gsl_matrix *out) {
    /*
    [ 1 - (((t*sin(t))/(2*cos(t) - 2) + 1)*(b^2 + c^2))/t^2,       (a*b*((t*sin(t))/(2*cos(t) - 2) + 1))/t^2 - c/2,       b/2 + (a*c*((t*sin(t))/(2*cos(t) - 2) + 1))/t^2]
    [       c/2 + (a*b*((t*sin(t))/(2*cos(t) - 2) + 1))/t^2, 1 - (((t*sin(t))/(2*cos(t) - 2) + 1)*(a^2 + c^2))/t^2,       (b*c*((t*sin(t))/(2*cos(t) - 2) + 1))/t^2 - a/2]
    [       (a*c*((t*sin(t))/(2*cos(t) - 2) + 1))/t^2 - b/2,       a/2 + (b*c*((t*sin(t))/(2*cos(t) - 2) + 1))/t^2, 1 - (((t*sin(t))/(2*cos(t) - 2) + 1)*(a^2 + b^2))/t^2]
    */
    double a = gsl_vector_get(in,0);
    double b = gsl_vector_get(in,1);
    double c = gsl_vector_get(in,2);
    double t = sqrt(a*a+b*b+c*c);
    if (t == 0) {
        gsl_matrix_set(out, 0, 0, 1);
        gsl_matrix_set(out, 1, 1, 1);
        gsl_matrix_set(out, 2, 2, 1);

        gsl_matrix_set(out, 0, 1, 0);
        gsl_matrix_set(out, 0, 2, 0);
        gsl_matrix_set(out, 1, 0, 0);
        gsl_matrix_set(out, 1, 2, 0);
        gsl_matrix_set(out, 2, 0, 0);
        gsl_matrix_set(out, 2, 1, 0);
    } else {
        gsl_matrix_set(out, 0, 0, 1 - (((t*sin(t))/(2*cos(t) - 2) + 1)*(pow(b,2) + pow(c,2)))/pow(t,2));
        gsl_matrix_set(out, 0, 1, (a*b*((t*sin(t))/(2*cos(t) - 2) + 1))/pow(t,2) - c/2);
        gsl_matrix_set(out, 0, 2, b/2 + (a*c*((t*sin(t))/(2*cos(t) - 2) + 1))/pow(t,2));
        gsl_matrix_set(out, 1, 0, c/2 + (a*b*((t*sin(t))/(2*cos(t) - 2) + 1))/pow(t,2));
        gsl_matrix_set(out, 1, 1, 1 - (((t*sin(t))/(2*cos(t) - 2) + 1)*(pow(a,2) + pow(c,2)))/pow(t,2));
        gsl_matrix_set(out, 1, 2, (b*c*((t*sin(t))/(2*cos(t) - 2) + 1))/pow(t,2) - a/2);
        gsl_matrix_set(out, 2, 0, (a*c*((t*sin(t))/(2*cos(t) - 2) + 1))/pow(t,2) - b/2);
        gsl_matrix_set(out, 2, 1, a/2 + (b*c*((t*sin(t))/(2*cos(t) - 2) + 1))/pow(t,2));
        gsl_matrix_set(out, 2, 2, 1 - (((t*sin(t))/(2*cos(t) - 2) + 1)*(pow(a,2) + pow(b,2)))/pow(t,2));
    }
}

void logSO3(const gsl_matrix *in, gsl_vector *out) {
    double t = acos((gsl_matrix_get(in,0,0)+gsl_matrix_get(in,1,1)+gsl_matrix_get(in,2,2)-1)/2);
    if (t==0) {
        gsl_vector_set(out,0,0);
        gsl_vector_set(out,1,0);
        gsl_vector_set(out,2,0);
    } else {
        gsl_vector_set(out,0,t*(gsl_matrix_get(in,2,1) - gsl_matrix_get(in,1,2))/(2*sin(t)));
        gsl_vector_set(out,1,t*(gsl_matrix_get(in,0,2) - gsl_matrix_get(in,2,0))/(2*sin(t)));
        gsl_vector_set(out,2,t*(gsl_matrix_get(in,1,0) - gsl_matrix_get(in,0,1))/(2*sin(t)));
    }
}

/*
SE(3) exponential
*/

void expSE3(const gsl_vector *in, gsl_matrix *out) {
    expSO3(in,out);
    gsl_vector *p = gsl_vector_alloc(3);
    gsl_matrix *T = gsl_matrix_alloc(3,3);

    gsl_vector_const_view h_w = gsl_vector_const_subvector(in,0,3);
    gsl_vector_const_view h_v = gsl_vector_const_subvector(in,3,3);

    tanSO3(&h_w.vector, T);
    gsl_blas_dgemv(CblasTrans, 1.0, T, &h_v.vector, 0.0, p);

    gsl_matrix_set(out, 0, 3, gsl_vector_get(p,0));
    gsl_matrix_set(out, 1, 3, gsl_vector_get(p,1));
    gsl_matrix_set(out, 2, 3, gsl_vector_get(p,2));

    gsl_matrix_set(out, 3, 0, 0);
    gsl_matrix_set(out, 3, 1, 0);
    gsl_matrix_set(out, 3, 2, 0);
    gsl_matrix_set(out, 3, 3, 1);

    gsl_vector_free(p);
    gsl_matrix_free(T);
}

void logSE3(const gsl_matrix *in, gsl_vector *out) {
    logSO3(in,out);

    gsl_matrix *T = gsl_matrix_alloc(3,3);
    gsl_vector *x = gsl_vector_alloc(3);

    gsl_vector_const_view p = gsl_matrix_const_subcolumn(in, 3, 0, 3);
    gsl_vector_const_view w = gsl_vector_const_subvector(out, 0, 3);
    itanSO3(&w.vector, T);
    gsl_blas_dgemv(CblasTrans, 1.0, T, &p.vector, 0.0, x);

    gsl_vector_set(out, 3, gsl_vector_get(x,0));
    gsl_vector_set(out, 4, gsl_vector_get(x,1));
    gsl_vector_set(out, 5, gsl_vector_get(x,2));

    gsl_matrix_free(T);
    gsl_vector_free(x);
}

void tanSE3(const gsl_vector *in, gsl_matrix *out) {
    /*
    [tanSO3(w), T]
    [0, tanSO3(w)]

    T=
    [                                                                                ((b^2 + c^2)*(2*t - 3*sin(t) + t*cos(t))*(a*d + b*e + c*f))/t^5 - ((t - sin(t))*(2*b*e + 2*c*f))/t^3, ((c*(sin(t)/t + (2*cos(t) - 2)/t^2) + a*b*(((3*sin(t))/t - 3)/t^2 - (2*cos(t) - 2)/(2*t^2)))*(a*d + b*e + c*f))/t^2 - ((sin(t)/t - 1)*(a*e + b*d))/t^2 - (f*(2*cos(t) - 2))/(2*t^2), (e*(2*cos(t) - 2))/(2*t^2) - ((sin(t)/t - 1)*(a*f + c*d))/t^2 - ((b*(sin(t)/t + (2*cos(t) - 2)/t^2) - a*c*(((3*sin(t))/t - 3)/t^2 - (2*cos(t) - 2)/(2*t^2)))*(a*d + b*e + c*f))/t^2]
    [ (f*(2*cos(t) - 2))/(2*t^2) - ((sin(t)/t - 1)*(a*e + b*d))/t^2 - ((c*(sin(t)/t + (2*cos(t) - 2)/t^2) - a*b*(((3*sin(t))/t - 3)/t^2 - (2*cos(t) - 2)/(2*t^2)))*(a*d + b*e + c*f))/t^2,                                                                                ((a^2 + c^2)*(2*t - 3*sin(t) + t*cos(t))*(a*d + b*e + c*f))/t^5 - ((t - sin(t))*(2*a*d + 2*c*f))/t^3, ((a*(sin(t)/t + (2*cos(t) - 2)/t^2) + b*c*(((3*sin(t))/t - 3)/t^2 - (2*cos(t) - 2)/(2*t^2)))*(a*d + b*e + c*f))/t^2 - ((sin(t)/t - 1)*(b*f + c*e))/t^2 - (d*(2*cos(t) - 2))/(2*t^2)]
    [ ((b*(sin(t)/t + (2*cos(t) - 2)/t^2) + a*c*(((3*sin(t))/t - 3)/t^2 - (2*cos(t) - 2)/(2*t^2)))*(a*d + b*e + c*f))/t^2 - ((sin(t)/t - 1)*(a*f + c*d))/t^2 - (e*(2*cos(t) - 2))/(2*t^2), (d*(2*cos(t) - 2))/(2*t^2) - ((sin(t)/t - 1)*(b*f + c*e))/t^2 - ((a*(sin(t)/t + (2*cos(t) - 2)/t^2) - b*c*(((3*sin(t))/t - 3)/t^2 - (2*cos(t) - 2)/(2*t^2)))*(a*d + b*e + c*f))/t^2,                                                                                ((a^2 + b^2)*(2*t - 3*sin(t) + t*cos(t))*(a*d + b*e + c*f))/t^5 - ((t - sin(t))*(2*a*d + 2*b*e))/t^3]

    */

    gsl_matrix *TSO = gsl_matrix_alloc(3,3);
    tanSO3(in, TSO);

    gsl_matrix *T = gsl_matrix_alloc(3,3);

    double a = gsl_vector_get(in,0);
    double b = gsl_vector_get(in,1);
    double c = gsl_vector_get(in,2);
    double d = gsl_vector_get(in,3);
    double e = gsl_vector_get(in,4);
    double f = gsl_vector_get(in,5);
    double t = sqrt(a*a+b*b+c*c);
    if (t == 0) {
        gsl_matrix_set(T, 0, 0, 0);
        gsl_matrix_set(T, 0, 1, f/2);
        gsl_matrix_set(T, 0, 2, -e/2);
        gsl_matrix_set(T, 1, 0, -f/2);
        gsl_matrix_set(T, 1, 1, 0);
        gsl_matrix_set(T, 1, 2, d/2);
        gsl_matrix_set(T, 2, 0, e/2);
        gsl_matrix_set(T, 2, 1, -d/2);
        gsl_matrix_set(T, 2, 2, 0);
    } else {
        gsl_matrix_set(T, 0, 0, ((pow(b,2) + pow(c,2))*(2*t - 3*sin(t) + t*cos(t))*(a*d + b*e + c*f))/pow(t,5) - ((t - sin(t))*(2*b*e + 2*c*f))/pow(t,3));
        gsl_matrix_set(T, 0, 1, ((c*(sin(t)/t + (2*cos(t) - 2)/pow(t,2)) + a*b*(((3*sin(t))/t - 3)/pow(t,2) - (2*cos(t) - 2)/(2*pow(t,2))))*(a*d + b*e + c*f))/pow(t,2) - ((sin(t)/t - 1)*(a*e + b*d))/pow(t,2) - (f*(2*cos(t) - 2))/(2*pow(t,2)));
        gsl_matrix_set(T, 0, 2, (e*(2*cos(t) - 2))/(2*pow(t,2)) - ((sin(t)/t - 1)*(a*f + c*d))/pow(t,2) - ((b*(sin(t)/t + (2*cos(t) - 2)/pow(t,2)) - a*c*(((3*sin(t))/t - 3)/pow(t,2) - (2*cos(t) - 2)/(2*pow(t,2))))*(a*d + b*e + c*f))/pow(t,2));
        gsl_matrix_set(T, 1, 0, (f*(2*cos(t) - 2))/(2*pow(t,2)) - ((sin(t)/t - 1)*(a*e + b*d))/pow(t,2) - ((c*(sin(t)/t + (2*cos(t) - 2)/pow(t,2)) - a*b*(((3*sin(t))/t - 3)/pow(t,2) - (2*cos(t) - 2)/(2*pow(t,2))))*(a*d + b*e + c*f))/pow(t,2));
        gsl_matrix_set(T, 1, 1, ((pow(a,2) + pow(c,2))*(2*t - 3*sin(t) + t*cos(t))*(a*d + b*e + c*f))/pow(t,5) - ((t - sin(t))*(2*a*d + 2*c*f))/pow(t,3));
        gsl_matrix_set(T, 1, 2, ((a*(sin(t)/t + (2*cos(t) - 2)/pow(t,2)) + b*c*(((3*sin(t))/t - 3)/pow(t,2) - (2*cos(t) - 2)/(2*pow(t,2))))*(a*d + b*e + c*f))/pow(t,2) - ((sin(t)/t - 1)*(b*f + c*e))/pow(t,2) - (d*(2*cos(t) - 2))/(2*pow(t,2)));
        gsl_matrix_set(T, 2, 0, ((b*(sin(t)/t + (2*cos(t) - 2)/pow(t,2)) + a*c*(((3*sin(t))/t - 3)/pow(t,2) - (2*cos(t) - 2)/(2*pow(t,2))))*(a*d + b*e + c*f))/pow(t,2) - ((sin(t)/t - 1)*(a*f + c*d))/pow(t,2) - (e*(2*cos(t) - 2))/(2*pow(t,2)));
        gsl_matrix_set(T, 2, 1, (d*(2*cos(t) - 2))/(2*pow(t,2)) - ((sin(t)/t - 1)*(b*f + c*e))/pow(t,2) - ((a*(sin(t)/t + (2*cos(t) - 2)/pow(t,2)) - b*c*(((3*sin(t))/t - 3)/pow(t,2) - (2*cos(t) - 2)/(2*pow(t,2))))*(a*d + b*e + c*f))/pow(t,2));
        gsl_matrix_set(T, 2, 2, ((pow(a,2) + pow(b,2))*(2*t - 3*sin(t) + t*cos(t))*(a*d + b*e + c*f))/pow(t,5) - ((t - sin(t))*(2*a*d + 2*b*e))/pow(t,3));
    }

    int i, j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            gsl_matrix_set(out,i,j,gsl_matrix_get(TSO,i,j));
            gsl_matrix_set(out,i+3,j+3,gsl_matrix_get(TSO,i,j));
            gsl_matrix_set(out,i+3,j,0);
            gsl_matrix_set(out,i,j+3,gsl_matrix_get(T,i,j));
        }
    }
    gsl_matrix_free(TSO);
    gsl_matrix_free(T);


}

void itanSE3(const gsl_vector *in, gsl_matrix *out) {
    /*
    [itanSO(w), T]
    [0, itanSO(w)]
    T =
    [ ((b^2 + c^2)*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + t^2 - 4))/(2*t^4*(cos(t) - 1)) - ((2*b*e + 2*c*f)*(2*cos(t) + t*sin(t) - 2))/(2*t^2*(cos(t) - 1)),       ((a*e + b*d)*(2*cos(t) + t*sin(t) - 2))/(2*t^2*(cos(t) - 1)) - f/2 - (a*b*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + t^2 - 4))/(2*t^4*(cos(t) - 1)),       e/2 + ((a*f + c*d)*(2*cos(t) + t*sin(t) - 2))/(2*t^2*(cos(t) - 1)) - (a*c*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + t^2 - 4))/(2*t^4*(cos(t) - 1))]
    [       f/2 + ((a*e + b*d)*(2*cos(t) + t*sin(t) - 2))/(2*t^2*(cos(t) - 1)) - (a*b*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + t^2 - 4))/(2*t^4*(cos(t) - 1)), ((a^2 + c^2)*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + t^2 - 4))/(2*t^4*(cos(t) - 1)) - ((2*a*d + 2*c*f)*(2*cos(t) + t*sin(t) - 2))/(2*t^2*(cos(t) - 1)),       ((b*f + c*e)*(2*cos(t) + t*sin(t) - 2))/(2*t^2*(cos(t) - 1)) - d/2 - (b*c*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + t^2 - 4))/(2*t^4*(cos(t) - 1))]
    [       ((a*f + c*d)*(2*cos(t) + t*sin(t) - 2))/(2*t^2*(cos(t) - 1)) - e/2 - (a*c*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + t^2 - 4))/(2*t^4*(cos(t) - 1)),       d/2 + ((b*f + c*e)*(2*cos(t) + t*sin(t) - 2))/(2*t^2*(cos(t) - 1)) - (b*c*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + t^2 - 4))/(2*t^4*(cos(t) - 1)), ((a^2 + b^2)*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + t^2 - 4))/(2*t^4*(cos(t) - 1)) - ((2*a*d + 2*b*e)*(2*cos(t) + t*sin(t) - 2))/(2*t^2*(cos(t) - 1))]
    */

    gsl_matrix *TSO = gsl_matrix_alloc(3,3);
    itanSO3(in, TSO);

    gsl_matrix *T = gsl_matrix_alloc(3,3);

    double a = gsl_vector_get(in,0);
    double b = gsl_vector_get(in,1);
    double c = gsl_vector_get(in,2);
    double d = gsl_vector_get(in,3);
    double e = gsl_vector_get(in,4);
    double f = gsl_vector_get(in,5);
    double t = sqrt(a*a+b*b+c*c);
    if (t == 0) {
        gsl_matrix_set(T, 0, 0, 0);
        gsl_matrix_set(T, 0, 1, -f/2);
        gsl_matrix_set(T, 0, 2, e/2);
        gsl_matrix_set(T, 1, 0, f/2);
        gsl_matrix_set(T, 1, 1, 0);
        gsl_matrix_set(T, 1, 2, -d/2);
        gsl_matrix_set(T, 2, 0, -e/2);
        gsl_matrix_set(T, 2, 1, d/2);
        gsl_matrix_set(T, 2, 2, 0);
    } else {
        gsl_matrix_set(T, 0, 0, ((pow(b,2) + pow(c,2))*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + pow(t,2) - 4))/(2*pow(t,4)*(cos(t) - 1)) - ((2*b*e + 2*c*f)*(2*cos(t) + t*sin(t) - 2))/(2*pow(t,2)*(cos(t) - 1)));
        gsl_matrix_set(T, 0, 1, ((a*e + b*d)*(2*cos(t) + t*sin(t) - 2))/(2*pow(t,2)*(cos(t) - 1)) - f/2 - (a*b*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + pow(t,2) - 4))/(2*pow(t,4)*(cos(t) - 1)));
        gsl_matrix_set(T, 0, 2, e/2 + ((a*f + c*d)*(2*cos(t) + t*sin(t) - 2))/(2*pow(t,2)*(cos(t) - 1)) - (a*c*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + pow(t,2) - 4))/(2*pow(t,4)*(cos(t) - 1)));
        gsl_matrix_set(T, 1, 0, f/2 + ((a*e + b*d)*(2*cos(t) + t*sin(t) - 2))/(2*pow(t,2)*(cos(t) - 1)) - (a*b*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + pow(t,2) - 4))/(2*pow(t,4)*(cos(t) - 1)));
        gsl_matrix_set(T, 1, 1, ((pow(a,2) + pow(c,2))*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + pow(t,2) - 4))/(2*pow(t,4)*(cos(t) - 1)) - ((2*a*d + 2*c*f)*(2*cos(t) + t*sin(t) - 2))/(2*pow(t,2)*(cos(t) - 1)));
        gsl_matrix_set(T, 1, 2, ((b*f + c*e)*(2*cos(t) + t*sin(t) - 2))/(2*pow(t,2)*(cos(t) - 1)) - d/2 - (b*c*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + pow(t,2) - 4))/(2*pow(t,4)*(cos(t) - 1)));
        gsl_matrix_set(T, 2, 0, ((a*f + c*d)*(2*cos(t) + t*sin(t) - 2))/(2*pow(t,2)*(cos(t) - 1)) - e/2 - (a*c*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + pow(t,2) - 4))/(2*pow(t,4)*(cos(t) - 1)));
        gsl_matrix_set(T, 2, 1, d/2 + ((b*f + c*e)*(2*cos(t) + t*sin(t) - 2))/(2*pow(t,2)*(cos(t) - 1)) - (b*c*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + pow(t,2) - 4))/(2*pow(t,4)*(cos(t) - 1)));
        gsl_matrix_set(T, 2, 2, ((pow(a,2) + pow(b,2))*(a*d + b*e + c*f)*(4*cos(t) + t*sin(t) + pow(t,2) - 4))/(2*pow(t,4)*(cos(t) - 1)) - ((2*a*d + 2*b*e)*(2*cos(t) + t*sin(t) - 2))/(2*pow(t,2)*(cos(t) - 1)));
    }

    int i, j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            gsl_matrix_set(out,i,j,gsl_matrix_get(TSO,i,j));
            gsl_matrix_set(out,i+3,j+3,gsl_matrix_get(TSO,i,j));
            gsl_matrix_set(out,i+3,j,0);
            gsl_matrix_set(out,i,j+3,gsl_matrix_get(T,i,j));
        }
    }
    gsl_matrix_free(TSO);
    gsl_matrix_free(T);
}


/*
SO(3) cayley transform
*/

void caySO3(const gsl_vector *in, gsl_matrix *out) {
    /*
    [ (a^2 - b^2 - c^2 + 4)/(a^2 + b^2 + c^2 + 4),       -(2*(2*c - a*b))/(a^2 + b^2 + c^2 + 4),        (2*(2*b + a*c))/(a^2 + b^2 + c^2 + 4)]
    [       (2*(2*c + a*b))/(a^2 + b^2 + c^2 + 4), -(a^2 - b^2 + c^2 - 4)/(a^2 + b^2 + c^2 + 4),       -(2*(2*a - b*c))/(a^2 + b^2 + c^2 + 4)]
    [      -(2*(2*b - a*c))/(a^2 + b^2 + c^2 + 4),        (2*(2*a + b*c))/(a^2 + b^2 + c^2 + 4), -(a^2 + b^2 - c^2 - 4)/(a^2 + b^2 + c^2 + 4)]
    */

    double a = gsl_vector_get(in,0);
    double b = gsl_vector_get(in,1);
    double c = gsl_vector_get(in,2);

    gsl_matrix_set(out, 0, 0, (pow(a,2) - pow(b,2) - pow(c,2) + 4)/(pow(a,2) + pow(b,2) + pow(c,2) + 4));
    gsl_matrix_set(out, 0, 1, -(2*(2*c - a*b))/(pow(a,2) + pow(b,2) + pow(c,2) + 4));
    gsl_matrix_set(out, 0, 2, (2*(2*b + a*c))/(pow(a,2) + pow(b,2) + pow(c,2) + 4));
    gsl_matrix_set(out, 1, 0, (2*(2*c + a*b))/(pow(a,2) + pow(b,2) + pow(c,2) + 4));
    gsl_matrix_set(out, 1, 1, -(pow(a,2) - pow(b,2) + pow(c,2) - 4)/(pow(a,2) + pow(b,2) + pow(c,2) + 4));
    gsl_matrix_set(out, 1, 2, -(2*(2*a - b*c))/(pow(a,2) + pow(b,2) + pow(c,2) + 4));
    gsl_matrix_set(out, 2, 0, -(2*(2*b - a*c))/(pow(a,2) + pow(b,2) + pow(c,2) + 4));
    gsl_matrix_set(out, 2, 1, (2*(2*a + b*c))/(pow(a,2) + pow(b,2) + pow(c,2) + 4));
    gsl_matrix_set(out, 2, 2, -(pow(a,2) + pow(b,2) - pow(c,2) - 4)/(pow(a,2) + pow(b,2) + pow(c,2) + 4));
}

void icaySO3(const gsl_matrix *in, gsl_vector *out) {
    /*
    (4*(r8 + r1*r8 - r2*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1)
    (4*(r3 - r2*r6 + r3*r5))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1)
    (4*(r4 + r4*r9 - r6*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1)
    */

    double r1 = gsl_matrix_get(in,0,0);
    double r2 = gsl_matrix_get(in,0,1);
    double r3 = gsl_matrix_get(in,0,2);
    double r4 = gsl_matrix_get(in,1,0);
    double r5 = gsl_matrix_get(in,1,1);
    double r6 = gsl_matrix_get(in,1,2);
    double r7 = gsl_matrix_get(in,2,0);
    double r8 = gsl_matrix_get(in,2,1);
    double r9 = gsl_matrix_get(in,2,2);

    gsl_vector_set(out,0,(4*(r8 + r1*r8 - r2*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1));
    gsl_vector_set(out,1,(4*(r3 - r2*r6 + r3*r5))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1));
    gsl_vector_set(out,2,(4*(r4 + r4*r9 - r6*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1));
}

void dcaySO3(const gsl_vector *in, gsl_vector *out) {

}

/*
SE(3) cayley
*/

void caySE3(const gsl_vector *in, gsl_matrix *out) {
    /*
    upper right
    4/(4+t^2) *
    f*(b/2 + (a*c)/4) - e*(c/2 - (a*b)/4) + d*(a^2/4 + 1)
    d*(c/2 + (a*b)/4) - f*(a/2 - (b*c)/4) + e*(b^2/4 + 1)
    e*(a/2 + (b*c)/4) - d*(b/2 - (a*c)/4) + f*(c^2/4 + 1)
    */
    caySO3(in,out);
    double a = gsl_vector_get(in, 0);
    double b = gsl_vector_get(in, 1);
    double c = gsl_vector_get(in, 2);
    double d = gsl_vector_get(in, 3);
    double e = gsl_vector_get(in, 4);
    double f = gsl_vector_get(in, 5);
    double t = a*a+b*b+c*c;

    gsl_matrix_set(out, 0, 3, 4/(4+t)*(f*(b/2 + (a*c)/4) - e*(c/2 - (a*b)/4) + d*(a*a/4 + 1)));
    gsl_matrix_set(out, 1, 3, 4/(4+t)*(d*(c/2 + (a*b)/4) - f*(a/2 - (b*c)/4) + e*(b*b/4 + 1)));
    gsl_matrix_set(out, 2, 3, 4/(4+t)*(e*(a/2 + (b*c)/4) - d*(b/2 - (a*c)/4) + f*(c*c/4 + 1)));

    gsl_matrix_set(out, 3, 0, 0);
    gsl_matrix_set(out, 3, 1, 0);
    gsl_matrix_set(out, 3, 2, 0);
    gsl_matrix_set(out, 3, 3, 1);
}

void icaySE3(const gsl_matrix *in, gsl_vector *out) {
    /*
    There could be terms that are being dropped, but I think that may just be an issue with

    (4*(r8 + r1*r8 - r2*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1)
    (4*(r3 - r2*r6 + r3*r5))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1)
    (4*(r4 + r4*r9 - r6*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1)
    (2*(p1 - p2*r2 + p1*r5 - p3*r3 + p1*r9 + p3*r2*r6 - p3*r3*r5 - p2*r2*r9 + p2*r3*r8 + p1*r5*r9 - p1*r6*r8))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1)
    (2*(p2 + p2*r1 - p1*r4 - p3*r6 + p2*r9 - p3*r1*r6 + p3*r3*r4 + p2*r1*r9 - p2*r3*r7 - p1*r4*r9 + p1*r6*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1)
    (2*(p3 + p3*r1 - p1*r7 + p3*r5 - p2*r8 + p3*r1*r5 - p3*r2*r4 - p2*r1*r8 + p2*r2*r7 + p1*r4*r8 - p1*r5*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1)
    */
    double r1 = gsl_matrix_get(in, 0, 0);
    double r2 = gsl_matrix_get(in, 0, 1);
    double r3 = gsl_matrix_get(in, 0, 2);
    double r4 = gsl_matrix_get(in, 1, 0);
    double r5 = gsl_matrix_get(in, 1, 1);
    double r6 = gsl_matrix_get(in, 1, 2);
    double r7 = gsl_matrix_get(in, 2, 0);
    double r8 = gsl_matrix_get(in, 2, 1);
    double r9 = gsl_matrix_get(in, 2, 2);

    double p1 = gsl_matrix_get(in, 0, 3);
    double p2 = gsl_matrix_get(in, 1, 3);
    double p3 = gsl_matrix_get(in, 2, 3);

    gsl_vector_set(out, 0, (4*(r8 + r1*r8 - r2*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1));
    gsl_vector_set(out, 1, (4*(r3 - r2*r6 + r3*r5))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1));
    gsl_vector_set(out, 2, (4*(r4 + r4*r9 - r6*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1));
    gsl_vector_set(out, 3, (2*(p1 - p2*r2 + p1*r5 - p3*r3 + p1*r9 + p3*r2*r6 - p3*r3*r5 - p2*r2*r9 + p2*r3*r8 + p1*r5*r9 - p1*r6*r8))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1));
    gsl_vector_set(out, 4, (2*(p2 + p2*r1 - p1*r4 - p3*r6 + p2*r9 - p3*r1*r6 + p3*r3*r4 + p2*r1*r9 - p2*r3*r7 - p1*r4*r9 + p1*r6*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1));
    gsl_vector_set(out, 5, (2*(p3 + p3*r1 - p1*r7 + p3*r5 - p2*r8 + p3*r1*r5 - p3*r2*r4 - p2*r1*r8 + p2*r2*r7 + p1*r4*r8 - p1*r5*r7))/(r1 + r5 + r9 + r1*r5 - r2*r4 + r1*r9 - r3*r7 + r5*r9 - r6*r8 + r1*r5*r9 - r1*r6*r8 - r2*r4*r9 + r2*r6*r7 + r3*r4*r8 - r3*r5*r7 + 1));

}

/*
Other useful stuff
*/

void adjointSE3(const gsl_vector *in, gsl_matrix *out) {
    // seems like some of the copying could be optimized away by the compiler, but do it manually anyways (definately a case of overly early optimization)
    double wx = gsl_vector_get(in, 0);
    double wy = gsl_vector_get(in, 1);
    double wz = gsl_vector_get(in, 2);
    double vx = gsl_vector_get(in, 3);
    double vy = gsl_vector_get(in, 4);
    double vz = gsl_vector_get(in, 5);

    //upper left
    gsl_matrix_set(out, 0, 0, 0);
    gsl_matrix_set(out, 1, 1, 0);
    gsl_matrix_set(out, 2, 2, 0);
    gsl_matrix_set(out, 0, 1, -wz);
    gsl_matrix_set(out, 0, 2, wy);
    gsl_matrix_set(out, 1, 0, wz);
    gsl_matrix_set(out, 1, 2, -wx);
    gsl_matrix_set(out, 2, 0, -wy);
    gsl_matrix_set(out, 2, 1, wx);

    //lower right
    gsl_matrix_set(out, 3, 3, 0);
    gsl_matrix_set(out, 4, 4, 0);
    gsl_matrix_set(out, 5, 5, 0);
    gsl_matrix_set(out, 3, 4, -wz);
    gsl_matrix_set(out, 3, 5, wy);
    gsl_matrix_set(out, 4, 3, wz);
    gsl_matrix_set(out, 4, 5, -wx);
    gsl_matrix_set(out, 5, 3, -wy);
    gsl_matrix_set(out, 5, 4, wx);

    //lower left
    gsl_matrix_set(out, 3, 0, 0);
    gsl_matrix_set(out, 4, 1, 0);
    gsl_matrix_set(out, 5, 2, 0);
    gsl_matrix_set(out, 3, 1, -vz);
    gsl_matrix_set(out, 3, 2, vy);
    gsl_matrix_set(out, 4, 0, vz);
    gsl_matrix_set(out, 4, 2, -vx);
    gsl_matrix_set(out, 5, 0, -vy);
    gsl_matrix_set(out, 5, 1, vx);

    //upper right
    int i, j;
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            gsl_matrix_set(out, i, j+3, 0);
        }
    }
}

void AdjointSE3(const gsl_matrix *in, gsl_matrix *out) {

    gsl_matrix *p_skew = gsl_matrix_alloc(3,3);
    gsl_matrix *corner = gsl_matrix_alloc(3,3);

    gsl_vector_const_view p = gsl_matrix_const_subcolumn(in, 3, 0, 3);
    gsl_matrix_const_view R = gsl_matrix_const_submatrix(in, 0, 0, 3, 3);
    skew(&p.vector, p_skew);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, p_skew, &R.matrix, 0.0, corner);

    int i, j;
    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            gsl_matrix_set(out, i, j, gsl_matrix_get(&R.matrix, i, j));
            gsl_matrix_set(out, i+3, j+3, gsl_matrix_get(&R.matrix, i, j));
            gsl_matrix_set(out, i+3, j, gsl_matrix_get(corner, i, j));
            gsl_matrix_set(out, i, j+3, 0);
        }
    }



}
