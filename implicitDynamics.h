#ifndef IMPLICITDYNAMICS_H
#define IMPLICITDYNAMICS_H

#include <gsl/gsl_matrix.h>

struct dynParams {
    double dt;
    const gsl_matrix *g;
    gsl_matrix *g_next;
    const gsl_matrix *y;
    gsl_matrix *y_next;
    const gsl_vector *q;
};

void rodDynamics(double s, double dt, const gsl_matrix *g, const gsl_vector *y, const gsl_vector *y_prev, const gsl_vector *q, gsl_vector *y_der);
void integrate(double dt, const gsl_matrix *g, const gsl_matrix *y, const gsl_vector *q, const gsl_vector *xi0, gsl_matrix *g_next, gsl_matrix *y_next);
void tipCondition(const gsl_matrix *g, const gsl_matrix *y, const gsl_vector *q, gsl_vector *condition);
int dynamics(const gsl_vector *xi0, void *params, gsl_vector *condition);
void solver(double dt, const gsl_matrix *g, const gsl_matrix *y, const gsl_vector *q, gsl_matrix *g_next, gsl_matrix *y_next);

#endif
