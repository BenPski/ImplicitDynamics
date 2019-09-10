#ifndef LIE_THEORY_H
#define LIE_THEORY_H

#include <gsl/gsl_matrix.h>

void flattenConfig(const gsl_matrix *in, gsl_vector *out);
void unflattenConfig(const gsl_vector *in, gsl_matrix *out);
void skew(const gsl_vector *in, gsl_matrix *out);
void unskew(const gsl_matrix *in, gsl_vector *out);
void se(const gsl_vector *in, gsl_matrix *out);
void unse(const gsl_matrix *in, gsl_vector *out);
void expSO3(const gsl_vector *in, gsl_matrix *out);
void logSO3(const gsl_matrix *in, gsl_vector *out);
void expSE3(const gsl_vector *in, gsl_matrix *out);
void logSE3(const gsl_matrix *in, gsl_vector *out);
void tanSE3(const gsl_vector *in, gsl_matrix *out);
void itanSE3(const gsl_vector *in, gsl_matrix *out);
void adjointSE3(const gsl_vector *in, gsl_matrix *out);
void AdjointSE3(const gsl_matrix *in, gsl_matrix *out);

void caySE3(const gsl_vector *in, gsl_matrix *out);
void icaySE3(const gsl_matrix *in, gsl_vector *out);

#endif
