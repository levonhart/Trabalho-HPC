#ifndef IVP_H
#define IVP_H

#include <stdlib.h>

#define MIN_DV 0.0000001
typedef double (ivp_function) (double, double*, size_t, double*);
typedef ivp_function *ivp_function_ptr;

void runge_kutta(const double t, const ivp_function_ptr y[], unsigned count[],
		unsigned size, double c[], int nconst[], int nump);

#endif /* end of include guard: IVP_H */
