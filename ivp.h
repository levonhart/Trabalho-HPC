#include <stdlib.h>

typedef double (ivp_function) (double, double*, size_t, double*);
typedef ivp_function *ivp_function_ptr;

void runge_kutta(ivp_function_ptr v[], unsigned count[], unsigned size,
		double constants[], unsigned nconst[]);
