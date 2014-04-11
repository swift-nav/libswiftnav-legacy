#include "f2c.h"
#include "signal1.h"
#ifdef __cplusplus
extern "C" {
#endif

size_t
#ifdef KR_headers
signal_(sigp, proc) integer *sigp; sig_pf proc;
#else
signal_(integer *sigp, sig_pf proc)
#endif
{
	int sig;
	sig = (int)*sigp;

	return (size_t)signal(sig, proc);
	}
#ifdef __cplusplus
}
#endif
