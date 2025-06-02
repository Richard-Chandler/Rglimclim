#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(closefiles)(void *, void *);
extern void F77_NAME(dailystats)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(glmfit)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(readdata2)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rglcfileconnect)(void *, void *, void *, void *);
extern void F77_NAME(scratchsync)(void *, void *, void *, void *);
extern void F77_NAME(scratchupdate)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(scropn)(void *, void *, void *, void *);
extern void F77_NAME(simulate)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(writedata)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(writeheader)(void *);
extern void F77_NAME(writelabels)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(zbqlini)(void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"closefiles",      (DL_FUNC) &F77_NAME(closefiles),       2},
    {"dailystats",      (DL_FUNC) &F77_NAME(dailystats),      11},
    {"glmfit",          (DL_FUNC) &F77_NAME(glmfit),          60},
    {"readdata2",       (DL_FUNC) &F77_NAME(readdata2),        9},
    {"rglcfileconnect", (DL_FUNC) &F77_NAME(rglcfileconnect),  4},
    {"scratchsync",     (DL_FUNC) &F77_NAME(scratchsync),      4},
    {"scratchupdate",   (DL_FUNC) &F77_NAME(scratchupdate),    6},
    {"scropn",          (DL_FUNC) &F77_NAME(scropn),           4},
    {"simulate",        (DL_FUNC) &F77_NAME(simulate),        45},
    {"writedata",       (DL_FUNC) &F77_NAME(writedata),       12},
    {"writeheader",     (DL_FUNC) &F77_NAME(writeheader),      1},
    {"writelabels",     (DL_FUNC) &F77_NAME(writelabels),     18},
    {"zbqlini",         (DL_FUNC) &F77_NAME(zbqlini),          1},
    {NULL, NULL, 0}
};

void R_init_Rglimclim(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
