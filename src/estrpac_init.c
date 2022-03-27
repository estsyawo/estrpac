#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .C calls */
extern void Kern_DL(void *, void *, void *, void *);
extern void Kern_Esc(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"Kern_DL",  (DL_FUNC) &Kern_DL,  4},
    {"Kern_Esc", (DL_FUNC) &Kern_Esc, 4},
    {NULL, NULL, 0}
};

void R_init_estrpac(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
