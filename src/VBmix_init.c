#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP constrClassif(SEXP, SEXP, SEXP, SEXP);
extern SEXP control(SEXP, SEXP);
extern SEXP couple(SEXP, SEXP);
extern SEXP dDirichlet(SEXP, SEXP, SEXP);
extern SEXP EM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP extractSimpleModel(SEXP, SEXP);
extern SEXP gdist(SEXP, SEXP, SEXP, SEXP);
extern SEXP getLabels(SEXP, SEXP);
extern SEXP getResp(SEXP, SEXP);
extern SEXP gmmdensity(SEXP, SEXP);
extern SEXP gmmgen(SEXP, SEXP);
extern SEXP gmmkmsock(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP jsmc(SEXP, SEXP, SEXP);
extern SEXP jsut(SEXP, SEXP);
extern SEXP klmc(SEXP, SEXP, SEXP);
extern SEXP klut(SEXP, SEXP);
extern SEXP mergeClassif(SEXP, SEXP, SEXP, SEXP);
extern SEXP mixKnn(SEXP, SEXP, SEXP, SEXP);
extern SEXP mkmeans(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mmppca(SEXP, SEXP, SEXP, SEXP);
extern SEXP mppca(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP multinomial(SEXP, SEXP);
extern SEXP mvndensity(SEXP, SEXP, SEXP, SEXP);
extern SEXP mvngen(SEXP, SEXP, SEXP);
extern SEXP mvnradiusdensity(SEXP, SEXP);
extern SEXP rDirichlet(SEXP, SEXP);
extern SEXP sampleClassif(SEXP, SEXP, SEXP, SEXP);
extern SEXP sort_index(SEXP, SEXP);
extern SEXP varbayes(SEXP, SEXP, SEXP, SEXP);
extern SEXP vbcomp(SEXP, SEXP, SEXP, SEXP);
extern SEXP vbconstr(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"constrClassif",      (DL_FUNC) &constrClassif,      4},
    {"control",            (DL_FUNC) &control,            2},
    {"couple",             (DL_FUNC) &couple,             2},
    {"dDirichlet",         (DL_FUNC) &dDirichlet,         3},
    {"EM",                 (DL_FUNC) &EM,                 9},
    {"extractSimpleModel", (DL_FUNC) &extractSimpleModel, 2},
    {"gdist",              (DL_FUNC) &gdist,              4},
    {"getLabels",          (DL_FUNC) &getLabels,          2},
    {"getResp",            (DL_FUNC) &getResp,            2},
    {"gmmdensity",         (DL_FUNC) &gmmdensity,         2},
    {"gmmgen",             (DL_FUNC) &gmmgen,             2},
    {"gmmkmsock",          (DL_FUNC) &gmmkmsock,          5},
    {"jsmc",               (DL_FUNC) &jsmc,               3},
    {"jsut",               (DL_FUNC) &jsut,               2},
    {"klmc",               (DL_FUNC) &klmc,               3},
    {"klut",               (DL_FUNC) &klut,               2},
    {"mergeClassif",       (DL_FUNC) &mergeClassif,       4},
    {"mixKnn",             (DL_FUNC) &mixKnn,             4},
    {"mkmeans",            (DL_FUNC) &mkmeans,            5},
    {"mmppca",             (DL_FUNC) &mmppca,             4},
    {"mppca",              (DL_FUNC) &mppca,              5},
    {"multinomial",        (DL_FUNC) &multinomial,        2},
    {"mvndensity",         (DL_FUNC) &mvndensity,         4},
    {"mvngen",             (DL_FUNC) &mvngen,             3},
    {"mvnradiusdensity",   (DL_FUNC) &mvnradiusdensity,   2},
    {"rDirichlet",         (DL_FUNC) &rDirichlet,         2},
    {"sampleClassif",      (DL_FUNC) &sampleClassif,      4},
    {"sort_index",         (DL_FUNC) &sort_index,         2},
    {"varbayes",           (DL_FUNC) &varbayes,           4},
    {"vbcomp",             (DL_FUNC) &vbcomp,             4},
    {"vbconstr",           (DL_FUNC) &vbconstr,           5},
    {NULL, NULL, 0}
};

void R_init_VBmix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

