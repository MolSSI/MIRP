#include "mirp/eri.h"
#include "mirp/math.h"
#include "mirp/mpfr_help.h"
#include "mirp/arb_help.h"
#include <Python.h>
#include <stdlib.h>

#define UNUSED(x) (void)(x)

PyObject * export_mirp_single_eri_double(PyObject * self, PyObject *args)
{
    UNUSED(self);

    double alpha[4];
    int lmn[4][3];
    double ABCD[4][3];

    if(!PyArg_ParseTuple(args, "(iii)d(ddd)(iii)d(ddd)(iii)d(ddd)(iii)d(ddd)",
                         &(lmn[0][0]), &(lmn[0][1]), &(lmn[0][2]), &(alpha[0]), &(ABCD[0][0]), &(ABCD[0][1]), &(ABCD[0][2]),
                         &(lmn[1][0]), &(lmn[1][1]), &(lmn[1][2]), &(alpha[1]), &(ABCD[1][0]), &(ABCD[1][1]), &(ABCD[1][2]),
                         &(lmn[2][0]), &(lmn[2][1]), &(lmn[2][2]), &(alpha[2]), &(ABCD[2][0]), &(ABCD[2][1]), &(ABCD[2][2]),
                         &(lmn[3][0]), &(lmn[3][1]), &(lmn[3][2]), &(alpha[3]), &(ABCD[3][0]), &(ABCD[3][1]), &(ABCD[3][2])))
        return NULL;



    double result;
    mirp_single_eri_double(&result,
                           lmn[0][0], lmn[0][1], lmn[0][2], alpha[0], ABCD[0],
                           lmn[1][0], lmn[1][1], lmn[1][2], alpha[1], ABCD[1],
                           lmn[2][0], lmn[2][1], lmn[2][2], alpha[2], ABCD[2],
                           lmn[3][0], lmn[3][1], lmn[3][2], alpha[3], ABCD[3]);

    PyObject * ret = PyFloat_FromDouble(result);
    return ret;
}


PyObject * export_mirp_prim_eri_double(PyObject * self, PyObject *args)
{
    UNUSED(self);

    double alpha[4];
    int am[4];
    double ABCD[4][3];

    if(!PyArg_ParseTuple(args, "id(ddd)id(ddd)id(ddd)id(ddd)",
                         &(am[0]), &(alpha[0]), &(ABCD[0][0]), &(ABCD[0][1]), &(ABCD[0][2]),
                         &(am[1]), &(alpha[1]), &(ABCD[1][0]), &(ABCD[1][1]), &(ABCD[1][2]),
                         &(am[2]), &(alpha[2]), &(ABCD[2][0]), &(ABCD[2][1]), &(ABCD[2][2]),
                         &(am[3]), &(alpha[3]), &(ABCD[3][0]), &(ABCD[3][1]), &(ABCD[3][2])))
        return NULL;



    size_t size = MIRP_NCART(am[0]) * MIRP_NCART(am[1]) * MIRP_NCART(am[2]) * MIRP_NCART(am[3]);
    double * result = (double *)malloc(size * sizeof(double));

    mirp_prim_eri_double(result,
                         am[0], alpha[0], ABCD[0],
                         am[1], alpha[1], ABCD[1],
                         am[2], alpha[2], ABCD[2],
                         am[3], alpha[3], ABCD[3]);

    PyObject * ret = PyList_New(size);
    for(size_t i = 0; i < size; i++)
        PyList_SetItem(ret, i, PyFloat_FromDouble(result[i]));

    free(result);
    return ret;
}


PyObject * export_mirp_eri_double(PyObject * self, PyObject *args)
{
    UNUSED(self);

    PyObject * alpha_p[4], * coeff_p[4];
    int am[4];
    double ABCD[4][3];

    if(!PyArg_ParseTuple(args, "i(ddd)O!O!i(ddd)O!O!i(ddd)O!O!i(ddd)O!O!",
                         &(am[0]), &(ABCD[0][0]), &(ABCD[0][1]), &(ABCD[0][2]), &PyList_Type, &(alpha_p[0]), &PyList_Type, &(coeff_p[0]),
                         &(am[1]), &(ABCD[1][0]), &(ABCD[1][1]), &(ABCD[1][2]), &PyList_Type, &(alpha_p[1]), &PyList_Type, &(coeff_p[1]),
                         &(am[2]), &(ABCD[2][0]), &(ABCD[2][1]), &(ABCD[2][2]), &PyList_Type, &(alpha_p[2]), &PyList_Type, &(coeff_p[2]),
                         &(am[3]), &(ABCD[3][0]), &(ABCD[3][1]), &(ABCD[3][2]), &PyList_Type, &(alpha_p[3]), &PyList_Type, &(coeff_p[3])))
        return NULL;


    int nprim[4];
    int ngeneral[4];
    nprim[0] = PyList_Size(alpha_p[0]);
    nprim[1] = PyList_Size(alpha_p[1]);
    nprim[2] = PyList_Size(alpha_p[2]);
    nprim[3] = PyList_Size(alpha_p[3]);

    ngeneral[0] = PyList_Size(coeff_p[0]) / nprim[0];
    ngeneral[1] = PyList_Size(coeff_p[1]) / nprim[1];
    ngeneral[2] = PyList_Size(coeff_p[2]) / nprim[2];
    ngeneral[3] = PyList_Size(coeff_p[3]) / nprim[3];

    size_t size = MIRP_NCART(am[0]) * MIRP_NCART(am[1]) * MIRP_NCART(am[2]) * MIRP_NCART(am[3]);
    size *= (ngeneral[0] * ngeneral[1] * ngeneral[2] * ngeneral[3]);
    double * result = (double *)malloc(size * sizeof(double));

    double * alpha[4];
    double * coeff[4];
    for(int i = 0; i < 4; i++)
    {
        alpha[i] = malloc(nprim[i] * sizeof(double));
        coeff[i] = malloc(nprim[i] * ngeneral[i] * sizeof(double));

        for(int j = 0; j < nprim[i]; j++)
            alpha[i][j] = PyFloat_AsDouble(PyList_GetItem(alpha_p[i], j));
        for(int j = 0; j < nprim[i]*ngeneral[i]; j++)
            coeff[i][j] = PyFloat_AsDouble(PyList_GetItem(coeff_p[i], j));
    }

    mirp_eri_double(result,
                    am[0], ABCD[0], nprim[0], ngeneral[0], alpha[0], coeff[0],
                    am[1], ABCD[1], nprim[1], ngeneral[1], alpha[1], coeff[1],
                    am[2], ABCD[2], nprim[2], ngeneral[2], alpha[2], coeff[2],
                    am[3], ABCD[3], nprim[3], ngeneral[3], alpha[3], coeff[3]);

    PyObject * ret = PyList_New(size);
    for(size_t i = 0; i < size; i++)
        PyList_SetItem(ret, i, PyFloat_FromDouble(result[i]));

    free(result);
    for(int i = 0; i < 4; i++)
    {
        free(alpha[i]);
        free(coeff[i]);
    }
    return ret;
}
