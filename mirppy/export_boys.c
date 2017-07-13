#include "mirp/kernels/boys.h"
#include "mirp/arb_help.h"
#include <Python.h>

#define UNUSED(x) (void)(x)

PyObject * export_mirp_boys_double(PyObject * self, PyObject *args)
{
    UNUSED(self);

    int m;
    double t;

    if(!PyArg_ParseTuple(args, "id", &m, &t))
        return NULL;

    double F[m+1];

    mirp_boys_double(F, m, t);

    PyObject * ret = PyList_New(m+1);

    for(int i = 0; i <= m; i++)
        PyList_SetItem(ret, i, PyFloat_FromDouble(F[i]));

    return ret;
}


PyObject * export_mirp_boys_interval(PyObject *self, PyObject *args)
{
    UNUSED(self);

    int m;
    long working_prec;
    const char * t = NULL;

    if(!PyArg_ParseTuple(args, "isl", &m, &t, &working_prec))
        return NULL;

    arb_t t_mp, F_mp[m+1];
    arb_init(t_mp);
    mirp_init_arb_arr(F_mp, m+1);
    arb_set_str(t_mp, t, working_prec);

    mirp_boys_interval(F_mp, m, t_mp, working_prec);

    PyObject * ret = PyList_New(m+1);

    for(int i = 0; i <= m; i++)
    {
        /* We need an (over)estimate on the number of digits we are going to print
         *
         * ndigits = log10(2) * precision (in bits) + a safety factor
         * If it is an integer, we use the working precision. The only
         * case where the result of the boys function is an integer is m = 0, t = 0,
         * and perhaps at very long range you might get 0.
         */
        long ndigits;
        slong f_prec = working_prec;
        if(!arb_is_int(F_mp[i]))
            f_prec = arb_rel_accuracy_bits(F_mp[i]);

        ndigits = 0.3010299956639812 * f_prec+4;

        char * strtmp = arb_get_str(F_mp[i], ndigits, ARB_STR_NO_RADIUS);

        PyObject * pystr = PyUnicode_FromString(strtmp);
        PyList_SetItem(ret, i, pystr);
        free(strtmp);
    }

    arb_clear(t_mp);
    mirp_clear_arb_arr(F_mp, m+1);

    return ret;
}

