// This file implements a fast C version of python functions written by Jesse Bloom, 2007-2009
// It is designed for fast implementation of phylogenetic likelihood calculations.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <numpy/arrayobject.h>
// USE_ACCELERATE_CBLAS is preprocessor command specifying that we perform matrix multiplications
// using the blas functions from the Mac OS X Accelerate cblas library.  Otherwise the matrix 
// multiplications are performed in pure C code.  If we do not import the Accelerate library, then
// we need to separately import the complex.h header for complex numbers
#ifdef USE_ACCELERATE_CBLAS
#  include <Accelerate/Accelerate.h> 
#else
#  include <complex.h>
#endif
//
// Documentation string for MList function
static char MList_docs[] = "This method computes a list of the mutation probabilities.\n\nThis method is a fast C implementation of the '_MList' method of ddg_inference.DDGSet objects.\nIts calculations parallel that of ddg_inference.DDGSet.MList.\nIt takes the following calling variables (in order):\n'uts' is a Python list of floats with meaning described in ddg_inference.DDGSet._MList doc string.\n'only_need' is a Python list of integers with meaning described in ddg_inference.DDGSet._MList doc string.\n\n'n_aa' is a python integer giving the number of amino acid types (typically 20).\n'length' is a python integer giving the length of the protein.\n'grs' is the _grs variables of the DDGSet.  It is a list of tuples giving the diagonalized gr matrices.\n\tEach tuple must contain matrices that are either all of the type 'float_' (C type double) or all of\n\tthe type 'complex_' (C type double complex).\n'mlist' is a numpy ndarray of type 'float_' (C type double) of the appropriate length, and in C-style\n\tcontiguous format.  This array is filled with the appropriate values, as described in DDGSet._MList.\n'residue_to_compute' is the number of the residue for which we compute the 'MList'.\nAll numpy ndarray arguments must be C-style contiguous.";

// This MList function creates a list of the M values giving the probabilities of mutations
static PyObject *MList(PyObject *self, PyObject *args) {
    // Calling variables are (in order): uts, only_need, n_aa, length, grs, mlist, residue_to_compute
    PyObject *uts, *only_need, *grs, *r_grs_tuple, *gr_diag, *p, *p_inv, *mlist;
    long n_aa, length, n_aa2, n_uts, residue_to_compute, i_ut, y, index, irowcolumn, only_need_index, only_need_i, only_need_i_n_aa;
    double *arr_mlist, *cp_inv, *cp, *cgr_diag, *cp_inv_exp;
    complex double *complex_cp_inv, *complex_cp, *complex_cgr_diag, *complex_cp_inv_exp, *complex_naa2_list, *complex_naa_list;
    double ut, exp_ut;
    complex double complex_exp_ut;
    int array_type;
#ifdef USE_ACCELERATE_CBLAS
    complex double complex_one = 1, complex_zero = 0;
    long index2;
#else
    complex double complex_mxy;
    double mxy;
    long x, yindex;
#endif
    // Parse the arguments.  
    if (! PyArg_ParseTuple( args, "O!O!llO!O!l", &PyList_Type, &uts, &PyList_Type, &only_need, &n_aa, &length, &PyList_Type, &grs, &PyArray_Type, &mlist, &residue_to_compute)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to MList.");
        return NULL;
    }
    // Error checking on arguments
    if (length < 1) { // length of the protein
        PyErr_SetString(PyExc_ValueError, "length is less than one.");
        return NULL;
    }
    if (n_aa < 1) { // number of amino acids.  Normally will be 20.
        PyErr_SetString(PyExc_ValueError, "n_aa is less than one.");
        return NULL;
    }
    if (PyList_GET_SIZE(grs) != length) { // make sure grs is of the same size as length
        // PyErr_SetString(PyExc_ValueError, "grs is not of the same size as length.");
        char errstring[200];
        sprintf(errstring, "grs is not of the same size as length: %ld, %ld", PyList_GET_SIZE(grs), length);
        PyErr_SetString(PyExc_ValueError, errstring);
        return NULL;
    }
    n_uts = PyList_GET_SIZE(uts); // number of entries in uts
    if (n_uts < 1) { // make sure there are entries in uts 
        PyErr_SetString(PyExc_ValueError, "uts has no entries.");
        return NULL;
    }
    if (PyList_GET_SIZE(only_need) != n_uts * length) { // make sure only_need is of correct size
        PyErr_SetString(PyExc_ValueError, "only_need is of wrong size.");
        return NULL;
    }
    if (! ((residue_to_compute >= 0) && (residue_to_compute < length))) {
        char errstring[200];
        sprintf(errstring, "Invalid value for residue_to_compute: %ld, %ld", residue_to_compute, length);
        PyErr_SetString(PyExc_ValueError, errstring);
        return NULL;
    }
    n_aa2 = n_aa * n_aa; // square of the number of amino acids
    // The results will be returned in a numpy ndarray 'float_' (C type double) array called mlist.
    // This array will be of size length * n_uts * n_aa2.
    arr_mlist = (double *) PyArray_DATA(mlist); // this is the data array of mlist
    long const sizeof_cp_inv = n_aa2 * sizeof(double);
    long const complex_sizeof_cp_inv = n_aa2 * sizeof(complex double);
    // Now begin filling arr_mlist with the appropriate values
    index = 0;
    only_need_index = residue_to_compute * n_uts;
    r_grs_tuple = PyList_GET_ITEM(grs, residue_to_compute);
    // gr_diag, p, and p_inv are the eigenvalues, left, and right diagonalizing matrices of gr
    gr_diag = PyTuple_GET_ITEM(r_grs_tuple, 0); 
    p = PyTuple_GET_ITEM(r_grs_tuple, 1);
    p_inv = PyTuple_GET_ITEM(r_grs_tuple, 2);
    // determine if these arrays are complex double or real doubles
    array_type = PyArray_TYPE(gr_diag);
    if (array_type == NPY_DOUBLE) { // array is of doubles, real not complex
        cp_inv_exp = (double *) malloc(sizeof_cp_inv); // allocate memory
        // Note that these next assignments assume that the arrays are C-style contiguous
        cp = PyArray_DATA(p);
        cp_inv = PyArray_DATA(p_inv);
        cgr_diag = PyArray_DATA(gr_diag);
        for (i_ut = 0; i_ut < n_uts; i_ut++) { // loop over ut values
            only_need_i = PyInt_AS_LONG(PyList_GET_ITEM(only_need, only_need_index)); // value of only_need
            only_need_index++;
            if (only_need_i == -1) { // we don't need to do anything for these entries in mlist
                index += n_aa2;
            } else { // we need to compute at least some entries in mlist
                ut = PyFloat_AS_DOUBLE(PyList_GET_ITEM(uts, i_ut));  // ut value
                for (irowcolumn = 0; irowcolumn < n_aa; irowcolumn++) {
                    exp_ut = exp(ut * cgr_diag[irowcolumn]);
                    for (y = irowcolumn * n_aa; y < (irowcolumn + 1) * n_aa; y++) {
                        cp_inv_exp[y] = cp_inv[y] * exp_ut;
                    }
                }
                if (only_need_i == -2) { // we need to compute all of these entries in mlist
#ifdef USE_ACCELERATE_CBLAS
                    // multiply the matrices using cblas_dgemm, and fill mlist with the results
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_aa, n_aa, n_aa, (double) 1.0, cp, n_aa, cp_inv_exp, n_aa, (double) 0.0, &arr_mlist[index], n_aa);
                    for (index2 = index; index2 < index + n_aa2; index2++) {
                        arr_mlist[index2] = fabs(arr_mlist[index2]);
                    }
                    index += n_aa2;
#else
                    // multiply the matrices in pure C code, and fill mlist with the results
                    for (x = 0; x < n_aa2; x += n_aa) {
                        for (y = 0; y < n_aa; y++) {
                            mxy = 0.0;
                            yindex = y;
                            for (irowcolumn = x; irowcolumn < x + n_aa; irowcolumn++) {
                                mxy += cp[irowcolumn] * cp_inv_exp[yindex];
                                yindex += n_aa;
                            }
                            arr_mlist[index++] = fabs(mxy);
                        }
                    }
#endif
                } else { // we need to compute entries in mlist only for x = only_need_i
                    only_need_i_n_aa = only_need_i * n_aa;
                    index += only_need_i_n_aa;
#ifdef USE_ACCELERATE_CBLAS
                    // do the matrix vector multiplication using cblas, and put results in mlist
                    cblas_dgemv(CblasRowMajor, CblasTrans, n_aa, n_aa, (double) 1.0, cp_inv_exp, n_aa, &cp[only_need_i_n_aa], 1, (double) 0.0, &arr_mlist[index], 1); 
                    for (index2 = index; index2 < index + n_aa2 - only_need_i_n_aa; index2++) {
                        arr_mlist[index2] = fabs(arr_mlist[index2]);
                    }
                    index += n_aa2 - only_need_i_n_aa;
#else
                    // do the matrix vector multiplication in pure C code, and put results in mlist
                    for (y = 0; y < n_aa; y++) {
                        mxy = 0.0;
                        yindex = y;
                        for (irowcolumn = only_need_i_n_aa; irowcolumn < only_need_i_n_aa + n_aa; irowcolumn++) {
                            mxy += cp[irowcolumn] * cp_inv_exp[yindex];
                            yindex += n_aa;
                        }
                        arr_mlist[index++] = fabs(mxy);
                    }
                    index += n_aa2 - only_need_i_n_aa - n_aa;
#endif
                }
            }
        }
        free(cp_inv_exp);
    } else if (array_type == NPY_CDOUBLE) { // array is complex doubles
        complex_cp_inv_exp = (complex double *) malloc(complex_sizeof_cp_inv); // allocate memory
        complex_naa2_list = (complex double *) malloc(n_aa2 * sizeof(complex double)); // allocate memory
        complex_naa_list = (complex double *) malloc(n_aa * sizeof(complex double)); // allocate memory
        // Note that these next assignments assume that the arrays are C-style contiguous
        complex_cp = PyArray_DATA(p);
        complex_cp_inv = PyArray_DATA(p_inv);
        complex_cgr_diag = PyArray_DATA(gr_diag);
        for (i_ut = 0; i_ut < n_uts; i_ut++) { // loop over ut values
            only_need_i = PyInt_AS_LONG(PyList_GET_ITEM(only_need, only_need_index)); // value of only_need
            only_need_index++;
            if (only_need_i == -1) { // we don't need to do anything for these entries in mlist
                index += n_aa2;
            } else { // we need to compute at least some entries in mlist
                ut = PyFloat_AS_DOUBLE(PyList_GET_ITEM(uts, i_ut));  // ut value
                for (irowcolumn = 0; irowcolumn < n_aa; irowcolumn++) {
                    complex_exp_ut = cexp(ut * complex_cgr_diag[irowcolumn]);
                    for (y = irowcolumn * n_aa; y < (irowcolumn + 1) * n_aa; y++) {
                        complex_cp_inv_exp[y] = complex_cp_inv[y] * complex_exp_ut;
                    }
                }
                if (only_need_i == -2) { // we need to compute all of these entries in mlist
#ifdef USE_ACCELERATE_CBLAS
                    // multiply the matrices using cblas_zgemm, and fill mlist with the results
                    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_aa, n_aa, n_aa, &complex_one, complex_cp, n_aa, complex_cp_inv_exp, n_aa, &complex_zero, complex_naa2_list, n_aa); 
                    for (index2 = 0; index2 < n_aa2; index2++) {
                        arr_mlist[index + index2] = cabs(complex_naa2_list[index2]);
                    }
                    index += n_aa2;
#else
                    // multiply the matrices in pure C code, and fill mlist with the results
                    for (x = 0; x < n_aa2; x += n_aa) {
                        for (y = 0; y < n_aa; y++) {
                            complex_mxy = (complex double) 0.0;
                            yindex = y;
                            for (irowcolumn = x; irowcolumn < x + n_aa; irowcolumn++) {
                                complex_mxy += complex_cp[irowcolumn] * complex_cp_inv_exp[yindex];
                                yindex += n_aa;
                            }
                            arr_mlist[index++] = cabs(complex_mxy);
                        }
                    }
#endif
                } else { // we need to compute entries in mlist only for x = only_need_i
                    only_need_i_n_aa = only_need_i * n_aa;
                    index += only_need_i_n_aa;
#ifdef USE_ACCELERATE_CBLAS
                    // do the matrix vector multiplication using cblas, and put results in mlist
                    cblas_zgemv(CblasRowMajor, CblasTrans, n_aa, n_aa, &complex_one, complex_cp_inv_exp, n_aa, &complex_cp[only_need_i_n_aa], 1, &complex_zero, complex_naa_list, 1);  
                    for (index2 = 0; index2 < n_aa; index2++) {
                        arr_mlist[index + index2] = cabs(complex_naa_list[index2]);
                    }
                    index += n_aa2 - only_need_i_n_aa;
#else
                    // do the matrix vector multiplication in pure C code, and put results in mlist
                    for (y = 0; y < n_aa; y++) {
                        complex_mxy = (complex double) 0.0;
                        yindex = y;
                        for (irowcolumn = only_need_i_n_aa; irowcolumn < only_need_i_n_aa + n_aa; irowcolumn++) {
                            complex_mxy += complex_cp[irowcolumn] * complex_cp_inv_exp[yindex];
                            yindex += n_aa;
                        }
                        arr_mlist[index++] = cabs(complex_mxy);
                    }
                    index += n_aa2 - only_need_i_n_aa - n_aa;
#endif
                }
            }
        }
        free(complex_cp_inv_exp);
        free(complex_naa2_list);
        free(complex_naa_list);
    } else { // array is of neither real nor complex doubles
        PyErr_SetString(PyExc_ValueError, "gr entries are neither double nor complex doubles.");
        return NULL;
    }
    return PyInt_FromLong((long) 1);
}
//
// Documentation string for dMList function
static char dMList_docs[] = "This method computes a list of the mutation probabilities.\n\nThis method is a fast C implementation of the '_dMList' method of ddg_inference.DDGSet objects.\nIts calculations parallel that of ddg_inference.DDGSet.dMList.\nIt takes the following calling variables (in order):\n'uts' is a Python list of floats with meaning described in ddg_inference.DDGSet._dMList doc string.\n'only_need' is a Python list of integers with meaning described in ddg_inference.DDGSet._dMList doc string.\n\n'n_aa' is a python integer giving the number of amino acid types (typically 20).\n'length' is a python integer giving the length of the protein.\n'grs' is the _grs variables of the DDGSet.  It is a list of tuples giving the diagonalized gr matrices.\n\tEach tuple must contain matrices that are either all of the type 'float_' (C type double) or all of\n\tthe type 'complex_' (C type double complex).\n'dmlist' is a numpy ndarray of type 'float_' (C type double) of the appropriate length, and in C-style\n\tcontiguous format.  This array is filled with the appropriate values, as described in DDGSet._dMList.\n'residue_to_compute' is the number of the residue for which we compute the 'dMList'.\n'iwt' is the integer code for the wildtype residue (no derivative is computed for this ddG)\n'brs' is a list of length 'n_aa', with entry 'xi' giving the Br matrix with respect to\n\tthe ddG value for mutating from the wildtype to 'xi'; undefined for 'xi == 'iwt'.  Brx is defined\n\tby Brx = P_inv dGr/dz P where Gr = P Gr_diag P_inv.\nAll numpy ndarray arguments must be C-style contiguous.";

// This dMList function creates a list of the dM values giving the probabilities of mutations
static PyObject *dMList(PyObject *self, PyObject *args) {
    // Calling variables are (in order): uts, only_need, n_aa, length, grs, dmlist, residue_to_compute, iwt, brs
    PyObject *uts, *only_need, *grs, *r_grs_tuple, *gr_diag, *p, *p_inv, *dmlist, *brs, *brz;
    long n_aa, length, n_aa2, n_uts, residue_to_compute, i_ut, x, y, index, only_need_index, only_need_i, only_need_i_n_aa, index2, iwt, n_aa3, z;
    double *arr_dmlist, *cp_inv, *cp, *cgr_diag, *cexpd, *cbrz, *cvrz, *cvrz_p_inv;
    complex double *complex_cp_inv, *complex_cp, *complex_cgr_diag, *complex_cexpd, *complex_naa2_list, *complex_naa_list, *complex_cvrz, *complex_cvrz_p_inv, *complex_cbrz;
    double ut, exp_utdx, exp_utdy, dx, dy;
    complex double complex_exp_utdx, complex_exp_utdy, complex_dx, complex_dy;
    int array_type;
#ifdef USE_ACCELERATE_CBLAS
    complex double complex_one = 1, complex_zero = 0;
#else
    complex double complex_dmxy, complex_v_p_inv_xy;
    double dmxy, v_p_inv_xy;
    long yindex, irowcolumn;
#endif
    // Parse the arguments.  
    if (! PyArg_ParseTuple( args, "O!O!llO!O!llO!", &PyList_Type, &uts, &PyList_Type, &only_need, &n_aa, &length, &PyList_Type, &grs, &PyArray_Type, &dmlist, &residue_to_compute, &iwt, &PyList_Type, &brs)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to dMList.");
        return NULL;
    }
    // Error checking on arguments
    if (length < 1) { // length of the protein
        PyErr_SetString(PyExc_ValueError, "length is less than one.");
        return NULL;
    }
    if (n_aa < 1) { // number of amino acids.  Normally will be 20.
        PyErr_SetString(PyExc_ValueError, "n_aa is less than one.");
        return NULL;
    }
    if (PyList_GET_SIZE(grs) != length) { // make sure grs is of the same size as length
        PyErr_SetString(PyExc_ValueError, "grs is not of the same size as length.");
        return NULL;
    }
    n_uts = PyList_GET_SIZE(uts); // number of entries in uts
    if (n_uts < 1) { // make sure there are entries in uts 
        PyErr_SetString(PyExc_ValueError, "uts has no entries.");
        return NULL;
    }
    if (PyList_GET_SIZE(only_need) != n_uts * length) { // make sure only_need is of correct size
        PyErr_SetString(PyExc_ValueError, "only_need is of wrong size.");
        return NULL;
    }
    if (! ((residue_to_compute >= 0) && (residue_to_compute < length))) {
        PyErr_SetString(PyExc_ValueError, "Invalid value for residue_to_compute.");
        return NULL;
    }
    if (! ((iwt >= 0) && (iwt < n_aa))) {
        PyErr_SetString(PyExc_ValueError, "Invalid value for iwt.");
        return NULL;
    }
    if (PyList_GET_SIZE(brs) != n_aa) { // make sure brs has one entry for each amino acid
        PyErr_SetString(PyExc_ValueError, "brs is not of length equal to n_aa");
        return NULL;
    }
    n_aa2 = n_aa * n_aa; // square of the number of amino acids
    n_aa3 = n_aa2 * n_aa; // cube of the number of amino acids
    // The results will be returned in a numpy ndarray 'float_' (C type double) array called dmlist.
    // This array will be of size length * n_uts * n_aa3.
    arr_dmlist = (double *) PyArray_DATA(dmlist); // this is the data array of dmlist
    long const sizeof_cexpd = n_aa2 * sizeof(double);
    long const complex_sizeof_cexpd = n_aa2 * sizeof(complex double);
    // gr_diag, p, and p_inv are the eigenvalues, left, and right diagonalizing matrices of gr
    r_grs_tuple = PyList_GET_ITEM(grs, residue_to_compute);
    gr_diag = PyTuple_GET_ITEM(r_grs_tuple, 0); 
    p = PyTuple_GET_ITEM(r_grs_tuple, 1);
    p_inv = PyTuple_GET_ITEM(r_grs_tuple, 2);
    // Now begin filling arr_dmlist with the appropriate values
    index = 0;
    only_need_index = residue_to_compute * n_uts;
    // determine if these arrays are complex double or real doubles
    array_type = PyArray_TYPE(gr_diag);
    if (array_type == NPY_DOUBLE) { // array is of doubles, real not complex
        // Note that these next assignments assume that the arrays are C-style contiguous
        cp = PyArray_DATA(p);
        cp_inv = PyArray_DATA(p_inv);
        cgr_diag = PyArray_DATA(gr_diag);
        cexpd = (double *) malloc(sizeof_cexpd); // allocate memory
        cvrz = (double *) malloc(sizeof_cexpd); // allocate memory
        cvrz_p_inv = (double *) malloc(sizeof_cexpd); // allocate memory
        for (i_ut = 0; i_ut < n_uts; i_ut++) { // loop over ut values
            only_need_i = PyInt_AS_LONG(PyList_GET_ITEM(only_need, only_need_index)); // value of only_need
            only_need_index++;
            if (only_need_i == -1) { // we don't need to do anything for these entries in dmlist
                index += n_aa3;
            } else { // we need to compute at least some entries in dmlist
                ut = PyFloat_AS_DOUBLE(PyList_GET_ITEM(uts, i_ut));  // ut value
                // Entries of cexpd are defined by D_xy = (exp(ut d_x) - exp(ut d_y)) / (d_x - d_y)
                // for x != y, and D_yy = ut exp(d_x ut)
                index2 = 0;
                for (x = 0; x < n_aa; x++) {
                    dx = cgr_diag[x];
                    exp_utdx = exp(ut * dx);
                    for (y = 0; y < n_aa; y++) {
                        if (y == x) {
                            cexpd[index2] = ut * exp_utdx;
                        } else {
                            dy = cgr_diag[y];
                            exp_utdy = exp(ut * dy);
                            cexpd[index2] = (exp_utdx - exp_utdy) / (dx - dy);
                        }
                        index2++;
                    }
                }
                for (z = 0; z < n_aa; z++) {
                    // compute derivative with respect to z
                    if (z == iwt) { // don't compute values with respect to wildtype
                        index += n_aa2;
                        continue;
                    }
                    brz = PyList_GET_ITEM(brs, z);
                    cbrz = PyArray_DATA(brz);
                    // cvrz is the element-by-element product of cbrz and cexpd
                    for (index2 = 0; index2 < n_aa2; index2++) {
                        cvrz[index2] = cbrz[index2] * cexpd[index2];
                    }
#ifdef USE_ACCELERATE_CBLAS
                    // multiply the matrices using cblas_dgemm
                    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_aa, n_aa, n_aa, (double) 1.0, cvrz, n_aa, cp_inv, n_aa, (double) 0.0, cvrz_p_inv, n_aa); // multiply cvrz and cp_inv into cvrz_p_inv
#else
                    // multiply the matrices in pure C code
                    // multiply cvrz and cp_inv into cvrz_p_inv
                    index2 = 0;
                    for (x = 0; x < n_aa2; x += n_aa) {
                        for (y = 0; y < n_aa; y++) {
                            v_p_inv_xy = 0.0;
                            yindex = y;
                            for (irowcolumn = x; irowcolumn < x + n_aa; irowcolumn++) {
                                v_p_inv_xy += cvrz[irowcolumn] * cp_inv[yindex];
                                yindex += n_aa;
                            }
                            cvrz_p_inv[index2++] = v_p_inv_xy;
                        }
                    }
#endif
                    if (only_need_i == -2) { // we need to compute all of these entries in dmlist
#ifdef USE_ACCELERATE_CBLAS
                        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_aa, n_aa, n_aa, (double) 1.0, cp, n_aa, cvrz_p_inv, n_aa, (double) 0.0, &arr_dmlist[index], n_aa); // multiply cp and cvrz_p_inv into arr_dmlist[index : index + n_aa2]
                        index += n_aa2;
#else
                        // multiply the matrices in pure C code, and fill dmlist with the results
                        // multiply cp and cvrz_p_inv into arr_dmlist[index : index + n_aa2]
                        for (x = 0; x < n_aa2; x += n_aa) {
                            for (y = 0; y < n_aa; y++) {
                                dmxy = 0.0;
                                yindex = y;
                                for (irowcolumn = x; irowcolumn < x + n_aa; irowcolumn++) {
                                    dmxy += cp[irowcolumn] * cvrz_p_inv[yindex];
                                    yindex += n_aa;
                                }
                                arr_dmlist[index++] = dmxy;
                            }
                        }
#endif
                    } else { // we need to compute entries in dmlist only for x = only_need_i
                        only_need_i_n_aa = only_need_i * n_aa;
                        index += only_need_i_n_aa;
#ifdef USE_ACCELERATE_CBLAS
                        // do the matrix vector multiplication using cblas, and put results in dmlist
                        cblas_dgemv(CblasRowMajor, CblasTrans, n_aa, n_aa, (double) 1.0, cvrz_p_inv, n_aa, &cp[only_need_i_n_aa], 1, (double) 0.0, &arr_dmlist[index], 1); 
                        index += n_aa2 - only_need_i_n_aa;
#else
                        // do the matrix vector multiplication in pure C code, and put results in mlist
                        for (y = 0; y < n_aa; y++) {
                            dmxy = 0.0;
                            yindex = y;
                            for (irowcolumn = only_need_i_n_aa; irowcolumn < only_need_i_n_aa + n_aa; irowcolumn++) {
                                dmxy += cp[irowcolumn] * cvrz_p_inv[yindex];
                                yindex += n_aa;
                            }
                            arr_dmlist[index++] = dmxy;
                        }
                        index += n_aa2 - only_need_i_n_aa - n_aa;
#endif
                    }
                }
            }
        }
        free(cexpd);
        free(cvrz);
        free(cvrz_p_inv);
    } else if (array_type == NPY_CDOUBLE) { // array is of complex doubles
        // Note that these next assignments assume that the arrays are C-style contiguous
        complex_cp = PyArray_DATA(p);
        complex_cp_inv = PyArray_DATA(p_inv);
        complex_cgr_diag = PyArray_DATA(gr_diag);
        complex_cexpd = (complex double *) malloc(complex_sizeof_cexpd); // allocate memory
        complex_cvrz = (complex double *) malloc(complex_sizeof_cexpd); // allocate memory
        complex_cvrz_p_inv = (complex double *) malloc(complex_sizeof_cexpd); // allocate memory
        complex_naa2_list = (complex double *) malloc(n_aa2 * sizeof(complex double)); // allocate memory
        complex_naa_list = (complex double *) malloc(n_aa * sizeof(complex double)); // allocate memory
        for (i_ut = 0; i_ut < n_uts; i_ut++) { // loop over ut values
            only_need_i = PyInt_AS_LONG(PyList_GET_ITEM(only_need, only_need_index)); // value of only_need
            only_need_index++;
            if (only_need_i == -1) { // we don't need to do anything for these entries in dmlist
                index += n_aa3;
            } else { // we need to compute at least some entries in dmlist
                ut = PyFloat_AS_DOUBLE(PyList_GET_ITEM(uts, i_ut));  // ut value
                // Entries of cexpd are defined by D_xy = (exp(ut d_x) - exp(ut d_y)) / (d_x - d_y)
                // for x != y, and D_yy = ut exp(d_x ut)
                index2 = 0;
                for (x = 0; x < n_aa; x++) {
                    complex_dx = complex_cgr_diag[x];
                    complex_exp_utdx = cexp(ut * complex_dx);
                    for (y = 0; y < n_aa; y++) {
                        if (y == x) {
                            complex_cexpd[index2] = ut * complex_exp_utdx;
                        } else {
                            complex_dy = complex_cgr_diag[y];
                            complex_exp_utdy = cexp(ut * complex_dy);
                            complex_cexpd[index2] = (complex_exp_utdx - complex_exp_utdy) / (complex_dx - complex_dy);
                        }
                        index2++;
                    }
                }
                for (z = 0; z < n_aa; z++) {
                    // compute derivative with respect to z
                    if (z == iwt) { // don't compute values with respect to wildtype
                        index += n_aa2;
                        continue;
                    }
                    brz = PyList_GET_ITEM(brs, z);
                    complex_cbrz = PyArray_DATA(brz);
                    // cvrz is the element-by-element product of cbrz and cexpd
                    for (index2 = 0; index2 < n_aa2; index2++) {
                        complex_cvrz[index2] = complex_cbrz[index2] * complex_cexpd[index2];
                    }
#ifdef USE_ACCELERATE_CBLAS
                    // multiply the matrices using cblas_zgemm
                    cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_aa, n_aa, n_aa, &complex_one, complex_cvrz, n_aa, complex_cp_inv, n_aa, &complex_zero, complex_cvrz_p_inv, n_aa); // multiply cvrz and cp_inv into cvrz_p_inv
#else
                    // multiply the matrices in pure C code
                    // multiply cvrz and cp_inv into cvrz_p_inv
                    index2 = 0;
                    for (x = 0; x < n_aa2; x += n_aa) {
                        for (y = 0; y < n_aa; y++) {
                            complex_v_p_inv_xy = 0.0;
                            yindex = y;
                            for (irowcolumn = x; irowcolumn < x + n_aa; irowcolumn++) {
                                complex_v_p_inv_xy += complex_cvrz[irowcolumn] * complex_cp_inv[yindex];
                                yindex += n_aa;
                            }
                            complex_cvrz_p_inv[index2++] = complex_v_p_inv_xy;
                        }
                    }
#endif
                    if (only_need_i == -2) { // we need to compute all of these entries in dmlist
#ifdef USE_ACCELERATE_CBLAS
                        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n_aa, n_aa, n_aa, &complex_one, complex_cp, n_aa, complex_cvrz_p_inv, n_aa, &complex_zero, complex_naa2_list, n_aa); // multiply cp and cvrz_p_inv into arr_dmlist[index : index + n_aa2]
                        for (index2 = 0; index2 < n_aa2; index2++) {
                            arr_dmlist[index + index2] = creal(complex_naa2_list[index2]);
                        }
                        index += n_aa2;
#else
                        // multiply the matrices in pure C code, and fill dmlist with the results
                        // multiply cp and cvrz_p_inv into arr_dmlist[index : index + n_aa2]
                        for (x = 0; x < n_aa2; x += n_aa) {
                            for (y = 0; y < n_aa; y++) {
                                complex_dmxy = 0.0;
                                yindex = y;
                                for (irowcolumn = x; irowcolumn < x + n_aa; irowcolumn++) {
                                    complex_dmxy += complex_cp[irowcolumn] * complex_cvrz_p_inv[yindex];
                                    yindex += n_aa;
                                }
                                arr_dmlist[index++] = creal(complex_dmxy);
                            }
                        }
#endif
                    } else { // we need to compute entries in dmlist only for x = only_need_i
                        only_need_i_n_aa = only_need_i * n_aa;
                        index += only_need_i_n_aa;
#ifdef USE_ACCELERATE_CBLAS
                        // do the matrix vector multiplication using cblas, and put results in dmlist
                        cblas_zgemv(CblasRowMajor, CblasTrans, n_aa, n_aa, &complex_one, complex_cvrz_p_inv, n_aa, &complex_cp[only_need_i_n_aa], 1, &complex_zero, complex_naa_list, 1); 
                        for (index2 = 0; index2 < n_aa; index2++) {
                            arr_dmlist[index + index2] = creal(complex_naa_list[index2]);
                        }
                        index += n_aa2 - only_need_i_n_aa;
#else
                        // do the matrix vector multiplication in pure C code, and put results in mlist
                        for (y = 0; y < n_aa; y++) {
                            complex_dmxy = 0.0;
                            yindex = y;
                            for (irowcolumn = only_need_i_n_aa; irowcolumn < only_need_i_n_aa + n_aa; irowcolumn++) {
                                complex_dmxy += complex_cp[irowcolumn] * complex_cvrz_p_inv[yindex];
                                yindex += n_aa;
                            }
                            arr_dmlist[index++] = creal(complex_dmxy);
                        }
                        index += n_aa2 - only_need_i_n_aa - n_aa;
#endif
                    }
                }
            }
        }
        free(complex_cexpd);
        free(complex_cvrz);
        free(complex_cvrz_p_inv);
        free(complex_naa2_list);
        free(complex_naa_list);
    } else { // array is of neither real nor complex doubles
        PyErr_SetString(PyExc_ValueError, "matrices are neither double nor complex doubles.");
        return NULL;
    }
    return PyInt_FromLong((long) 1);
}
//
// Documentation string for ComputeGr function
static char ComputeGr_docs[] = "This is C extension to compute the Gr matrices.\n\nThis method implements code expressed in python in the ddg_inference.DDGSet._ComputeGr method.\nIt returns the numpy ndarray 'gr'.\nCalling variables are (in order):\n'n_aa' a python integer giving the number of amino acids (typically 20).\n'ddgs' the python list indexed_ddgs defined in the DDGSet._ComputeGr method.\n'p0sums' the python list of this name defined in the DDGSet._ComputeGr method.\n'indexed_c_set' the python list of this name defined in the DDGSet._ComputeGr method.";

// This ComputeGr function creates the Gr matrix for a residue
static PyArrayObject *ComputeGr(PyObject *self, PyObject *args) {
    // Calling variables are (in order): n_aa, ddgs, p0sums, indexed_c_set
    long n_aa, xi, yi, index, ddgrxy;
    PyObject *ddgs, *p0sums, *indexed_c_set, *gr;
    double *arr_gr;
    double gryy, grxy, frxy;
    if (! PyArg_ParseTuple( args, "lO!O!O!", &n_aa, &PyList_Type, &ddgs, &PyList_Type, &p0sums, &PyList_Type, &indexed_c_set)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to ComputeGr.");
        return NULL;
    }
    npy_intp dims[] = {n_aa, n_aa}; // dimensions of gr
    gr = PyArray_SimpleNew(2, dims, NPY_DOUBLE);   // create gr as a numpy ndarray
    arr_gr = (double *) PyArray_DATA(gr); // this is the data array of gr
    index = 0; // initial index to ddgs and indexed_c_set
    for (yi = 0; yi < n_aa; yi++) {
        gryy = 0.0;
        for (xi = 0; xi < n_aa; xi++) {
            if (xi == yi) {
                index++;
            }
            else {
                ddgrxy = PyInt_AS_LONG(PyList_GET_ITEM(ddgs, index));
                if (ddgrxy <= 0) {
                    frxy = 1.0;
                }
                else {
                    ddgrxy--;
                    frxy = PyFloat_AS_DOUBLE(PyList_GET_ITEM(p0sums, ddgrxy));
                }
                grxy = frxy * PyFloat_AS_DOUBLE(PyList_GET_ITEM(indexed_c_set, index));
                index++;
                gryy -= grxy;
                arr_gr[xi * n_aa + yi] = grxy;
            }
        }
        arr_gr[yi * n_aa + yi] = gryy;
    }
    return (PyArrayObject *) gr; // return the computed gr matrix
}

// Module documentation string
static char cddg_inference_docs[] = "This module contains fast implementations of bottleneck calculations of the ddg_inference.DDGSet object.\n\nAll of the methods implemented here closely parallel the pure Python implementations offered by the\n\tDDGSet object.  However, these methods are substantially faster.  The methods are:\n'MList' which does the work of DDGSet.MList\n'dMList' which does the work of DDGSet.dMList'\n'ComputeGr' which does some of the work of DDGSet._ComputeGr";

// The module functions
static PyMethodDef cddg_inference_funcs[] = {
    {"MList", (PyCFunction) MList, METH_VARARGS, MList_docs},
    {"dMList", (PyCFunction) dMList, METH_VARARGS, dMList_docs},
    {"ComputeGr", (PyCFunction) ComputeGr, METH_VARARGS, ComputeGr_docs},
    {NULL}
};

// Initialize the module
void initcddg_inference(void) {
    Py_InitModule3("cddg_inference", cddg_inference_funcs, cddg_inference_docs);
    import_array(); // this command is necessary to use the numarray C API
}
