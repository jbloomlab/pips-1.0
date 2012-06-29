// This file implements a fast C version of a python function written by Jesse Bloom, 2007
// It is designed for fast implementation of phylogenetic likelihood calculations.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <numpy/arrayobject.h>
//
// GLOBAL VARIABLES
// underflow is an integer specifying how often we scale the conditional probabilities
//    to avoid numerical underflow.
// n_internal is the integer number of internal nodes in the tree (usual index is i).
// n_tips is the integer number of tips in the tree (usual index is t).
// n_aa is the integer number of amino acids (usual index is x or y).  Each amino
//    acid type is encoded by an integer x, 0 <= x < n_aa.  
// length is the integer sequence length.  Residues are numbered by the usual index
//    r, where 0 <= r < length.
// tips is an integer array of length * n_tips, with tips[r * n_tips + t] (0 <= t < n_tips,
//    0 <= r < length) giving the identity of amino acid r at tip t as x, where 0 <= x < n_aa.
//    If an amino acid is a gap, it is encoded as -1.
// descendents and descendent_is_tip are both integer arrays of length 2 * n_internal.         
//    They give information about the descendents of internal node i, where 0 <= i < n_aa.
//    If descendent_is_tip[2 * i] is 1, then the right descendent of i is a tip node.  In
//    this case, the identity of the amino acid at residue r of this descendent is given by
//    element r * n_tips + descendents[2 * i] of tips (so 0 <= descendents[2 * i] < n_tips).
//    If descendent_is_tip[2 * i] is 0, then the right descendent of i is not a tip node.  In
//    this case, descendents[2 * i] gives the number j of the internal node (0 <= j < n_internal)
//    that is the right descendent of node j.  Crucially, the nodes are numbered so that 
//    j will always be less than i.  The information about the left descendent is stored
//    similarly, but in element 2 * i + 1 of descendents and descendent_is_tip.  Because of
//    the ordering where j < i, this guarantees that the last entry (i = n_aa - 1) corresponds
//    to the root node of the tree.
// n_unique is the integer number of unique branch lengths, which is <= n_tips + n_internal
// globals_initialized is a boolean switch that is 0 if the global variables are not
//    initialized, and 1 otherwise.
long globals_initialized = 0;
long n_internal, n_tips, n_aa, length, n_unique, underflow;
long *tips, *descendents, *descendent_is_tip, *unique_uts_index_n_aa2;

// Documentation string for Initialize method
static char Initialize_docs[] = "Initializes the module with variables about the phylogenetic tree.\n\nThis function must be called before using the ResidueLogLikelihood or ResidueLogLokelihoodAndDeriv function.  It specifies information about \n\tthe phylogenetic tree data that is used by the ResidueLogLikelihood and ResidueLogLikelihoodAndDeriv function.\n\tTo reset this data, call this function again.\nThe arguments are, in order:\n\tlength - an integer specifying the length of the protein sequences.\n\tn_aa - an integer giving the number of different amino acids, typically 20.\n\tn_tips - an integer giving the number of tip nodes in the phylogenetic tree.\n\tn_internal - an integer giving the number of internal nodes in the phylogenetic tree.\n\tdescendents - see explanation below for descendent_is_tip\n\ttips - a list of length * n_tips integer entries.  This list gives information about the sequences\n\t\tof the tip nodes.  Element r * n_tips + t gives the amino acid at  residue r (0 <= r < length) of tip\n\t\tnode t (0 <= t < n_tips).  The amino acids are encoded as integers x with 0 <= x < n_aa.\n\t\tEach amino acid is indicated by a unique integer in this range.  Gaps are indicated by -1.\n\tdescendent_is_tip - both descendents and descendent_is_tip are lists of length 2 * n_internal.\n\t\tThey give information about the descendents of the internal nodes.  Information about internal node i\n\t\t(where 0 <= i < n_internal).  If the right descendent of node i is a tip node, then\n\t\tdescendent_is_tip[2 * i] is the integer 1.  In this case, descendents[2 * i] gives the number j\n\t\t(where 0 <= j < i) of the internal node that is the right descendent of i.  Note that the nodes\n\t\tmust be numbered so that j < i if j is a descendent of i.  If descendent_is_tip[2 * i] is the\n\t\tinteger 0, then descendents[2 * i] is the number of the tip node t (0 <= t < n_tips) that is the right\n\t\tdescendent of i.  The information about the left descendent of node i is stored in elements 2 * i + 1\n\t\tof descendent_is_tip and descendents.\n\t\n\t\tprobabilities of the amino acids at these positions.\n\n\tn_unique - the integer number of unique branch lengths, which will always be <= n_tips + n_internal.\n\tunique_uts_index_n_aa2 - a list of n_tips + n_internal integers.  This specifies the indices that will be\n\t\tused for different tips in the mlist calling variable.\n\tunderflow - an integer specifying how often we scale the conditional probabilities to avoid underflow.\n\t\tA larger value specifies that we scale more often to avoid underflow.  Must be > 0.\nIf the function executes successfully, it returns 1.  Otherwise it raises an exception.";

// This Initialize function initializes the global variables accessed by the module
static PyObject *Initialize(PyObject *self, PyObject *args) {
    // Calling variables are (in order): length, n_aa, n_tips, n_internal, tips, descendents, descendent_is_tip, n_unique, unique_uts_index_n_aa2, underflow
    PyObject *call_tips, *call_descendents, *call_descendent_is_tip, *call_unique_uts_index_n_aa2, *entry;
    long i, ientry, n_aa2;
    // Parse the arguments.  This assigns values to length, n_aa, n_tips, and n_internal
    if (! PyArg_ParseTuple( args, "llllO!O!O!lO!l", &length, &n_aa, &n_tips, &n_internal, &PyList_Type, &call_tips, &PyList_Type, &call_descendents, &PyList_Type, &call_descendent_is_tip, &n_unique, &PyList_Type, &call_unique_uts_index_n_aa2, &underflow)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to Initialize.");
        return NULL;
    }
    // Error checking on non-list arguments
    if (length < 1) {
        PyErr_SetString(PyExc_ValueError, "length is less than one.");
        return NULL;
    }
    if (n_aa < 1) {
        PyErr_SetString(PyExc_ValueError, "n_aa is less than one.");
        return NULL;
    }
    n_aa2 = n_aa * n_aa;
    if (n_tips < 1) {
        PyErr_SetString(PyExc_ValueError, "n_tips is less than one.");
        return NULL;
    }
    if (n_internal < 1) {
        PyErr_SetString(PyExc_ValueError, "n_internal is less than one.");
        return NULL;
    }
    if ((n_unique < 1) || (n_unique > n_internal + n_tips)) {
        PyErr_SetString(PyExc_ValueError, "n_unique is not greater than zero and less than n_tips + n_internal.");
        return NULL;
    }
    if (underflow < 1) {
        PyErr_SetString(PyExc_ValueError, "underflow is not greater than zero.");
        return NULL;
    }
    if (globals_initialized) {
        // free existing allocated memory for global arrays
        if (tips == NULL || descendents == NULL || descendent_is_tip == NULL) {
            PyErr_SetString(PyExc_ValueError, "Some arguments not initialized.");
            return NULL;
        }
        free(tips);
        free(descendents);
        free(descendent_is_tip);
        free(unique_uts_index_n_aa2);
        globals_initialized = 0; 
    }
    // allocate memory for global arrays
    globals_initialized = 1;
    tips = (long *) malloc(length * n_tips * sizeof(long));
    descendents = (long *) malloc(2 * n_internal * sizeof(long));
    descendent_is_tip = (long *) malloc(2 * n_internal * sizeof(long));
    unique_uts_index_n_aa2 = (long *) malloc((n_tips + n_internal) * sizeof(long));
    // assign global variable array values
    // First, assign tips from call_tips
    if (PyList_Size(call_tips) != length * n_tips) {
        PyErr_SetString(PyExc_ValueError, "tips is of the wrong length.");
        return NULL;
    }
    for (i = 0; i < length * n_tips; i++) {
        entry = PyList_GetItem(call_tips, i);
        if (! PyInt_Check(entry)) {
            PyErr_SetString(PyExc_ValueError, "tips contains a non-integer entry.");
            return NULL;
        }
        ientry = PyInt_AsLong(entry);
        if (! (-1 <= ientry && ientry < n_aa)) {
            PyErr_SetString(PyExc_ValueError, "tips contains an integer entry out of the valid range.");
            return NULL;
        }
        tips[i] = ientry;
    }
    // Now assign descendent_is_tip from call_descendent_is_tip
    if (PyList_Size(call_descendent_is_tip) != 2 * n_internal) {
        PyErr_SetString(PyExc_ValueError, "descendent_is_tip is of the wrong length.");
        return NULL;
    }
    for (i = 0; i < 2 * n_internal; i++) {
        entry = PyList_GetItem(call_descendent_is_tip, i);
        if (! PyInt_Check(entry)) {
            PyErr_SetString(PyExc_ValueError, "descendent_is_tip contains a non-integer entry.");
            return NULL;
        }
        ientry = PyInt_AsLong(entry);
        if (! (ientry == 0 || ientry == 1)) {
            PyErr_SetString(PyExc_ValueError, "descendent_is_tip contains an integer other than 0 or 1.");
            return NULL;
        }
        descendent_is_tip[i] = ientry;
    }
    // Now assign descendents from call_descendents
    if (PyList_Size(call_descendents) != 2 * n_internal) {
        PyErr_SetString(PyExc_ValueError, "descendents is of the wrong length.");
        return NULL;
    }
    for (i = 0; i < 2 * n_internal; i++) {
        entry = PyList_GetItem(call_descendents, i);
        if (! PyInt_Check(entry)) {
            PyErr_SetString(PyExc_ValueError, "descendents contains a non-integer entry.");
            return NULL;
        }
        ientry = PyInt_AsLong(entry);
        if (! ((descendent_is_tip[i] && 0 <= ientry && ientry < n_tips) || (! descendent_is_tip[i] && 0 <= ientry && ientry < n_internal))) {
            PyErr_SetString(PyExc_ValueError, "descendents contains an integer entry out of the valid range.");
            return NULL;
        }
        descendents[i] = ientry;
    }
    // Now assign unique_uts_index_n_aa2 from call_unique_uts_index_n_aa2
    if (PyList_Size(call_unique_uts_index_n_aa2) != n_tips + n_internal) {
        PyErr_SetString(PyExc_ValueError, "unique_uts_index_n_aa2 is of the wrong length.");
        return NULL;
    }
    for (i = 0; i < n_tips + n_internal; i++) {
        entry = PyList_GetItem(call_unique_uts_index_n_aa2, i);
        if (! PyInt_Check(entry)) {
            PyErr_SetString(PyExc_ValueError, "unique_uts_index_n_aa2 contains a non-integer entry.");
            return NULL;
        }
        ientry = PyInt_AsLong(entry);
        if (! ((ientry >= 0) && (ientry < n_unique * n_aa2))) {
            PyErr_SetString(PyExc_ValueError, "unique_uts_index_n_aa2 contains an entry outside the valid range.");
            return NULL;
        }
        unique_uts_index_n_aa2[i] = ientry;
    }
    // Return 1
    return PyInt_FromLong((long) 1);
}

// Documentation string for IsInitialized function
static char IsInitialized_docs[] = "Returns True if function Initialize has already been called; False otherwise.";

// IsInitialized function
static PyObject *IsInitialized(PyObject *self, PyObject *args) {
    // There are no calling variables
    if (globals_initialized) { // Initialize has already been called
        Py_RETURN_TRUE;
    } else { // Initialize has not yet been called
        Py_RETURN_FALSE;
    }
}

// Documentation string for ResidueLogLikelihood function
static char ResidueLogLikelihood_docs[] = "Computes the log likelihood for a residue.\n\nThis is a C function that implements the same calculation that is performed by the '_ResidueLogLikelihood\n\tmethod of the 'DDGSet' class defined in the 'ddg_inference.py' module.\n\tThe calling variables are (in order) mlist, pi, and r. The information about the structure of the phylogenetic\n\ttree is stored in the global variables initialized by calling the 'Initialize' function.";

// ResidueLogLikelihood function
static PyObject *ResidueLogLikelihood(PyObject *self, PyObject *args) {
    // Calling variables are (in order): mlist, r
    PyObject *mlist, *pi;
    double *p;
    double rloglikelihood, piy, piyr, piyl, pr, r_log_scaling_sum, scale;
    long r, index, i, right_index, left_index, tr, tl, xr, xl, indexr, indexl, jl_times_n_aa, jr_times_n_aa, r_times_n_tips, x, y, i_times_n_aa, jl, jr, root, r_times_n_aa, root_times_n_aa, indexry, indexly;
    long const n_aa2 = n_aa * n_aa;
    // Make sure globals are initialized
    if (! globals_initialized) {
        PyErr_SetString(PyExc_ValueError, "Call to ResidueLogLikelihood without a previous call to Initialize.");
        return NULL;
    }
    // Parse the arguments. 
    if (! PyArg_ParseTuple( args, "O!O!l", &PyArray_Type, &mlist, &PyList_Type, &pi, &r)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to ResidueLogLikelihood.");
        return NULL;
    }
    // Make sure the array is of the type double -- this is just error checking
    if (PyArray_TYPE(mlist) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "mlist is not an array of doubles (type float_).");
        return NULL;
    }
    // Make sure the array has the proper number of elements -- this is just error checking
    // Note that we never check that they are correctly dimensioned (i.e. nd = 1);  I'm not sure how to do this
    if (PyArray_SIZE(mlist) != n_aa2 * n_unique) {
        PyErr_SetString(PyExc_ValueError, "mlist does not have the correct number of entries.");
        return NULL;
    }
    // Make sure that r specifies a valid residue
    if (! ((r >= 0) && (r < length))) {
        PyErr_SetString(PyExc_ValueError, "Invalid value of r.");
        return NULL;
    }
    // The way that we access elements of mlist is as follows
    //    double m;
    //    long i;
    //    m = * ((double *) PyArray_GETPTR1(mlist, i)); 
    // assigns m to be the value of element i of mlist.
    //
    // Now loop over all residues to compute the log likelihood
    // For each residue, the calculation basically consists of computing the elements of the double
    //    array p.  It is recomputed for reach residue.  p[i * n_aa + y] is the cumulative probability
    //    that internal node i is residue y given the information in the subtree below i.  Because of
    //    the ordering of the internal node numbering, the entries for i = n_aa - 1 correspond to the
    //    root node of the tree.  Note that if node j is part of the subtree of node i, j < i.  This
    //    is ensured by the node numbering that should have been passed to Initialize.
    // Allocate p
    p = (double *) malloc(n_internal * n_aa * sizeof(double)); // allocate p
    r_log_scaling_sum = 0.0; // sum of scaling factor corrections
    r_times_n_tips = r * n_tips;
    for (i = 0; i < n_internal; i++) { // i loops over all internal nodes
        right_index = 2 * i; // index of right descendent in descendents, descendent_is_tip
        left_index = right_index + 1; // index of left descendent in descendents, descendent_is_tip
        i_times_n_aa = i * n_aa;
        if (descendent_is_tip[right_index] && descendent_is_tip[left_index]) { // both descendents are tips
            tr = descendents[right_index];
            tl = descendents[left_index];
            xr = tips[r_times_n_tips + tr];
            xl = tips[r_times_n_tips + tl];
            if (xr == -1) { // a gap on right descendent
                if (xl == -1) { // both left and right descendents are gaps
                    for (y = 0; y < n_aa; y++) {
                        p[i_times_n_aa + y] = 1.0;
                    }
                }
                else { // gap is only on right descendent
                    index = unique_uts_index_n_aa2[tl] + xl * n_aa;
                    for (y = 0; y < n_aa; y++) {
                        p[i_times_n_aa + y] = * ((double *) PyArray_GETPTR1(mlist, index));
                        index++;
                    }
                }
            }
            else if (xl == -1) { // gap only on left descendent {
                index = unique_uts_index_n_aa2[tr] + xr * n_aa;
                for (y = 0; y < n_aa; y++) {
                    p[i_times_n_aa + y] = * ((double *) PyArray_GETPTR1(mlist, index));
                    index++;
                }
            }
            else { // no gap on either descendent
                indexr = unique_uts_index_n_aa2[tr] + xr * n_aa;
                indexl = unique_uts_index_n_aa2[tl] + xl * n_aa;
                for (y = 0; y < n_aa; y++) {
                    p[i_times_n_aa + y] = (* ((double *) PyArray_GETPTR1(mlist, indexr))) * (* ((double *) PyArray_GETPTR1(mlist, indexl)));
                    indexr++;
                    indexl++;
                }
            }
        }
        else if (descendent_is_tip[right_index] && ! descendent_is_tip[left_index]) { // only right descendent is tip
            tr = descendents[right_index];
            xr = tips[r_times_n_tips + tr];
            jl = descendents[left_index];
            jl_times_n_aa = jl * n_aa;
            if (xr == -1) { // right descendent is a gap
                indexl = unique_uts_index_n_aa2[n_tips + jl];
                for (y = 0; y < n_aa; y++) {
                    indexly = indexl + y;
                    piy = 0.0;
                    for (x = 0; x < n_aa; x++) {
                        piy += (* ((double *) PyArray_GETPTR1(mlist, indexly))) * p[jl_times_n_aa + x];
                        indexly += n_aa;
                    }
                    p[i_times_n_aa + y] = piy;
                }
            }
            else { // not a gap
                indexl = unique_uts_index_n_aa2[n_tips + jl];
                indexr = unique_uts_index_n_aa2[tr] + xr * n_aa;
                for (y = 0; y < n_aa; y++) {
                    indexly = indexl + y;
                    piy = 0.0;
                    for (x = 0; x < n_aa; x++) {
                        piy += (* ((double *) PyArray_GETPTR1(mlist, indexly))) * p[jl_times_n_aa + x];
                        indexly += n_aa;
                    }
                    p[i_times_n_aa + y] = piy * (* ((double *) PyArray_GETPTR1(mlist, indexr)));
                    indexr++;
                }
            }
        }
        else if (! descendent_is_tip[right_index] && descendent_is_tip[left_index]) { // only left descendent is tip
            jr = descendents[right_index];
            jr_times_n_aa = jr * n_aa;
            tl = descendents[left_index];
            xl = tips[r_times_n_tips + tl];
            if (xl == -1) { // left descendent is gap
                indexr = unique_uts_index_n_aa2[n_tips + jr];
                for (y = 0; y < n_aa; y++) {
                    indexry = indexr + y;
                    piy = 0.0;
                    for (x = 0; x < n_aa; x++) {
                        piy += (* ((double *) PyArray_GETPTR1(mlist, indexry))) * p[jr_times_n_aa + x];
                        indexry += n_aa;
                    }
                    p[i_times_n_aa + y] = piy;
                }
            }
            else { // not a gap
                indexr = unique_uts_index_n_aa2[n_tips + jr];
                indexl = unique_uts_index_n_aa2[tl] + xl * n_aa;
                for (y = 0; y < n_aa; y++) {
                    indexry = indexr + y;
                    piy = 0.0;
                    for (x = 0; x < n_aa; x++) {
                        piy += (* ((double *) PyArray_GETPTR1(mlist, indexry))) * p[jr_times_n_aa + x];
                        indexry += n_aa;
                    }
                    p[i_times_n_aa + y] = piy * (* ((double *) PyArray_GETPTR1(mlist, indexl)));
                    indexl++;
                }
            }
        }
        else { // both internal nodes
            jr = descendents[right_index];
            jr_times_n_aa = jr * n_aa;
            jl = descendents[left_index];
            jl_times_n_aa = jl * n_aa;
            indexr = unique_uts_index_n_aa2[n_tips + jr];
            indexl = unique_uts_index_n_aa2[n_tips + jl];
            for (y = 0; y < n_aa; y++) {
                piyr = 0.0;
                piyl = 0.0;
                indexry = indexr + y;
                indexly = indexl + y;
                for (x = 0; x < n_aa; x++) {
                    piyr += (* ((double *) PyArray_GETPTR1(mlist, indexry))) * p[jr_times_n_aa + x];
                    piyl += (* ((double *) PyArray_GETPTR1(mlist, indexly))) * p[jl_times_n_aa + x];
                    indexry += n_aa;
                    indexly += n_aa;
                }
                p[i_times_n_aa + y] = piyr * piyl;
            }
        }
        if (!((i + 1) % underflow)) { // scale the probabilities for the node to avoid underflow
            scale = p[i_times_n_aa]; // find the maximum value over all y, starting with y = 0
            for (y = 1; y < n_aa; y++) {
                piy = p[i_times_n_aa + y];
                if (piy > scale) {
                    scale = piy;
                }
            }
            if (scale != 0.0) { // we cannot scale by a scaling factor of 0
                for (y = 0; y < n_aa; y++) {
                    p[i_times_n_aa + y] /= scale;
                }
                r_log_scaling_sum += log(scale); // add to the log scaling factor correction
            }
        }
    }
    // we have now computed p for the root node.  Compute overall likelihood for this residue
    root = n_internal - 1; // index of the root node
    pr = 0.0;
    r_times_n_aa = r * n_aa;
    root_times_n_aa = root * n_aa;
    for (y = 0; y < n_aa; y++) {
        pr += p[root_times_n_aa + y] * PyFloat_AsDouble(PyList_GetItem(pi, r_times_n_aa + y));
    }
    // perform a primitive check for underflow.  Note that this not a very good check -- there could still be a lot of underflow not detected by this check.
    if (pr == 0.0) {
        char error_string[1000];
        sprintf(error_string, "Underflow occurred in the computation of the likelihood for residue %d.", (int) r);
        PyErr_SetString(PyExc_ValueError, error_string);
        return NULL;
    }
    rloglikelihood = log(pr);
    rloglikelihood += r_log_scaling_sum; // correct for scaling
    free(p); // free memory allocated for p
    // error checking to see if there is a problem associated with rloglikelihood being nan (not a number)
    if isnan(rloglikelihood) {
        fprintf(stderr, "\nError associated with residue log likelihood being nan for residue %d.\n", (int) r);
        fprintf(stderr, "The value of pr (the unscaled likelihood) is %f.\n", pr);
        fprintf(stderr, "The value of r_log_scaling_sum (the log scaling factor) is %f.\n", r_log_scaling_sum);
        // We do not raise an exception, since we assume this will be done by the calling function which
        // will also print additional diagnostic information.
    }
    // construct the return value from rloglikelihood, turning it into a Python float
    return PyFloat_FromDouble(rloglikelihood);
}

// Documentation string for ResidueLogLikelihoodAndDeriv function
static char ResidueLogLikelihoodAndDeriv_docs[] = "Computes the log likelihood and its derivative for a residue.\n\nThis is a C function that implements the same calculation that is performed by the '_ResidueLogLikelihoodAndDeriv\n\tmethod of the 'DDGSet' class defined in the 'ddg_inference.py' module.\n\tThe calling variables are (in order) mlist, dmlist, pi, dpi, r, iwt. The information about the structure of the phylogenetic\n\ttree is stored in the global variables initialized by calling the 'Initialize' function.\nReturned variable is 2-tuple (rloglikelihood, drloglikelihood) where drloglikelihood is list\n\twith element indexed by z giving the derivative of rloglikelihood with respect to ddG for residue z.";

// ResidueLogLikelihoodAndDeriv function
static PyObject *ResidueLogLikelihoodAndDeriv(PyObject *self, PyObject *args) {
    // Calling variables are (in order): mlist, dmlist, pi, dpi, r, iwt
    PyObject *mlist, *dmlist, *pi, *dpi, *drloglikelihood;
    double *p, *dp, *dpiyr, *dpiyl, *dpr;
    double rloglikelihood, piy, piyr, piyl, pr, r_log_scaling_sum, scale, mr, ml, mly, mry, pjlx, pjrx, piry, pry;
    long r, index, i, right_index, left_index, tr, tl, xr, xl, indexr, indexl, jl_times_n_aa, jr_times_n_aa, r_times_n_tips, x, y, z, i_times_n_aa, i_times_n_aa2, jl, jr, root, r_times_n_aa, root_times_n_aa, indexry, indexly, iwt, dp_index, dm_index, dm_indexr, dm_indexl, z_n_aa2, jl_times_n_aa2, jr_times_n_aa2, dm_indexly, dm_indexry, dpjlx_index, dpjrx_index, root_times_n_aa2, dpi_index;
    long const n_aa2 = n_aa * n_aa;
    // Make sure globals are initialized
    if (! globals_initialized) {
        PyErr_SetString(PyExc_ValueError, "Call to ResidueLogLikelihoodAndDeriv without a previous call to Initialize.");
        return NULL;
    }
    // Parse the arguments. 
    if (! PyArg_ParseTuple( args, "O!O!O!O!ll", &PyArray_Type, &mlist, &PyArray_Type, &dmlist, &PyList_Type, &pi, &PyList_Type, &dpi, &r, &iwt)) {
        PyErr_SetString(PyExc_TypeError, "Invalid calling arguments to ResidueLogLikelihoodAndDeriv.");
        return NULL;
    }
    // Make sure the arrays are of the type double -- this is just error checking
    if (PyArray_TYPE(mlist) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "mlist is not an array of doubles (type float_).");
        return NULL;
    }
    if (PyArray_TYPE(dmlist) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "mlist is not an array of doubles (type float_).");
        return NULL;
    }
    // Make sure the arrays have the proper number of elements -- this is just error checking
    // Note that we never check that they are correctly dimensioned (i.e. nd = 1);  I'm not sure how to do this
    if (PyArray_SIZE(mlist) != n_aa2 * n_unique) {
        PyErr_SetString(PyExc_ValueError, "mlist does not have the correct number of entries.");
        return NULL;
    }
    if (PyArray_SIZE(dmlist) != n_aa * n_aa2 * n_unique) {
        PyErr_SetString(PyExc_ValueError, "dmlist does not have the correct number of entries.");
        return NULL;
    }
    // Make sure that r specifies a valid residue
    if (! ((r >= 0) && (r < length))) {
        PyErr_SetString(PyExc_ValueError, "Invalid value of r.");
        return NULL;
    }
    // Make sure that iwt specifies a valid amino acid index
    if (! ((iwt >= 0) && (iwt < n_aa))) {
        PyErr_SetString(PyExc_ValueError, "Invalid value of iwt.");
        return NULL;
    }
    // The way that we access elements of mlist is as follows
    //    double m;
    //    long i;
    //    m = * ((double *) PyArray_GETPTR1(mlist, i)); 
    // assigns m to be the value of element i of mlist.
    //
    // Now loop over all residues to compute the log likelihood
    // For each residue, the calculation basically consists of computing the elements of the double
    //    array p.  It is recomputed for reach residue.  p[i * n_aa + y] is the cumulative probability
    //    that internal node i is residue y given the information in the subtree below i.  Because of
    //    the ordering of the internal node numbering, the entries for i = n_aa - 1 correspond to the
    //    root node of the tree.  Note that if node j is part of the subtree of node i, j < i.  This
    //    is ensured by the node numbering that should have been passed to Initialize.
    // Likewise, dp is a list of length n_internal * n_aa2.  dp[i * n_aa2 + y * n_aa + z] is the
    //    derivative of p[i * n_aa + y] with respect to the ddG for mutating the residue from
    //    its wildtype identity to z, with entries for z equal to the wildtype being undefined.
    //
    // Allocate p and dp
    p = (double *) malloc(n_internal * n_aa * sizeof(double)); // allocate p
    dp = (double *) malloc(n_internal * n_aa2 * sizeof(double)); // allocate dp
    dpiyr = (double *) malloc(n_aa * sizeof(double)); // allocate dpiyr
    dpiyl = (double *) malloc(n_aa * sizeof(double)); // allocate dpiyl
    dpr = (double *) malloc(n_aa * sizeof(double)); // allocate dpr
    for (z = 0; z < n_internal * n_aa2; z++) {
        dp[z] = 0.0;
    }
    r_log_scaling_sum = 0.0; // sum of scaling factor corrections
    r_times_n_tips = r * n_tips;
    for (i = 0; i < n_internal; i++) { // i loops over all internal nodes
        right_index = 2 * i; // index of right descendent in descendents, descendent_is_tip
        left_index = right_index + 1; // index of left descendent in descendents, descendent_is_tip
        i_times_n_aa = i * n_aa;
        i_times_n_aa2 = i * n_aa2;
        if (descendent_is_tip[right_index] && descendent_is_tip[left_index]) { // both descendents are tips
            tr = descendents[right_index];
            tl = descendents[left_index];
            xr = tips[r_times_n_tips + tr];
            xl = tips[r_times_n_tips + tl];
            if (xr == -1) { // a gap on right descendent
                if (xl == -1) { // both left and right descendents are gaps
                    for (y = 0; y < n_aa; y++) {
                        p[i_times_n_aa + y] = 1.0;
                        dp_index = i_times_n_aa2 + y * n_aa;
                        for (z = 0; z < n_aa; z++) {
                            if (z != iwt) {
                                dp[dp_index + z] = 0.0;
                            }
                        }
                    }
                }
                else { // gap is only on right descendent
                    index = unique_uts_index_n_aa2[tl] + xl * n_aa;
                    dm_index = n_aa * unique_uts_index_n_aa2[tl] + xl * n_aa;
                    for (y = 0; y < n_aa; y++) {
                        p[i_times_n_aa + y] = * ((double *) PyArray_GETPTR1(mlist, index));
                        dp_index = i_times_n_aa2 + y * n_aa;
                        for (z = 0; z < n_aa; z++) {
                            if (z != iwt) {
                                dp[dp_index + z] = * ((double *) PyArray_GETPTR1(dmlist, dm_index + n_aa2 * z));
                            }
                        }
                        index++;
                        dm_index++;
                    }
                }
            }
            else if (xl == -1) { // gap only on left descendent {
                index = unique_uts_index_n_aa2[tr] + xr * n_aa;
                dm_index = n_aa * unique_uts_index_n_aa2[tr] + xr * n_aa;
                for (y = 0; y < n_aa; y++) {
                    p[i_times_n_aa + y] = * ((double *) PyArray_GETPTR1(mlist, index));
                    dp_index = i_times_n_aa2 + y * n_aa;
                    for (z = 0; z < n_aa; z++) {
                        if (z != iwt) {
                            dp[dp_index + z] = * ((double *) PyArray_GETPTR1(dmlist, dm_index + n_aa2 * z));
                        }
                    }
                    index++;
                    dm_index++;
                }
            }
            else { // no gap on either descendent
                indexr = unique_uts_index_n_aa2[tr] + xr * n_aa;
                indexl = unique_uts_index_n_aa2[tl] + xl * n_aa;
                dm_indexr = n_aa * unique_uts_index_n_aa2[tr] + xr * n_aa;
                dm_indexl = n_aa * unique_uts_index_n_aa2[tl] + xl * n_aa;
                for (y = 0; y < n_aa; y++) {
                    mr = * ((double *) PyArray_GETPTR1(mlist, indexr));
                    ml = * ((double *) PyArray_GETPTR1(mlist, indexl));
                    p[i_times_n_aa + y] = mr * ml;
                    dp_index = i_times_n_aa2 + y * n_aa;
                    for (z = 0; z < n_aa; z++) {
                        if (z != iwt) {
                            z_n_aa2 = z * n_aa2;
                            dp[dp_index + z] = (* ((double *) PyArray_GETPTR1(dmlist, dm_indexr + z_n_aa2))) * ml + mr * (* ((double *) PyArray_GETPTR1(dmlist, dm_indexl + z_n_aa2)));
                        }
                    }
                    indexr++;
                    indexl++;
                    dm_indexr++;
                    dm_indexl++;
                }
            }
        }
        else if (descendent_is_tip[right_index] && ! descendent_is_tip[left_index]) { // only right descendent is tip
            tr = descendents[right_index];
            xr = tips[r_times_n_tips + tr];
            jl = descendents[left_index];
            jl_times_n_aa = jl * n_aa;
            jl_times_n_aa2 = jl_times_n_aa * n_aa;
            if (xr == -1) { // right descendent is a gap
                indexl = unique_uts_index_n_aa2[n_tips + jl];
                dm_indexl = n_aa * indexl;
                for (y = 0; y < n_aa; y++) {
                    indexly = indexl + y;
                    dm_indexly = dm_indexl + y;
                    piy = 0.0;
                    dp_index = i_times_n_aa2 + y * n_aa;
                    for (x = 0; x < n_aa; x++) {
                        mly = * ((double *) PyArray_GETPTR1(mlist, indexly));
                        pjlx = p[jl_times_n_aa + x];
                        piy += mly * pjlx;
                        dpjlx_index = jl_times_n_aa2 + x * n_aa;
                        for (z = 0; z < n_aa; z++) {
                            if (z != iwt) {
                                dp[dp_index + z] += (* ((double *) PyArray_GETPTR1(dmlist, dm_indexly + z * n_aa2))) * pjlx + mly * dp[dpjlx_index + z];
                            }
                        }
                        indexly += n_aa;
                        dm_indexly += n_aa;
                    }
                    p[i_times_n_aa + y] = piy;
                }
            }
            else { // not a gap
                indexl = unique_uts_index_n_aa2[n_tips + jl];
                indexr = unique_uts_index_n_aa2[tr] + xr * n_aa;
                dm_indexl = n_aa * indexl;
                dm_indexr = unique_uts_index_n_aa2[tr] * n_aa + xr * n_aa;
                for (y = 0; y < n_aa; y++) {
                    indexly = indexl + y;
                    dm_indexly = dm_indexl + y;
                    piy = 0.0;
                    dp_index = i_times_n_aa2 + y * n_aa;
                    for (x = 0; x < n_aa; x++) {
                        mly = * ((double *) PyArray_GETPTR1(mlist, indexly));
                        pjlx = p[jl_times_n_aa + x];
                        piy += mly * pjlx;
                        dpjlx_index = jl_times_n_aa2 + x * n_aa;
                        for (z = 0; z < n_aa; z++) {
                            if (z != iwt) {
                                dp[dp_index + z] += (* ((double *) PyArray_GETPTR1(dmlist, dm_indexly + z * n_aa2))) * pjlx + mly * dp[dpjlx_index + z];
                            }
                        }
                        indexly += n_aa;
                        dm_indexly += n_aa;
                    }
                    mr = * ((double *) PyArray_GETPTR1(mlist, indexr));
                    p[i_times_n_aa + y] = piy * mr;
                    for (z = 0; z < n_aa; z++) {
                        if (z != iwt) {
                            dp[dp_index + z] *= mr;
                            dp[dp_index + z] += (* ((double *) PyArray_GETPTR1(dmlist, dm_indexr + z * n_aa2))) * piy;
                        }
                    }
                    indexr++;
                    dm_indexr++;
                }
            }
        }
        else if (! descendent_is_tip[right_index] && descendent_is_tip[left_index]) { // only left descendent is tip
            jr = descendents[right_index];
            jr_times_n_aa = jr * n_aa;
            jr_times_n_aa2 = jr_times_n_aa * n_aa;
            tl = descendents[left_index];
            xl = tips[r_times_n_tips + tl];
            if (xl == -1) { // left descendent is gap
                indexr = unique_uts_index_n_aa2[n_tips + jr];
                dm_indexr = n_aa * indexr;
                for (y = 0; y < n_aa; y++) {
                    indexry = indexr + y;
                    dm_indexry = dm_indexr + y;
                    piy = 0.0;
                    dp_index = i_times_n_aa2 + y * n_aa;
                    for (x = 0; x < n_aa; x++) {
                        mry = * ((double *) PyArray_GETPTR1(mlist, indexry));
                        pjrx = p[jr_times_n_aa + x];
                        piy += mry * pjrx;
                        dpjrx_index = jr_times_n_aa2 + x * n_aa;
                        for (z = 0; z < n_aa; z++) {
                            if (z != iwt) {
                                dp[dp_index + z] += (* ((double *) PyArray_GETPTR1(dmlist, dm_indexry + z * n_aa2))) * pjrx + mry * dp[dpjrx_index + z];
                            }
                        }
                        indexry += n_aa;
                        dm_indexry += n_aa;
                    }
                    p[i_times_n_aa + y] = piy;
                }
            }
            else { // not a gap
                indexr = unique_uts_index_n_aa2[n_tips + jr];
                indexl = unique_uts_index_n_aa2[tl] + xl * n_aa;
                dm_indexr = n_aa * indexr;
                dm_indexl = unique_uts_index_n_aa2[tl] * n_aa + xl * n_aa;
                for (y = 0; y < n_aa; y++) {
                    indexry = indexr + y;
                    dm_indexry = dm_indexr + y;
                    piy = 0.0;
                    dp_index = i_times_n_aa2 + y * n_aa;
                    for (x = 0; x < n_aa; x++) {
                        mry = * ((double *) PyArray_GETPTR1(mlist, indexry));
                        pjrx = p[jr_times_n_aa + x];
                        piy += mry * pjrx;
                        dpjrx_index = jr_times_n_aa2 + x * n_aa;
                        for (z = 0; z < n_aa; z++) {
                            if (z != iwt) {
                                dp[dp_index + z] += (* ((double *) PyArray_GETPTR1(dmlist, dm_indexry + z * n_aa2))) * pjrx + mry * dp[dpjrx_index + z];
                            }
                        }
                        indexry += n_aa;
                        dm_indexry += n_aa;
                    }
                    ml = * ((double *) PyArray_GETPTR1(mlist, indexl));
                    p[i_times_n_aa + y] = piy * ml;
                    for (z = 0; z < n_aa; z++) {
                        dp[dp_index + z] *= ml;
                        dp[dp_index + z] += (* ((double *) PyArray_GETPTR1(dmlist, dm_indexl + z * n_aa2))) * piy;
                    }
                    indexl++;
                    dm_indexl++;
                }
            }
        }
        else { // both internal nodes
            jr = descendents[right_index];
            jr_times_n_aa = jr * n_aa;
            jr_times_n_aa2 = jr_times_n_aa * n_aa;
            jl = descendents[left_index];
            jl_times_n_aa = jl * n_aa;
            jl_times_n_aa2 = jl_times_n_aa * n_aa;
            indexr = unique_uts_index_n_aa2[n_tips + jr];
            indexl = unique_uts_index_n_aa2[n_tips + jl];
            dm_indexr = n_aa * indexr;
            dm_indexl = n_aa * indexl;
            for (y = 0; y < n_aa; y++) {
                dp_index = i_times_n_aa2 + y * n_aa;
                piyr = 0.0;
                piyl = 0.0;
                indexry = indexr + y;
                indexly = indexl + y;
                dm_indexry = dm_indexr + y;
                dm_indexly = dm_indexl + y;
                for (z = 0; z < n_aa; z++) {
                    dpiyr[z] = 0.0;
                    dpiyl[z] = 0.0;
                }
                for (x = 0; x < n_aa; x++) {
                    pjrx = p[jr_times_n_aa + x];
                    mry = * ((double *) PyArray_GETPTR1(mlist, indexry));
                    piyr += mry * pjrx;
                    pjlx = p[jl_times_n_aa + x];
                    mly = * ((double *) PyArray_GETPTR1(mlist, indexly));
                    piyl += mly * pjlx;
                    dpjrx_index = jr_times_n_aa2 + x * n_aa;
                    dpjlx_index = jl_times_n_aa2 + x * n_aa;
                    for (z = 0; z < n_aa; z++) {
                        if (z != iwt) {
                            z_n_aa2 = z * n_aa2;
                            dpiyr[z] += (* ((double *) PyArray_GETPTR1(dmlist, dm_indexry + z_n_aa2))) * pjrx + mry * dp[dpjrx_index + z];
                            dpiyl[z] += (* ((double *) PyArray_GETPTR1(dmlist, dm_indexly + z_n_aa2))) * pjlx + mly * dp[dpjlx_index + z];
                        }
                    }
                    indexry += n_aa;
                    indexly += n_aa;
                    dm_indexry += n_aa;
                    dm_indexly += n_aa;
                }
                p[i_times_n_aa + y] = piyr * piyl;
                for (z = 0; z < n_aa; z++) {
                    if (z != iwt) {
                        dp[dp_index + z] = dpiyr[z] * piyl + piyr * dpiyl[z];
                    }
                }
            }
        }
        if (!((i + 1) % underflow)) { // scale the probabilities for the node to avoid underflow
            scale = p[i_times_n_aa]; // find the maximum value over all y, starting with y = 0
            for (y = 1; y < n_aa; y++) {
                piy = p[i_times_n_aa + y];
                if (piy > scale) {
                    scale = piy;
                }
            }
            if (scale != 0.0) { // we cannot scale by a scaling factor of 0
                for (y = 0; y < n_aa; y++) {
                    p[i_times_n_aa + y] /= scale;
                    dp_index = i_times_n_aa2 + y * n_aa;
                    for (z = 0; z < n_aa; z++) {
                        if (z != iwt) {
                            dp[dp_index + z] /= scale;
                        }
                    }
                }
                r_log_scaling_sum += log(scale); // add to the log scaling factor correction
            }
        }
    }
    // we have now computed p for the root node.  Compute overall likelihood for this residue
    root = n_internal - 1; // index of the root node
    pr = 0.0;
    for (z = 0; z < n_aa; z++){
        dpr[z] = 0.0;
    }
    r_times_n_aa = r * n_aa;
    root_times_n_aa = root * n_aa;
    root_times_n_aa2 = root * n_aa2;
    for (y = 0; y < n_aa; y++) {
        piry = PyFloat_AsDouble(PyList_GetItem(pi, r_times_n_aa + y));
        pry = p[root_times_n_aa + y];
        pr += pry * piry;
        dp_index = root_times_n_aa2 + n_aa * y;
        dpi_index = n_aa2 * r + y;
        for (z = 0; z < n_aa; z++) {
            if (z != iwt) {
                dpr[z] += PyFloat_AsDouble(PyList_GetItem(dpi, dpi_index + z * n_aa)) * pry + piry * dp[dp_index + z];
            }
        }
    }
    // perform a primitive check for underflow.  Note that this not a very good check -- there could still be a lot of underflow not detected by this check.
    if (pr == 0.0) {
        char error_string[1000];
        sprintf(error_string, "Underflow occurred in the computation of the likelihood for residue %d.", (int) r);
        PyErr_SetString(PyExc_ValueError, error_string);
        return NULL;
    }
    rloglikelihood = log(pr); // logarithm of scaled probability
    drloglikelihood = PyList_New((Py_ssize_t) n_aa);
    for (z = 0; z < n_aa; z++) {
        PyList_SetItem(drloglikelihood, z, PyFloat_FromDouble(dpr[z] / pr));
    }
    rloglikelihood += r_log_scaling_sum; // correct for scaling
    free(p); // free memory allocated for p
    free(dp);
    free(dpiyr);
    free(dpiyl);
    free(dpr);
    // error checking to see if rloglikelihood is nan (not a number)
    if isnan(rloglikelihood) {
        fprintf(stderr, "\nError associated with residue log likelihood being nan for residue %d.\n", (int) r);
        fprintf(stderr, "The value of pr (the unscaled likelihood) is %f.\n", pr);
        fprintf(stderr, "The value of r_log_scaling_sum (the log scaling factor) is %f.\n", r_log_scaling_sum);
        // We do not raise an exception; we assume this is done by the calling function which
        // will also print additional diagnostic information.
    }
    // construct the return value as a Python tuple
    return PyTuple_Pack(2, PyFloat_FromDouble(rloglikelihood), drloglikelihood);
}

// Module documentation string
static char clog_likelihood_docs[] = "C extension for performing core calculations of the Python log likelihood calculations.\n\nThis module defines four functions, 'Initialize', 'IsInitialized', 'ResidueLogLikelihood', and ResidueLogLikelihoodAndDeriv'.\n\tThe 'Initialize' function should be called at the\n\tinitialization of an instance of the object, and stores information about the data and phylogenetic tree\n\tin global variables.  The 'ResidueLogLikelihood' and 'ResidueLogLikelihoodAndDeriv' function can be called repeatedly to compute the\n\tlog likelihood.\n\t'IsInitialized' reports whether 'Initialize' has already been called.\nUsing these C extensions is over 100 times faster than using the pure Python implementation\n\tof the calculations.  The C extensions parallel the pure Python code in the flow of calculation.";

// The module functions
static PyMethodDef clog_likelihood_funcs[] = {
    {"Initialize", (PyCFunction) Initialize, METH_VARARGS, Initialize_docs},
    {"IsInitialized", (PyCFunction) IsInitialized, METH_VARARGS, IsInitialized_docs},
    {"ResidueLogLikelihood", (PyCFunction) ResidueLogLikelihood, METH_VARARGS, ResidueLogLikelihood_docs},
    {"ResidueLogLikelihoodAndDeriv", (PyCFunction) ResidueLogLikelihoodAndDeriv, METH_VARARGS, ResidueLogLikelihoodAndDeriv_docs},
    {NULL}
};

// Initialize the module
void initclog_likelihood(void) {
    Py_InitModule3("clog_likelihood", clog_likelihood_funcs, clog_likelihood_docs);
    import_array(); // this command is necessary to use the numarray C API
}
