from distutils.core import setup, Extension
import sys

sys.path.append('/Library/Python/2.6/site-packages') # hack needed to find right package on my Mac

# test for required modules
try:
    import numpy
except ImportError:
    raise ImportError, "This package requires numpy.  Currently, numpy cannot be imported."
try:
    import transcendental
except ImportError:
    raise ImportError, "This package requires transcendental.  Currently, transcendental cannot be imported."

# list of C extensions
clog_likelihood = Extension('pips.clog_likelihood', sources = ['src/clog_likelihood.c'], include_dirs = [numpy.get_include()])
if sys.platform == 'darwin': # we will perform matrix multiplications with Accelerate cblas
    link_args = ['-framework Accelerate'] 
    macros = [('USE_ACCELERATE_CBLAS', None)] 
else:  # we will perform matrix multplications in pure C code.
    link_args = []
    macros = []
cddg_inference = Extension('pips.cddg_inference', sources = ['src/cddg_inference.c'], include_dirs = [numpy.get_include()], extra_compile_args = ['-O3'], extra_link_args = link_args, define_macros = macros)

# main setup command
setup(
    name = 'pips', 
    fullname = 'Phylogenetic Inference of Protein Stability', 
    version = '1.0', 
    author = 'Jesse D. Bloom', 
    author_email = 'jesse.bloom@gmail.com', 
    url = 'None yet created.', 
    description = 'Infers protein ddG values from sequences of known phylogeny using Bayesian inference.',  
    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'License :: Free for non-commercial use',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    platforms = 'Tested on Mac OS X.',
    packages = ['pips'],
    package_dir = {'pips' : 'src'},
    scripts = [
        'scripts/pips_analysis.py', 
        'scripts/pips_build_tree_and_alignment.py',
        'scripts/pips_run_cupsat.py',
        'scripts/pips_run_foldx.py',
        'scripts/pips_consensus.py',
        'scripts/pips_analyze_selected_mutations.py',
        'scripts/pips_correlate_selected_mutations.py',
        ],
    ext_modules = [clog_likelihood, cddg_inference]
)
