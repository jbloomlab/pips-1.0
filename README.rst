=====================================================
PIPS (Phylogenetic Inference of Protein Stability)
=====================================================

.. contents:: Contents
   :depth: 2

Overview
-----------
`pips`_ is a `Python`_ software package written by `Jesse Bloom`_. `pips`_ source code is available `on GitHub`_. 

This program is designed to analyze phylogenies to identify potentially stabilizing mutations.

References
-----------
The references for `pips`_ are:

* Jesse D. Bloom, Jagannath S. Nayak, and David Baltimore. "A computational-experimental approach identifies mutations that enhance surface expression of an oseltamivir-resistant influenza neuraminidase." PLoS One. 6:e22201 (2011)  http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0022201

* Jesse D. Bloom and Matthew J. Glassman. “Inferring stabilizing mutations from protein phylogenies: application to influenza hemagglutinin.” PLoS Comput. Biol. 5:e1000349 (2009) http://www.ncbi.nlm.nih.gov/pubmed/19381264

Requirements
-------------
`pips`_ is a `Python`_ package. It requires the following:

    * `Python`_: `pips`_ has been tested with `Python`_ version 2.6* and 2.7*. It is not known if `pips`_ will work with other versions of `Python`_.

    * `numpy`_

    * `gcc`_ or some similar C compiler is required to build the C extensions.

    * Either `scipy`_ or the `Python`_ ``transcendental`` package are required. Originally, `pips`_ used ``transcendental``, but this no longer appears to be available on the web, so you might need to use `scipy`_.

Installation
-------------
Download `pips`_ from the repository `on GitHub`_. Build and install `pips`_ with the following commands::

    python setup.py build
    python setup.py install

The last command might need to be replaced with::

    sudo python setup.py install

if you want to install globally and do not have superuser privileges by default, or by::

    python setup.py install --user

if you want to install locally.

Using `pips`_
------------------
Installing `pips`_ will install all the modules, which you can then import in `Python`_.

It also installs the following scripts, which are found in the ``./scripts/`` subdirectory of the main `pips`_ package. These scripts are::

    pips_analysis.py                     
    pips_correlate_selected_mutations.py
    pips_analyze_selected_mutations.py   
    pips_run_cupsat.py
    pips_build_tree_and_alignment.py     
    pips_run_foldx.py
    pips_consensus.py

Unfortunately, these scripts are not currently documented externally. However, if you look in the documentation string at the beginning of each script, it will explain how to create an input file appropriate for running the script.

The `pips`_ package also contains an ``./examples/`` subdirectory that shows example analyses -- these may also be helpful.

.. _`pips`: https://github.com/jbloom/pips-1.0
.. _`on GitHub`: https://github.com/jbloom/pips-1.0
.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Python`: https://www.python.org/
.. _`numpy`: http://www.numpy.org/
.. _`gcc`: http://gcc.gnu.org/
.. _`scipy`: http://www.scipy.org/
