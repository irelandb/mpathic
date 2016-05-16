.. _simulate_library:

==========================================
``simulate_library``
==========================================

.. contents::

Overview
-------------
``simulate_library`` is a program within the sortseq_tools package which creates a library of
random mutants from an initial wild type sequence and mutation rate.

After you install `sortseq_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: sortseq_tools.sortseq_for_doc
   :func: parser
   :prog: sortseq_tools
   :path: simulate_library

   

   
Examples Inputs and outputs
-----------

To generate a library of mutated sequences you could use the command::
    
    sortseq_tools simulate_library -w ACAGGGTTAC -n 50000 -m .2

This will output a simulated library of the form:: 

    ct    seq
    100   ACAGGGTTAC
    85    ACGGGGTTAC
    ...


     




.. include:: weblinks.txt
