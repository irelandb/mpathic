.. _simulate_sort:

==========================================
``simulate_sort``
==========================================

.. contents::

Overview
-------------
``simulate_sort`` is a program within the sortseq_tools package which simulates
performing a Sort Seq experiment.

After you install `sortseq_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: sortseq_tools.sortseq_for_doc
   :func: parser
   :prog: sortseq_tools
   :path: simulate_sort

   

   
Example Input and Output
-----------

The input table to this function must contain a counts (ct) column and a sequences
column. If the sequences are DNA the
column should be labeled 'seq', if they are RNA, the column should be labeled
'seq_rna', and if they are amino acid sequences they should be labeled 'seq_pro'.
By default, sequencing the library will also be simulated, and these results
will be contained in bin 0. The total number of simulated counts will be equal
to the number of library counts.

Example Input Table::

   ct    seq    
   5   AGGTA
   1   AGTTA
   ...

Example Output Table::

   seq    ct    ct_0     ct_1     ct_2 ...
   AGGTA  5     1        2        1
   AGTTA  1     0        1        0
   ...

The output table will contain the original columns, along with the sorted columns (ct_0, ct_1, ct_2, ...)

An example command to execute this analysis::

    sortseq_tools simulate_sort -i my_library.txt -nm LogNormal -o my_sorted.txt


.. include:: weblinks.txt
