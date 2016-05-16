.. _profile_mut:

==========================================
``profile_mut``
==========================================

.. contents::

Overview
-------------
``profile_mut`` is a program within the sortseq package which calculates the
mutation rate at each position in the library.

After you install `sortseq_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: sortseq_tools.sortseq_for_doc
   :func: parser
   :prog: sortseq_tools
   :path: profile_mut

   
If a wild type sequence is passed, the program will calculate the mutation rate
with respect to this sequence. Otherwise it will assume that the most prevalent
base at each position is the wild type.
   
Examples
-----------

The input format for this module must be in the standard form. 
The input table must contain a column labeled ct which is the sum
of counts for each of the sequencing bins (bin_0 through bin_K). The following 
columns should be labeled ct_0 ... ct_K, and should contain the sequencing counts
for each bin 0 through K. Bin 0 must be the unsorted or unselected bin. If you did
not sequence the unsorted bin, you should omit column ct_0 and start with ct_1.
The last column should contain sequence identity. If the sequences are DNA the
column should be labeled 'seq', if they are RNA, the column should be labeled
'seq_rna', and if they are amino acid sequences they should be labeled 'seq_pro'.
Order of columns does matter.

Example Input Table::

   ct    ct_0     ct_1     ct_2    seq
   13    1        5        7      ACG
   28    8        5        5      GGT
   ...

The output from this function will display a column with the position, a
column with mutation rate, and a column with error in calculated mutation rate.

Example Output Table::

   pos    mut    mut_err
   0      .4       .03      
   1      .2       .01    
   ...

Example Command to Run the analysis::
   
    sortseq profile_mut --err -i data.txt -o mutrate.txt

.. include:: weblinks.txt
