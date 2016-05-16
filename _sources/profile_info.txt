.. _profile_info:

==========================================
``profile_info``
==========================================

.. contents::

Overview
-------------
``profile_info`` is a program within the sortseq_tools package which calculates
the mutual information between base identity at a given position and expression
for each position in the given data set. This analysis is best used to identify regions which are important to transcription, such as transcription
factor binding sites.

After you install `sortseq_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: sortseq_tools.sortseq_for_doc
   :func: parser
   :prog: sortseq_tools
   :path: profile_info

   

   
Example Input and Output
-----------

The input must be in the standard form for the sortseq package. It must have
 a column for sequences and 
columns of counts for each bin. For selection experiments, ct_0 should label the
pre-selection library and ct_1 should be the post selection library. For MPRA
experiments, ct_0 should label the sequence library counts, and ct_1 should
label the mRNA counts.

Example Input Table::

   ct    ct_0     ct_1     ct_2    seq
   13    1        5        7      ACG
   28    8        5        5      GGT
   ...

Example Output Table::

    pos    info    info_err
    0      .02     .004
    1      .04     .004
    ...

The mutual information is given in bits.

An example command to run this analysis is::

    sortseq profile_info -i sorted_library.txt --err --start 20 --end 80 -o info_profile.txt --type DNA
    



.. include:: weblinks.txt
