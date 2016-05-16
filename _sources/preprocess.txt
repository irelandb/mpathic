.. _preprocess:

==========================================
``preprocess``
==========================================

.. contents::

Overview
-------------
``preprocess`` is a program within the sortseq_tools package which combines
individual tabular data files,FASTA, or FASTQ files and combines them into a
data format that sortseq_tools can use. All other preprocessing, such as stitching
paired end reads together, or filtering based on barcode or quality score, must
be handled before running this module.

After you install `sortseq_tools`_, this program will be available to run at the command line. 

Command-line usage
---------------------
.. argparse::
   :module: sortseq_tools.sortseq_for_doc
   :func: parser
   :prog: sortseq_tools
   :path: preprocess

   

   
Examples
-----------

The core input is a tabular data file which specifies which data files contain
the sequences from which bin. The file must contain a 'bin' column and a 'file'
column.

Example Input Table::

   bin  file
   0	my_bin_0.fasta
   1	my_bin_1.fasta

The output of this module will be a data set in a tabular format that can
be fed to the rest of the modules in the sortseq program.

Example Output Table::

   ct    ct_0     ct_1     ct_2    seq
   13    1        5        7      ACG
   28    8        5        5      GGT
   ...


Example command to run the analysis::

   sortseq preprocess -i my_fileslist.txt -o mydataset.txt --type DNA



.. include:: weblinks.txt
