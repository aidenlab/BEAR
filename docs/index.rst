.. Polar documentation master file, created by
   sphinx-quickstart on Tue Apr  7 10:31:01 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Polar: pipeline for viral diagnostic
====================================

Polar can be run on a SLURM cluster or on a single machine.
In either case, you must install the following required software:

* `BWA <https://github.com/lh3/bwa>`__
* `Samtools <http://www.htslib.org/download/>`__
* `Minimap2 <https://github.com/lh3/minimap2>`__
* `Megahit <https://github.com/voutcn/megahit>`__

Single Machine Quick Start
==========================

#. Install required software
#. Clone repository

   .. code-block:: bash

      git clone https://github.com/aidenlab/Polar.git

#. Run test

   .. code-block:: bash

      cd Polar/test
      ../align_serial.sh


SLURM Quick Start
=================

#. Ensure you have installed required software.
#. Clone repository

   .. code-block:: bash

      git clone https://github.com/aidenlab/Polar.git

#. Modify the variables at the top of align_slurm.sh to
   correspond to your system's load, commands, and queues.
#. Run test

   .. code-block:: bash

      cd Polar/test
      ../align_serial.sh

Detailed Guide
==============

Polar is the viral diagnostic pipeline, currently designed
for SARS-Cov-2. For more information on the protocol, see
our paper here:

The pipeline takes as input paired-end sequencing reads
and creates as output a PDF with the result of the test for the virus
(positive or negative). The PDF also includes other qualitative and 
quantitative measures, detailed below.

Repository
^^^^^^^^^^

The Polar code can be found at https://github.com/aidenlab/Polar.git
The repository contains the viral genomes


.. toctree::
   :maxdepth: 2
   :caption: Contents:


* :ref:`search`