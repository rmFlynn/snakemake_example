.. RMNP_Pipeline documentation pipeline file, created by
.. Rory on whatever day this is.

.. _get_it_working-configure_the_pipline:

=========================================
Configure the Pipeline
=========================================

Before you can run the pipeline, you must configure it. One of the awesome things about Snakemake is that it does not depend on an ephemeral command with a hidden default to run. All the setting for a Snakemake pipeline are explicitly stated in a config file, `assuming you are following best practices <https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html>`_.

The pipeline comes as a package deal, the scripts that are run, the tools that are used, and the results that are created are all frozen in time once you run your pipeline. This makes reproducibility a snap provided you don't re-organized the data after the fact. In snake make there is a place for everything and everything has its place!

Clone/Download the Repository:
-------------------------

You will need to get a copy of the empty pipeline, you can don’t all the data from this repository by clicking :download:`HERE <https://github.com/rmFlynn/snakemake_example/archive/refs/heads/main.zip>`. Or, better yet, with the command

.. code-block:: bash

 git clone https://github.com/rmFlynn/snakemake_example.git

Possibly you will want to name your pipeline folder, if so use:

.. code-block:: bash

  git clone https://github.com/rmFlynn/snakemake_example.git name_of_pipeline

The fact that the full pipeline is under git version control is a bonus, although you don’t need to commit your changes they are automatically being tracked locally and can be used to find problems. It will for example track the changes to the config file.

Change directories to the new folder you just made with the git command ``cd RMNP_pipline`` this is where the rest of this tutorial will assume your working directory is.

If you can't use git in your instance, or it is more trouble than it is worth, you can download the pipeline without any problems.

Edit the Config File
-------------------------

You need to edit the config file to fit your use case. There is a crazy number of options in this file, but most are self-explanatory, and you will not need to change all of them

 - The config lives in the config directory, and is a yaml file which is a fancy sort of text file, so open it up with your favorite text editor, be that the one on your desktop or nano or vi/vim/neovim. Then you can tweak these settings.

  .. code-block:: bash

    vi config/config.yaml

  When you are editing the config file, it is advised that you leave anything that you are not using as it is.

  The config file contains comments and  example code, to differentiate between the two the comments start with a traditional ``#`` and the code with ``#>``  to use the code delete the ``#>`` and nothing else.

 - First set the type of analyze you want to do there are two options meta_t and meta_g if you don't know the difference then you will want to review :ref:`explain_the_pipeline`

   here is an example on lines 1-3

   .. literalinclude:: ../../config/config.yaml
        :language: YAML
        :lines: 1-3


   The choice here should be clear to you.

 - Set the paths to FASTAs. There are many options depending on how your data is organized, FASTAs are so fun like that. You can use interleaved or non interleaved reads, impute or ascribe names, and fastas can be zipped or unzipped files.

   By default there is some example code on lines 4-20.

   .. literalinclude:: ../../config/config.yaml
        :language: YAML
        :lines: 4-20

   You will Notice that there are 3 options for how you enter your reads. Interleaved, paired and named. All of these are simple to use, but you may need to be aware of additional options for each.

    - Interleaved reads are simple to enter, the format is that of a list, in YAML lists start with ``-`. Each item in this list is a wildcard. So for example, you could type ``- path/*.fasta`` to enter an entire folder of interleaved FASTQ files.

      For one such file you would type:

      .. code-block:: YAML

          interleaved_reads:
             - path/*.fasta

      For 2 such files you would type

      .. code-block:: YAML

          interleaved_reads:
             - path1/*.fasta
             - path2/*.fasta

     And you can do this for as many wild cards as you want. you can also specify each file individually

      .. code-block:: YAML

          interleaved_reads:
             - path/file1.fasta
             - path/file2.fasta

      Note that the named of your samples will be imputed from the fasta names, so for above your samples would be named file1 and file2 in all results.

    - Paired reads are similar but more nuanced to enter, the format is that of a list again, but each item must be a wild card. The names of the sample are imputed like above. So if you have only one sample and that sample made forward reads  ``name_R1.fasta`` and reverse reads ``name_R2.fasta`` then the results would reference the sample as ``name``, and you would enter that record as``- path/name_*.fasta``, or similar.

      Part of the name imputing process for paired reads is identifying the forward and reverse fastaq for each sample. For example, if you have a folder with many samples of this type with each pair having the same name with the additions of ``_R1`` and ``_R2`` for forward and reversed reads, then they would be entered like so.

      .. code-block:: YAML

          paired_reads:
             - path/*.fasta

      provided the labeling of forward and reverse reads are as described, the program will find all the samples.

      However, if the labels are not ``_R1``, and ``_R2`` for example they could be ``.a``, and ``.b`` then you must edit some lower level seating. On lines 34-37 you will find the binning section, it looks like this by default.

       .. literalinclude:: ../../config/config.yaml
            :language: YAML
            :lines: 34-37

       So in this situation, you would change this portion to:

      .. code-block:: YAML

         binning:
           forward_id: '.a'
           backward_id: '.b'
           sickle_quality_type: 'sanger' # options include 'sanger' and 'illumina'


    - The named reads are where to go if none of the above works for you. The sample names are not part of the file names, things are inconsistent or otherwise complex. In this section, you must list out your files exhaustively. For example, let's say that we have 2 samples one paired one interleave and all the sample names disconnected from the files. We would change the placeholder from:

       .. literalinclude:: ../../config/config.yaml
            :language: YAML
            :lines: 16-20

      to:

      .. code-block:: YAML

          my_sample_one:
             R1: idk.fastq.gz
             R2: idk_rev.fastq.gz
          my_sample_two:
             inter: idk_again.fastq.gz

     The results will use these assigned names instead of the file name.

   At no time should you delete parts of the config you are not using, this includes this section.

 - You also need to set two locations for bowtie1 if you are making a database for this run. You can of course not set a location for the raw if you are using a pre-built db. However, if you do make the database, you should definitely save it in a good place so it can be used later.

   ``bowtie1_database_raw: an/fa/file/for/bowtie2.fa
     bowtie1_database_built: resources/bowtie2_databases/full_db # The permanent home for the database``

 - Then set a few more locations.

   ``gff: the/genes/file.gff # only for meta-t
     fasta_dir: the/dir/of/fastas # only for meta-g``
 - Finally, you can pick your outputs. These are mostly self-explanatory but I will still expand on them later.
   Simply change the Boolean values next to the output that you want to ``True``.

       .. literalinclude:: ../../config/config.yaml
            :language: YAML
            :lines: 28-33

  get_coverm: False

