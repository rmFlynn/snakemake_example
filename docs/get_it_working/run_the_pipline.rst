.. RMNP_Pipeline documentation pipeline file, created by
.. Rory on whatever day this is.

.. _get_it_working-run_the_pipline:

=========================================
Let's run This Pipeline!
=========================================


This part is straightforward, but still should be followed closely. Recall that a comment is just a ``#`` and example code is ``#>``.

Because we use Slurm as part of the pipeline, it only makes sense to run this pipeline in a slurm session itself.
The session can be very small, it just needs to monitor the pipeline, you can of course run this in a screen if that is your preference.
The code below is a jumping off point to get you started using Snakemake

.. literalinclude:: ../../slurm.sh
     :language: YAML
     :linenos:


What Makes a Snakemake Command
______

Snakemake is run with the command ``snakemake``, using this command from the root of the pipeline. It will then find the Snakefile in the workflow directory.
When using Snakemake with Slurm, a command needs 3 things to be valid.

 -  A slurm profile, this describes default arguments for Slurm, and the user specific settings for Slurm to run. You should have this by following the steps in :ref:`get_it_working-before_you_start`. use ``--profile`` to call it, if it is named Slurm like in this tutorial then it will be ``--profile slurm``.
 - The number of jobs that Slurm can run simultaneously. There may be a limit to how many jobs a person can have in the queue at a time. Most people have paid too much for their HPC for it to be half empty, so you may want to set this high, even if that means jobs in the queue.
   Use ``-j`` for this, for example ``-j 10``
 - The number of cores the jobs can demand. Note that the Snakemake profile will limit the number of cores in a single threaded job. But some jobs will try to run as many threads as possible.
   This limits the number of cores a max size job can have. The person running your HPC most likely gave you an upper limit, but you can run more jobs at a time if this number is lower. An example of this argument would be ``-c 10``

Step by Step
____________

 - first customize the Slurm settings of these files, at minimum you must put your email in.
 - the first step is to activate an environment that has Snakemake installed, this can be done with the source command.
   On line 15 there is an example of such a string.
 - Next do a dry run like on line 22 this will show you how many jobs will run and if there is a problem with the pipeline it will give you the error without wasting time.
   There is no need at all
 - If you are doing a test run then you may now want to use the ``--keep-incomplete`` and ``--notemp`` flags like on line 24. These will prevent the erasure of temporary files on success or failure. This will let you pickup where you left off if the pipeline failed, or you need to add something to the output.
 - If this is not a test, you can use a command like that on line 20 to start the pipeline.

