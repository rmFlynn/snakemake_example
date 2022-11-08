.. RMNP_Pipeline documentation pipeline file, created by
.. Rory on whatever day this is.

.. _get_it_working-before_you_start:

=====================
Before You Start
=====================

Hopefully, you have looked at the home page and have an idea about why this is a Snakemake pipeline, but you may not know how to use Snakemake in your environment. All environments are good, but this pipeline was made to work with Slurm. To work with Surm you need to set up a Slurm profile, this is fast, and you hopefully only need to do it once.

Get to the Point
________________

For Slurm scripts to work for you in Snakemake, a config file must live in your home directory AKA '~' to tell Slurm what to do. The location of this file will be `~/.config/snakemake/slurm/config.yaml`. Making a file in this location will create a profile called Slurm. You can make other profiles if you want to have different Slurm settings, but we will skip that nuance for now.

**This is not RMNP specific,** you will use the same profile for any Snakemake pipeline. That said, if you have already set up a profile, you may want to compare it to this one as the changes could have unforeseen consequences. For example, if you have your profile set to email you, then it may overload your email.


Set a profile
_____________

If you are pulling this repository from GitHub, you can copy this profile from the repo its self. If you want, you can download this file :download:`HERE <../../example_slurm_profile.yaml>` or find it here in the git repo `HERE <https://github.com/rmFlynn/snakemake_example/example_slurm_profile.yaml>`_. The content of the file itself can also be copied below.


**Before you continue, you must edit the partition of this file**, it should match the partition you would normally use in a standalone Slurm script.

.. literalinclude:: ../../example_slurm_profile.yaml
     :language: YAML
     :emphasize-lines: 15
     :linenos:


The file above should document its self, so we will not harp on here. If you want to work with Snakemake more often, you will find it useful to understand it fully.

