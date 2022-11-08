.. RMNP_Pipeline documentation pipline file file, created by
.. Rory on what ever day this is.

.. get_it_working-before_you_start:

=====================
Before You Start
=====================

Hopefully you have looked at the home page and have an idea about why this is a snakemake pipeline but you may not know how to use Snakemake in your environment. All environments are good but this pipeline was made to work with Slurm. To work with slurm you need to set up a slurm profile this is fast and you hopefully only need to do it once.

Get to the Point
________________

For slurm scripts to work for you in snakemake, a config file must live in your home directory AKA '~' to tell slurm what to do. The location of this file will be `~/.config/snakemake/slurm/config.yaml` if you are on zenith run
``cp /home/projects-wrighton-2/projects-flynn/Snakemake_Slurm_RMNP_Example/slurm_and_snake_make/example_config.yaml ~/.config/snakemake/slurm/config.yaml`` where a template is stored. You can also download a profile here. This will make a profile called slurm. You can make other profiles if you want to have different slurm setings but we are going a head to far.


