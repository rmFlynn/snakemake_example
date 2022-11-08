.. RMNP_Pipeline documentation pipline file file, created by
.. Rory on what ever day this is.

.. get_it_working-configure_the_pipline:

=========================================
Lets run This Pipeline!
=========================================


If you are here then hopefuly:

First clone/download the repository this git repository:
-------------------------

```
git clone https://github.com/WrightonLabCSU/RMNP_pipline.git
```

The fact that the full pipeline is under git version control is key, although you don’t need to commit your changes they are automatically being tracked locally and can be used to find problems. It will for example track the changes to the config file in the next section

Change directories to the new folder you just made with the git command ``cd RMNP_pipline`` this is where the rest of this tutorial will assume your working directory is.

If you can't use git in your instance, or it is more troble than it is worth you can download the data from here

Edit the Config File
-------------------------

You need to edit the config file to fit your use case. There is a crazy number of options in this file but most are self explanatory, and you will not need to change all of them

 - The config lives in the config directory, and is a yaml file wich is a fancy sort of text file so open it up with your favorite text editor be that the one on you desktop or nano or vi/vim/neovim. Then you can tweak these settings.

   ``vi config/config.yaml``

 - First set the type of analyse you want to do there are two options

   ``type: 'meta_t' # type can be 'meta_t' or 'meta_g'``

 - Set the paths to FASTAs. There are many options depending on how your data is organized, FASTAs are so fun like that. You can use interleaved or non interleaved reads, impute ore ascribe names to the FASTAs, and you can add them as zipped or unzipped files.
   ::
      inputs:
        interleaved_reads:
           #>- read0.fastq
           #>- read1.fastq
           # alternately
           #>- folder/path/*.fastq
        paired_reads:
          #- '/path/one/name0*.fastq*'
          #- '/path/two/name1*.fastq*'
        named_reads:
        # Here you can put nothing, or asmany reads as you want with the format:
        any_name_you_want:
          #>R0: something_R1.fastq.gz
          #>R1: something_R2.fastq.gz
          # Alternately
          #>inter: something_interleaved.fastq.gz
        # This is one database file should it be more?
        # You set both the input and output for the database, that bowtie1 will use. These files may be used for many``

   One side note is that if you use the paired_reads as your input, you need to also check that the forward and reversed reads are identified constantly with the strings in the binning section. We discuss this in the next section.

 - You also need to set two locations for bowtie1 if you are making a database for this run. You can of coarse not set a location for the raw if you are using a pre-built db. However, if you do make the database you should definitely save it in a good place so it can be used latter.

   ``bowtie1_database_raw: an/fa/file/for/bowtie2.fa
     bowtie1_database_built: resources/bowtie2_databases/full_db # The permanent home for the database``

 - Then set a few more location.

   ``gff: the/genes/file.gff # only for meta-t
     fasta_dir: the/dir/of/fastas # only for meta-g``

