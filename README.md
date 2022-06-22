# Snakemake workflow: Read Mapping, No Problem


A Snakemake workflow for the Wrighton Lab read processing pipeline.


## Usage

Using a Snakemake pipeline is different from a simple script, but it is not difficult and in many ways it is easier. The key difference is going to be that you will fill out a config file instead of providing commands. Also, instead of just using a script, you will use an entire structure of folders that will organize and log the running of your pipeline. 
### Lets Go Step by Step!
#### First clone this git repository:
```
git clone https://github.com/WrightonLabCSU/RMNP_pipline.git
```
The fact that the full pipeline is under git version control is key, although you don’t need to commit your changes they are automatically being tracked locally and can be used to find problems. It will for example track the changes to the config file in the next section

Change directories to the new folder you just made with the git command `cd RMNP_pipline` this is where the rest of this tutorial will assume your working directory is. 

#### Edit the Config File

Change directories to 
# TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning the user or organization.
* Replace `<name>` with the workflow name (can be the same as `<repo>`).
* Replace `<description>` with a description of what the workflow does.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.


