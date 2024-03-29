# README snRNA & HiTi projects
This is the README file for the project "Application of molecular barcoding and
engineered viral vectors for the determination of causative transcriptomic
changes in Parkinson’s disease and generation of a fluorescent reporter model
for studying cell-cell communication in the brain"

The link to the code files are found at
https://github.com/Atte-Oskari-Rasanen/snRNA_and_HiTi_project

################################################################################

The project consists of two parts: snRNA-seq and HiTi project.

The analysis pipeline for snRNA-seq can be found inside the folder snRNA project
containing separate folders based on the technology used (Scifi and SPLiT) as
well as helper scripts inside the helper_scripts folder.

Creation of separate virtual conda environments for the projects is recommended.
To this end, the required packages are found in the environment yml-files
found in both project folder.

Anaconda: https://docs.anaconda.com/anaconda/install/index.html


Furthermore, the following packages are needed:
Starcode (1.4)- https://github.com/gui11aume/starcode/releases/tag/1.4
Mview (1.6.7) - https://sourceforge.net/projects/bio-mview/files/bio-mview/
Cutadapt (4.2) - https://cutadapt.readthedocs.io/en/stable/installation.html
bbmap (latest) - https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/


The SPLiT folder contains the following files:
**Dockerfile** - contains the Parse biosciences analysis pipeline and the
required dependencies
installed in it. Can be built using the command:
> docker build . -t parsepipe

**parse_analyse_rat.sh** - contains the script which executes the pipeline on
the sample files

**analysis_env.yml** - required packages for the analysis. To install, run
> conda env create -f analysis_env.yml

3 Parse_analyse R files:
**1** - Preprocessing, analysis, subsetting and addition of single cell-level
annotations
**2**  - Grouping the DA cells into bins of a-syn and tagBFP. Performing DGE
between various conditions, visualising singificantly upregulated/downregulated
genes
**3**  - DGE between conditions, Pathway analysis

The helper_scripts folder contains the following scripts:
**1_linnarson_ref_processing.py** - Performs subsetting the reference to the
relevant cell types and samples each cell type so that all contain a maximum of
1000 cell columns from the cell x gene count matrix. Generation of a reference
file from Linnarson's dataset (file L5.Loom) found at their
site: http://loom.linnarssonlab.org/ .
**2_da_subtype_annotation.py** - Annotation script for annotating DA
cells into different subtypes.
**3_create_final_annotation_file.py** - Takes in the manually annotated
glutamatergic subset data and the annotations of the Seurat object containing the annotations of the SingleR and the more specific DA subtype annotations generated by
2_da_subtype_annotation.py script. Based on these, the final annotation file
is created by merging the data. This data is reintegrated back into the Seurat
object in the 1_Parse_analyse.R

Inside Scifi folder the following files are found:
**scifi5_preanalysis.py** - Preanalysis of Scifi5 data. Generation of the count matrix
and merging the rat (oDT) data with the barcoded AAVs (WP)
**scifi6_preanalysis.py** -  Preanalysis of Scifi6 data. Generation of the count matrix and merging the rat (oDT) data with the barcoded AAVs (WP) and categorising
based on the given groups
**scifi5_main_analysis.R** - Main analysis run using Seurat

Inside Arc-HiTi folder one can find samples characterised by folder name starting
with HITI. Each HITI sample subfolder contains main.py python script for the analysis workflow of the given sample. The script calls on functions from scripts_main.py found in the root.The folder amplicon_files contains DNA sequences of the reference files with the whole fluorescent protein and the Arc gene along with information about the primer binding sites and primer templates. In the root folder there is also the following files:

**scripts_main.py** - Contains all the functions and classes used for the
analysis. The functions are called in the analysis workflows performed
individually in the main.py script found within each sample subfolder.
**conda_env.yml** - The required python packages and their versions. to install,
run the following command:
> conda env create -f conda_env.yml
