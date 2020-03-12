# Analysing Peptide Expression in Proteomic Cancer data. 

A collection of python modules used to analyse and filter proteomic cancer data.  
They filter PEAKS output files until the peptides of interest remain.
The peptide lists get compared to one another and can be visualised in venn diagrams.
The resulting files are utilized in the other [repository](https://bitbucket.org/kcduong/proteogenomics_stats/)

## Getting Started

The lst_to_fasta_converter folder contains small scripts that were used to convert coding region identifier tool 
outputs to databases suitable for PEAKS input. These can be ignored.

Peaks_peptide_comparison.py compares PEAKS output csv files for similarity and differences in peptides.  
These comparisons can be visualised with peptide_venn.py as Venn diagrams.  
Unknown_peptide_seeker.py filters out peptides that are present in reference protein databases.   
Human_only_db.py filters out every non-human record from a reference protein database.

## Run
The modules can be run separately but there's also a pipeline, which combines a bunch of the aforementioned modules 
and generates the needed results in one go. 

### Built with
Python 3.7

And the following external libraries:
```
| Library        	| Version 	|
|-----------------	|---------	|
| BioPython       	| 1.75    	|
| matplotlib      	| 3.1.2   	|
| matplotlib-venn 	| 0.11.5  	|
| NumPy           	| 1.17.4  	|
| pandas          	| 0.25.3  	|
| SciPy           	| 1.3.3   	|
| statsmodels     	| 0.11.0  	|
```