## Project: Analysis of microarray data from a LNCap/SETD7 KD experiment
### Study ID: 
### Scientist: Chao Wang
### Data Analysis: Chao Wang, Davit Sargsyan 
### Created: 11/14/2017 

---

## From Chao's Email (11/13/2017):
Values of the 4 groups in the Excel files.    
May not care about the comparisons of these files.    

## Original data
https://cloudbox.dls.rutgers.edu/?ShareToken=3DE3F0B9ADC4D74CD788C3B132556598C5C817B4
From Christian Bixby, RUCDR, 11/17/2017

## Daily Logs
### 01/18/2018
* Hitmap of most differentially expressed genes.

### 01/17/2018
* Added normalized expression plot for all genes; corrected all labels.

### 01/13/2018
* Split the repository: moved all documents related to the no-replica manuscript to a new repo *single.repl.stat*.        
* Added hitmaps for expressions in each group corresponding to the hitmaps of differences.    
* For the manuscript, use the following version of the scripts:    
a. *lncap_setd7kd_data_preprocessing_v2.R*    
b. *lncap_setd7kd_venn_diagram_v6.R*    
c. *lncap_setd7kd_analysis_v6.R*    
d. *chaos_genes_values.R*

### 01/06/2018
* Added SD estimators and log2 changes (analysis v7)   
* Added db/db data from NCBI for the method paper

### 01/06/2018
* Saved selected genes' table and values (analysis v6)

### 12/19/2107
* Added script for subsetting data  keeping only differentially methylated genes (list of genes by Chao, based on hitmaps)

### 12/02/2017
* Added diff vs. mean plot and SD estimation, test stats and p-values for single-replica experiments (**lncap_setd7kd_analysis_v6.R**)

### 11/21/2017
* Modified hitmaps and Venn diagrams

### 11/20/2017
* Finished data preprocessing    
* Added hitmaps and Venn diagrams

### 11/19/2017
* Downloaded original files from RUCDR    
* Started preprocessing

### 11/14/2017
* Downlowded data from Chao's emails and converted to CSV (***C:\git_local\data*** folder)    
* Processed and merged LNCaP and Setd7 data sets    
* Heatmaps of values and differences

### 11/15/2017
* Made hitmaps and Venn diagram    
* ToDo: hitmap of (LNCaP - Setd7) vs. (Setd7 PEITC - Setd7)
