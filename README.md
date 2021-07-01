## prediction and pattern discovery for drug combination based on multimodel friendship features

Our paper is at https:

![ABC](.//Flowchart.jpg) 

## Installation
First install the dependencies via conda:
 * sklearn
 * numpy
 * Python >= 3.7

 
### Quick start
We provide an example script to run experiments on our dataset: 

- Run `python FSM.py`: predict drug-drug combinations, and evaluate the results with cross-validation. 

**Note:** For a detailed description of the data, please refer to the section below.
 
 
## Data Description  
This drug information is in the `data/` folder.  

**drugs.txt**: Drugbank id for all drugs.   
**cancer_combination.txt**: The gold standard for cancer combinations.    
**hypertension_combination.txt**: The gold standard for hypertension combinations.   
**pathway_tanimoto.txt**: Similarity of drug pairs calculated from drug pathway feature.   
**atc_tanimotol.txt**: Similarity of drug pairs calculated from drug atc feature.   
**chemical_tanimoto.txt**: Similarity of drug pairs calculated from drug chemical structure feature.   
**offside_tanimoto.txt**: Similarity of drug pairs calculated from drug adverse feature.   
**target_tanimoto.txt**: Similarity of drug pairs calculated from drug target feature.   
   


## Citation


### Contacts
If you have any questions,contact the maintainer of this package: HG_Chen 13247702278@163.com


