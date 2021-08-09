# 1. Introduction

BUFET2 is the new version of BUFET and it is an open-source software under the GPL v.3 licence,.

The BUFET algorithm generates the empirical distribution of the genes in a class, that are targeted by miRNAs using two statistical measures (left- and two-sided overlap). Then it calculates two p-values for the respective gene class.

# 2. Using BUFET2
Users have the option to either download and compile the source code for BUFET2 or use the Docker image [https://hub.docker.com/r/diwis/bufet2](diwis/bufet2).


## 2.1. Download and compile source code

### 2.1.1. Prerequisites
1. A system with at least 4GB of RAM
2. Linux and MacOS
3. Python interpreter (>= version 2.7) that can run from the command line.
4. g++ 4.8 and above.

### 2.1.2. Compiling BUFET

1. Download the source code from GitHub:
``` bash
git clone https://github.com/athenarc/BUFET2
```
2. Navigate inside the BUFET2 folder:
```bash
cd bufet2
```

3. Compile the code:
```bash
make
``` 

This will compile the code and create a .bin file. The .bin file must be in the same folder as the .py file at all times, in order for the program to execute correctly.

## 2.2. Executing BUFET2

### 2.2.1 Prerequisites

1. Input miRNA file, which is a text file containing only the names of differentially expressed miRNAs, each on a separate line. For example:
<pre>
hsa-miR-132-5p
hsa-miR-132-3p
</pre>

2. Gene synonym data file from NCBI, http://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz. Decompress the file with:
```bash
gzip -d All_Mammalia.gene_info.gz
```
    
3. Gene annotation data file retrieved by GO, KEGG, PANTHER, DisGeNET, etc. The file must contain the following CSV format for each line:
<pre>
gene_name|pathway_id|pathway_name
</pre>

4. miRNA-gene interactions file, which has the following CSV format for each line:
<pre>
miRNA_name|gene_name
</pre>
Note that all files listed above can contain header lines starting with the "#" character.

### 2.2.2. Script Execution
#### 2.2.2.1. Using the compiled source code

Navigate inside the folder containing the .py and .bin files and run the following command:
```bash
python3 bufet2.py [OPTIONS]
```

#### 2.2.2.2. Using the Docker image
1. Create a folder containing the files:
```bash
mkdir data
```
2. Run the script
```bash
docker run --rm -v data:/data diwis/bufet2 python3 /home/bufet2/bufet2.py [OPTIONS]
```

#### 2.2.2.3. Options

By default, the python script verifies that all input files exist, that they are not empty and that they have the correct format. Since the file check leads to increased execution times, it can be disabled by using the "--disable-file-check". However, we recommend that the file check remains enabled, since non-existing or empty files can crash the C++ core.

The script options are listed below:
<pre>
"-miRNA [filename]": path to the input miRNA file.
"-interactions [filename]": path to miRNA-gene interactions file.
"-synonyms [filename]": path to the gene synonym data file.
"-annotations [filename]": path to gene annotation data file.
"-output [filename]": path to output file. Created if it doesn't exist. Default: "output.txt"
"-iterations [value]": number of random miRNA groups to test against. Default value: 1000000
"-processors [value]": the number of cores to be used in parallel. Default value: system cores-1.
"-species [species_name]": specify either "human" or "mouse". Default species: "human". Further species can be added by following the instructions in section xx.
"--left-sided-only": calculate only left-sided p-values.
"--two-sided-only": calculate only two-sided p-values.
"--no-synonyms: disables use of synonyms (synonyms file not required by the script)"
"--print-involved-genes": Run a script after BUFET2 completion that calculates the genes involved in each gene class that are targeted by the miRNA group"
"-involved-genes-filename [filename]: path to the involved genes output file; works only when --print-involved-genes is invoked. Default: involved-genes.txt"
"--disable-interactions-check": disables the validation for the interactions file (not recommended).
"--disable-annotations-check": disables the validation for the annotations file (not recommended).
"--disable-synonyms-check": disables the validation for the synonyms file (not recommended).
"--disable-file-check": disables the validation for all files (not recommended).
'-h" or "--help": print help message and exit
</pre>

#### 2.2.2.4. Adding more species for synonym matching
Open bufet.py and in line 212 add the label you want to use for the species and the taxonomy ID inside the dictionary, enclosed by single or double quotes (', "). The script can now be run using the new label as an argument to the "-species" option.

#### 2.2.2.5. Print target genes involved in each annotation class
Use the argument
<pre>
--print-involved-genes
</pre>
to print the common genes between each annotation class (GO category) and each miRNA in the sample under examination in a file. You can also use the

<pre>
-involved-genes-filename <fileName>
</pre>
argument to specify an output file name other than the default (involved-genes.txt). More specifically, the structure of the output is:
<pre>
>OntologyClass1
miRNA1     Gene1,Gene2,Gene3
miRNA2     
miRNAn     Gene5,Gene6,Gene1
</pre>

An example of the output be seen bellow:
<pre>
>GO:0034199
hsa-miR-6834-3p 
hsa-miR-3142    
hsa-miR-4306    ADCY4,PRKAR2A,ADCY2
hsa-miR-3613-3p PRKAR2A,PRKAR2B,ADCY8,PRKAR1A
hsa-miR-30d-3p  
</pre>

#### 2.2.2.6. Disable gene synonym matching
Use the argument
<pre>
--no-synonyms
</pre>
to disable matching for gene synonyms. In this case a synonyms file is not required and if one is provided, it will be ignored.

#### 2.2.2.7. Example
1. Collect all input files inside the "data" folder.
2. Run the following:
    1. Using compiled script:
        ```bash
        python3 bufet2.py -miRNA data/<mirna_filename> -annotations data/<gene_class_filename> -output data/output.txt -iterations 1000000 -synonyms data/<synonyms_file_name>    
        ```
    with synonyms enabled or
        ```bash
        python3 bufet2.py -miRNA data/<mirna_filename> -annotations data/<gene_class_filename> -output data/output.txt -iterations 1000000 --no-synonyms   
        ```
    to run with synonyms disabled. This should produce a file output.txt in the data folder.
    
    2. Using Docker image:
        ```bash
        run --rm -v data:/data diwis/bufet2 python3 /home/bufet2/bufet2.py -miRNA /data/<mirna_filename> -ontology /data/<gene_class_filename> -iterations 1000000 -synonyms /data/<synonyms_file_name> 
        ```
    with synonyms enabled or 
        ```bash
        run --rm -v data:/data diwis/bufet2 python3 /home/bufet2/bufet2.py -miRNA /data/<mirna_filename> -ontology /data/<gene_class_filename> -iterations 1000000 --no-synonyms
        ```
    to run with synonyms disabled. This should produce a file output.txt in the data folder.

## Using the BUFET2 reverse search module:
TODO

## 4. Reproduction of the BUFET paper's experiments
The dataset.zip file contains all files required to reproduce the p-values experiments contained in the paper. The input miRNA lists are contained in the `inputs/` folder, the gene class data files are contained in the `annotations/` folder and the miRNA-to-gene datasets are contained in the `interactions/` folder.

1. Unzip the dataset.zip file:
```bash
unzip dataset.zip
```
2. Navigate inside the `datasets` folder:
```bash
cd datasets
```
3. Run the software (using Docker) using KEGG and miRTarBase:
    1. Alzheimer's disease:
    ```bash
    docker run --rm -v $(pwd):/data diwis/bufet2 python3 /home/bufet2/bufet2.py -miRNA /data/inputs/alzheimers_mirnas.txt -interactions /data/interactions/mirtarbase.txt -annotations /data/annotations/kegg.csv -iterations 1000000 -output /data/output_alz.txt --no-synonyms

    ```
    2. Breast cancer:
    ```bash
    docker run --rm -v $(pwd):/data diwis/bufet2 python3 /home/bufet2/bufet2.py -miRNA /data/inputs/breast_cancer_mirnas.txt -interactions /data/interactions/mirtarbase.txt -annotations /data/annotations/kegg.csv -iterations 1000000 -output /data/output_br.txt --no-synonyms

    ```
    3. Gastric cancer:
    ```bash
    docker run --rm -v $(pwd):/data diwis/bufet2 python3 /home/bufet2/bufet2.py -miRNA /data/inputs/gastric_cancer_mirnas.txt -interactions /data/interactions/mirtarbase.txt -annotations /data/annotations/kegg.csv -iterations 1000000 -output /data/output_gas.txt --no-synonyms

    ```
    4. Non-small cell lung cancer:
    ```bash
    docker run --rm -v $(pwd):/data diwis/bufet2 python3 /home/bufet2/bufet2.py -miRNA /data/inputs/non_small_cell_lung_cancer_mirnas.txt -interactions /data/interactions/mirtarbase.txt -annotations /data/annotations/kegg.csv -iterations 1000000 -output /data/output_nscl.txt --no-synonyms

    ```
