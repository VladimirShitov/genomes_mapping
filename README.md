# Genomes mapping

Project for creating a total list of genes for a bacteria.

Pipeline can be run by downloading repository and executing *main.py* file. Before that, please download all the required 
python packages:
* **numpy** – pipeline uses some arrays from numpy
* **matplotlib** – plots
* **seaborn** – plots
* **pandas** – working with tables
* **biopython** – for running BLAST and using NCBI api
* **tqdm** – progress bar

Also, please open file constants.py and set your email adress to the EMAIL constant.

File main.py executes following scripts from *./functions*:
1. **set_enviroment.py**
2. **download_genomes.py**
3. **create_db.py**
4. **align_to_db.py**

All the steps, that take a lot of time have **tqdm** inside, so you can always see progress bar in your python console. 
*main.py* writes result of it's work to a file *./logs/pipeline.log*. All the steps off the pipeline described below in details.

## set_enviroment.py

Creates subdirectories needed for the other steps of the pipeline.

Directory tree after work of the pipeline will look like this (except for .py files):
```bash
    .
    ├── data – folder where all the genomes and output files will be saved
    │   ├── assembly_summary.txt – table with all organism assemblies from refseq
    │   ├── GCF_000005845.2_ASM584v2 – folder with files for that assembly
    │   │   ├── GCF_000005845.2_ASM584v2_db_alignment.csv – generated by align_to_db()
    │   │   ├── GCF_000005845.2_ASM584v2_feature_table.txt – from NCBI
    │   │   ├── GCF_000005845.2_ASM584v2_genomic.fna – from NCBI
    │   │   ├── GCF_000005845.2_ASM584v2_protein.faa – from NCBI
    │   │   └── GCF_000005845.2_ASM584v2_protein.gpff – from NCBI
    │   ├── GCF_000006665.1_ASM666v1
    │   ...
    ├── logs
    │   ├── alignment.log – log of align_to_database()
    │   ├── database.log – log of create_database()
    │   ├── downloading.log – log of download_genomes()
    │   └── pipeline.log – general log
    └── plots – folder with graphs
```

## download_genomes.py
Downloads *assembly_summary.txt* from NCBI ftp-site to *./data* folder. For example, for Escherichia coli the file will be downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Escherichia_coli/ . 
If you want to use another organism, please set it's name to **ORGANISM** constant in *./constants.py*. (Note: in that case you will also have to set assembly id of the reference organism to a **REFERENCE** constant)

Then, only assemblies with **assembly_level** "Complete Genome" or "Chromosome" are selected. For each of them 
the following files will be downloaded:
* **protein.faa.gz** – contains amino acid sequences.
* **feature_table.txt.gz** – contains information about genes, i.e, *start* and *end* in a genome.
* **genomic.fna.gz** – complete nucleotide sequence of a genome.
* **protein.gpff.gz** – protein sequence and annotation for each gene.
Then all the files are extracted from their archives.

Information about the work process of a function is being written to a *./logs/downloading.log*.

## create_db.py
From all of the protein sequences of every genome creates a list of such sequences, that any pair of sequences from that list
has **identity < 50%**. Following algorithm is used for each genome (by "genome" I understand folder in *./data*)
starting from reference:
1. If genome doens't have *feature_table.txt* or *protein.faa* files, delete the whole folder.
2. Else perform a BLAST search for each gene against the genes in database (for the first run with reference database is empty). If a gene doesn't have BLAST hit with **E value > 0.001** and **identity > 50%**, add this gene to the database.

In the end we will get blast protein database with *cur_db* prefix and list of all genes *current_db.faa* in *./blastdb* folder.

Note, that number of genes in database doesn't fully represent number of unique genes of the organism. Imagine following situation: 3 genes *A*, *B* and *C*, such that *A* identity with *B* and with *C* > 50%, but identity of *B* and *C* is less than 50%.

![Illustration 1](/images/1.png)

Imagine, that algorithm met the gene *A* first. Then later genes *B* and *C* will align on it well enough, and **won't be added** to the database.

![Illustration 2](/images/2.png)

Now imagine, that algorithm met the gene *C* first. Then later the gene *A* will align on it well and **will not be added** to the database. That means, that gene *B* **will be added** to the database too. The same is true for situation when algorithm meets gene *B* first.

![Illustration 3](/images/3.png)

This means, that we will need to align all the genes on database genes later. Then we could define **gene** as a linked component in obtained graph. But we will have to be careful with situations like this:

![Illustration 4](/images/4.png)

Do *A* and *X* represent the same gene?
