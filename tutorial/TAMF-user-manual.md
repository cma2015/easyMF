### Brief introduction
- TAMF is designed to provide an easily accessible large-scale transcriptomic data analysis platform to help botanists even with little bioinformatics background completing complicated processing with terminal-based applications that are not user-friendly. It smoothly integrated the TAMF tools into Galaxy scientific analysis platform to provide web-based, easy-to-use and thoroughly tested tools enable users perfom comparative analysis. By integrating three major matrix factorization algorithmes (PCA [principle component analysis]; ICA [independent component analysis]; NMF [non-negative matrix factorization]), TAMF enables users perform a series of analysis based on pattern matrix (PM) and amplitude matrix (AM). For PM-based deep mining, four sub-modules are implemented to perform **Clustering analysis**, **single-cell analysis** (identify cell-types in single-cell RNA-Seq), **Spatial-course analysis** (illuminate the spatial dependency of sampled regions/voxels), and **Time-course analysis** (track the gene expression dynamics across temporal variations along the developmental stages); For AM-based deep mining, **Functional gene discovery** and **Pathway activity analysis** are provided to perform gene function prediction, gain biological insights into GWAS result, as well as active pathways, respectively. The TAMF project is hosted on GitHub (https://github.com/cma2015/TAMF) and can be accessed from http://bioinfo.nwafu.edu.cn:4002. In addition, in order to enable large-scale analysis, we also provided a standardized Docker image: [TAMF Docker image](https://hub.docker.com/r/malab/TAMF/).


### TAMF Docker image installation
- Step 1: [Docker installation](https://github.com/cma2015/TAMF/blob/master/tutorial/Docker-installation.md)
- Step 2: TAMF installation from Docker Hub
  ```bash
  # Pull latest version of TAMF from Docker Hub
  $ docker pull malab/tamf
  ```
- Step 3: Qucikly start
  ```bash
  $ docker run -it -p 4002:4002 malab/tamf bash
  $ cd galaxy
  $ bash run.sh
  ```
  Then you can access TAMF instance via http://localhost:4002


### Upload data

#### Download TAMF test data
- Download test data from [TAMF GitHub project](https://github.com/cma2015/TAMF). Click **Clone or download** (see figure below), and download the ZIP compressed files into your local device and then uncompress it. 

  ![Download data](https://github.com/cma2015/TAMF/blob/master/TAMF_images/1.png)


- For users who installed [Git](https://gist.github.com/derhuerst/1b15ff4652a867391f03), the following command can be used to download TAMF project to local device.
  ```bash
  git clone https://github.com/cma2015/TAMF.git
  ```
#### Upload regular file
- Click **Get Data** in the homepage (see figure below) of TAMF to upload files.


  ![Upload files](https://github.com/cma2015/TAMF/blob/master/TAMF_images/2.png)


  And then you will see the following interface:

  ![](https://github.com/cma2015/TAMF/blob/master/TAMF_images/3.png)


  Next, Click the button **Choose local file** and select a file you would like to upload (e.g. upload the file ```default_ara_gene_expression.txt``` in the directory ```/your directory/TAMF/test_data/Matrix preparation/```), you will see the following interface:


  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/4.png)


  Then click **Start** to upload file.


### MATRIX PREPARATION
  In this module, four sub-functions are provided to obtain the high quality gene expression matrix. 

  In this tutorial, we will show examples about how to use this module to obtain high quality (whose expression structure consistance with the whole expression profile: the Pearson correlation coefficient between samples with the first component [generated from sample-based PCA] higher than 0.75) gene expression for subsequent analysis.

- **User-specific** 
  Gene expression matrix (row represents the genes and the column represents the samples) are required to run this function correctly. Here, we take the file ```default_ara_gene_expression.txt``` located in ```/your directory/TAMF/test_data/Matrix preparation/``` as an example (see figure below). 

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/6.png)

   Then three outputs will be returned:

   **Output 1**: ```statistics analysis on gene expression matrix```, the explained (black dot) and cumulative explained (red dot) variance along the summary gene set patterns.

    ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/7.png)

    **Output 2**: ```prepared gene expression matrix```, the remained high-quality gene expression matrix.

    ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/8.png)

    **Output 3**: ```process information```

    ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/9.png)

- **Fetch from GSM**
  In this function, gene expression matrix will be automatically downloaded according to the GSM provided by user. Here, we use the file ```default_ara_GSM.txt``` located in directory ```/your directory/TAMF/test_data/Matrix preparation/``` to show its usage (see figure below): 

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/10.png)

  Then three outputs will be returned:

  **Output 1**: ```statistics analysis on gene expression matrix```, the explained (black dot) and cumulative explained (red dot) variance along the summary gene set patterns.

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/11.png)

  **Output 2**: ```prepared gene expression matrix```

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/12.png)

   **Output 3**: ```process information```

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/13.png)

- **Arabidopsis thaliana (Default)**
  
  For this function, we provided a built-in gene expression with 1096 samples and 20356 genes in Arabidopsis thaliana obtained from [He et *al*., 2016](https://onlinelibrary.wiley.com/doi/abs/10.1111/tpj.13175), see following figure to see how to use this function:

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/5.png)  


- **Fetch from Key Words** 
  In this function, users can input **Key words** they are interested, then TAMF will fetch the GSM accessions through **Key Words** input by users (see following pictures to see corresponding input and output). 

    **Input:** **words**: input the key words in the box.

    ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/14.png)

    **Output 1**: ```GSM information```, GSM accessions.

    ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/15.png)


###  MATRIX FACTORIZATION
Matrix factorization (MF) is a class of unsurpervised learning methods to decomposing a high-dimensional matrix into a low-dimensional matrix. In this module, we implemented three major MF algorithms (PCA, Principle Component Analysis; ICA, Independent component analysis; NMF, non-negative matrix factorization) to implement matrix factorization, however, to maximize TAMF's usage, user-specific methods are also supported.

- **For PCA, ICA or NMF**
To run this module, the only required input is the prepated gene expression matrix, here, we use the file ```default_ara_prepared_gene_expression.txt``` located in the fold ```/your directory/TAMF/test_data/Matrix factorization/``` to show how to use this module (see following figure). 

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/16.png)

  Then three outputs will be returned: 

  **Output 1**: ```amplitude matrix```

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/17.png)

  **Output 2**: ```pattern matrix```

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/18.png)

  **Output 3**: ```statistics analysis of the decomposition```

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/19.png)

- **User-specific**

  In this function, the decomposed amplitude and pattern matrix with the sample format with test data (```default_ara_PCA_amplitude_matrix.txt``` and ```default_ara_PCA_pattern_matrix.txt``` located in directory ```/your directory/TAMF/test_data/Matrix factorization/```), see following figure:

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/20.png)


### AM-BASED DEEP MINING

MF will generate amplitude matrix and pattern matrix. In AM-based deep mining, two functions including **Functional gene discovery** and **Pathway activity analysis** are implemented.

#### Knowledge data for functional gene discovery

  This module is used to generate knowledge data (**functional gene probabilities construction** and **gene description**) required for sub-function (**Biological interpretation of GWAS result** ) in next module (**Functional gene discovery**). To run this module, the following parameters are required:

  **Input 1:** **amplitude matrix** generated from the prepared gene expression matrix through Matrix factorization tool. Here, we use the file ```default_ara_PCA_amplitude_matrix.txt``` located in directory ```/your directory/TAMF/test_data/Knowledge data/```, and users can select them in the **History** panel.

  **Input 2:** **cpu number** to parallelly compute

  **Input 3:** **functional gene annotation**. Currently, we provide two options: i): choose the **species** in **fetch the functional gene annotation file from EnsemblPlants** window (automaticly fetching from EnsemblPlants); ii): upload the **functional gene annotation** ```default_ara_functional_annota.txt``` and **gene description** ```default_ara_gene_description.txt``` located in the directory ```/your directory/TAMF/test_data/Knowledge data/``` and select them in the **History** panel. 

  **Input 4:** **evidence code**. Golden standard evidence codes need to be selected (default are IDA, IEP, IGI, IMP, IPI, TAS).  
  
  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/21.png)


  Then three files will be returned:


  **Output 1**: ```gene2function construction```, functional gene probabilities construction (the numeric value represents each gene's functionally probability of each unit).

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/22.png)

  **Output 2**: ```golden standard functional terms information```

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/23.png)

  **Output 3**: ```gene description```

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/24.png)


#### Functional gene discovery
  
This module is used to perform functional gene discovery. Currently, we provide two sub-functions: prioritizing the new candidates of the seed genes (**a set of genes for a function/pathway/phenotype**), and gaining biological insight into GWAS results.

For this module, the following inputs are required:

**Input 1:** **amplitude matrix** generated from the prepared gene expression matrix through Matrix factorization tool. Here, we use the file ```default_ara_PCA_amplitude_matrix.txt``` located in directory ```/your directory/TAMF/test_data/Functional gene discovery/```, and users can select them in the **History** panel.

**Input 2:** **cpu number** to parallelly compute

Then, select the functional gene dicsovery options to run. Currently, we provide two sub-functions:

- **new candidates of a gene set**

  For this function, the amplitude matrix and seed genes are required. Here, we use the file ```default_ara_PCA_amplitude_matrix.txt``` located in the directory ```/your directory/TAMF/test_data/Functional gene discovery/``` to show its usage (see figure below). 

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/25.png)

  Then you will obtain two outputs

  **Output 1**: ```AUSR```, area under the curve of the self-ranked genes plot of the gene set through leave-one-out cross-validation

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/26.png)

  **Output 2**: ```new candidates```, the new candidate genes.

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/27.png)


- **Biological interpretation of GWAS result**
  
  For this function, three inputs are required including:

  **Input 3:** **SNPMapGene GWAS result**, a transformed GWAS result (the SNPs map to genes) data generated from a csv format (internal by comma) file (the raw GWAS data) through an R function [SNPMapGene.R](https://github.com/cma2015/RAP2/blob/master/SNPMapGene.R). Upload files (```default_ara_SNPMapGene.txt```) in directory ```/your directory/TAMF/test_data/Functional gene discovery/``` and select them in the **History** panel. 

  **Input 4:** **significant level** default is the 1e-05

  **Input 5:** **functional gene probabilities construction**, select the way to obtain the functional gene probabilities construction. Currently, we provide two options: **Default (Arabidopsis thaliana genes (TAIR10))** or **User-specific**  generated from the **Knowledge data** module ((```default_ara_gene2geneset.txt```) and (```default_ara_gene_description.txt```) in directory ```/your directory/TAMF/test_data/Functional gene discovery/``` )).

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/28.png)

  **Output 1**: ```manhattan plot ``` of the GWAS result

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/29.png)

  **Output 2**: ```visulization of the filter```, heatmap of the significant genes above the significant level

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/30.png)

  **Output 3**: ```enrichment analysis```, enrichment analysis of the significant genes above the significant level 

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/31.png)

  **Output 4**: ```prioritize the potential genes```, prioritize the potential related genes below the significant level 

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/32.png)

  **Output 5**: ```filter false related genes```, the high quality relevant genes after filtering the false related genes above the significant level

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/33.png) 


#### Pathway activity analysis

  This module is used to find the actived pathways and their interactive network. 
  To run this function, four inputs are required including:

  **Input 1**: **prepared gene expression matrix**, generated through the **Matrix preparation** module. Upload file (```default_ara_prepared_gene_expression.txt```) in directory ```/your directory/TAMF/test_data/Pathway activity analysis/```, and select it in the **History** panel. 

  **Input 2**: **amplitude matrix**, generated through the **Matrix factorization** module. Upload file (```default_ara_NMF_amplitude_matrix.txt```) in directory ```/your directory/TAMF/test_data/Pathway activity analysis/```, and select it in the **History** panel. 

  **Input 3**: **pattern matrix** generated through the **Matrix factorization** module. Upload file (```default_ara_NMF_pattern_matrix.txt```) in directory ```/your directory/TAMF/test_data/Pathway activity analysis/```, and select it in the **History** panel. 

  **Input 4**: **pathway annotation** file. Upload file (```default_ara_pathway_annota.txt```) in directory ```/your directory/TAMF/test_data/Pathway activity analysis/```, and select it in the **History** panel. 

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/34.png) 

  
   Then two outputs will be returned:

  **Output 1**: ```actived pathways```, the differentially expressed pathways

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/35.png) 

  **Output 2**: ```interactive network```, the interactive network among the actived pathways

  ![genome-level](https://github.com/cma2015/TAMF/blob/master/TAMF_images/36.png) 


### PM-BASED DEEP MINING
#### Clustering analysis

  This module is used to automatically cluster the samples and return the clustered sample information.
  
  To run this function, two inputs are required including: 

  **Input 1**: the **pattern matrix** generated through the **Matrix factorization** module. Upload file (```default_ara_PCA_pattern_matrix.txt```) in directory ```/your directory/TAMF/test_data/Clustering analysis/```, and select it in the **History** panel. 

  **Input 2**: select the cluster model
  
  ![assembled transcripts](https://github.com/cma2015/TAMF/blob/master/TAMF_images/37.png)

  Click **Execute** to run this function, then four outputs will be returned:

  **Output 1**: ```cluster visulization```

  ![assembled transcripts](https://github.com/cma2015/TAMF/blob/master/TAMF_images/38.png)

  **Output 2**: ```cluster information```, the cluster of the samples in each selected method

  ![assembled transcripts](https://github.com/cma2015/TAMF/blob/master/TAMF_images/39.png)

#### Knowledge data for single cell analysis

  This function generates the **Spec score** for **Single-cell analysis** module. The following parameters need to be input:

  **Input 1:** **cpu number** to parallelly compute

  **Input 2:** **known cell type prepared gene expression matrix** generated from **Matrix preparation**. Upload file (```default_ara_gene_expression.txt```) in directory ```/your directory/TAMF/test_data/Knowledge data/```, and select it in the **History** panel.    
  
  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/43.png)

  Then **Spec score** files will be returned

  ![Upload regular file](https://github.com/cma2015/TAMF/blob/master/TAMF_images/44.png)


#### Single cell analysis

  This module is used to identify cell types of the unknown single cell samples.

  To run this function, five inputs are required including:

  **Input 1**: The **single cell gene expression matrix**. Upload file (```default_ara_gene_expression.txt```) in directory ```/your directory/TAMF/test_data/Single cell analysis/```, and select it in the **History** panel.   
  
  **Input 2**: Select the **decomposition algorithm**

  **Input 3**: The **Spec score** generated from **known cell type prepared gene expression matrix** through **Knowledge data for single cell analysis** module. Upload file (```default_ara_SpecScore.txt```) in directory ```/your directory/TAMF/test_data/Single cell analysis/```, and select it in the **History** panel.  

  **Input 4**: **cpu number** to parallelly compute  

  **Input 5**: Select whether there are **disturb genes**. If yes, then upload them. Upload file (```default_ara_disturb_gene.txt```) in directory ```/your directory/TAMF/test_data/Single cell analysis/```, and select it in the **History** panel.  

  ![nucleic-acid feature](https://github.com/cma2015/TAMF/blob/master/TAMF_images/40.png)

  Then two outputs will be returned:

  **Output 1**: ```single cell type visulization```

  ![assembled transcripts](https://github.com/cma2015/TAMF/blob/master/TAMF_images/41.png)

  **Output 2**: ```identity cell type```, the identification of the unknown single cell samples

  ![assembled transcripts](https://github.com/cma2015/TAMF/blob/master/TAMF_images/42.png)


#### Spatial-course analysis

  This module is used to illuminate the spatial dependency of sampled regions/voxels.

  To run this function, four inputs are required including:

  **Input 1**: The **prepared gene expression matrix** generated through **Matrix preparation** module. Upload file (```default_maize_prepared_gene_expression.txt```) in directory ```/your directory/TAMF/test_data/Spatial-course analysis/```, and select it in the **History** panel.   
  
  **Input 2**: The **amplitude matrix** generated through **Matrix factorization** module. Upload file (```default_maize_NMF_amplitude_matrix.txt```) in directory ```/your directory/TAMF/test_data/Spatial-course analysis/```, and select it in the **History** panel.

  **Input 3**: The **pattern matrix** generated through **Matrix factorization** module. Upload file (```default_maize_NMF_pattern_matrix.txt```) in directory ```/your directory/TAMF/test_data/Spatial-course analysis/```, and select it in the **History** panel.  

  **Input 4**: Select the **species** to **fetch the annotation inforamtion from EnsemblPlants**  

  ![nucleic-acid feature](https://github.com/cma2015/TAMF/blob/master/TAMF_images/45.png)

  Then two outputs will be returned:

  **Output 1**: ```enrichment analysis for each tissue```

  ![assembled transcripts](https://github.com/cma2015/TAMF/blob/master/TAMF_images/46.png)

  **Output 2**: ```spatial specific genes expression module```

  ![assembled transcripts](https://github.com/cma2015/TAMF/blob/master/TAMF_images/47.png)


#### Time-course analysis

  This module is used to track the gene expression dynamics across temporal variations along the developmental stages.


  To run this function, six inputs are required including:

  **Input 1**: The **prepared gene expression matrix** generated through **Matrix preparation** module. Upload file (```default_maize_prepared_gene_expression.txt```) in directory ```/your directory/TAMF/test_data/Time-course analysis/```, and select it in the **History** panel.   
  
  **Input 2**: The **amplitude matrix** generated through **Matrix factorization** module. Upload file (```default_maize_NMF_amplitude_matrix.txt```) in directory ```/your directory/TAMF/test_data/Time-course analysis/```, and select it in the **History** panel.

  **Input 3**: The **pattern matrix** generated through **Matrix factorization** module. Upload file (```default_maize_NMF_pattern_matrix.txt```) in directory ```/your directory/TAMF/test_data/Time cell analysis/```, and select it in the **History** panel.  

  **Input 4**: The **developmental information of the samples**. Upload file (```default_maize_sample_stage.txt```) in directory ```/your directory/TAMF/test_data/Time cell analysis/```, and select it in the **History** panel.   

  **Input 5**: **cpu number** to parallelly compute   

  **Input 6**: Select the **species** to **fetch the annotation inforamtion from EnsemblPlants**  

  ![nucleic-acid feature](https://github.com/cma2015/TAMF/blob/master/TAMF_images/48.png)

  Then two outputs will be returned:

  **Output 1**: ```significant genes of each module```

  ![assembled transcripts](https://github.com/cma2015/TAMF/blob/master/TAMF_images/49.png)

  **Output 2**: ```visulization of the speficic module```, gene expression dynamic patterns.

  ![assembled transcripts](https://github.com/cma2015/TAMF/blob/master/TAMF_images/50.png)

