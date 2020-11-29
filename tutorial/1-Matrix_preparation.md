<div align="center"><h1><b>easyMF User Mannual</b></h1></div>

<div align="center">(version 1.0)</div>

easyMF is a user-friendly web platform that aims to facilitate biological discovery from large-scale transcriptome data through matrix factorization (MF). It offers several functional tools for gene expression matrix generation, and metagene-based exploratory analysis including sample clustering, signature gene identification, functional gene discovery, subtype cell detection, and pathway activity inference.

- easyMF project is hosted on https://github.com/cma2015/easyMF.
- easyMF docker image is available in https://hub.docker.com/r/malab/easymf.
- easyMF demo server can be accessed via [http://easyMF.omicstudio.cloud](http://deepea.omicstudio.cloud/).
- The following part shows installation of easyMF docker image and detailed documentation for each function in easyMF.



## 0. Matrix Preparation

Matrix factorization is typically started with an input of a gene expression matrix (genes in rows and individual samples in columns), which prompt us to design this module including three functions to prepare a high-quality gene expression matrix for downstream analysis.

The gene expression matrix can be automatically generated from raw reads using a bioinformatics pipeline (see following figure, Figure S1).

![Figure S1](/images/Figure S1.png)

This module consists of three functions: **Data Retrieval**, **Expression Matrix Generation**. and **Expression Quality Control**.

  <table class="fl-table">
  <thead>
    <tr>
      <th width="15%">Functions/Tools</th>
      <th width="15%">Description</th>
      <th width="15%">Inputs</th>
      <th width="15%">Outputs</th>
      <th width="15%">Time (test data)</th>
      <th width="15%">Program</th>
      <th width="25%">References</th>
    </tr>
  </thead>
  <tbody>
      <tr>
          <td>Data Retrieval</td>
          <td>Retrieve genome sequences, genome annotation, and RNA-Seq data automatically from public databases</td>
          <td>Select a species; Database version; Data type</td>
          <td>Genome sequences (in terms of Reference genome sequence); Genome annotation (in terms of Reference genome annotation); RNA-Seq data (in terms of Raw RNA-Seq data)</td>
          <td>Depends on network speed</td>
          <td>In-house scripts</td>
          <td>This study</td>
      </tr>
      <tr>
          <td rowspan="4">Expression Matrix Generation</td>
          <td rowspan="4">Generate a gene expression matrix (genes in rows and individual samples in columns) through raw RNA-Seq quality control, read-genome alignment, and gene expression abundance calculation</td>
          <td rowspan="4">Genome sequence and annotation; RNA-Seq data</td>
          <td rowspan="4">Gene expression matrix</td>
          <td rowspan="4">~ 4 mins</td>
          <td>fastp (Raw RNA-Seq quality control)</td>
          <td>Chen <I>et al</I>., 2018</td>
      </tr>
      <tr>
          <td>HISAT2 (Read-genome alignment)</td>
          <td>Kim <I>et al</I>., 2015</td>
      </tr>
      <tr>
          <td>StringTie (Gene expression abundance calculation)</td>
          <td>Pertea <I>et al</I>., 2015</td>
      </tr>
      <tr>
          <td>In-house scripts</td>
          <td>This study</td>
      </tr>
      <tr>
          <td rowspan="2">Expression Quality Control</td>
          <td rowspan="2">Generate a high-quality gene expression matrix through removing lowly expressed genes, outlier samples, or batch effects</td>
          <td rowspan="2">Raw gene expression matrix</td>
          <td rowspan="2">High-quality gene expression matrix</td>
          <td rowspan="2">~ 10s</td>
          <td>In-house scripts (Removing lowly expressed genes and outlier samples)</td>
          <td>This study</td>
      </tr>
      <tr>
          <td>sva (Removing batch effects)</td>
          <td>Leek <I>et al</I>., 2012</td>
      </tr>

## 1. Data Retrieval

Data Retrieval can be used to retrieve genome sequences, genome annotation, and RNA-Seq data automatically from public databases.

#### Inputs

For retrieving **genome sequences and annotation**, users need to select the option "Obtain Genome Sequences and Annotation" firstly.

- **Select a species**: This option provides the Latin name of 61 species to be selected.

- **Database version**: Ensembl releases from 25 to 47 are listed.

- **Data type**: Genome sequences (.fasta), or annotation (.gtf).

For retrieving **RNA-Seq** data, users need to select the option "Obtain RNA-Seq data".

- **Fecth data through data ID or ftp address**: easyMF provides two ways to download RNA-Seq data. 

  If users select **Fetch data through data ID**, easyMF download RNA-Seq data by NCBI's tool *sratoolkit* (version 2.3.5) using RNA-Seq IDs (such as SRR1765337) user provided. 

  If users select **Fetch data through data address**, easyMF download RNA-Seq data by wget using HTTP/FTP addresses user provided.


#### Outputs

For "Obtain Genome Sequences and Annotation"

- Sequences in FASTA format or annotation in GTF format.

For "Obtain RNA-Seq data"

- A compressed file containing RNA-Seq data in tar.gz format.


#### How to use this function

- The following screenshot shows us how to download genome sequences and annoatation using easyMF.

	![1-1](../assets/img/1-1.png)
	
- The following screenshot shows us how to download RNA-Seq data using easyMF.

	![1-1](../assets/img/1-1.png)
	

## 2. Expression Matrix Generation

This function can be used to generate a gene expression matrix (genes in rows and individual samples in columns) through raw RNA-Seq quality control, read-genome alignment, and gene expression abundance calculation.

#### Inputs

In **Data** section

- **Reference genome sequence**: Reference genome sequence in FASTA format used for read-genome alignment.

- **Reference genome annotation**: Reference genome annotation in GTF format used to estimate gene expression abundance.

- **Raw RNA-Seq data**: A compressed file containing RNA-Seq data in tar.gz format.

In **Parameters** section, easyMF need users set parameters used for "RNA-Seq quality control" and "Read-genome alignment".

For "RNA-Seq quality control"

- **Minimum read length**: A threshold of read length that reads shorter than the length will be discarded.

- **The quality value that a base is qualified**: A threshold of base quality value to trim low-quality reads.

For "Read-genome alignment"

- **Minimum intron length for RNA-Seq alignment**

- **Maximum intron length for RNA-Seq alignment**

#### Output

- **Raw Gene Expression Matrix**: A gene expression matrix (genes in rows and individual samples in columns).


#### How to use this function

- The following screenshot shows us how to generate a gene expression matrix using easyMF.

	![1-1](../assets/img/1-1.png)
		

## 3. Expression Quality Control

Once gene expression matrix generated, to accurately implement MF-based analysis, quality of the gene expression matrix need to be improved, which can be operated through three different dimensions including **removing lowly expressed genes**, **removing outlier samples**, **removing batch effects**.

#### Inputs

For **Removing lowly expressed genes**

- **Expression value of expressed genes**: Expression value of genes regarded as expressed.

- **Minimum sample number**: The number of samples of expressed genes.

easyMF provides default values for these two parameters: **Expression value of expressed genes** (default as 1) and **Minimum sample number** (default as 3), which means genes regarded as expressed with expression value greater than 1 in at least 3 samples.

For **Removing outlier samples**

- **Threshold of potential repeat samples**: Expression values between two samples are almost identical.

- **Threshold of low-quality samples**: The structure of the specific sample is away from the others.

For **Removing batch effects**

- **Sample information**: RNA-Seq samples with batch information. In the text, the first column presents sample IDs, and the second column presents batch information distinguished with Arabic numerals.

	```
		SRR1765379	1
		SRR1765380	1
		SRR1765337	2
		SRR1765338	2
	```
	
#### Output

- **High-quality gene expression matrix**: A high-quality gene expression matrix (genes in rows and individual samples in columns).

#### How to use this function

- The following screenshot shows us how to generate a high-quality gene expression matrix using easyMF.

	![1-1](../assets/img/1-1.png)






