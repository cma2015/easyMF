<div align="center"><h1><b>easyMF User Mannual</b></h1></div>

<div align="center">(version 1.0)</div>

easyMF is a user-friendly web platform that aims to facilitate biological discovery from large-scale transcriptome data through matrix factorization (MF). It offers several functional tools for gene expression matrix generation, expression matrix factorization, and metagene-based exploratory analysis including sample clustering, signature gene identification, functional gene discovery, subtype cell detection, and pathway activity inference.

- easyMF project is hosted on https://github.com/cma2015/easyMF.
- easyMF docker image is available in https://hub.docker.com/r/malab/easymf.
- easyMF demo server can be accessed via [http://easyMF.omicstudio.cloud](http://deepea.omicstudio.cloud/).
- The following part shows installation of easyMF docker image and detailed documentation for each function in easyMF.



## 0. Metagene-based Deep Mining Using AM

Amplitude matrix (AM), a matrix with genes in rows and metagenes in columns, describes gene-level relationships. In current version of easyMF, users can make use of AM for functional gene discovery and pathway activity inference.

This module consists of two functions: **Functional Gene Discovery**, and **Pathway Activity Inference**.

  <table class="fl-table">
  <thead>
    <tr>
      <th width="15%">Functions/Tools</th>
      <th width="15%">Description</th>
      <th width="15%">Inputs</th>
      <th width="15%">Outputs</th>
      <th width="10%">Time (test data)</th>
      <th width="15%">Program</th>
      <th width="15%">References</th>
    </tr>
  </thead>
  <tbody>
      <tr>
          <td>Functional Gene Discovery</td>
          <td>Calculate gene score and rank genes based on the probability of their association with a specific biology function</td>
          <td>Amplitude matrix; A set of genes with a specific characteristic</td>
          <td>Gene score and rank; Area under the self-ranked curve (AUSR) plot</td>
          <td>~ 10s</td>
          <td>In-house scripts</td>
          <td>Fehrmann <I>et al</I>., 2015</td>
      </tr>
      <tr>
          <td>Pathway Activity Inference</td>
          <td>Examine the pathway activity for any gene set of interest</td>
          <td>Amplitude matrix; Pathway annotation; A set of genes with a specific characteristic</td>
          <td>Actived pathways</td>
          <td>~ 1 mins</td>
          <td>In-house scripts</td>
          <td>This study</td>
      </tr>


##  1. Functional Gene Discovery 

Functional gene discovery can be used to calculate gene score and rank genes based on the probability of their association with a specific biology function.

#### Inputs

- **Amplitude matrix**:  An amplitude matrix of AM coefficients with genes in rows and metagenes in columns. Here is an example:

<table class="fl-table">
    <tr>
        <td></td>
        <td width="20%">Metagene 1</td>
        <td width="20%">Metagene 2</td>
        <td width="20%">Metagene 3</td>
        <td width="20%">...</td>
        <td width="20%">Metagene n</td>
    </tr>
    <tr>
        <td>Zm00001d053636</td>
        <td>0.080</td>
        <td>-0.889</td>
        <td>1.504</td>
        <td>...</td>
        <td>2.029</td>
    </tr>
    <tr>
        <td>Zm00001d053632</td>
        <td>1.338</td>
        <td>0.729</td>
        <td>-0.113</td>
        <td>...</td>
        <td>-0.049</td>
    </tr>
    <tr>
        <td>...</td>
        <td>...</td>
        <td>...</td>
        <td>...</td>
        <td>...</td>
        <td>...</td>
    </tr>
    <tr>
        <td>Zm00001d053635</td>
        <td>-1.674</td>
        <td>0.036</td>
        <td>-0.047</td>
        <td>...</td>
        <td>-0.494</td>
    </tr>
*  **Functional genes**: A set of genes associated with a specific biology function, such as enriched in a phenotype of interest. If users select **Upload a file with functional gene IDs from local disk**, a newline-delimited file containing gene IDs needs to be provided; if users select **Enter functional gene IDs**, gene IDs need to be separated by comma. Here are two examples:

​       A newline-delimited file containing gene IDs for **Upload a file with functional gene IDs from local disk**:

```
Zm00001d053636
Zm00001d053632
Zm00001d053630
...
Zm00001d053635
```

​     Comma-separated gene IDs for **Enter functional gene IDs**:

```
Zm00001d053636,Zm00001d053632,Zm00001d053630,...,Zm00001d053635
```


#### Outputs

- **Gene score and rank**: Summary of gene prioritization results. Each column shows **Gene ID**, **Score**, **Rank**, and **Annotation**. The higher ranking of a gene, the more related to the biological function.  Here is an example:

<table class="fl-table">
    <tr>
        <td width="25%">Gene ID</td>
        <td width="25%">Score</td>
        <td width="25%">Rank</td>
        <td width="25%">Annotation</td>
    </tr>
    <tr>
        <td>Zm00001d053636</td>
        <td>1</td>
        <td>1</td>
        <td>Label</td>
    </tr>
    <tr>
        <td>Zm00001d053632</td>
        <td>0.888</td>
        <td>3</td>
        <td>Label</td>
     </tr>
    <tr>
        <td>Zm00001d004839</td>
        <td>1</td>
        <td>1</td>
        <td>Unlabel</td>
    </tr>
    <tr>
        <td>...</td>
        <td>...</td>
        <td>...</td>
        <td>...</td>
    </tr>
    <tr>
        <td>Zm00001d053635</td>
        <td>0.92</td>
        <td>2</td>
        <td>Unlabel</td>
    </tr>


- **Area under the self-ranked curve (AUSR) plot**: A plot of ratio (*Ra*) along the *y* axis versus self-rank along the *x* axis, where rank represents the ranks of all positive genes, *Ra*(*l*) represents the ratio of ranks lower than a pre-defined level of *l*.

  ![06_00](/easyMF_images/06_00_Functional_gene_discovery.png)


#### How to use this function

- Test data for this function are in directory `Test_data/03_Metagene-based_Deep_Mining_Using_AM` including.

- The following screenshot shows us how to implement functional gene discovery using easyMF.

  **Step 1**: upload test data in directory `Test_data/03_Metagene-based_Deep_Mining_Using_AM` to history panel;

  ![06-01](/easyMF_images/06_01_Functional_gene_discovery.png)

  **Step 2** input the corresponding files and appropriate parameters, then run the function.

  ![06_02](/easyMF_images/06_02_Functional_gene_discovery.png)

#### Running time

This step will cost ~ 10s for the test data.



## 2. Pathway Activity Inference

Pathway activity inference can be used to examine the pathway activity for any gene set of interest.

#### Inputs

- **Amplitude matrix**: An amplitude matrix of AM coefficients with genes in rows and metagenes in columns. Here is an example:

  <table class="fl-table">
      <tr>
          <td></td>
          <td width="20%">Metagene 1</td>
          <td width="20%">Metagene 2</td>
          <td width="20%">Metagene 3</td>
          <td width="20%">...</td>
          <td width="20%">Metagene n</td>
      </tr>
      <tr>
          <td>Zm00001d053636</td>
          <td>0.080</td>
          <td>-0.889</td>
          <td>1.504</td>
          <td>...</td>
          <td>2.029</td>
      </tr>
      <tr>
          <td>Zm00001d053632</td>
          <td>1.338</td>
          <td>0.729</td>
          <td>-0.113</td>
          <td>...</td>
          <td>-0.049</td>
      </tr>
      <tr>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
          <td>...</td>
      </tr>
      <tr>
          <td>Zm00001d053635</td>
          <td>-1.674</td>
          <td>0.036</td>
          <td>-0.047</td>
          <td>...</td>
          <td>-0.494</td>
      </tr>

- **Pathway annotation**: A pathway annotation file, which contains **Gene ID**, **Pathway ID**, and **Pathway name** separated by a tab character. Here is an example:

  | Gene ID        | Pathway ID | Pathway name                 |
  | :------------- | :--------- | :--------------------------- |
  | Zm00001d042869 | zma00010   | Glycolysis / Gluconeogenesis |
  | Zm00001d025586 | zma00010   | Glycolysis / Gluconeogenesis |
  | Zm00001d039089 | zma00020   | Citrate cycle (TCA cycle)    |
  | Zm00001d037278 | zma00030   | Pentose phosphate pathway    |

- **Gene set**: A list of genes used to estimate pathway activity.


#### Outputs

- **Pathway activity**: The activity of each pathway. Each column shows **Pathway**, ***P*-value**, **FDR**, **Term**, **Significant**, **Annotate** and AM coefficient in each metagene. In the result, active pathway can be obtained through a *p*-value filtration of each pathway information.

| Pathway  | P-value | FDR   | Term                            | Significant | Annotate | Metagene1 | Metagenen2 |
| :------- | :------ | :---- | :------------------------------ | :---------- | :------- | :-------- | :--------- |
| zma00592 | 0.077   | 0.977 | alpha-Linolenic acid metabolism | 3           | 33       | 0.001     | 0.022      |
| zma04146 | 0.101   | 0.934 | Peroxisome                      | 2           | 67       | 0.016     | 0.035      |
| zma00010 | 0.177   | 0.864 | Glycolysis / Gluconeogenesis    | 2           | 102      | 0.031     | 0.112      |
| zma00906 | 0.227   | 0.981 | Carotenoid biosynthesis         | 2           | 27       | 0         | 0.012      |

#### How to use this function

- Test data for this function are in directory `Test_data/03_Metagene-based_Deep_Mining_Using_AM` including `01_Amplitude_matrix` and `01_Functional_gene_list`.

- The following screenshot shows us how to examine the pathway activity using easyMF.

  **Step 1**: upload test data in directory `Test_data/03_Metagene-based_Deep_Mining_Using_AM` to history panel;

  ![07_01](/easyMF_images/07_01_Pathway_activity_inference.png)

  **Step 2** input the corresponding files and appropriate parameters, then run the function.

  ![07_02](/easyMF_images/07_02_Pathway_activity_inference.png)

#### Running time

This step will cost ~ 1 mins for the test data.