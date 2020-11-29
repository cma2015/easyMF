<div align="center"><h1><b>easyMF User Mannual</b></h1></div>

<div align="center">(version 1.0)</div>

easyMF is a user-friendly web platform that aims to facilitate biological discovery from large-scale transcriptome data through matrix factorization (MF). It offers several functional tools for gene expression matrix generation, and metagene-based exploratory analysis including sample clustering, signature gene identification, functional gene discovery, subtype cell detection, and pathway activity inference.

- easyMF project is hosted on https://github.com/cma2015/easyMF.
- easyMF docker image is available in https://hub.docker.com/r/malab/easymf.
- easyMF demo server can be accessed via [http://easyMF.omicstudio.cloud](http://deepea.omicstudio.cloud/).
- The following part shows installation of easyMF docker image and detailed documentation for each function in easyMF.



## 0. Metagene-based Deep Mining Using AM

Amplitude matrix (AM), a matrix with genes in rows and metagenes in columns, describes gene-level relationships. In current version of easyMF, users can make use of AM for functional gene discovery and pathway activity inference.

This module consists of three functions: **Functional Gene Discovery**, and **Pathway Activity Inference**.

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

- **Amplitude matrix**:  An amplitude matrix with genes in rows and metagenes in columns.
- **Seed genes**: A list of genes which are associated with a specific biology function.


#### Outputs

- **Gene score and rank**: Summary of gene prioritization result containing score and rank information of each gene.

- **Area under the self-ranked curve (AUSR) plot**: A plot of ratio (*Ra*) along the *y* axis versus self-rank along the *x* axis, where rank represents the ranks of all positive genes, *Ra*(*l*) represents the ratio of ranks lower than a pre-defined level of *l*.


#### How to use this function

- The following screenshot shows us how to implement functional gene discovery using easyMF.

	![1-1](../assets/img/1-1.png)
	
	

## 2. Pathway Activity Inference

Pathway activity inference can be used to examine the pathway activity for any gene set of interest.

#### Inputs

- **Amplitude matrix**: An amplitude matrix with genes in rows and metagenes in columns. 

- **Pathway annotation**：A pathway annotation file, which contains **Gene ID**, **Pathway ID**, and **Pathway name** separated by a tab character. Here is an example:

  | Gene ID        | Pathway ID | Pathway name                 |
  | :------------- | :--------- | :--------------------------- |
  | Zm00001d042869 | zma00010   | Glycolysis / Gluconeogenesis |
  | Zm00001d025586 | zma00010   | Glycolysis / Gluconeogenesis |
  | Zm00001d039089 | zma00020   | Citrate cycle (TCA cycle)    |
  | Zm00001d037278 | zma00030   | Pentose phosphate pathway    |

- **Gene set**: A list of genes used to estimate pathway activity.


#### Outputs

- **Pathway activity**: The activity of each pathway. In the result, active pathway can be obtained through a *p*-value filtration of each pathway information.


#### How to use this function

- The following screenshot shows us how to implement pathway activity inference  using easyMF.

  ![1-1](../assets/img/1-1.png)

