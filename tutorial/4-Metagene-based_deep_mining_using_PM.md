<div align="center"><h1><b>easyMF User Mannual</b></h1></div>

<div align="center">(version 1.0)</div>

easyMF is a user-friendly web platform that aims to facilitate biological discovery from large-scale transcriptome data through matrix factorization (MF). It offers several functional tools for gene expression matrix generation, and metagene-based exploratory analysis including sample clustering, signature gene identification, functional gene discovery, subtype cell detection, and pathway activity inference.

- easyMF project is hosted on https://github.com/cma2015/easyMF.
- easyMF docker image is available in https://hub.docker.com/r/malab/easymf.
- easyMF demo server can be accessed via [http://easyMF.omicstudio.cloud](http://deepea.omicstudio.cloud/).
- The following part shows installation of easyMF docker image and detailed documentation for each function in easyMF.



## 0. Metagene-based Deep Mining Using PM

Pattern matrix (AM), a matrix with metagenes in rows and samples in columns, describes sample-level relationships. In current version of easyMF, users can make use of PM for sample clustering, temporal and spatial transcriptome analysis, and subtype cell detection.

This module consists of three functions: **Sample Clustering Analysis**, **Temporal-spatial Analysis**, and **Subtype Cell Detection**.

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
          <td rowspan="6">Sample Clustering Analysis</td>
          <td rowspan="6">Cluter samples automatically based on the pattern matrix</td>
          <td rowspan="6">Pattern matrix</td>
          <td rowspan="6">Cluster information; Cluster visualization; </td>
          <td rowspan="6">~ 15s</td>
          <td>mclust</td>
          <td>Scrucca <I>et al</I>., 2016</td>
      </tr>
      <tr>
          <td>apcluster</td>
          <td>Bodenhofer <I>et al</I>., 2011</td>
      </tr>
      <tr>
          <td>SSE</td>
          <td>This study</td>
      </tr>
      <tr>
          <td>fpc</td>
          <td>Hennig, 2013</td>
      </tr>
      <tr>
          <td>vegan</td>
          <td>Dixon, 2003</td>
      </tr>
      <tr>
          <td>gap</td>
          <td>Maechler <I>et al</I>., 2012</td>
      </tr>
      <tr>
          <td rowspan="3">Temporal-spatial Analysis</td>
          <td rowspan="3">Detect functional modules and identity signature genes </td>
          <td rowspan="3">Gene expression matrix; Amplitude matrix; Pattern matrix; Sample information</td>
          <td rowspan="3">Signature genes of each module; GO enrichment analysis of each module; Visualization of specific module</td>
          <td rowspan="3">~ 5 mins</td>
          <td>In-house scripts</td>
          <td>This study</td>
      </tr>
      <tr>
          <td>cogaps</td>
          <td>Stein-O'Brien <I>et al</I>., 2017</td>
      </tr>
      <tr>
          <td>topGO</td>
          <td>Alexa and Rahnenführer, 2009</td>
      </tr>
      <tr>
          <td rowspan="4">Subtype Cell Detection</td>
          <td rowspan="4">Identity a cell type of unknown single-cell RNA-Seq samples</td>
          <td rowspan="4">Gene expression matrix; Spec score</td>
          <td rowspan="4">Identification results of cell type; Single cell type visulization</td>
          <td rowspan="4">~ 30s</td>
          <td>prcomp</td>
          <td>This study</td>
      </tr>
      <tr>
          <td>ica</td>
          <td>Helwig, 2015</td>
      </tr>
      <tr>
          <td>bignmf</td>
          <td>Helwig, 2015</td>
      </tr>
      <tr>
          <td>In-house scripts</td>
          <td>This study</td>
      </tr>   

## 1. Sample Clustering Analysis

In the current version, easyMF provides six optional algorithms (mclust, apcluster, SSE, fpc, vegan , and gap) to cluster samples using PM coefficients. The clusters are visualized in plots, as well as tables, providing a quick overview of the relationships between samples.

#### Inputs

- **Pattern matrix**:  A pattern matrix with metagenes in rows and samples in columns.

- **Cluster algorithms**: easyMF provides six cluster algorithms to be cluster samples including mclust, apcluster, SSE, fpc, vegan, and gap.

  **Note**: easyMF supports multi-selection for different algorithms.

#### Outputs

- **cluster visulization**:  A dot plot of the clustering result.
- **cluster information**: A cluster of the samples in each selected method.

#### How to use this function

- The following screenshot shows us how to implement functional gene discovery using easyMF.

  ![1-1](../assets/img/1-1.png)



## 2. Temporal-spatial Analysis

easyMF can be used to determine the extent to which genes change over time in response to perturbations (e.g., developmental time), and does so by integrating gene expression values, and gene- and sample-level relationships. It can also be used to identify signature genes dominated at specific compartments of the transcriptomes with spatial resolution in individual tissue samples (spatial transcriptomes).

#### Inputs

In **Data** section

- **Gene expression matrix**: A gene expression matrix generated by the module **Matrix Preparation**.
- **Amplitude matrix**: An amplitude matrix decomposed by the gene expression matrix input.
- **Pattern matrix**: A pattern matrix decomposed by the gene expression matrix input.
- **Sample information**: Sample information containing development stages or spatial compartments.

In **Parameters** section

- **Threshold of Pearson Correlation Coefficient**: A Pearson Correlation Coefficient value used for signature gene identification.
- **Threshold of P-value**: A *P*-value used for signature gene identification.
- **Select a species**: Species name which can be selected from drop-down menu is used to annotate gene information.

#### Outputs

- **Signature genes of each metagene**: Summary of signature genes in each metagene, which contains **Metagene**, **Gene ID**, and **Gene description**.
- **GO enrichment analysis of each metagene**: Summary of GO enrichment results of each metagene.
- **Visualization of metagenes**: Hierarchical clustering analysis of pattern matrix.

#### How to use this function

- The following screenshot shows us how to implement temporal-spatial analysis using easyMF.

  ![1-1](../assets/img/1-1.png)

  

## 3. Subtype Cell Detection

Single-cell RNA-Seq, which measures gene expressions at the level of a single cell, has been developed as a powerful technique to investigate the function of individual cells. easyMF can be used to identify cell types of unknown cells, which is one of key steps in the process of single-cell transcriptome analysis.

#### Inputs

- **Gene expression matrix**: Gene count matrix of the single cell.

- **Decomposition options**: Currently, easyMF provides three algorithms to decompose gene expression matrix including **PCA**, **ICA**, and **NMF**.

- **Spec score**: Spec score of the known type of tissues. Here's an example:

  | Atrichoblast | Columella | Cortex | Endodermis |       |
  | :----------- | :-------- | :----- | :--------- | ----- |
  | ORF25        | 0         | 0.127  | 0          | 0     |
  | AT2G26550    | 0         | 0      | 0          | 0     |
  | AT2G26570    | 0         | 0      | 0.102      | 0.591 |
  | PSBA         | 0         | 0      | 0.105      | 0.105 |

In addition, **Spec score** can also be obtained through a known tissue expression matrix.

#### Outputs

- **Single cell type visualization**
- **Cell type detection result**: Summary of cell type detection results.