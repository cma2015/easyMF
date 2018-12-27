### Brief introduction
- CAFU is a Galaxy-based bioinformatics framework for comprehensive assembly and functional annotation of unmapped RNA-seq data from single- and mixed-species samples which integrates plenty of existing NGS analytical tools and our developed programs, and features an easy-to-use interface to manage, manipulate and most importantly, explore large-scale unmapped reads. Besides the common process of reads cleansing, reads mapping, unmapped reads generation and novel transcription assembly, CAFU optionally offers the multiple-level evidence analysis of assembled transcripts, the sequence and expression characteristics of assembled transcripts, and the functional exploration of assembled transcripts through gene co-expression analysis and genome-wide association analysis. Taking the advantages of machine learning (ML) technologies, CAFU also effectively addresses the challenge of classifying species-specific transcript assembled using unmapped reads from mixed-species samples. The CAFU project is hosted on GitHub(https://github.com/cma2015/CAFU) and can be accessed from http://bioinfo.nwafu.edu.cn:4001. In addition, in order to enable large-scale analysis, we also provided a standardized Docker image: [CAFU Docker image](https://hub.docker.com/r/malab/cafu/).


### CAFU Docker image installation
- Step 1: [Docker installation](https://github.com/cma2015/CAFU/blob/master/tutorial/Docker_installation.md)
- Step 2: CAFU installation from Docker Hub
  ```bash
  # Pull latest version of CAFU from Docker Hub
  $ docker pull malab/cafu
  ```
- Step 3: Qucikly start
  ```bash
  $ docker run -it -p 80:80 malab/cafu bash
  $ cd /home/galaxy
  $ bash run.sh
  ```
  Then you can access CAFU instance via http://localhost:80


### Upload data

#### Download CAFU test data
- Download test data from [CAFU GitHub project](https://github.com/cma2015/CAFU). Click **Clone or download** (see figure below), and download the ZIP compressed files into your local device and then uncompress it. 

  ![Download data](https://github.com/cma2015/CAFU/blob/master/CAFU_images/1.png)


- For users who installed [Git](https://gist.github.com/derhuerst/1b15ff4652a867391f03), the following command can be used to download CAFU project to local device.
  ```bash
  git clone https://github.com/cma2015/CAFU.git
  ```
#### Upload regular file
- Click **Get Data** in the homepage (see figure below) of CAFU to upload files.


  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/2.png)


  And then you will see the following interface:

  ![](https://github.com/cma2015/CAFU/blob/master/CAFU_images/3.png)


  Next, Click the button **Choose local file** and select a file you would like to upload (e.g. upload the file ```mapping_info``` in the directory ```/your directory/CAFU/test_data/SE RNA-Seq/```), you will see the following interface:

  
  ![Upload regular file](https://github.com/cma2015/CAFU/blob/master/CAFU_images/Picture7.png)

  
  Then click **Start** to upload file.

#### Upload a collection
- Similar to **Upload regular file**, click **Get Data** first (see figure below):


  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/2.png)

  
  And then you will see the following interface:
  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/5.png)

  
  Afterwards, select a list files to upload as a collection (e.g. upload all files with ZIP suffix in the folder ```/your directory/CAFU/test_data/SE RNA-Seq/```):
  
  **Note:** A collcetion also permits one file.

  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/6.png)


  Click **Start** to upload, after finishing uploading, click **Build** (see figure below):

  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/7.png)

  
  Then enter a name for your collection and click **Create list** to finish.

  
  ![Upload files](https://github.com/cma2015/CAFU/blob/master/CAFU_images/8.png)


### UNMAPPED READ EXTRACTION
In this module, we provide an example for each module to show how to use to perform unmapped reads extraction.
- **Quality control**

  In this function, we implemented FastQC (Andrews *et al*., 2010) to enable users to perform quality control. In this tutorial, we will use a list of single-end RNA-Seq collection file (located in the folder ```/your directory/CAFU/test_data/SE RNA-Seq/```) to perform quality control.
  
  To run this module correctly, upload the files in the folder ```/your directory/CAFU/test_data/SE RNA-Seq/``` with **ZIP** suffix as a collection named as ```SE-RNA-Seq```（see section **Upload a collection** to see how to upload a list of RNA-Seq datasets as a collection).

  
  ![Quality control](https://github.com/cma2015/CAFU/blob/master/CAFU_images/9.png)

  
  Finally, click **Execute** to start performing quality control.
  Once the analysis is finished, a basic text (RawData) and a HTML output file (Webpage) will be returned. The HTML output contains:

  - Basic Statistics
  - Per base sequence quality
  - Per sequence quality scores
  - Per base sequence content
  - Per base GC content
  - Per sequence GC content
  - Per base N content
  - Sequence Length Distribution
  - Sequence Duplication Levels
  - Overrepresented sequences
  - Kmer Content

  You can access the output by click results in the right history bar (see figure below).


  ![Quality control](https://github.com/cma2015/CAFU/blob/master/CAFU_images/10.png)


- **Trim raw reads**

  In this function, poly-A/T is firstly trimmmed using fqtrim (Pertea, 2015), and then high-quality reads (e.g. score > 20) are retained by Trimmomatic (Bolger et al., 2014). The reads less than 20bp (in default) are discarded. 
  
  Here, we used the same collection (```SE-RNA-Seq```) with last step (**Quality control**) to trim raw reads.


  ![Quality control](https://github.com/cma2015/CAFU/blob/master/CAFU_images/11.png)


  For each FASTQ formatted file, CAFU returns two outputs including:

  - ```fqtrim_DATA_ID.fastq```: poly-A/T trimmed RNA-Seq using fqtrim.

  - ```Trimmed_DATA_ID.fastq```: High-quality RNA-Seq generated by Trimmomatic.

- **Extract Unmapped Reads**

  This function integrates several bioinformatic softwares including HISAT2 (Kim et al., 2015), SAMTools (Li et al., 2009), and BEDTools (Quinlan et al., 2010) to extract unmapped reads. Firstly, HISAT2 is used to align user-defined high-quality reads to reference genome. Secondly, SAMTools is used to extract unmapped alignments in SAM format from alignments generated by HISAT2 with parameter "-f 4". Finally, BEDTools is used to convert the unmapped alignments in SAM format to FASTQ format. In addition, CAFU also supports dual RNA-Seq analysis by aligning RNA-Seq data to multiple reference genomes.

  To run this function, three inputs (see figure below) are required. Here, we still use single-end RNA-Seq files (upload these files as a collection named as ```SE-RNA-Seq```) in directory ```/your directory/CAFU/test_data/SE RNA-Seq/``` to show how to use this function.
 
  **Inputs description:**
  
  **Input 1:** A **collection** of reference genome. Upload files (```maize.fa.zip```) in directory ```/your directory/CAFU/test_data/genomes/``` as a collection named as ```Ref_Genome``` (user-defined name). See section **Upload collection file** to learn how to upload a collection.
  
  **Input 2:** A **collection** of trimmed single-end RNA-Seq files (```Trimmed_DATA_ID.fastq```) generated from **Trim Raw Reads**.
  
  **Input 3:** A **regular** file containing mapping-operation information. Upload files in directory ```/your directory/CAFU/test_data/SE RNA-Seq/mapping_info``` as a regular file.
  
  **Note:** **Input 3** is a semicolon seperated matrix which contains two columns. The first column contains RNA-Seq ID in each experiment. The second column is the corresponding reference genome ID for each experiment.
  ```bash
  SRR2144382;maize
  SRR2144383;maize
  SRR2144384;maize
  SRR2144385;maize
  SRR2144410;maize
  SRR2144411;maize
  SRR2144442;maize
  SRR2144443;maize
  ```
 
  ![Extrat unmapped reads](https://github.com/cma2015/CAFU/blob/master/CAFU_images/Picture8.png)


  Then click **Execute** to start extracting unmapped reads using CAFU. After finishing this process, the final unmapped reads named as ```all_unmapped_reads.fastq``` and the number of unmapped reads per sample named as ```Status_of_number_of_unmapped_reads``` (see figure below) will be returned.


  ![unmapped reads output](https://github.com/cma2015/CAFU/blob/master/CAFU_images/13.png)


- **Remove Contamination**

  Considering that sequences obtained from impure nucleic acid preparations may contain DNA from sources other than the sample. Those sequence contaminations are a serious concern to the quality of the data used for downstream analysis. So in this function, potential contamination sequences are removed using Deconseq (Schmieder *et al*., 2011) with user-defined coverage and identity (e.g., 0.95) by aligning input reads to a contamination database. In the current version of CAFU, 3,529 bacterial and 81 viral reference genomes (downloaded from NCBI on 2018/11/05) are provided as default, however, user-defined contamination database is also supported. 

  Here we use the ```Default contamination and single-end unmapped reads``` as an example (see figure below):

  ![remove contamination](./review-2/Picture1.png)


  Then the clean reads with FASTQ format will be returned.


###  *DE NOVO* TRANSCRIPT ASSEMBLY OF UNMAPPED READS
- **Assemble Unmapped Reads**

  In this function, three steps including ***de novo* assembly of unmapped reads**, **remove redundancy of transcript fragments**, and **re-assembly transcript fragments** will be sequentially performed to assemble unmapped reads. In this tutorial, we will use the unmapped reads generated by module ```Extract Unmapped Reads``` or ```Remove Contamination``` as the input (see figure below).

  
  ![unmapped reads output](https://github.com/cma2015/CAFU/blob/master/CAFU_images/14.png)
  

  Then assembled transcripts named as ```Unmapped_reads_de_novo_assembled_transcripts```  will be returned (see figure below).


  ![assembled transcripts](https://github.com/cma2015/CAFU/blob/master/CAFU_images/15.png)


### EVIDENCE SUPPORT OF ASSEMBLED TRANSCRIPTS
- **Expression-level Evidence**

  This function allows users to eliminate assembled transcripts with low read coverage and/or low expression abundance, which are likely assembly artifacts. RNA-Seq reads used for **Assemble Unmapped Reads** are mapped with newly assembled transcripts or/and reference transcripts by using bowtie2 (Langmead et al., 2012). CAFU outputs the read coverage of assembled transcripts at single-base resolution using BEDTools (Quinlan et al., 2010), and estimates the expression abundance of all transcripts in terms of FPKM (Fragments Per Kilobase Million) using RSEM (Li et al., 2011). Assembled transcripts with low read coverage (e.g., less than 10) and/or low expression (e.g., FPKM less than 1) in the majority of samples (e.g., 80%) are discarded.

  NOTE: RNA-Seq for calculating expression abundance and read coverage of transcripts are used the data used in *de novo* assembly. Thus, users only require to input the newly assembled transcripts from unmapped reads (generated from the function **Assemble Unmapped Reads**) or/and reference transcripts with FASTA format. 

  ![expression-level](https://github.com/cma2015/CAFU/blob/master/CAFU_images/18.png)
  
  
  Then three files will be returned:

  
  **Output 1**: ```Read coverage of each transcript in each sample```, a matrix whose rows represent transcripts, and columns represent read coverage.

  **Output 2**: ```Expression abundance of each transcript in each sample```, a matrix whose rows represent transcripts, and columns represent expression abundance.

  **Output 3**: ```Confident transcript ID by Expression-level Evidence```, a confident transcript ID filtered by both read coverage results and expression abundance.


- **Genome-level Evidence**
  
  This function can be used to identify *de novo*-assembled transcripts missing from the existing genome annotation. All the assembled transcripts are aligned to the reference genome sequences of the species of interests using GMAP (Wu *et al*., 2005), and the best genomic matches with high identity (e.g., ≥ 95%) and coverage (e.g., ≥ 95%) are selected. Users can also eliminate assembled transcripts with no introns, which could represent either noise or pseudogenes.
  
  To run this function, at least three inputs are required including:

  **Input 1**: Reference genome sequences of the species of interest (input as a data collection). Upload files (```maize.fa.zip, Oryza.fa.zip, Sorghum.fa.zip```) in directory ```/your directory/CAFU/test_data/genomes/``` as a collection named as ```Ref_Genome``` (user-defined name).

  **Input 2**: Assembled transcript sequences generated from function **Assemble Unmapped Reads**. Test data (```assembled_transcript.fasta```) is in directory ```/your directory/CAFU/test_data/others/```.

  **Input 3**: Genome annotation (GFF/GTF format) of corresponding species with RNA-Seq samples. Upload files (```maize.gff3.zip```) in directory ```/your directory/CAFU/test_data/genomes/```.

  **Input 4**: A character indicating the reference genome name (e.g. maize).


  ![genome-level](https://github.com/cma2015/CAFU/blob/master/CAFU_images/19.png)


  Then four outputs will be returned:

  **Output 1**: ```Integrated GMAP results of newly assembled transcripts against all reference genome sequences```: GMAP alignment results (coverage and identity) of each assembled transcript against all reference genome sequences.

  **Output 2**: ```Confitent transcript information```: Confident transcript information filtered by high coverage and identity.

  **Output 3**: ```The same/similar-intron transcript ID```: Assembled transcript ID which possess the same/similar intron with corresponding species reference transcripts.

  **Output 4**: ```Novel transcript ID```: Novel transcripts missing in the existing genome annotation.

- **Transcript-level Evidence**

  This function can be used to select assembled transcripts with high similarity comparing to other well-annotated transcripts, such as full-length transcripts generated from single-molecule real-time sequencing and/or high-quality transcripts annotated in closely related species. After aligning assembled transcripts with other well-annotated transcripts with GAMP (Wu *et al*., 2005), CAFU outputs the best transcript alignments with high identity (e.g., ≥ 95%) and coverage (e.g., ≥ 95%).
  
  To run this function, at least two inputs are required including:

  **Input 1**: Reference sequence of well-annotated transcripts, such as full-length transcripts generated from single-molecule real-time sequencing and/or high-quality transcripts annotated in closely related species. Upload files (```ref_trans.fasta.zip, ref_trans_1.fasta.zip```) in directory ```/your directory/CAFU/test_data/Transcripts/``` as a collection named as ```Well-annotated transcript sequences``` (user-defined name).

  **Input 2**: Assembled transcript sequences generated from the function **Assemble Unmapped Reads**. Test data (```assembled_transcript.fasta```) is in directory ```/your directory/CAFU/test_data/others/```.


   ![genome-level](https://github.com/cma2015/CAFU/blob/master/CAFU_images/22.png) 

   
   Then three outputs will be returned:

   **Output 1**: ```Integrated GMAP results of newly assembled transcripts against all reference transcript sequences```, GMAP alignment results (coverage and identity) of each assembled transcript against all other well-annotated transcripts.

   **Output 2**: ```Confident transcript information```, Confident transcript information filtered by high coverage and identity.

   **Output 3**: ```Confident transcript ID```, IDs of transcripts that prossess high similarity comparing to other well-annotated transcripts.

- **Protein-level Evidence**

  In this function, coding potential evidence of transcripts is fistly evaluated using CPC2 (Kang *et al*., 2017). Then for coding transcripts, Pfam (Finn *et al*., 2014) will be used to identify the protein families. 
  
  Here, we use the file ```/your directory/CAFU/test_data/others/assembled_transcript.fasta``` to execute this function.

  Then three outputs will be returned:

  **Output 1**: CPC2 output. A tab seperated matrix contains seven columns. Each column shows the sequence ID, putative peptide length, Fickett score, isoelectric point, the integrity of the orf, coding probability and the coding/noncoding classification label. More details about this output can be seen from [CPC2 official website](http://cpc2.cbi.pku.edu.cn/help.php). 

  **Output 2**: Confident transcript ID. 

  **Output 3**: A tab seperated matrix contains transcript ID, alignment start, alignment end, envelope start, envelope end, Hmm access, Hmm name, Type of domain, Hmm start, Hmm end, Hmm length, Bit score, E-value, Significance, Clan, etc.


### SPECIES ASSIGNMENT OF ASSEMBLED TRANSCRIPTS
- **Species Assignment of Transcripts**

  SAT (Species Assignment of Transcripts) is a machine learning-based toolkit used for species assignment of transcritps assembled using unmapped reads from mixed-species (eg., pathogen-host) samples. 
  
  In this tutorial, we will show how to use SAT using all the FASTA format files in directory ```/your directory/CAFU/test_data/SAT/```. To run SAT, at least three files are required including:

  **Input 1**: The positive coding sequences (CDS) with FASTA format.  

  **Input 2**: The negative coding sequences (CDS) with FASTA format.

  **Input 3**: The coding sequences (CDS) with FASTA format assembled using unmapped reads from mixed-species (e.g., pathogen-host) 

  ![assembled transcripts](https://github.com/cma2015/CAFU/blob/master/CAFU_images/16.png)

  Then click **Execute** to run this function, then four outputs will be returned:

  **Output 1**: ```Probabilistic score of each transcript```, A probabilistic score of each transcript, the higher, the more probable being positive.

  **Output 2**: ```Eight commonly used measures under specidied threshold```, A bar plot evaluating the measures (Sn, Sp, Pr, Acc, MCC, Fscore, AUC and AUPR) under specified threshold. 

  **Output 3**: ```The presicion recall curves in k-fold cross validation```, The PR curve of k-fold cross-validation. 

  **Output 4**: ```The receiver operating curves in k-fold cross validation```, The ROC curve of k-fold cross-validation. 

### SEQUENCE CHARACTERIZATION OF ASSEMBLED TRANSCRIPTS

- **Characterize Nucleic-acid Feature**

  This function allows users to explore the nucleic-acid similarity between assembled and reference transcripts in terms of the distribution of transcript length and G+C content.

  In this function, two features can be analyzed (see figure below). To run this function, two inputs are required including:

  **Input 1**: The sequences of assembled transcripts (or confident/novel transcripts) derived from unmapped reads. Test data (```assembled_transcript.fasta```) is in directory ```/your directory/CAFU/test_data/others/```.
  
  **Input 2**: The sequences of transcripts from the existing genome annotation. Test data (```ref_trans.fasta.zip```) is in directory ```/your directory/CAFU/test_data/Transcripts/```.

  ![nucleic-acid feature](https://github.com/cma2015/CAFU/blob/master/CAFU_images/21.png)

  Then three outputs will be returned:

  - **For Transcript length**

    **Assembled transcript length**: Length of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    **Reference transcript length**: Length of transcripts from the existing genome annotation.

    **Length distribution comparison**: Length distribution comparison between assembled transcripts and reference transcript.

  - **For GC content**


    **Assembled transcript GC content**: GC content of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    **Reference transcript GC content**: GC content of transcripts from the existing genome annotation.

    **Length distribution comparison**: GC content distribution comparison between assembled transcripts and reference transcript.


- **Characterize Amino-acid Feature**

  This function allows users to explore the amino-acid features similarity used in SAT between assembled and reference transcripts.

  Here, we take an example (see figure below) to show how to use this function to compare k-mer frequency of assembled and reference transcripts. 

  ![nucleic-acid feature](https://github.com/cma2015/CAFU/blob/master/CAFU_images/23.png)
  
  The outputs contain:
  
  **Note**: The detailed feature descriptions are available at http://bioinformatics.hitsz.edu.cn/BioSeq-Analysis/
  
  - **For Kmer**

    ```Assembled transcript K-mer (k = 1)```: K-mer frequency of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript K-mer (k = 1)```: K-mer frequency of transcripts from the existing genome annotation.

    ```Assembled transcript K-mer (k = 2)```: K-mer frequency of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.
    ```Reference transcript K-mer (k = 2)```: K-mer frequency of transcripts from the existing genome annotation.
  - **For DR**

    ```Assembled transcript DR```: Distance-based residues encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript DR```: Distance-based residues encoding of transcripts from the existing genome annotation.
  - **For AC**

    ```Assembled transcript AC```: Auto-covariance encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript AC```: Auto-covariance encoding of transcripts from the existing genome annotation.
  - **For CC**

    ```Assembled transcript CC```: Cross-covariance encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads
    
    ```Reference transcript CC```: Cross-covariance encoding of transcripts from the existing genome annotation.
  - **For ACC**

    ```Assembled transcript ACC```: Auto-cross-covariance encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript ACC```: Auto-cross-covariance encoding of transcripts from the existing genome annotation.
  - **PDT**

    ```Assembled transcript PDT```: Kmer of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript PDT```: Kmer of transcripts from the existing genome annotation.

  - **For PC-PseAAC**

    ```Assembled transcript PC-PseAAC```: Physicochemical distance transformation encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript PC-PseAAC```: Physicochemical distance transformation encoding of transcripts from the existing genome annotation.
  - **For SC-PseAAC**

    ```Assembled transcript SC-PseAAC```: General series correlation pseudo amino acid composition encoding of assembled transcripts (or confident/novel transcripts) derived from unmapped reads.

    ```Reference transcript SC-PseAAC```: General series correlation pseudo amino acid composition encoding of transcripts from the existing genome annotation.

- **Detect Alternative Splicing Events**

  This function is used to detect alternative splicing events in assembled transcripts using an R package SGSeq (Goldstein *et al*., 2016).

  The only input of this function is GFF/GTF annotation file, here we use the test data located in ```/your directory/CAFU/test_data/others/AS_test.gtf.zip``` to run this function.

  ![nucleic-acid feature](https://github.com/cma2015/CAFU/blob/master/CAFU_images/28.png)

  Then two outputs will be returned:

  **Output 1**: ```Alternative splicing results```, Alternative splicing events of all transcripts.

  **Output 2**: ```Read counts of all genome features in GTF file```, Read counts supported for all genome features in GTF file provided by users. 

### EXPRESSION PROFILES OF ASSEMBLED TRANSCRIPTS
- **Analyze Condition Specificity**

  This function identifies a set of transcripts highly expressed under different conditions. The condition specificity of a transcript for condition type T is defined using the formula described in (Ma et al., 2014).

  Here, we use the test data ```assembled_transcript_expression, RNA-Seq_sample_information``` in the folder ```/your directory/CAFU/test_data/others/``` to show its u