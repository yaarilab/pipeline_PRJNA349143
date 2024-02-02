
A pipeline for BCR repertoire libraries from the form of - UMI Barcoded Illumina MiSeq 325+275 paired-end 5’RACE BCR mRNA. that were produced in the same fashion as those in [VanderHeiden et al. 2017](https://journals.aai.org/jimmunol/article/198/4/1460/109668/Dysregulation-of-B-Cell-Repertoire-Formation-in)

Library preperation and sequencing method:

The sequences were amplified specific primers  1. AbSeq R1 Human IG Primers.fasta , 2. AbSeq R2 TS.fasta 3. AbSeq Human IG InternalCRegion.fasta.
The generated libraries were then sequenced with Illumins MiSeq 325+275. 



Each read was sequenced from one end of the target cDNA, so that the togather the two reads cover the entire the full variable region of the Ig heavy chain. The sequencing used an uneven number of cycles for the two reads, with a 2x300 kit. The V(D)J reading frame progressed from the start of read 2 to the start of read 1. Read 1 was in the reverse complement orientation, had a partial C-region, and was 325 nucleotides long. Read 2 was 275 nucleotides long, containing the 5’RACE template switch site with a 17-nucleotide UMI barcode preceding it.


Input files:

1. Two fastq file of paired-end sequencing
2. primer files
3. Assemble pairs reference file

To test the pipeline:

We recommend to test the pipeline using a small example from the original reads that can download using the fastq-dump command.

```bash
fastq-dump --split-files -X 25000 SRR4026043
```

And upload directly to dolphinnext. 


Output files:

1. {sampleName}_collapse-unique.fastq
2. {sampleName}_atleast-2.fastq
3. log tab file for each steps
4. report for some of the steps


Pipeline container:

* Docker: immcantation/suite:4.3.0


Sequence processing steps:

* Quality control, UMI annotation and primer masking
	1. FilterSeq quality
	2. MaskPrimer score

	
* Generation of UMI consensus sequences
	3. PairSeq
	4. BuildConsensus


* Paired-end assembly of UMI consensus sequences
	5. PairSeq	
	6. AssemblePairs sequential 

	
* Deduplication and filtering
	7. MaskPrimer align
	8. ParseHeaders collapse
	9. CollapseSeq
	10. SplitSeq group


Primers used:


* [AbSeq R1 Human IG Primers](https://bitbucket.org/kleinstein/presto/src/master/examples/VanderHeiden2017/AbSeq_R1_Human_IG_Primers.fasta)

* [AbSeq R2 TS](https://bitbucket.org/kleinstein/presto/src/master/examples/VanderHeiden2017/AbSeq_R2_TS.fasta)

* [AbSeq Human IG InternalCRegion](https://bitbucket.org/kleinstein/presto/src/master/examples/VanderHeiden2017/AbSeq_Human_IG_InternalCRegion.fasta)

Reference used:

* [IMGT_Human_IG_V](https://bitbucket.org/kleinstein/presto/src/master/examples/VanderHeiden2017/IMGT_Human_IG_V.fasta)