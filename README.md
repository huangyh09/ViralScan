# STARsolo-ViralScan

This repository summarise a pipeline to efficiently scan viral transcripts in 
single-cell RNA-seq data with 
[STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)

Thanks to the flexible support from STARsolo, the pipeline has two main 
features:

1. Default input: the unmapped reads in the bam file from CellRanger or STARsolo
   by mapping reads to the host (e.g., human or mouse)
2. Default output: UMI counts per virus per cell


## Step 1a: Construct viral genomes (provided)

One can generate any customised list of viral genomes. 
Here, we compiled 833 viral sequences, including 762 viruses from 
[VirTect](https://github.com/WGLab/VirTect),
7 other viruses and 64 consensus sequences of human endogenous retroviruses 
(HERVs) compiled by 
[Vargiu et al](https://retrovirology.biomedcentral.com/articles/10.1186/s12977-015-0232-y)

The reference viral genomes are stored in the 
[viral_reference/viruses_833.fasta](./viral_reference/viruses_833.fasta)

Note, in order to avoid false positives, we only used the CSD for the following
viruses, as their remaining regions have a poly-A or poly-T region that may 
give false positives:
* virus169: NC_004102.1_Hepatitis_C_virus_genotype_1,_complete_genome 
* virus170: NC_009823.1_Hepatitis_C_virus_genotype_2,_complete_genome
* virus174: NC_009827.1_Hepatitis_C_virus_genotype_6,_complete_genome
* virus194: NC_022518.1_Human_endogenous_retrovirus_K113_complete_genome
* virus515: NC_001672.1_Tick-borne_encephalitis_virus,_complete_genome


## Step 1b: Make gene annotation (provided)

In order to use STARsolo, we have to obtain a gene annotation file in GTF 
format. Many databases, e.g., 
[NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore), have good annotations 
of common virues. Here, as focusinig on scan viral existence, we introduce a 
pseudo definition of genes by using the whole genome as a gene, hence only 
relying on the input viral genomes.

You can use the [Viral_GTF_maker.py](./Viral_GTF_maker.py) as follows,

```bat
python Viral_GTF_maker.py -f viral_reference/viruses_833.fasta -o YOUR_OUTPUT_DIRECTORY
```

Then you will see three output files in the output folder:
* viruses_833.numbered.fasta: same sequences with number virus as simplier id
* viruses_833.numbered.gtf: the gene annotation file
* viruses_833.id_notes.tsv: the notes for viruses


## Step 1c: Build STAR reference

With the gene annotation at hand, now we can build the reference for running 
STARsolo.

For 10x Genomics data, you can use the following settings, but you can also 
refer to 
[STAR's documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)


```bat
STAR --runThreadN 3 --runMode genomeGenerate \
    --genomeDir YOUR_Ref_DIR --genomeFastaFiles viruses_833.numbered.fasta \
    --sjdbGTFfile viruses_833.numbered.gtf \
    --sjdbOverhang 97 --genomeSAindexNbases 8
```


## Step 2: Align unmapped reads
Usually, the unmapped reads are kept in the bam file when aligning it to the 
host (e.g., human or mouse), for example, CellRanger for 10x Genomics data.
You can obtain the unmapped reads by using flag 4 (``-f 4``) in the bam file; 
see more discussions 
[here](https://kb.10xgenomics.com/hc/en-us/articles/360004689632).

Once again, STAR has nice support to take bam as input and additional command 
line, so you can re-align the unmapped reads to viral genomes by a single 
command line, for example,


```bat
STAR --runThreadN 20 --genomeDir YOUR_Ref_DIR \
    --soloType Droplet --soloCBwhitelist YOUR_cell_list \
    --readFilesIn YOUR_CellRanger_outs_Bam_file --readFilesType SAM SE \
    --readFilesCommand samtools view -f 4 \
    --soloInputSAMattrBarcodeSeq CR UR \
    --soloInputSAMattrBarcodeQual CY UY \
    --outSAMtype BAM Unsorted 
```

**Note 1**, here you can use the whitelist provided by cellranger (see 
[STARsolo's note](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)
). However, you could also directly use the called cells from cellranger, e.g.,
in the filtered_feature_bc_matrix folder. The command line can do the trick:

```bat
zcat filtered_feature_bc_matrix/barcodes.tsv.gz | awk '{print substr($0, 1, length($0)-2)}' > YOUR_cell_list
```

**Note 2**, You may want to quickly view the total reads mapped to each virus. 
You can use the versatile samtools by its ``idxstat`` (you may need to sort the 
bam file first).

```bat
samtools sort Aligned.out.bam -o Aligned.sortedByCoord.out.bam
samtools index Aligned.sortedByCoord.out.bam

samtools idxstat Aligned.sortedByCoord.out.bam
```

