---
title: "• Sequence repositories"
subtitle: "Wonders of the scientific world!"

date: 2022-11-03T00:00:00+01:00

fontawesome: true
linkToMarkdown: true

toc:
  enable: true
  keepStatic: false
  auto: false
code:
  copy: true
math:
  enable: true
share:
  enable: false
comment:
  enable: false
---

## Sequence repositories
While basic methods of sequencing have existed since the 1970s, it is in the 21st century that the field exploded. This was due to key technological advancements, collectively known as "high-throughput" or "next-generation" sequencing. Generating high-quality sequences, at a fraction of the price, made genomics accessible to even modest research groups. 

### The genetic revolution
To illustrate this, [sequencing a human genome in 2000 required many specialist research groups, several years, and $100 million dollars](https://www.genome.gov/about-genomics/fact-sheets/Sequencing-Human-Genome-cost); by 2010 this had dropped to a matter of weeks and $100,000. In 2022 a human genome (or at least a full exome) can be sequenced for a few hundred dollars, in a matter of days - as such, aggregated bioinformatic data roughly follow [Moore's Law](https://www.wikiwand.com/en/Moore's_law), doubling every 18 months.

Given the vast amount of data being produced, international consortia were established, with the goal of making any and all bioinformatic research data available to the world, for free (well, ultimately funded by the tax-payer, as with all public research). Most journals will require that sequence information (produced in the research process leading to publication) is deposited with one of the major repositories. These repos should be seen as wonders of the modern scientific world! The openly-accessible data they contain drive forward research in medicine, biology, mathematics, information theory, even nanotechnology and biocomputing.

{{< admonition type="info" title="Some big repositories" open=false >}}
- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/): this is the USA-based national repo administered by the National Centre for Biotechnology Information, containing a wide-variety of databases. The biggest of these are "genome", "gene", and "protein", as well as raw (pre-assembly) reads, for example in "SRA", the sequence read archive. As of 2022, it contains around 2.5 billion sequences, of around 17 trillion nucleotides.
- [Ensembl](https://www.ensembl.org/index.html): the repository of the European Bioinformatics Institute. It is accompanied by a taxonomy-focused repository, [EnsembleGenomes](https://ensemblgenomes.org/), which is particularly useful for microbial research.
- [ExPASy](https://www.expasy.org/): this is the archive for the Swiss Bioinformatics Institute. Originally is was focussed on proteomics and protein sequences, and today is a major hub for genomics focussed on cell-lines, gene function, drug-design and pathogen classification.
- And there are more, which we encourage you to look for.
{{< /admonition >}}

## Acquiring bioinformatic data
Using a browser, we will explore a repository, and get hold of some data to work with in BioPython. A key point here is that, as a starting point for answering a research question, we should become familiar with searching "as a human", before automating such searches.

Let's visit [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). It is worth getting familiar with at least the major sections of GenBank.

### Making your own FASTA files
We are going to search for a well-characterised, ultra-conserved gene, that we are going to use for some common bioinformatics tasks. By the dropdown, select "Protein" - this is the name of the database we will be searching. Note the others in the menu - some are major DBs like _gene_ and _genome_, others are more niche. In the search box, enter `signal recognition particle 72`. This is a ubiquitous gene in eukaryotes, and arguably in prokaryotes too.

The search results will be shown - above the results, change the results per page from `20 per page` to `200 per page`. Have a look down the results, and see what details you can recognise. Open one by clicking the link, and have a look over the details that typically accompany an accession. 

### Exercise 1
{{< admonition type="note" title="Exercise 1" open=true >}}
Sequences come in a variety of formats. Look at the following data, and identify what format it is in
- Format 1
```
>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGT
```
- Format 2
```
5 36
Molossus_molossus          MASGGGVAVASLWTEVNRCGQSGDYTRALKALTKI
Rousettus_aegyptiacus      MSDVAKWEKQLERLLAEGGESAKVIAAIDKILSAS
Pipistrellus_kuhlii        MKILQKAEDLCRHSLSEDSDGTEEDPQAELAIIHG
Rhinolophus_ferrumequinum  MASGGSGGVSVPALWSEVNRYGQNGDFTRALKTVN
Phyllostomus_discolor      MKILQKAEDLCRHSLSEDSDGTEEDPQAELAIIHG
```
- Format 3
```
@MN00537:51:000H2K25G:1:11101:2213:1092 1:N:0:9
CTCCAGTCCTTACTCCCATATCTAACCTCTTACCCCTACNTCATAGGTANACATTTTA
+
ZTTHJFFFFAAAABBBABBACCCAAAADFFFFFABDDDFFAAHHJIILPPPLUUUVZZ
```
- Format 4
```
LOCUS       AB000263             368 bp    mRNA    linear   PRI 05-FEB-1999
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   AB000263
ORIGIN      
        1 acaagatgcc attgtccccc ggcctcctgc tgctgctgct ctccggggcc acggccaccg
       61 ctgccctgcc cctggagggt ggccccaccg gccgagacag cgagcatatg caggaagcgg
```
- Format 5
```
ID   AB000263 standard; RNA; PRI; 368 BP.
XX
AC   AB000263;
XX
DE   Homo sapiens mRNA for prepro cortistatin like peptide, complete cds.
XX
SQ   Sequence 368 BP;
     acaagatgcc attgtccccc ggcctcctgc tgctgctgct ctccggggcc acggccaccg        60
     ctgccctgcc cctggagggt ggccccaccg gccgagacag cgagcatatg caggaagcgg       120
```
- Format 6
```
gi|6273285|gb|AF191659.1|AF191      TATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273284|gb|AF191658.1|AF191      TATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273290|gb|AF191664.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273289|gb|AF191663.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273291|gb|AF191665.1|AF191      TATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
                                    ******          ********  **** ********* *********
```
{{< /admonition >}}

{{< admonition type="tip" title="Solutions" open=false >}}
- Format 1: [FASTA](https://www.wikiwand.com/en/FASTA_format)
- Format 2: [PHYLIP](https://www.wikiwand.com/en/PHYLIP)
- Format 3: [FASTQ](https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html)
- Format 4: [GenBank](https://www.genomatix.de/online_help/help/sequence_formats.html)
- Format 5: [EMBL](https://www.genomatix.de/online_help/help/sequence_formats.html)
- Format 6: [CLUSTAL](https://bioperl.org/formats/alignment_formats/ClustalW_multiple_alignment_format)
{{< /admonition >}}

### Exercise 2
{{< admonition type="note" title="Exercise 2" open=true >}}
**1.** Go back to your results list. We are going to make a `FASTA` file of one of these sequences.
- return to your list of results, searching the `Protein` database for `signal recognition particle 72`
- scroll down until you find a species you like the sound of
- make sure your record has an amino acid length of **between 600 and 690 AA**
- avoid records that say, `isoform`, but don't worry if it says `partial`
- click on `FASTA` to see the sequence for the record
- select and copy _just the FASTA data_, ie, starting with the `>` character, and ending with the last AA position
- go to your Jupyter Lab session, create a new file, and paste this data in
- save this file `srp_single.fas`

**2.** We also want to create a `multiple FASTA` file - which is just a `FASTA` file with more than one sequence in it.
- return to your list of results and if you haven't already, ask it to show `200 per page` 
- have a browse down this list, and tick the select box on **ten** records:
  - make sure each of these ten are from a different species (eg, don't get more than one from _Homo sapiens_)
  - as above, make sure each one has a record with an amino acid length of **between 600 and 690 AA** and avoids `isoform`s
- Once you have selected ten that you like, scroll to the very bottom (or the very top) of the records, and click on the `Summary`. From the menu it shows, choose `FASTA (text)`. This will open a new page, with the protein sequences for the gene SRP72, for the species you chose. Select all of this text (`ctrl-a`), and copy it. Now, go to your Jupyter Lab tab, and create a new text file. Paste your FASTA sequences into the new file, and save it as `srp_multi.fas`
{{< /admonition >}}


