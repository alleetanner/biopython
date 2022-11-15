---
title: "• Sequence repositories"
subtitle: "Python as a sequence search tool"

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
To illustrate this, [sequencing a human genome in 2000 required many specialist research groups, several years, and $100 million dollars](https://www.genome.gov/about-genomics/fact-sheets/Sequencing-Human-Genome-cost); by 2010 this had dropped to a matter of weeks and $100,000. In 2022 a a human genome (or at least a full exome) can be sequenced for a few hundred dollars, in a matter of days - as such, aggregated bioinformatic data follows a parallel of [Moore's Law](https://www.wikiwand.com/en/Moore's_law), roughly doubling every 18 months.

Given the vast amount of data being produced, international consortia were established, with the goal of making any and all bioinformatic research data available to the world, for free (well, ultimately funded by the tax-payer, as with all public research). Most journals will require that sequence information (produced in the research process leading to publication) is deposited with one of the major repositories. These repos should be seen as wonders of the modern scientific world! The openly-accessible data they contain drive forward research in medicine, biology, mathematics, information theory, even nanotechnology and biocomputing.

{{< admonition type="info" title="Some big repositories" open=false >}}
- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/): this is the USA-based national repo administered by the National Centre for Biotechnology Information, containing a wide-variety of databases. The biggest of these are "genome", "gene", and "protein", as well as raw (pre-assembly) reads, for example in "SRA", the sequence read archive. As of 2022, it contains around 2.5 billion sequences, of around 17 trillion nucleotides.
- [Ensembl](https://www.ensembl.org/index.html): the repository of the European Bioinformatics Institute. It is accompanied by a taxonomy-focused repository, [EnsembleGenomes](https://ensemblgenomes.org/), which is particularly useful for microbial research.
- [ExPASy](https://www.expasy.org/): this is the archive for the Swiss Bioinformatics Institute. Originally is was focussed on proteomics and protein sequences, and today is a major hub for genomics focussed on cell-lines, gene function, drug-design and pathogen classification.
- And there are more, which we encourage you to look for.
{{< /admonition >}}

## Using repositories manually
Before we use BioPython to carry out a query, we will get familiar with a sequence archive using a browser. The key point here is that, as a starting point for answering a research question, we should become familiar with searching "as a human", before automating such searches.

Let's visit [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). It is worth getting familiar with at least the major sections of GenBank.

### Making your own FASTA files
We are going to search for a well-characterised, ultra-conserved gene, that we are going to use for some common bioinformatics tasks. By the dropdown, select "Protein" - this is the name of the database we will be searching. Note the others in the menu - some are major DBs like gene and genome, others are more niche. In the search box, enter `signal recognition particle 72`. This is a ubiquitous gene in eukaryotes, and arguably in prokaryotes too.

The search results will be shown - above the results, change the results per page from `20 per page` to `100 per page`. Have a look down the results, and see what details you can recognise. Open one by clicking the link, and have a look over the details that typically accompany an accession. 

Go back to your results list. Now, we are going to select a few of these, and get the protein sequences as a FASTA file. Select by ticking the box next to a record:
- choose ten records
- make sure each of then ten are from a different species (eg, don't get more than one from _Homo sapiens_)
- make sure each one has a record with an amino acid length of between 600 and 690 aa,
- don't worry if the record says `partial`, as long as it meets the above requirements.
Once you have selected 10 (or so) that you like, scroll to the very bottom (or the very top) of the records, and click on the "Summary". From the menu it shows, choose "FASTA (text)". This will open a new page, with the protein sequences for the gene SRP72, for the species you chose. Select all of this text (ctrl-a), and cut it. Now, go to your Jupyter Lab tab, and create a new text file. Paste your FASTA sequences into the new file, and save it (ideally in the folder you are working in today) as `srp72_multi.fas`.

For the exercises in the following section, we want to make one more FASTA file. This will be similar to `srp72_multi.fas` that you just made, but contain a single sequence. Go back to NCBI, and your list of SRP72 matches. Choose one more record (again, meeting the above criteria), from a species you have not already chosen. Click through and look at the main record, then find the FASTA for that gene. As before, copy and paste it into a new file in Jupyter, with the name `srp72_single.fas`.

### Exercise
{{< admonition type="Exercise" title="Sequence formats" open=true >}}
Sequences come in a variety of formats. Look at the following data, and identify what format it is in
- Format 1
```
>gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
LCLYTHIGRNIYYGSYLYSETWNTGIMLLLITMATAFMGYVLPWGQMSFWGATVITNLFSAIPYIGT
```
- Format 2
```
5 34
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
{{< /admonition >}}

{{< admonition type="warning" title="Solutions" open=false >}}
- Format 1: [FASTA](https://www.wikiwand.com/en/FASTA_format)
- Format 2: [PHYLIP](https://www.wikiwand.com/en/PHYLIP)
- Format 3: [FASTQ](https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html)
- Format 4: [GenBank](https://www.genomatix.de/online_help/help/sequence_formats.html)
- Format 5: [EMBL](https://www.genomatix.de/online_help/help/sequence_formats.html)
{{< /admonition >}}
