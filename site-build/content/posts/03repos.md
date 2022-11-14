---
title: "• Genetic data repositories"
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

### Sequence repositories
While basic methods of sequencing have existed since the 1970s, it is in the 21st century that the field exploded. This was due to key technological advancements, collectively known as "high-throughput" or "next-generation" sequencing. Generating high-quality sequences, at a fraction of the price, made genomics accessible to even modest research groups. To illustrate this, [sequencing a human genome in 2000 required many specialist research groups, several years, and $100 million dollars](https://www.genome.gov/about-genomics/fact-sheets/Sequencing-Human-Genome-cost); by 2010 this had dropped to a matter of weeks and $100,000, and in 2022 a human genome (or at least a full exome) can be sequenced for around $800, in a matter of days. As such, aggregated bioinformatic data has followed a parallel of [Moore's Law](https://www.wikiwand.com/en/Moore's_law), roughly doubling every 18 months.

Given the vast amount of data being produced, international consortia were established, with the goal of making any and all bioinformatic research data available to the world, for free (well, ultimately funded by the tax-payer, as with all public research). Most journals will require that sequence information (produced in the research process leading to publication) is deposited with one of the major repositories. These repos should be seen as wonders of the modern scientific world! The openly-accessible data they contain drive forward research in medicine, biology, mathematics, information theory, even nanotechnology and biocomputing.

{{< admonition type="info" title="Some big repositories" open=true >}}
- [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/): this is the USA-based national repo administered by the National Centre for Biotechnology Information, containing a wide-variety of databases. The biggest of these are "genome", "gene", and "protein", as well as raw (pre-assembly) reads, for example in "SRA", the sequence read archive. As of 2022, it contains around 2.5 billion sequences, of around 17 trillion nucleotides.
- [Ensembl](https://www.ensembl.org/index.html): the repository of the European Bioinformatics Institute. It is accompanied by a taxonomy-focused repository, [EnsembleGenomes](https://ensemblgenomes.org/), which is particularly useful for microbial research.
- [ExPASy](https://www.expasy.org/): this is the archive for the Swiss Bioinformatics Institute. Originally is was focussed on proteomics and protein sequences, and today is a major hub for genomics focussed on cell-lines, gene function, drug-design and pathogen classification.
- And there are more, which we encourage you to look for.
{{< /admonition >}}

### Using BioPython to query a repository

Before we use BioPython to carry out a query, we will get familiar with a sequence archive using a browser. The key point here is that, as a starting point for answering a research question, we should become familiar with searching "as a human", before automating such searches.

Let's visit [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). It is worth getting familiar with at least the major sections of GenBank. If the interface feels a bit out-of-date, you'll just have to get used to that! Research never moves as fast a commerce in terms of presentation, but the data remains the same. Consider the clunkiness as endearing!

We are going to search for a well-characterised, ultra-conserved gene, that we are going to use for some common bioinformatics tasks. Our goals here are to find homologous genes from different species, compare and align these, and infer a phylogeny.



```
code goes herer
```
next setion blakfahuir wrafguighbarg

#### Exercise
{{< admonition type="question" title="Exercise" open=true >}}
a question box
- q1
- q2
{{< /admonition >}}

#### Exercise
{{< admonition type="note" title="Solution" open=false >}}
a note box
- dfdsafg
- q2
{{< /admonition >}}
