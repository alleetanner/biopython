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
While basic methods of sequencing have existed since the 1970s, it is in the 21st century that the field exploded. This was due to key technological advancements, collectively known as "high-throughput" or "next-generation" sequencing. Generating high-quality sequences, at a fraction of the price, made genomics accessible to even modest research groups. To illustrate this, [sequencing a human genome in 2000 required many specialist research groups, several years, and $100 million dollars](https://www.genome.gov/about-genomics/fact-sheets/Sequencing-Human-Genome-cost); by 2010 this had dropped to a matter of weeks and $100,000. In 2022 a a human genome (or at least a full exome) can be sequenced for a few hundred dollars, in a matter of days - as such, aggregated bioinformatic data follows a parallel of [Moore's Law](https://www.wikiwand.com/en/Moore's_law), roughly doubling every 18 months.

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

We are going to search for a well-characterised, ultra-conserved gene, that we are going to use for some common bioinformatics tasks. By the dropdown, select "Protein" - this is the name of the database we will be searching. Note the others in the menu - some are major DBs like gene and genome, others are more niche. In the search box, enter `signal recognition particle 72`. This is a ubiquitous gene in eukaryotes, and arguably in prokaryotes too.

The search results will be shown - above the results, change the results per page from `20 per page` to `100 per page`. Have a look down the results, and see what details you can recognise. Open one by clicking the link, and have a look over the details that typically accompany an accession. 

Go back to your results list. Now, we are going to select a few of these, and get the protein sequences as a FASTA file. Select by ticking the box next to a record:
- choose ten records
- make sure each of then ten are from a different species (eg, don't get more than one from _Homo sapiens_)
- make sure each one has a record with an amino acid length of between 600 and 690 aa,
- don't worry if the record says `partial`, as long as it meets the above requirements.
Once you have selected 10 (or so) that you like, scroll to the very bottom (or the very top) of the records, and click on the "Summary". From the menu it shows, choose "FASTA (text)". This will open a new page, with the protein sequences for the gene SRP72, for the species you chose. Select all of this text (ctrl-a), and cut it. Now, go to your Jupyter Lab tab, and create a new text file. Paste your FASTA sequences into the new file, and save it (ideally in the folder you are working in today).

W





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
