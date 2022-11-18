---
title: "• Further tools and resources"
subtitle: "Bringing shell scripts and python scripts together"

date: 2022-11-06T00:00:00+01:00

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

## The `bio` package
We used this command line package [`bio`](https://www.bioinfo.help/). As of late 2022, this is clean and streamlined set of command line tools for bioinformatics tasks. I appreciate that this is accompanied by a book which is paywalled but I will look into whether I can provide access for students. However, just the command line glossary is very useful, which you can install with
```shell
bio --download
```
after this is downloaded, you can get (very succinct) definitions of bioinformatics terms, for example
```shell
bio explain nucleotide_binding_site
```
also, see that asking for an ambiguous term, for example
```shell
bio explain nucleotide
```
will return a list of topics relavent to the query (in this case, nucleotides), that you can then ask for the entry on. I find this a useful way of clarifying terms, and a way of starting more in-depth research of a particular topic.

The `bio` package also has some nice [https://www.bioinfo.help/bio-tips.html](one-liners) which are well worth getting familiar with.

## The `diytranscriptomics.com` community and workspace
In [DIYtranscriptomics](https://diytranscriptomics.com/) you can find courses, lectures and guides on transcriptomic topics. This work does use `R` as their processing language, so that might be useful to see how things are done there.
