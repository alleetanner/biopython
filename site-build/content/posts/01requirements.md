---
title: "• Working with sequences"
subtitle: "The ultimate biological microscope"

date: 2022-11-02T00:00:00+01:00

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
  enable: true
comment:
  enable: true
---

{{< admonition type="note" title="Glossary" open=false >}}
The words **library**, **module** and **package** roughly mean the same thing; it is entirely natural to be confused by this!
- **Library** : a collection of python files that expand the ability of python, using the `import` command. These are like accessories or modifications, typically giving you access to powerful, professional functions and classes written by collaborations of expert coders.
  - A library can be [standard](https://docs.python.org/3/library/index.html) (it comes built into python), for example `time` or `math`, or
  - 3rd party, so, not a core part of python itself and will need to be installed.
- **Package** : a library that is available for delivery to a package manager, such as `pip` or `conda`. [PyPI](https://pypi.org/) is the main package storehouse for python.
- **Module** : anything that is imported to a main running script. Libraries and packages are made of modules.
{{< /admonition >}}


### Python experience
To get the most out of this session, you should have some experience with Python. In particular, you should
* have experience to beginner or intermediate level Python - for example, having attended our sessions on either of these would be an ideal introduction
* be familiar with the biological background to genetics or genomics. For example, knowing the relationships between nucleotides and amino acids and having a general knowledge of what these are and what they do is sufficient.
* bdbkj
* some experience of using the [command line](https://alleetanner.github.io/intro-to-command-line/) is also expected.

We will recap some important coding concepts like `import`ing, creating functions using `def`, and how to work with dataframes. Streamlit encourages a "linear" coding style, and uses the Pandas dataframe as a foundation.sdfio no

### Installing required software
You will need Python version 3.10 or higher installed on your machine. You can check what version you have installed with the command `python3 --version`. If you need to upgrade, search for a guide to follow on how to do that for your operating system.

We also need to install some packages

Normal text body.

{{< admonition type="warning" open=true >}}
- Warning1
- Warning2
{{< /admonition >}}
{{< admonition type="info" open=true >}}
- Info1
- Info2
{{< /admonition >}}

### Another section
blahfdajkfhu aefuibahj rwefyugr

```
code goes herer
```
next setion blakfahuir wrafguighbarg

### Exercise
{{< admonition type="question" title="Exercise" open=true >}}
Let's put together some of the functions we have used so far to manipulate some sequence information. Look at the script below and complete the gaps. When run on the command line, your script should neatly output the manipulations to the incoming raw DNA.

Answer these questions once you are finished:
- Which of these two reading frames is probably correct?
- Why?
- 
```
from Bio.Seq import Seq


raw_dna = Seq("TTACCAAAAACCCCTTTGGGAAAGCAT")

print(f"My incoming raw DNA:       {____________________}")

complement = ____________________

print(f"Complementary strand:      {____________________}")

transcribed_forward = ____________________

print(f"Complement transcription:  {____________________}")

amino_acids_forward = ____________________

print(f"Amino acid translation:    {____________________}")



dna_reverse = ____________________

print(f"As a reverse complement:   {____________________}")

transcribed_backward = ____________________

print(f"Reverse transcription:     {____________________}")

amino_acids_backward = ____________________

print(f"Amino acids reading <-     {____________________}")
```

{{< admonition type="warning" title="Solution" open=false >}}
```
from Bio.Seq import Seq


raw_dna = Seq("TTACCAAAAACCCCTTTGGGAAAGCAT")

print(f"My incoming DNA:           {raw_dna}")

complement = raw_dna.complement()

print(f"Complementary strand:      {complement}")

transcribed_forward = complement.transcribe()

print(f"Complement transcription:  {transcribed_forward}")

amino_acids_forward = transcribed_forward.translate()

print(f"Amino acid translation:    {amino_acids_forward}")


# From here we will work with reverse complements

dna_reverse = raw_dna.reverse_complement()

print(f"As a reverse complement:   {dna_reverse}")

transcribed_backward = dna_reverse.transcribe()

print(f"Reverse transcription:     {transcribed_backward}")

amino_acids_backward = transcribed_backward.translate()

print(f"Amino acids reading <-     {amino_acids_backward}")
```
{{< /admonition >}}
{{< /admonition >}}
