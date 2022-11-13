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
  enable: false
comment:
  enable: false
---


### `Seq`

The main module to get to grips with in Biopython is Seq. This is the primary interface to sequence data that you will work with. It is imported as:

```
from Bio.Seq import Seq
```

Using `Seq` we can create sequence objects from standard Python strings:

```
my_dna = Seq("AGTACACTGGTT")
```

#### RNA

You can the corresponding RNA sequence from a DNA sequence by using the transcribe method:

```
my_rna = my_dna.transcribe()
```

Once you have an RNA sequence, you can again do standard operations on it, such as getting the complement:

```
my_rna.complement_rna()
```

Seq('UCAUGUGACCAA')

It is also possible to convert back from an RNA sequence to a DNA sequence:

my_dna_from_rna = my_rna.back_transcribe()

Which, if it's working correctly should give us back the original data:

my_dna == my_dna_from_rna

True

Translation

Once we have an RNA sequence, you can get the expressed protein with translate:

my_protein = my_rna.translate()

my_protein

Seq('STLV')




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
{{< /admonition >}}

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

