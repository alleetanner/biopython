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

`Seq` is a class in the BioPython library - which means that we have both data (in this case, sequence information) and functions attributed to the object. Since those functions are part of a class, they are also known as "methods". If this is confusing, don't worry too much! It just means that we can make special strings with `Seq`, and each also has special actions that can be done to it. It makes manipulating sequencing information easy, with a clean syntax. `Seq` is imported like this:

```
from Bio.Seq import Seq
```

Let's start by creating a short string of DNA. We can then call `Seq` on it to create a new object, called `dna_seq`. 

```
dna_string = "TTACCAAAAACCCCTTTGGGAAAGCAT"
dna_seq = Seq(raw_dna)
```

This command _could_ be done with a single command of `dna_seq = Seq("TTACCAAAAACCCCTTTGGGAAAGCAT")`, but we'll keep it as two steps for clarity. You might be reading sequence information in from a file, so it is worth being explicit about where the data is coming from. We can examine these, so lets do that with

```
print(type(dna_seq))
print(dir(dna_seq))
```
These will print out what kind of object our `dna_seq` is, as well as the methods associated with it. Don't worry too much about all the dunders there, but notice the methods, at the end of the list.


### Transcription

You'll see one of the methods is `transcribe()`, so let's do that and turn our DNA into RNA:

```
rna_seq = dna_seq.transcribe()
print(rna_seq)
```

We might also want to get the complement of this RNA strand. Since `rna_seq` is also a `Seq` object, we can call one of the methods, `complement_rna()`:

```
rna_complement = rna_seq.complement_rna()
print(rna_complement)
```

It is also possible to convert back from an RNA sequence to a DNA sequence:

```
dna_from_rna = rna_seq.back_transcribe()
```

Which, if it's working correctly should give us back the original data:

```
print(dna_seq == dna_from_rna)
```


### Translation
The next thing we might want to do with this nucleotide sequence is see what amino acid sequence it would translate to. 
```
protein_seq = rna_seq.translate()
print(protein_seq)
```
Note that this will always translate forward using the first nucleotide as the open reading frame. Try altering the raw DNA sequence to break this kind of translation. There are two primary ways of doing this - what are they?




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

