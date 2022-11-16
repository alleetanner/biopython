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


## BioPython's `Seq` class

`Seq` is a class in the BioPython library - which means that we have both data (in this case, sequence information) and functions attributed to the object. Since those functions are part of a class, they are also known as "methods". If this is confusing, don't worry too much! It just means that we can make special strings with `Seq`, and each also has special actions that can be done to it. It makes manipulating sequencing information easy, with a clean syntax. 

In your Jupyter Lab text editor, start a new script by importing `Seq` like this:

```python
from Bio.Seq import Seq
```

Let's start by creating a short string of DNA. We can then call `Seq` on it to create a new object, called `dna_seq`. 

```python
dna_string = "CATTACTTTGGCGGAAAAATATTT"
dna_seq = Seq(dna_string)
print(dna_seq)
```

(This command _could_ be done with a single command of `dna_seq = Seq("CATTACTTTGGCGGAAAAATATTT")`, but we'll keep it as two steps for clarity. You might be reading sequence information in from a file, so it is worth being explicit about where the data is coming from.)

Save your file, and in the terminal run the file using `python`. 

We can examine our `dna_seq` object, so let's do that with:

```python
print(type(dna_seq))
print(dir(dna_seq))
```

These will print out what kind of object our `dna_seq` is, as well as the methods associated with it. Don't worry too much about all the `__dunders__` there, but notice the methods, at the end of the list. You can learn more about how to use these methods in [the BioPython documentation for `Seq`](https://biopython.org/wiki/Seq). Remove these `type` and `dir` lines from your script, they will make your terminal a bit messy and hard to read!

{{< admonition type="info" title="Remember to save your script!" open=false >}}
It is very easy to forget to save your script before running it! Jupyter Lab will tell you if a file is modified without saving, with the circle in the tab title. A `ctrl-S` will save it quickly. Also remember to be in the right folder in the terminal, when you are running your script.
{{< /admonition >}}

### Transcription

You'll see one of the methods of `Seq` is `transcribe()`, so let's do that and turn our DNA into RNA:

```python
rna_seq = dna_seq.transcribe()
print(rna_seq)
```

We might also want to get the complement of this RNA strand. Since `rna_seq` is also a `Seq` object, we can call one of the methods, `complement_rna()`:

```python
rna_complement = rna_seq.complement_rna()
print(rna_complement)
```

It is also possible to convert back from an RNA sequence to a DNA sequence:

```python
dna_from_rna = rna_seq.back_transcribe()
```

Which, if it's working correctly should give us back the original data - we can check that by asking "are these two objects identical?", like this:

```python
if dna_seq == dna_from_rna:
    print("Yes, they are identical.")
else:
    print("No, they differ.")
```

### Translation
The next thing we might want to do with this nucleotide sequence is see what amino acid sequence it would translate to (let's assume our sequence is protein-coding for a tiny polypeptide!). 

```python
protein_seq = rna_seq.translate()
print(protein_seq)
```

Note that this will always translate forward using the first nucleotide as the open reading frame. Try altering the raw DNA sequence to break this kind of translation (you should be getting a `BiopythonWarning`, or a `Bio.Data.CodonTable.TranslationError` - read these tracebacks and pick out the important information it is giving you). There are two primary ways of doing this - what are they?

### Exercise
{{< admonition type="note" title="Exercise" open=true >}}
Let's put together some of the methods we have used so far, and arrange them in functions. Create a new script in Jupyter Lab, and build a DNA translator, using the script below as a guide (ie, fill the `___` gaps!). When run on the command line, your script should neatly output the manipulations to the incoming raw DNA.

Once you have a working script, use it to answer these questions:
- Which of the two reading directions is probably correct?
- Why?
- Is this 5' to 3', or 3' to 5'?

```python
from Bio.Seq import Seq


def complement_transcribe_translate(___):
    """
    Very good students fill out the docstring too! 
    """
    print(f"My incoming DNA:           {___}")
    
    complement = ___.___()
    print(f"Complementary strand:      {___}")
    
    transcribed_forward = ___.___()
    print(f"Complement transcription:  {___}")
    
    amino_acids_forward = ___.___()
    print(f"Amino acid translation:    {___}")


def reverse_translate(___):
    """
    Very good students fill out the docstring too! 
    """
    dna_reverse = ___.___()
    print(f"As a reverse complement:   {___}")
    
    transcribed_backward = ___.___()
    print(f"Reverse transcription:     {___}")
    
    amino_acids_backward = ___.___()
    print(f"Amino acids reading <-     {___}")


# here is my incoming DNA sequence, made into a Seq object
dna_seq = Seq("TTACCAAAAACCCCTTTGGGAAAGCAT")

# call my function "complement_transcribe_translate" from above
complement_transcribe_translate(___)

# call my function "reverse_translate" from above
reverse_translate(___)
```
{{< /admonition >}}

{{< admonition type="tip" title="Solution" open=false >}}
```python
from Bio.Seq import Seq


def complement_transcribe_translate(dna):
    """
    Take some DNA, and print out some conversions:
    the complementary strand,
    the transcription of the complementary strand,
    the amino acid translation of this complement.
    
    Args: DNA as a Seq object.
    
    Returns: nothing, this just prints to the terminal.
    """
    
    print(f"My incoming DNA:           {dna}")
    
    complement = dna.complement()
    print(f"Complementary strand:      {complement}")
    
    transcribed_forward = complement.transcribe()
    print(f"Complement transcription:  {transcribed_forward}")
    
    amino_acids_forward = transcribed_forward.translate()
    print(f"Amino acid translation:    {amino_acids_forward}")
    
    
def reverse_translate(dna):
    """
    Take some DNA, and print out some conversions:
    the reverse complement,
    the reverse transcription,
    the amino acid translation of reverse transcription.
    
    Args: DNA as a Seq object.
    
    Returns: nothing, this just prints to the terminal.
    """
    
    dna_reverse = dna.reverse_complement()
    print(f"As a reverse complement:   {dna_reverse}")
    
    transcribed_backward = dna_reverse.transcribe()
    print(f"Reverse transcription:     {transcribed_backward}")
    
    amino_acids_backward = transcribed_backward.translate()
    print(f"Amino acids reading <-     {amino_acids_backward}")


# here is my incoming DNA sequence, made into a Seq object
dna_seq = Seq("TTACCAAAAACCCCTTTGGGAAAGCAT")

# call my function "complement_transcribe_translate" from above
complement_transcribe_translate(dna_seq)

# call my function "reverse_translate" from above
reverse_translate(dna_seq)
```
{{< /admonition >}}

