---
title: "• Alignments"
subtitle: "Homology at the molecular level."

date: 2022-11-05T00:00:00+01:00

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

## Multiple Sequence Alignments

Multiple Sequence Alignments are a collection of sequences which have been arranged to form a *matrix* with every sequence being the same length. This arrangement is usually done by insertion of gap. Each position along an alignment is *an assertion of homology*: we are expressing an evolutionary relationship between each nucleotide, codon, or amino acid (characteristics of the *genotype*). This is analagous to a trait character matrix, where our observations are on morphological, behavioural or other outcomes of gene expression (characteristics of the *phenotype*).

Alignment can be regarded as a matrix of letters, where each row is held as a `SeqRecord` object internally. The `MultipleSeqAlignment` object holds this kind of data, and the `AlignIO` module is used for reading and writing these as various file formats.

For this section we need to download some files which we will be reading in. Create a small Python script to download these from a URL:

```python
import urllib
for i in ["PF05371_seed.sth", "dummy_aln.phy"]:
    urllib.request.urlretrieve(f"https://milliams.com/courses/biopython/{i}", i)
```

### Parsing or Reading Sequence Alignments

Much like `SeqIO`, `AlignIO` contains two functions for reading in sequence alignments:
- `read()` will return a single MultipleSeqAlignment object
- `parse()` will return an iterator which gives MultipleSeqAlignment objects

Both functions require two arguments:
- A string specifying the name of the file to open.
- A string specifying the alignment format. [See here for a full listing of supported formats](http://biopython.org/wiki/AlignIO)

### Single alignments
Let's start with a single alignments file which contains the seed alignment for the `Phage_Coat_Gp8 (PF05371)` PFAM entry. The file contains a lot of annotation information but let's just go ahead and load it in to see how it looks:

```python
from Bio import AlignIO
aln_seed = AlignIO.read("PF05371_seed.sth", "stockholm")
print(aln_seed)
```

```
Alignment with 7 rows and 52 columns
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73
```

Note in the above output the sequences have been elided in the middle (`...`). We could instead write our own code to format this as we please by iterating over the rows as `SeqRecord` objects and printing the first 50 values of each sequence:

```python
for record in aln_seed:
    print(f"{record.seq[:50]} - {record.id}")
```

```
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSS - COATB_BPIKE/30-81
AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVS - Q9T0Q8_BPIKE/1-52
DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSS - COATB_BPI22/32-83
AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS - COATB_BPM13/24-72
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFAS - COATB_BPZJ2/1-49
AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTS - Q9T0Q9_BPFD/1-49
FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVS - COATB_BPIF1/22-73
```

With any supported file format, we can load an alignment in the same way just by changing the format string. For example, use `phylip` for PHYLIP files, `nexus` for NEXUS files or `emboss` for the alignments output by the EMBOSS tools.

### Multiple Alignments

In general, Alignment files often contain multiples alignments, and to read these files we must use the AlignIO.parse function.

We have previously downloaded a file called `dummy_aln.phy` which contains some dummy alignment information in PHYLIP format. If we wanted to read this in using `AlignIO` we could use:

```python
aln_dummy = AlignIO.parse("dummy_aln.phy", "phylip")
for alignment in aln_dummy:
    print(alignment)
    print("---")
```

```
Alignment with 5 rows and 6 columns
AAACCA Alpha
AAACCC Beta
ACCCCA Gamma
CCCAAC Delta
CCCAAA Epsilon
---
Alignment with 5 rows and 6 columns
AAACAA Alpha
AAACCC Beta
ACCCAA Gamma
CCCACC Delta
CCCAAA Epsilon
---
Alignment with 5 rows and 6 columns
AAAAAC Alpha
AAACCC Beta
AACAAC Gamma
CCCCCA Delta
CCCAAC Epsilon
---
```

The `.parse()` function returns an iterator. If we want to keep all the alignments in memory at once, then we need to turn the iterator into a list, just as we did with `SeqIO.parse`:

```python
alignments = list(AlignIO.parse("dummy_aln.phy", "phylip"))
second_aln = alignments[1]
print(second_aln)
```

```
Alignment with 5 rows and 6 columns
AAACAA Alpha
AAACCC Beta
ACCCAA Gamma
CCCACC Delta
CCCAAA Epsilon
```

### Writing Alignments

Now we’ll look at `AlignIO.write()` which is used for alignments output (writing files).

This function takes three arguments:

- Some MultipleSeqAlignment objects
- A string specifying a handle or a filename to write to
- A lower case string specifying the sequence format.

We start by creating a `MultipleSeqAlignment` object the hard way (by hand!). Note we create some `SeqRecord` objects to construct the alignment from.

```python
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

align1 = MultipleSeqAlignment([
    SeqRecord(Seq("ACTGCTAGCTAG"), id="toto"),
    SeqRecord(Seq("ACT-CTAGCTAG"), id="titi"),
    SeqRecord(Seq("ACTGCTAGDTAG"), id="tata"),])

print(align1)
```

```
Alignment with 3 rows and 12 columns
ACTGCTAGCTAG toto
ACT-CTAGCTAG titi
ACTGCTAGDTAG tata
```

Now let's try to output, in PHYLIP format, these alignments in a file with the `Phage_Coat_Gp8` alignments.

```python
my_alignments = [align1, aln_seed]
AlignIO.write(my_alignments, "mixed.phy", "phylip")
```

## Exercise
{{< admonition type="note" title="Exercise" open=true >}}
Read in the alignment in `PF05371_seed.sth` and write it out in PHYLIP format.
{{< /admonition >}}
