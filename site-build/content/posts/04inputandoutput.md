---
title: "• Working with FASTA"
subtitle: "A standard format for bioinformatic data"

date: 2022-11-04T00:00:00+01:00

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

## FASTA files
### `SeqIO`
`SeqIO` is used for managing sequences as _inputs_ and _outputs_. We are going to learn three methods from `SeqIO`:
- `SeqIO.read` for bringing in single sequences
- `SeqIO.parse` for working with multiple sequences, and
- `SeqIO.write` for creating new files.
In the previous section, we created two `FASTA` files: one a single sequence, the other a multi sequence. Let's start be reading in a single sequence

#### `SeqIO.read`
Start a fresh Python script. Let's read in the single `FASTA` file:
```python
from Bio import SeqIO

fasta_file = "srp_single.fas"
srp72_record = SeqIO.read(fasta_file, "fasta")

print(srp72_record)
```
Have a look at the output; you'll see that it has been turned into a `Seq` object. We can also just look at the sequence by using the attribute `.seq` on the object:

```python
print(srp72_record.seq)
```

### Working with FASTAs containing multiple sequences
FASTA files often contain many sequences. These can represent a range of genomic data, for example transcriptomic fragments, a set of genes all from the same organism, or a set of the same gene from a range of organisms. 

To work with these kinds of files, BioPython provides `SeqIO.parse`. This gives you an iterable object, which can, for example, can be looped with a `for` loop. Here we are working with the `srp72.fas` file that you made in the previous section, which should contain multiple sequences.

```python
fasta_batch = SeqIO.parse("srp72.fas", "fasta")

for record in fasta_batch:
    print(record.description)
for record in fasta_batch:
    print(record.description)
```

Note that the type of the object returned by `SeqIO.parse` is known as a *generator*. One of the features of a generator in Python is that you can only loop over them once before they are exhausted (this is very useful for memory and speed-critical processes). In the above code snippet, if there was a second loop taking in `fasta_batch`, it would do nothing!

Typically to resolve this we *cast* the generator to a list, which can be iterated over as many times as you like:

```python
fasta_batch = list(SeqIO.parse("srp72.fas", "fasta"))

for record in fasta_batch:
    print(record.description)
# this repeat will now work!
for record in fasta_batch:
    print(record.description)
```

### Exercise
{{< admonition type="note" title="Exercise" open=true >}}
Create a Python script:
- Create a function which takes two arguments: [1] an iterable object made by `SeqIO.parse`, [2] a target sequence motif to find.
  - In the function, loop over each record in the iterable.
  - If the target sequence is found in a record, `print` something to the terminal to say it was found, and which record it was.
- Outside of the function, create a variable with your batch of fastas in it (see the `SeqIO.parse` section above).
- Outside of the function, create a variable containing a string to search for (keep it short and simple (say, 3 AAs), and demonstrate a positive match and a negative match).
- Call your function, passing in the iterable and the target motif.
- *Bonus!*  Have the function report if no results were found.
{{< /admonition >}}

{{< admonition type="tip" title="Solution" open=false >}}
```python
from Bio import SeqIO


def fasta_search(fastas, target):
    """
    Search for a sequence motif in multiple FASTA records.
    
    Args:    fastas - SeqIO iterable
             target - the sequence to find (string)
    Returns: Nothing.
    """
    # create hit success boolean
    match_success = False
    
    for record in fasta_batch:
        
        if target in record.seq:
            print(f"MATCH FOUND FOR TARGET {target}")
            print(record.description)
            print(len(record.seq))
            # HIT! toggle my success boolean
            match_success = True
    
    # report if nothing was found.
    if match_success == False:
        print(f"No matches found for target {target}")


# parse in my fasta file
fasta_batch = list(SeqIO.parse("srp72.fas", "fasta"))
# declare what I am looking for, as a string
motif = "RDA"

# call the function!
fasta_search(fasta_batch, motif)
```
{{< /admonition >}}

## Output
We have seen how `Bio.SeqIO.parse` is used for sequence input (reading files), and now we’ll look at `Bio.SeqIO.write` which is for sequence output (writing files). This function takes three arguments: [1] some `SeqRecord` object(s), [2] a filename to write to, and [3] a sequence format.

Here is an example, where we start by creating a few `SeqRecord` objects the hard way (by hand, rather than by loading them from a file):

```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

rec1 = SeqRecord(
    Seq(
        "MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
        "SSACKSLKEAFTPLGISQQLLWWLQ",
    ),
    id="gi|14150838|gb|AAK54648.1|AF376133_1",
    description="chalcone synthase [Cucumis sativus]",
)

rec2 = SeqRecord(
    Seq(
        "YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ"
        "DMVVVEIPKLGKEAAVKAIKEWGQAA",
    ),
    id="gi|13919613|gb|AAK7873142.1|",
    description="chalcone synthase [Fragaria vesca subsp. bracteata]",
)

rec3 = SeqRecord(
    Seq(
        "QQGCFAGGTVLRLAKDLAENNRGAKRYMYLTEEAAVKALLILKENPSMCEYMAPSLDARQ"
        "VVCSEITAVTFRGPEAAVKAIKEWGQ",
    ),
    id="gi|48373667|gb|AAK33312.1|",
    description="chalcone synthase [Quercus robur]",
)

# create a list of SeqRecords
my_records = [rec1, rec2, rec3]

# run SeqIO.write to create an output file
SeqIO.write(my_records, "my_example.fas", "fasta")
```
{{< admonition type="warning" title="Be careful!" open=false >}}
Writing out files like this *will permanently overwrite any existing files of the same name!* Always be careful with names and paths that you are writing to. 
{{< /admonition >}}

### Exercise
{{< admonition type="note" title="Exercise" open=true >}}
Have a little walk, or even do some stretches. Maybe just in the room if it is raining. Outside is better, but admittedly it is often winter. Sometimes just looking out the window will do.
{{< /admonition >}}



