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

To work with these kinds of files, BioPython provides `SeqIO.parse`. Unlike `SeqIO.read`, this creates an iterable object (it is actually a _single use_ iterable, called a _generator_): so, this can be used as an input for a loop. Here we are working with the `srp_multi.fas` file that you made in the previous section, which contains multiple sequences.

```python
from Bio import SeqIO

fasta_file = "srp_multi.fas"
fasta_batch = SeqIO.parse(fasta_file, "fasta")

for record in fasta_batch:
    print(record.description)
```

We can only loop over a generator once before they are exhausted (this is very useful for memory and speed-critical processes). In the above code snippet, if there was a second loop taking in `fasta_batch`, it would do nothing!

If we wish to loop again over a generator, we usually *cast* the generator to a list, which can be iterated over as many times as you like:

```python
from Bio import SeqIO

fasta_file = "srp_multi.fas"
fasta_batch = list(SeqIO.parse(fasta_file, "fasta"))

for record in fasta_batch:
    print(record.description)
# this repeated loop will now work!
for record in fasta_batch:
    print(record.description)
```

## Output
We have seen how `Bio.SeqIO.parse` is used for sequence input (reading files), and now we’ll look at `Bio.SeqIO.write` for creating files. This function takes three arguments: [1] some `SeqRecord` object(s), [2] a filename to write to, and [3] a sequence format.

Here is an example, where we start by creating a few `SeqRecord` objects (by hand, rather than by loading them from a file):

```python
from Bio import SeqIO
from Bio.Seq import Seq

fasta_file = "srp_multi.fas"

srp72_parsed = SeqIO.parse(fasta_file, "fasta")

# run SeqIO.write to create an output file
SeqIO.write(srp72_parsed, "my_example.fas", "fasta")
```

### Exercises
{{< admonition type="note" title="Exercise 1" open=true >}}
Create a new Python file, containing this:
```python
from Bio import SeqIO

fasta_file = "srp_multi.fas"
srp72_parsed = SeqIO.parse(fasta_file, "fasta")

for record in srp72_parsed:
    print(record.id)

SeqIO.write(srp72_parsed, "output_file.fas", "fasta")
```
- Run this script, and examine the contents of `output_file.fas`... which will be empty! Why?
- Edit the script so that it correctly outputs a fasta file.
{{< /admonition >}}
{{< admonition type="tip" title="Solution" open=false >}}
You could fix this in a couple of ways; firstly we could cast to a list:
```python
from Bio import SeqIO

fasta_file = "srp_multi.fas"
srp72_parsed = list(SeqIO.parse(fasta_file, "fasta"))

for record in srp72_parsed:
    print(record.id)

SeqIO.write(srp72_parsed, "output_file.fas", "fasta")
```
(or we could avoid exhausting the generator by skipping the `print` loop)
```python
from Bio import SeqIO

fasta_file = "srp_multi.fas"
srp72_parsed = SeqIO.parse(fasta_file, "fasta")

SeqIO.write(srp72_parsed, "output_file.fas", "fasta")
```
{{< /admonition >}}

{{< admonition type="note" title="Exercise 2" open=true >}}
- Ensure your script from exercise one uses `SeqIO` to parse your incoming FASTA file to a list.
- Add the single sequence in your file `srp_single.fas` to that list.
- Output a new FASTA file which includes this sequences (so, it should have 11 sequences in it)
{{< /admonition >}}

{{< admonition type="tip" title="Solution" open=false >}}
 ```python
from Bio import SeqIO

fasta_file = "srp_multi.fas"

srp72_parsed = list(SeqIO.parse(fasta_file, "fasta"))

fasta_to_add = "srp_single.fas"
new_fasta = SeqIO.read(fasta_to_add, "fasta")

srp72_parsed.append(new_fasta)

SeqIO.write(srp72_parsed, "output_file.fas", "fasta")
```
{{< /admonition >}}

{{< admonition type="note" title="Exercise 3" open=true >}}
- Have a look at the [documentation for `SeqIO`](https://biopython.org/wiki/SeqIO). 
- Modify your script from exercise 2 so that it outputs in tabular format.
{{< /admonition >}}

{{< admonition type="tip" title="Solution" open=false >}}
The third argument of `SeqIO.write` is the output format, so we can change the outputting line to 
```python
SeqIO.write(srp72_parsed, "output_file.fas", "tab")
```
{{< /admonition >}}

{{< admonition type="warning" title="Be careful!" open=false >}}
Writing out files like this *will permanently overwrite any existing files of the same name!* Always be careful with names and paths that you are writing to. 
{{< /admonition >}}

### Exercise [bonus]
{{< admonition type="note" title="Bonus exercise" open=false >}}
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

### Exercise
{{< admonition type="note" title="Exercise" open=false >}}
Have a little walk, or even do some stretches. Maybe just in the room if it is raining. Outside is better, but admittedly it is often winter. Sometimes just looking out the window will do.
{{< /admonition >}}



