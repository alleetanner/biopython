---
title: "• Input and Output"
subtitle: "Running sequence information through BioPython."

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

## GenBank files
We'll start by looking at the contents of a GenBank file. Through this you will also learn how to use `urllib`, a standard Python library for getting raw data from internet pages. Create a new file, and add these lines

```
import urllib

url = "https://raw.githubusercontent.com/alleetanner/biopython/main/resources/sonichh_gopherus.gbk"

urllib.request.urlretrieve(url, "sonic-hh-gopherus.gbk")
```

`urllib` is being asked to go to an address (first argument), and save it in a file called `sonic-hh-gopherus.gbk` (the second argument). Now we are going to use BioPython to read it into an object, with `SeqIO.parse`. We need to specify the filename (the one we just downloaded using `urllib`), and the format it is in. If we pass "genbank" (or "gb") as the second argument then it will read it as a GenBank file:

```
record_iterator = SeqIO.parse("sonic-hh-gopherus.gbk", "genbank")
```

We are calling our object variable `record_iterator`, because it could contain many records, but for now it just contains one. We can access this single record by calling [`next`](https://docs.python.org/3/library/functions.html#next):

```
gbk_record = next(record_iterator)
```

GenBank files usually contain a lot more information than a FASTA so more of the fields of the SeqRecord will be filled in. This means that we can, for example, see the record contents:

```
print(gbk_record.annotations)
```
So, have the above five lines of code saved in a file, and run it in the Jupyter terminal. Inspect the output to get familiar with what information a GenBank record typically contains. (Some records are comprehensively annotated, others can be very sparse! This one is quite good.) We could also isolate just the sequence data with

```
print(gbk_record.seq)
```

## FASTA files
### `SeqRecord` objects
If you have a FASTA file with a single sequence in it, we can read it in with `SeqIO.read`. This takes two arguments, the name of the file you want to open, and the format that file is in.

```
srp72_record = SeqIO.read("srp72_single.fas", "fasta")
print(srp72_record)
```

As with a GenBank file, you can also isolate just the sequence information by using `.seq` on the object:

```
print(srp72_record.seq)
```

### Working with FASTAs containing multiple sequences
FASTA files often contain many sequences. These can represent a range of genomic data, for example transcriptomic fragments, a set of genes all from the same organism, or a set of the same gene from a range of organisms. 

To work with these kinds of files, BioPython provides `SeqIO.parse`. This gives you an iterable object, which can, for example, can be looped with a `for` loop. Here we are working with the `srp72.fas` file that you made in the previous section, which should contain multiple sequences.

```
fasta_batch = SeqIO.parse("srp72.fas", "fasta")

for record in fasta_batch:
    print(record.description)
for record in fasta_batch:
    print(record.description)
```

Note that the type of the object returned by `SeqIO.parse` is known as a *generator*. One of the features of a generator in Python is that you can only loop over them once before they are exhausted (this is very useful for memory and speed-critical processes). In the above code snippet, if there was a second loop taking in `fasta_batch`, it would do nothing!

Typically to resolve this we *cast* the generator to a list, which can be iterated over as many times as you like:

```
fasta_batch = list(SeqIO.parse("srp72.fas", "fasta"))

for record in fasta_batch:
    print(record.description)
# this repeat will now work!
for record in fasta_batch:
    print(record.description)
```


Normal text body.
### Exercise
{{< admonition type="exercise" title="Exercise" open=true >}}
Create a Python script:
- Create a function which takes two arguments: [1] an iterable object made by `SeqIO.parse`, [2] a target sequence motif to find.
  - In the function, loop over each record in the iterable.
  - If the target sequence is found in a record, `print` something to the terminal to say it was found, and which record it was.
- Outside of the function, create a variable with your batch of fastas in it (see the `SeqIO.parse` section above).
- Outside of the function, create a variable containing a string to search for (keep it short and simple (say, 3 AAs), and demonstrate a positive match and a negative match).
- Call your function, passing in the iterable and the target motif.
- *Bonus* Have the function report if no results were found.
{{< /admonition >}}

{{< admonition type="warning" title="Solution" open=false >}}
```Python3
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
