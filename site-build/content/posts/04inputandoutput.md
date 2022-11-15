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

# GenBank files
We'll start by looking at the contents of a GenBank file. Through this you will also learn how to use `urllib`, a standard Python library for getting raw data from internet pages. Create a new file, and add these lines

```
import urllib
urllib.request.urlretrieve("https://raw.githubusercontent.com/alleetanner/biopython/main/resources/sonichh_gopherus.gbk", "sonic-hh-gopherus.gbk")
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








# FASTA files
## `SeqRecord` objects

If you have a FASTA file with a single sequence in it, the simplest way to read it in is with SeqIO.read. This takes two arguments, the name of the file you want to open, and the format that file is in.

### Section title
* bullet1
* bullt2

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
{{< admonition type="question" title="Questions" open=true >}}
a question box
- q1
- q2
{{< /admonition >}}
