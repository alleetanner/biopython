from Bio import SeqIO
from Bio.Blast import NCBIWWW


blast_query = SeqIO.read("homo_srp72.fas", "fasta")

print(blast_query)

results = NCBIWWW.qblast("blastp", "nt", blast_query.seq)

print(results)
