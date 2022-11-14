from Bio import SeqIO


fasta_file = "srp72.fas"


for record in SeqIO.parse(fasta_file, "fasta"):
    print(record.description)
                        
