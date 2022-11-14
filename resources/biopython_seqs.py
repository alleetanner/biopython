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