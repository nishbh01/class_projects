import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import  Entrez, Medline, SeqIO


# NCBI Eutilites # data extraction for mRNA transcripts
from Bio import Entrez
Entrez.email = 'nbhanda1@ramapo.edu'
handle = Entrez.esearch(db='nucleotide', term='Tumor Suppressors mRNA AND human[orgn]', retmax = 2000)
record = Entrez.read(handle)
handle.close()
Id_record = record['IdList']

# fetching the data 
fasta = Entrez.efetch(db='nucleotide', rettype = 'fasta', retmode = 'text', id=Id_record)
fasta_read = fasta.read()
fasta.close()
with open("suppressors_transcript.fasta", 'w') as f:
    f.write(fasta_read)
f.close()



from Bio import SeqIO

sequences = SeqIO.parse("suppressors_transcript.fasta", "fasta")
sfiltered_sequences = [seq for seq in sequences if "NM" in seq.id]
SeqIO.write(sfiltered_sequences, "snucleotide_seq.fasta", "fasta")
print(len(sfiltered_sequences))

isoleucine_codons = ['ATA', 'ATC', 'ATT']

def codon_frequency(sequences, codons):
    frequency = {}
    for i in range(len(sequences)):
        sequence = sequences[i]
        codon = str()
        
        for nucleotide in sequence:
            codon += nucleotide
        
            if len(codon) % 3 == 0:
                if codon in isoleucine_codons:
                    try:
                        frequency[codon] += 1
                    except:
                        frequency[codon] = 0
                    codon = str()
                else:
                    codon = str()
    return frequency

input_file = "snucleotide_seq.fasta"
sequences1 = [str(record.seq) for record in SeqIO.parse(input_file, "fasta")]
frequency = codon_frequency(sequences1, isoleucine_codons)

print("Isoleucine codon frequency:")
for codon, count in frequency.items():
    print(f"{codon}: {count}")
codons = frequency.keys()
count = frequency.values()
plt.bar(codons, count, color = ['r', 'b', 'g'])
plt.show()