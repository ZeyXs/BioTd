from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# exo 1

maSequence1 = Seq("AATGACGTTTACGTTTCATAGTTAA")

print("Exo 1:")
print(len(maSequence1))
print(maSequence1[:3])
print(maSequence1[-5:])
print(maSequence1.find("CCG"))
print(maSequence1.count("TTT"))
        
print([i for i in range(len(maSequence1)) if maSequence1.startswith("TT", i)])

# exo 2 

maSequence2 = Seq("ATGACGTTTACGTTTCATAGTTAA")

print("Exo 2:")
maSequenceComplement = maSequence2.reverse_complement()
print(maSequenceComplement)
print(maSequenceComplement.translate())

nom = maSequence2.complement().count("GC")*2
denom = len(maSequence2)

print(str(round((nom/denom)*100, 2)) + "%")

# exo 3

print("Exo 3:")

uneSequence = list(SeqIO.parse("resources/uneSequence.gb", "genbank"))
print(uneSequence[0].id)
#print(uneSequence[0].seq)

seqFasta = SeqIO.convert("resources/uneSequence.gb", "genbank", "output/seqFasta.fasta", "fasta")
seqFasta = SeqIO.convert("resources/sequences.gb", "genbank", "output/sequencesFasta.fasta", "fasta")

def convert_from_genbank(file_path: str):
    for i in range(len(seq[0])):
        seq = list(SeqIO.parse(file_path, "genbank"))
        print(f"L'identifiant de la séquence est {seq[i].id}")
        print("Sa suite de lettres est :")
        print(seq[i].seq)
    
convert_from_genbank("resources/uneSequence.gb")