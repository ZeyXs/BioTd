from Bio import SeqIO

# exo 1

seq = "TCTGTTAACCATCCACTTCG"

def occurences_count(seq: str) -> str:
    a = seq.count("A")
    c = seq.count("C")
    g = seq.count("G")
    t = seq.count("T")
    return f"A: {a}, C: {c}, G: {g}, T: {t}"

def occurences(seq: str) -> dict:
    occ = {"A": 0, "C": 0, "G": 0, "T": 0}
    for c in seq:
        occ[c] += 1
    return occ

print("Exo 1.1:", occurences_count(seq))
print("Exo 1.2:", occurences(seq))

# exo 2

seq = "TCTGTTATGACCATCCACTTCG"

def seq_valide(seq: str) -> bool:
    return ("ATG" in seq, seq.find("ATG"))
    
print("Exo 2:", seq_valide(seq))

# exo 3

seq1= "ATGCGTAGTCGT"
seq2= "AGGTTCGTATG" 
seq_list = [seq1,seq2]

def occurences_list(seq: list) -> dict:
    occ = {"A": 0, "C": 0, "G": 0, "T": 0}
    for i in range(len(seq)):
        for c in seq[i]:
            occ[c] += 1
    return occ

def seq_valide_list(seq: list) -> bool:
    full_seq = ""
    for i in range(len(seq)):
        full_seq += seq[i]
    return ("ATG" in full_seq, full_seq.find("ATG"))

print("Exo 3:", occurences_list(seq_list))
print("Exo 3:", seq_valide_list(seq_list))

# exo 4

seq1 = "TACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGTCGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCCCTTCCCGTGGGGGGAACGCGTGCATTTATCAGATCAAAACCAACCCGGTCAGCCCCCTCCCGGCTCCGGCCGGGGGTCGGGCGCCGGCGGCTTTGGTGACTCTAGATAACCTCGGGCCGATCGCACGTCCCCGTGGCGGCGACGACCCATTCGA"
seq2 = "TACCTGGTTGATCCTGCCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCATGTCTAAGTACGCACGGCCGGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGGTCGCTCGCTCCTCTCCTACTTGGATAACTGTGGTAATTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCCCTTCCCGTGGGGGGAACGCGTGCATTTATCAGATCAAAACCAACCCGGTCAGCCCCCTCCCGGCTCCGGCCGGGGGTCGGGCGCCGGCGGCTTTGGTGACTCTAGATAACCTCGGGCCGATCGCACGTCCCCGTGGCGGCGACGACCCATTCGA"

def diff(seq1: str, seq2: str) -> str:
    diff = 0
    length = min(len(seq1), len(seq2))
    for i in range(length):
        if seq1[i] != seq2[i]:
            diff += 1
    return f"{diff} résidus différents sur {length}, soit {round(((length-diff)/length) * 100, 2)}% d'identité"

print("Exo 4:", diff(seq1, seq2))

# exo 5

def open_fasta(file_path: str):
    with open("resources/" + file_path) as file:
        for line in file:
            if line[0] != ">":
                print(line.strip())

print("Exo 5:")
open_fasta("maSequence.fasta")

# exo 6

def open_fasta_lower(file_path: str):
    with open("resources/" + file_path) as file:
        for line in file:
            if line[0] != ">":
                print(line.strip().lower())

print("Exo 6:")
open_fasta_lower("maSequence.fasta")

# exo 7

def entete(file_path: str):
    entetes = []
    with open("resources/" + file_path) as file:
        for line in file:
            if line[0] == ">":
                entetes.append(line)
    with open("output/maSeqResultat.txt", 'w') as file:
        for entete in entetes:
            file.write(str(entete))

print("Exo 7.1: (Fichier créer dans ./output/maSeqResultat.txt)")
entete("maSequence.fasta")

def num_accession(file_path: str):
    with open("output/" + file_path) as file:
        for line in file:
            header_list = line.split("|")
            print(header_list[3])

print("Exo 7.2:")
num_accession("maSeqResultat.txt")

# exo 8

def purinique_pyrimidique():
    pp_dict = {'G': 'purinique', 'A': 'pyrimidinique', 'C': 'purinique', 'T': 'pyrimidinique'}
    u_input = input("Renseigner une base : ")
    if u_input in pp_dict:
        print(pp_dict[u_input])

purinique_pyrimidique()

# exo 9

def convert_to_genbank(file_path: str):
    input_handle = open("resources/" + file_path)
    output_handle = open("output/maSequence.gb", "w")

    sequences = list(SeqIO.parse(input_handle, "fasta"))

    for seq in sequences:
        seq.annotations["molecule_type"] = 'DNA'
    
    SeqIO.write(sequences, output_handle, "genbank")

    output_handle.close()
    input_handle.close()

print("Exo 9 : Fichier créer dans ./output/maSequence.gb")
convert_to_genbank("mesSequences.fasta")

# exo 10

class ORF():
    
    def __init__(self, orf: str):
        self.orf = orf
        self.c_start = "ATG"
        self.c_stop = ["TAA", "TGA", "TAG"]
    
    def codons_start(self):
        return self.orf.count(self.c_start)
    
    def pos_codons_start(self):
        n = []
        for i in range(len(self.orf)):
            if c[i:i+3] == self.c_start:
                n.append(i)
        return n
    
    def pos_codons_stop(self):
        n = []
        for i in range(len(self.orf)):
             if c[i:i+3] in self.c_stop:
                n.append(i)
        return n

    def analyse_genome():
        pass

    def __str__(self) -> str:
        return self.orf
    