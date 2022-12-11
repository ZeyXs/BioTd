from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

# Exercice 1

def exo1():
    mes_seq = list(SeqIO.parse("resources/sequences_VIH.gb", "genbank"))
    print("------------------")
    for seq in mes_seq:
        print(f"La séquence {seq.name} a une longueur de {len(seq)} nucléotides.\nSa description est : {seq.description}.\nSon type moléculaire est : {seq.annotations['molecule_type']}\nSa date de dernière mise à jour dans la banque de données est : {seq.annotations['date']}.")
        print("------------------")

#exo1()
   
# Exercice 2

def exo2():
    mes_seq = list(SeqIO.parse("resources/sequence_Dengue.gb", "genbank"))
    for seq in mes_seq:
        SeqIO.write(seq, f"output/{seq.name}.gb", 'genbank')
        SeqIO.write(seq, f"output/{seq.name}.fasta", 'fasta')

#exo2()

# Exercice 3

def exo3(file: str):
    seq = list(SeqIO.parse(f"output/{file}", "genbank"))
    gene = len([elm for elm in seq[0].features if elm.type == "gene"])
    seq = list(SeqIO.parse(f"output/{file}", "genbank"))
    gene_list = [elm.qualifiers["gene"][0] + f" ({elm.location.start}..{elm.location.end})" for elm in seq[0].features if elm.type == "gene"]
    print(f"Organisme : {seq[0].annotations['organism']}\nNombre de gènes : {gene} \nListe de gènes : {' '.join(gene_list)}")
    
#exo3("NC_001802.gb")

# Exercice 4

def exo4(file: str):
    seq_file = list(SeqIO.parse(f"output/{file}", "genbank"))
    gene_list = [elm.qualifiers["gene"][0] for elm in seq_file[0].features if elm.type == "gene"]
    seq_list = [elm.qualifiers["translation"][0] for elm in seq_file[0].features if elm.type == "CDS"]
    id_list = [elm.qualifiers["protein_id"][0] for elm in seq_file[0].features if elm.type == "CDS"]
    with open(f"output/seq_proteiques.fasta", "w") as f:
        for i in range(len(gene_list)):
            f.write(f">{id_list[i]} ; {gene_list[i]} ; {len(seq_list[i])}\n")
            f.write('\n'.join([seq_list[i][i:i + 60] for i in range(0, len(seq_list[i]), 60)]) + "\n")
    
exo4("NC_001802.gb")
    

