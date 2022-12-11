from Bio.Seq import Seq
from Bio.Align.Applications import MafftCommandline
from Bio import SeqIO

# 2
def exo_2a(file: str):
    SeqIO.convert(f"resources/{file}", "genbank", "output/APBA1-nucleique.fasta", "fasta")

exo_2a("APBA1.gb")

def exo_2b(file: str):
    seq_file = list(SeqIO.parse(f"resources/{file}", "genbank"))
    output_file = open("output/APBA1-proteique.fasta", "w")
    i = 0
    for seq in seq_file:
        for elm in seq.features:
            if elm.type == "CDS":
                entete = "Séquence de base" if i == 0 else f"Séquence variant {i}"
                translation = elm.qualifiers["translation"][0]
                translation_formatted = '\n'.join([translation[i:i+60] for i in range(0, len(translation), 60)])
                output_file.write(f">APBA1 - {entete} - {seq.id} - " + elm.qualifiers["protein_id"][0] + "\n" + translation_formatted + "\n")
                i += 1
             
exo_2b("APBA1.gb")

def exo_2c(file: str):
    command = MafftCommandline(input="output/APBA1-nucleique.fasta")
    stdout, stderr = command()
    with open(f"output/aln-{file}.fasta", "w") as fd:
        fd.write(stdout)
    seq_file = list(SeqIO.parse(f"output/aln-{file}.fasta", "fasta"))
    fd = open(f"output/comparaison-{file}.fasta", "w")
    fd.write("position   lettre_base   lettre_variant1   lettre_variant2   lettre_variant3   lettre_variant4\n")
    
    for i in range(len(seq_file[0].seq)):
        fd.write(f"   {i+1} " + " "*(12-len(str(i+1))) + f"{seq_file[0].seq[i]}              {seq_file[1].seq[i]}                 {seq_file[2].seq[i]}                 {seq_file[3].seq[i]}                 {seq_file[4].seq[i]}\n")
        
    fd.close()
            
exo_2c("APBA1-proteique")

