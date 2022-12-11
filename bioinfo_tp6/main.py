from Bio import SeqIO
from Bio import Entrez
from Bio import pairwise2
from Bio.Align.Applications import MafftCommandline

# 1
def exo1():    
    # Requestgene_list
    print("Request...")
    Entrez.email = "basile.gauttron@etu.umontpellier.fr"
    request_hiv = Entrez.esearch(db="nucleotide", term="HIV1 [orgn] AND RefSeq OR HIV2 [orgn] AND RefSeq")
    result_hiv = Entrez.read(request_hiv)
    request_hiv.close()
    
    # Parsing
    fic_hiv = Entrez.efetch(db="nucleotide", id=result_hiv["IdList"], rettype="gb")
    
    # Output
    with open("output/HIV.gb", "w") as fd:
        fd.write(fic_hiv.read())


exo1()

# 2 
def exo2():
    # Retrieving protein
    seq_file = list(SeqIO.parse("output/HIV.gb", "genbank"))
    protein_id = [elm.qualifiers["protein_id"][0] for elm in seq_file[0].features if elm.type == "CDS" and elm.qualifiers["gene"][0] == "asp"][0]
    
    print(protein_id)
    # Request
    print("Request...")
    Entrez.email = "basile.gauttron@etu.umontpellier.fr"
    request_protein = Entrez.esearch(db="Protein", term=protein_id)
    result_protein = Entrez.read(request_protein)
    request_protein.close()
    
    # Parsing
    fic_protein = Entrez.efetch(db="Protein", id=result_protein["IdList"], rettype="gb")
    seq_protein = SeqIO.parse(fic_protein, "genbank")
    
    # Output
    SeqIO.write(seq_protein, "output/prot-asp.gb", "genbank")
    
#exo2()

# 3
def exo3():
    seq_file = list(SeqIO.parse("output/HIV.gb", "genbank"))
    protein_list = [elm.qualifiers["protein_id"][0] for elm in seq_file[1].features if elm.type == "CDS"]
    
    # Request
    print("Request...")
    Entrez.email = "basile.gauttron@etu.umontpellier.fr"
    search_proteins = ' '.join(protein_list)
    request_proteins = Entrez.esearch(db="Protein", term=search_proteins)
    result_proteins = Entrez.read(request_proteins)
    request_proteins.close()
    
    # Parsing
    fic_protein = Entrez.efetch(db="Protein", id=result_proteins["IdList"], rettype="gb")
    seq_protein = SeqIO.parse(fic_protein, "genbank")
    
    # Output
    SeqIO.write(seq_protein, "output/prot-HIV2.gb", "genbank")
    SeqIO.convert("output/prot-HIV2.gb", "genbank", "output/prot-HIV2.fasta", "fasta")
    
#exo3()

# 4
def exo4():
    dict_ = {}
    hiv1 = list(SeqIO.parse("output/prot-HIV1.fasta", "fasta"))
    hiv2 = list(SeqIO.parse("output/prot-HIV2.fasta", "fasta"))
    
    for i in range(len(hiv1)):
        dict_[i] = (0, 0)
        for j in range(len(hiv2)):
            pw = pairwise2.align.globalxx(hiv1[i].seq, hiv2[j].seq)
            print(pw[0].score)
            if dict_[i][0] < pw[0].score:
                dict_[i] = (pw[0].score, j)
        SeqIO.write([hiv1[i], hiv2[dict_[i][1]]], f"output/comparaison/comparaison_{i+1}.fasta", format="fasta")
                
#exo4()

# 5
def exo5():
    # Request
    print("Request...")
    Entrez.email = "basile.gauttron@etu.umontpellier.fr"
    request_apba = Entrez.esearch(db="Nucleotide", term="NM_001163.4 XM_005251968.4 XM_011518617.3 XM_017014670.2 XM_047423300.1")
    result_apba = Entrez.read(request_apba)
    request_apba.close()
        
    # Parsing
    fic_apba = Entrez.efetch(db="Nucleotide", id=result_apba["IdList"], rettype="gb")
    seq_apba = SeqIO.parse(fic_apba, "genbank")
    
    # Output
    SeqIO.write(seq_apba, "output/APBA1.gb", "genbank")
    SeqIO.convert("output/APBA1.gb", "genbank", "output/APBA1.fasta", "fasta")
    
    # Alignement
    mcl = MafftCommandline(input="output/APBA1.fasta")
    stdout, stderr = mcl()
    with open("output/aln-APBA1.fasta", "w") as fd:
        fd.write(stdout)
    
#exo5()