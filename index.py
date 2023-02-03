from Bio.Seq import Seq

seq1 = Seq("ACGTAGCTACGATCACAGCTA")

# print(len(seq1)) TAMANHO DA SEQUENCIA

# print(seq1.transcribe()) TRANSCREVE||RESULTADO = RNA

# print(seq1.translate()) TRADUZ||RESULTADO = PROTEINA

# print(seq1.reverse_complement()) REVERSO COMPLEMENTAR


# SEQUENCIAS SÃO É DIFERENTE DE UMA STRING MAS ELE PODE SER TRADADA COMO UMA
# print(seq1[0])

# for i,n in enumerate(seq1):
#     print(i,n)

# print(seq1.count('CTACGA')) BUSCA SUBSTRINGS

# print(seq1[::-1]) REVERTER SEQUENCIA

from Bio.SeqUtils import GC

# print(GC(seq1))

# seq2 = str(seq1) CONVERTENDO SEQUENCIAS EM STRING

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
# for i in SeqIO.parse("seq.fasta","fasta"):
#     print(i.id,":",i.seq)

# gravar = []

# for i in SeqIO.parse("novoseq.fasta","fasta"):
#     # print(i.id,":",i.description,i.seq,"\n")
#     nome = i.id
#     descricao = i.description
#     seq = i.seq

#     # coletando id do uniprot
#     id_uniprot = re.findall('\|.*\|',descricao)
#     id_uniprot = id_uniprot[0].replace("|","")
#     # print(id_uniprot)

#     # objeto de gravação
#     aux = SeqRecord(seq,id_uniprot,description='')
#     gravar.append(aux)

# # criando um novo arquivo fasta
# SeqIO.write(gravar,"seq2.fasta","fasta")

from Bio.PDB import *

pdb = PDBList()
pdb.retrieve_pdb_file('4MDP')

parser = MMCIFParser()

estrutura = parser.get_structure("4MDP",'md/4mdp.cif')

# estrutura -> modelos -> cadeias -> residuos -> átomos

# for modelo in estrutura:
    # print(modelo)
    # for cadeia in modelo:
        # print(cadeia)
        # for residuo in cadeia:
            # nome = residuo.get_resname()
            # if nome != "HOH":
                # print(nome,residuo.id[1])

            # for atomo in residuo:
                # print(atomo.id)
                # print(atomo.coord)
                # print(atomo)
                

# distancia euclidiana entra 1ys 475 e leu 468 - ca
R1 = estrutura[0]['A'][475]['CA']
R2 = estrutura[0]['A'][476]['CA']
# R2 = estrutura[0]['A'][468]['CA']

distancia = R1 - R2
print(distancia,'angstrons')