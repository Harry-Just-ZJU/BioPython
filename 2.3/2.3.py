
# 处理成‘字符串’
from Bio.Seq import Seq
my_seq = Seq("GATCG")
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))

# 计数
from Bio.Seq import Seq
print("AAAA".count("AA"))
print(Seq("AAAA").count("AA"))

# 算GC含量
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
my_seq= Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC')
print(gc_fraction(my_seq))

# 切片
from Bio.Seq import Seq
my_seq= Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
print(my_seq[4:12])
print(my_seq[0::3])
print(my_seq[1::3])
print(my_seq[2::3])
print(my_seq[::-1])

# fasta格式
fasta_format_string = ">Name\n%s\n" % my_seq
print(fasta_format_string)

# 连接序列
from Bio.Seq import Seq
protein_seq= Seq("EVRNAK")
dna_seq =Seq("ACGT")
print(protein_seq + dna_seq)

from Bio.Seq import Seq
list_of_seqs =[Seq("ACGT"), Seq("AACC"), Seq("GGTT")]
concatenated = Seq("")
for s in list_of_seqs:
    concatenated += s
print(concatenated)

from Bio.Seq import Seq
contigs =[Seq("ATG"),Seq("ATCCCG"),Seq("TTGCA")]
spacer = Seq("N"*10)
print(spacer.join(contigs))

# 改大小写
from Bio.Seq import Seq
dna_seq =Seq("acgtACGT")
print(dna_seq.upper())

print("GTAC" in dna_seq)
print("GTAC" in dna_seq.upper())

# 转录
from Bio.Seq import Seq
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
messenger_rna = coding_dna.transcribe()
print(messenger_rna)
print(messenger_rna.back_transcribe())

# 翻译
from Bio.Seq import Seq
messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG")
print(messenger_rna.translate())
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print(coding_dna.translate())

print(coding_dna.translate(table="Vertebrate Mitochondrial"))
print(coding_dna.translate(table=2))
print(coding_dna.translate(to_stop=True))
print(coding_dna.translate(table=2, to_stop=True))
print(coding_dna.translate(table=2, stop_symbol="@"))

from Bio.Seq import Seq
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA" + \
           "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT" + \
           "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT" + \
           "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT" + \
           "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")

print(gene.translate(table="Bacterial"))
print(gene.translate(table="Bacterial", to_stop=True))
print(gene.translate(table="Bacterial", cds=True))

# 翻译表
from Bio.Data import CodonTable

standard_table= CodonTable.unambiguous_dna_by_name["Standard"]
mito_table= CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]

standard_table= CodonTable.unambiguous_dna_by_id[1]
mito_table= CodonTable.unambiguous_dna_by_id[2]

print(standard_table)
print(mito_table)

print(mito_table.stop_codons)
print(mito_table.start_codons)
print(mito_table.forward_table["ACG"])

