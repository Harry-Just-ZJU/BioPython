
# 处理数据
from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
print(my_seq)   

print(my_seq.complement()) 

print(my_seq.reverse_complement())   


# 解析文件
from Bio import SeqIO

# FASTA
for seq_record in SeqIO.parse('ls_orchid.fasta', 'fasta'):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

# GenBank
for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))