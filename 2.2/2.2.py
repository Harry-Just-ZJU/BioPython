
from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
print(my_seq)               # AGTACACTGGT

print(my_seq.complement())  # TCATGTGACCA

print(my_seq.reverse_complement())   # ACCAGTGTACT

from Bio import SeqIO
for seq_record in SeqIO.parse('ls_orchid.fasta', 'fasta'):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

for seq_record in SeqIO.parse("ls_orchid.gbk", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))