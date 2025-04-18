from Bio.Seq import Seq
my_seq = Seq("GATCG")
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))

from Bio.Seq import Seq
print("AAAA".count("AA"))
print(Seq("AAAA").count("AA"))

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
my_seq= Seq('GATCGATGGGCCTATATAGGATCGAAAATCGC')
print(gc_fraction(my_seq))

from Bio.Seq import Seq
my_seq= Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
print(my_seq[4:12])
print(my_seq[0::3])
print(my_seq[1::3])
print(my_seq[2::3])
print(my_seq[::-1])

fasta_format_string = ">Name\n%s\n" % my_seq
print(fasta_format_string)

from Bio.Seq import Seq
protein_seq= Seq("EVRNAK")
dna_seq =Seq("ACGT")
print(protein_seq+ dna_seq)

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

from Bio.Seq import Seq
dna_seq =Seq("acgtACGT")
print(dna_seq.upper())

print("GTAC" in dna_seq)
print("GTAC" in dna_seq.upper())

