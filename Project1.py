# Step 1: Biological Seq Selection


# Step 2: Sequence Quality & Basic Analysis
from Bio import SeqIO
record = SeqIO.read("mouse1.fasta", "fasta")
print ("The id of the sequence is:", record.id)
print ("The name of the sequence is:", record.name)
print ("The description of the sequence is:", record.description)
print ("The length of the sequence is:", len(record))


# Step 3: Sequence Filtering & Validation
for aa in record:
    count = record.count(aa)
    print(aa, ":", count)
length = len(record)
if length < 100:
    print("Sequence rejected: too short for functional analysis.")
else:
    print("Sequence accepted based on length. Yes, this sequence is suitable for downstream analysis.")


# Step 4: Homology Search
from Bio.Blast import NCBIWWW
from Bio import SeqIO
record = SeqIO.read("mouse1.fasta", "fasta")
result_handle = NCBIWWW.qblast(
    program="blastp",
    database="nr",
    sequence=record.seq
)
with open("blast_result1.xml", "w") as out:
    out.write(result_handle.read())
print("BLAST search completed")

# This sequence analysis shows a highly significant hit to pre-mRNA processing protein PRP39. 
# The top hit shows 100% similairty with an E-value of 0.
# This shows that the hypothetical protein is highly reliable and biologically meaningful match.
# This protein is conserved across species and indicates it's role at a cellular level.

# Step 5: Functional Annotation
from Bio.Blast import NCBIXML
with open ("blast_result1.xml") as b:
    blast_record = NCBIXML.read(b)
first_alignment = blast_record.alignments[0]
print (first_alignment.title)
print (first_alignment.length)