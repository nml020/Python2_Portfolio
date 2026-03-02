# NatalieLeBlanc_Python2_Portfolio
This is the portfolio of the Python 2 codes that I learned through Winter 2025-2026 in the Python 2 class (BISC 4503).

## Sequence Objects
In this analysis, we used Biopython for DNA sequence creation, modification, transcription and translation, codon table visualization, etc. 

```python
# Import Seq type
from Bio.Seq import Seq
```


```python
# Create DNA sequence
my_seq = Seq("GATCG")
```


```python
# Sequence Loop
# Print index and base
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# Print sequence length
print(len(my_seq))
```

    5



```python
# Print first base
print(my_seq[0])
```

    G



```python
# Print last base
print(my_seq[4])
```

    G



```python
# Print third base
print(my_seq[2])
```

    T



```python
# Count occurrence of "AA"
Seq("AAAA").count("AA")
```




    2




```python
# Create longer DNA sequence
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
# Find sequence length
len(my_seq)
```




    32




```python
# Coung G bases
my_seq.count("G")
```




    9




```python
# Count percentage of "GC"
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    46.875




```python
# Import GC function
from Bio.SeqUtils import gc_fraction
```


```python
# Load DNA sequence again
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
# Calculate GC fraction on sequence
gc_fraction(my_seq)
```




    0.46875




```python
# Find segment of sequence
my_seq[4:12]
```




    Seq('GATGGGCC')




```python
# Start with first character and grab every third nucleotide at the end 
my_seq[0::3]
```




    Seq('GCTGTAGTAAG')




```python
# Reverse previous sequence
my_seq[1::3]
```




    Seq('AGGCATGCATC')




```python
# Find the nucleotide at position 2 
my_seq[2:3]
```




    Seq('T')




```python
# Reverse original sequence
my_seq[::-1]
```




    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')




```python
# State original sequence
str(my_seq)
```




    'GATCGATGGGCCTATATAGGATCGAAAATCGC'




```python
# Create FASTA string
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
# Print FASTA sequence
print(fasta_format_string)
```

    >Name
    GATCGATGGGCCTATATAGGATCGAAAATCGC
    



```python
# Create two sequences 
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
# Join sequences together (1 with 2)
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
# Join two sequences together (2 with 1)
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
# Create list with specific sequences
contigs = [Seq("ATG"), Seq("ATCCG"), Seq("TTGCA")]
```


```python
# Create spacer for sequence
spacer = Seq("N" *10)
```


```python
# Join sequences with spacer
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCGNNNNNNNNNNTTGCA')




```python
# Create sequence object
dna_seq = Seq("acgtACGT")


```


```python
# Sequence stored and ran 
dna_seq
```




    Seq('acgtACGT')




```python
# Convert sequence to uppercase
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# Convert sequence to lowercase
dna_seq.lower()
```




    Seq('acgtacgt')




```python
# Convert back to uppercase
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
# Try to find segment in sequence
"gtac" in dna_seq
```




    False




```python
# Try to find segment in sequence
"GTAC" in dna_seq
```




    False




```python
# Print sequence
dna_seq
```




    Seq('acgtACGT')




```python
# Convert sequence to uppercase
dna_seq = dna_seq.upper()
```


```python
# Find segment in new uppercase sequence
"GTAC" in dna_seq
```




    True




```python
# Re-print sequence
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
# Find sequence complement
my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')




```python
# Find reverse complement of sequence
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')




```python
# Protein sequence created
# Print sequence complement
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
# Create coding DNA sequence
coding_dna = Seq("ATGGCCATTGTAATGGCCGCTCAAAGGGTGCCCGATAG")
```


```python
# Print coding sequence
coding_dna
```




    Seq('ATGGCCATTGTAATGGCCGCTCAAAGGGTGCCCGATAG')




```python
# Find reverse complement of sequence
template_dna = coding_dna.reverse_complement()
```


```python
# Print reverse complememnt 
template_dna
```




    Seq('CTATCGGGCACCCTTTGAGCGGCCATTACAATGGCCAT')




```python
# Print coding sequence again
coding_dna
```




    Seq('ATGGCCATTGTAATGGCCGCTCAAAGGGTGCCCGATAG')




```python
# Transcription of DNA to RNA
messenger_rna = coding_dna.transcribe()
```


```python
# Print new RNA sequence
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGCCGCUCAAAGGGUGCCCGAUAG')




```python
# Find reverse complement od DNA and convert to RNA
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGCCGCUCAAAGGGUGCCCGAUAG')




```python
# Transcribe RNA back to DNA
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGCCGCTCAAAGGGTGCCCGATAG')




```python
# Transcribe sequence to RNA
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGCCGCUCAAAGGGUGCCCGAUAG')




```python
# Translate RNA sequence to protein sequence
messenger_rna.translate()
```





    Seq('MAIVMAAQRVPD')




```python
# Translate DNA using "Vertebrate Mitchondrial" code
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMAAQ*VPD')




```python
# Translate DNA using Table name
coding_dna.translate(table = 2)
```




    Seq('MAIVMAAQ*VPD')




```python
# Translation ends at stop codon
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMAAQRVPD')




```python
# Apply ending at stop codon to Table 2
coding_dna.translate(table = 2, to_stop=True)
```




    Seq('MAIVMAAQ')




```python
# Stop symbol changed
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMAAQ!VPD')




```python
# Insert bacterial sequence
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTCCCGTCACAATTACAGAGCGATCGTGATAATTCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
# Sequence given title
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLHFPWFWSLPWQHRLRKLRSRHNYRAIVIIRGYYWDGGHWRDHGWW...HR*')




```python
# Translation ends at stop codon 
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLHFPWFWSLPWQHRLRKLRSRHNYRAIVIIRGYYWDGGHWRDHGWW...HHR')




```python
# Check sequence
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLHFPWFWSLPWQHRLRKLRSRHNYRAIVIIRGYYWDGGHWRDHGWW...HHR')




```python
# Import codon tables
from Bio.Data import CodonTable
```


```python
# Load standard table
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
# Load "Vertebrate Mitochondrial" table
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
# Display standard table
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# Display mitochondrial table
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
# Print stop codons for mitoochondrial table
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
# Print start codons for mitochodnrial tables
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
# Create sequence 
seq = Seq("ACGT")
```


```python
# Confirm sequence
"ACGT" == seq1
```




    True




```python
# Confirm sequence
seq1 =="ACGT"
```




    True




```python
# Create undefined sequence
unknown_seq = Seq(None, 10)
```


```python
# Print undefined sequence
unknown_seq
```




    Seq(None, length=10)




```python
# Find sequence length
len(unknown_seq)
```




    10




```python
# Create new sequence with given length
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
# Find segment at specific location
seq[1000:1020]
```




    Seq(None, length=20)




```python
# Find segment at specific location
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
# Find segment and length at specific location
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
# Create sequence
seq = Seq("ACGT")
```


```python
# Create undefined sequence
undefined_seq = Seq(None, length = 10)
```


```python
# Combine sequence with undefined sequence
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
# Create sequence
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
# Import MutableSeq
from Bio.Seq import MutableSeq
```


```python
# Make sequence editable
mutable_seq = MutableSeq(my_seq)
```


```python
# Print sequence
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# Change bse
mutable_seq[5] = "C"
```


```python
# Print sequence
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# Remove base
mutable_seq.remove("T")
```


```python
# Print sequence
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
# Reverse sequence
mutable_seq.reverse()
```


```python
# Print sequence
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# Convert sequence to original
new_seq = Seq(mutable_seq)
```


```python
# Print sequence
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
# Import functions
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
# State string of sequence
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
# Find reverse complement of string
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
# Transcribe string
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
# Back transcribe string
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
# Translation of string to proteins
translate(my_string)
```




    'AVMGRWKGGRAAG*'


## Sequence Annotations
In this analysis, we used Biopython and learned how to take sequences with records, look at their information, see their annotations, rotate sequences, and generate reverse complements of them.

```python
# Import SeqRecord
from Bio.SeqRecord import SeqRecord
```


```python
# Import Seq
from Bio.Seq import Seq
```


```python
# Create sequence 
simple_seq = Seq("GATC")
```


```python
# Create SeqRecord
simple_seq_r = SeqRecord(simple_seq)
```


```python
# Print sequence
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# Set record ID
simple_seq_r.id = "AC12345"
```


```python
# Make sequence description
simple_seq_r.description = "Made up sequence for VDB Computational Biology Class"
```


```python
# Print description
print(simple_seq_r.description)
```

    Made up sequence for VDB Computational Biology Class



```python
# Print saved sequence
simple_seq_r.seq
```




    Seq('GATC')




```python
# Display record of sequence information
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for VDB Computational Biology Class', dbxrefs=[])




```python
# Add anotation
simple_seq_r.annotations["evidence"] = "None. This is just an example"
```


```python
# Print annotation
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example



```python
# Display record again
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for VDB Computational Biology Class', dbxrefs=[])




```python
# Add quality score
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
# Print letter annotations
print(simple_seq_r.letter_annotations)
```

    {'phred_quality': [40, 40, 38, 30]}



```python
# FASTA source website
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
# Import SeqIO
from Bio import SeqIO
```


```python
# Read FASTA file
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
# Print sequence record
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# Print sequence
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# Show record ID
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
# Show record description 
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# Show references
record.dbxrefs
```




    []




```python
# Show annotations
record.annotations
```




    {}




```python
# Show features
record.features
```




    []




```python
# Link to GenBank source 
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
# Read Genbank file
record = SeqIO.read("NC_005816.GB.txt", "genbank")
```


```python
# Print sequence record
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# Print sequence
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
# Show record ID
record.id
```




    'NC_005816.1'




```python
# Show record name
record.name
```




    'NC_005816'




```python
# Show record description
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# Show annotations
record.letter_annotations
```




    {}




```python
# Find annotation length 
len(record.annotations)
```




    13




```python
# Give source of annotations
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
# Show database cross references
record.dbxrefs
```




    ['Project:58037']




```python
# Show record length
len(record.features)
```




    41




```python
# Import SeqFeature
from Bio import SeqFeature
```


```python
# Create start position after base 5
start_pos = SeqFeature.AfterPosition(5)
```


```python
# Create end position between bases 9 and 10
end_pos = SeqFeature.BetweenPosition(10, left=9, right=10)
```


```python
# Create location with start and end positions
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
# Print location
print(my_location)
```

    [>5:(9^10)]



```python
# Give start position
my_location.start
```




    AfterPosition(5)




```python
# Give end location
my_location.end
```




    BetweenPosition(10, left=9, right=10)




```python
# Convert end position to integer
int(my_location.end)
```




    10




```python
# Convert start location to integer
int(my_location.start)
```




    5




```python
# Create location from 5 to 9 
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
# Print location
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
# Import SeqRecord
from Bio.SeqRecord import SeqRecord
```


```python
# Print record information
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGSPSETHLDSMVGQALFGD"
...         "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIGDHLREVGLTFHLLKDVPGLISK"
...         "NIEKSLKEDFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
...         "SSAC"),
                  id="gi|14150838|gn|AAK54648.1|AF376133_1",
            description="chalcone synthase [Cucumis sativus]",
...
                  )
```


```python
# Print record in FASTA format 
print(record.format("fasta"))
```

    >gi|14150838|gn|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGSPSETHLDSMVGQALFG
    DGAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIGDHLREVGLTFHLLKDVPGLIS
    KNIEKSLKEDFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGN
    MSSAC
    



```python
# Import SeqIO
from Bio import SeqIO
```


```python
# Read FASTA file
record = SeqIO.read("NC_005816.GB.txt", "genbank")
```


```python
# Print sequence record
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# Print record length 
len(record)
```




    9609




```python
# Print number of annotated features
len(record.features)
```




    41




```python
# Print feature 20 
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
# Print feature 21
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# Print subsequence from position 4300-4800
sub_record = record[4300:4800]
```


```python
# Find length of subsequence
len(sub_record)
```




    500




```python
# Give length of record features
len(sub_record.features)
```




    2




```python
# Print sub-record feature 0
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
# Print sub-record feature 1
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
# Print all record feature 0
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
# Print all record feature 1
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
# Give record annotations 
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
# Show database cross references
sub_record.dbxrefs
```




    []




```python
# Add typology information
sub_record.annotations["topology"] = "linear"
```


```python
# Show typology annotations
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
# Give record ID
sub_record.id
```




    'NC_005816.1'




```python
# Give record name
sub_record.name
```




    'NC_005816'




```python
# Print description of record
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
# Edit record description
sub_record.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'
```


```python
# Print new description
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
# Print first 200 characters of genbank
print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...



```python
# Print record from Genbank
record = SeqIO.read("NC_005816.GB.txt", "genbank")
```


```python
# Print record
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# Give record length
len(record)
```




    9609




```python
# Give length of annotations
len(record.features)
```




    41




```python
# Show database cross references
record.dbxrefs
```




    ['Project:58037']




```python
# Print annotation names
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# Shift record putting first 2000 bases at the end 
shifted = record[2000:] + record[:2000]
```


```python
# Print shifted sequence
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
# Give length of shifted sequence 
len(shifted)
```




    9609




```python
# Count features in shifted sequence
len(shifted.features)
```




    40




```python
# Print annotation names
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
# Show database cross references
shifted.dbxrefs
```




    []




```python
# Show database cross references for shifted sequences
shifted.dbxrefs = record.dbxrefs[:]
```


```python
# Show database cross references
shifted.dbxrefs
```




    ['Project:58037']




```python
# Copy annotations from original sequence
shifted.annotations = record.annotations.copy()
```


```python
# Print updated annotation key
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
# Print original record
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
# Print: ID, length, features, dbxrefs, annotations
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
# Create sequence reverse complement 
rc = record.reverse_complement(id = "Testing")
```


```python
# Print reverse complement
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
# Print summary of reverse complement record
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0



## Sequence I/O
In this analysis, we imported code from GenBank and FASTA files while creating reverse complements, converting, and writing resulting sequences.
```python
# Import SeqIO
from Bio import SeqIO
```


```python
# Create loop for FASTA file
# Print sequence ID
# Print DNA sequence and length
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# Create loop for GenBank file
# Print record ID
# Print sequence and length
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    Z78533.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    Z78532.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    Z78531.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    Z78530.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    Z78529.1
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    Z78527.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    Z78526.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    Z78525.1
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    Z78524.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    Z78523.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    Z78522.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    Z78521.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    Z78520.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    Z78519.1
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    Z78518.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    Z78517.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    Z78516.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    Z78515.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    Z78514.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    Z78513.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    Z78512.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    Z78511.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    Z78510.1
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    Z78509.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    Z78508.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    Z78507.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    Z78506.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    Z78505.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    Z78504.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    Z78503.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    Z78502.1
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    Z78501.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    Z78500.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    Z78499.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    Z78498.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    Z78497.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    Z78496.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    Z78495.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    Z78494.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    Z78493.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    Z78492.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    Z78491.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    Z78490.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    Z78489.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    Z78488.1
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    Z78487.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78486.1
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    Z78485.1
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    Z78484.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    Z78483.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    Z78482.1
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    Z78481.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    Z78480.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    Z78479.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    Z78478.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    Z78477.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    Z78476.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    Z78475.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    Z78474.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    Z78473.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    Z78472.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    Z78471.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78470.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    Z78469.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    Z78468.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    Z78467.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    Z78466.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    Z78465.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    Z78464.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    Z78463.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    Z78462.1
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    Z78461.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    Z78460.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    Z78459.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    Z78458.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    Z78457.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    Z78456.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    Z78455.1
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    Z78454.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    Z78453.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    Z78452.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    Z78451.1
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    Z78450.1
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    Z78449.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    Z78448.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    Z78447.1
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    Z78446.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    Z78445.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    Z78444.1
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    Z78443.1
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    Z78442.1
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    Z78440.1
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# Create list of identifiers for GenBank file
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
# List identifiers
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
# Create iterator for FASTA file
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
# Retrieve first sequence
first_record = next(record_iterator)
```


```python
# Print first sequence ID
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
# Print sequence description
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
# Retrieve second sequence
second_record = next(record_iterator)
```


```python
# Print second sequence ID
print(second_record.id)
```

    gi|2765657|emb|Z78532.1|CCZ78532



```python
# Print sequence description
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
# Store GenBank files in a list
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
# Print amount of records
print("found %i records" % len(records))
```

    found 94 records



```python
# Print record label
print("The last record")
last_record = records[-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
# Print record label
print("The first record")
first_record = records[-3]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703



```python
# Create new GenBank iterator
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
# Create first record
first_record = next(record_iterator)
```


```python
# Print record
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
# Print record with annotations
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
# Print annotation keys
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
# Print annotation field values
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
# Print "source" annotations
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
# Print organism name
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
# Create list to store species names
all_species = []
```


```python
# Add all organisms to species storange from GenBank
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
# Print all species names
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
# Same as code above
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
]
```


```python
# Print all species
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
# Create empty list for species storage
all_species = []
```


```python
# Add all organisms to species storange from FASTA
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
# Print species names
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
# Create iterator for FASTA 
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
# Retrieve first sequence record
first_record = next(record_iterator)
```


```python
# Print first record ID
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
# Edit ID
first_record.id = "new_id"
```


```python
# Confirm ID
first_record.id
```




    'new_id'




```python
# Print first record
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new_id', name='gi|2765658|emb|Z78533.1|CIZ78533', description='gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
# Update record description
first_record.description = first_record.id + " " + "Mutations induced randomly"
```


```python
# Print edited record
print(first_record.format("fasta"[:200]))
```

    >new_id Mutations induced randomly
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    



```python
# Import Seq
from Bio.Seq import Seq
```


```python
# Import SeqRecord
from Bio.SeqRecord import SeqRecord
```


```python
# Create sequence of proteins
# Assign sequence identifier
# Add sequence desccription
rec1 = SeqRecord(
Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFDG"
    "GAGAVIVGSDPDLSVERPLYELVWTGATLLPSDEGAIDGHLREVGLTFHLLKDVPGLISK"
    "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGKEEKMRATREVLSEYGNM"
    "SSAC",
   ),
    id = "gi|14150838|gb|AAK54658.1|AF376133_1",
    description = "chalcone synthase [Cucumis sativus]",
)
```


```python
# Create sequence of proteins
# Assign sequence identifier
# Add sequence desccription
rec2 = SeqRecord(
Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQ"
       "DMVVEIPKLGKEAAVAIKEQGQ",
   ),
    id = "gi | 13919613 |gb| AAK33142.1",
    description = "chalcone synthase [Fragaria vesca subsp. bracteata]",

)
```


```python
# Create sequence of proteins
# Assign sequence identifier
# Add sequence desccription
rec3 = SeqRecord(
Seq("MVTVEEFRRAQCAEGPATVMMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFLRMC"
    "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
    "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN"
    "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV"
    "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW"
    "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT"
    "TGEGLEWGVLFGFGPGLTVETVVLHSVAT",
   ),
    id="gi|13925890|gb|AAK49457.1|",
    description="chalcone synthase [Nicotiana tabacum]",
)
```


```python
# Create list of sequences
my_records = [rec1, rec2, rec3]
```


```python
# Import SeqIO
# FASTA format
from Bio import SeqIO
SeqIO.write(my_records, "my_example.fasta", "fasta")
```




    3




```python
# Retrieve GenBank file
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
# Files converted to FASTA
count = SeqIO.write(records, "my_example.fasta", "fasta")
print("Coverted %i records" % count)
```

    Coverted 94 records



```python
# Create loop through Genbank files
# Print record ID
# Print reverse complement
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement())
```

    Z78533.1
    GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGTCAACCGCCCTCCATTGGGGTTCGGAAGGGTTCTCCTGCCTGTCCGACAAGCACGACAACATGGGGGATTCGCACGGCAGCTGCTGCCAGCATCCGTCCACCTCTTGCCGGGTTCCGGGCCATCAAAACACCAGCTCTTGGACCCGCCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACCACGCCGGCCTGGCTGTATGCCGGGCAAGCATTGGCAGGAGAGAGACGAAGCGCGACGCCCAAGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCAAGATTCCCGTTTGCGAGAGTCATCAAAATTCATTGGGCACGCGACAGCACGCCGCCGCTCCGGGTTTTGGGGAAGACAATGCCATTCGCCGGTGATGCTTTCATATGGCTTGGCGCCCAAACTGCGCCGGGCTAGAGGTTCAAACCCGCCATGGACGCTCCCGAGGCGGCCCAACAACAAATCAGGGTCACCACGGGAGCAATGCCCCCGGTGAGCTGAGTACACCGGTCCTCCGGATTCACTCGATCGTTTATTCCACGGTCTCATCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78532.1
    GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAATCAACTGCCCAATGGTTATTTTAGCTCCATTGGGGTTCAATTAGGTTCTTGTGTAGGTTCGAAAAAATACAACAACATGGGGGATTCAAATAGCAGCCTTATGACTGTTAGCATTCTCCACCTCGTGCCACATTCCTACCCATCAAAGCAACAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCATTCACATCCGTATAATGCCAGCTTAGCGATATGCCAAGCAAGCATTGGTAGGAGAGACGCAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTAGCTGCGTTCTTCATCGATGTTAGAGCCAAGATATCCGTTGCGAGAGTCATAAAAATTCACTGGCATGCTCAGTAGCATACTGCCCCTCTGATTTTTTCTGACAATAATGTCATTCATCAGTGATGCTTTATATGACTTGGCGCAAACTGCACCGTACTAAAGTTCAAACCTGCCATGAAAGCTCTTGAGGAGGCCCAACAACAAAGCAGGGTCACGACAAAAGCAGTGCCACGACGAGCTGAGTTACCACAGGTCCTCCAGATTCACTCGATCATATATTCTGTTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78531.1
    TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCCCGAGGGCGGTTTTGCGTCCACTGGGGGTTCAAACAGGTTCTTCTGTAGGCCTGACGAGTACGACAACACGGGGGATTCGAACGGCAGCCTTGCGGCTGCCAGCATCCTCCACCTACTGCCAAGTTCCGGGCCATCAGAACACCGATGCTTAGACCCGCCGCACCTAGGCGCAAGGGGCCAATCTTTCACGTCCTCACAATGACAGACTAGCCGCATGCCAATCAAGCATTATCAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAGATATCCCGTTGCCGAGAGTCGTCATAAATTCATTGGCACGCGCAACAGACGCCGCCCCTCCGAGTTTTTCTTGACAAAAATGCCATCCATCGGTGACGCTCCATATGACTTGGCGCAAACTGCGCCGCGCTAGACGTTCAAACCCGCCATGAAAGTTCCCGAGGCGGCCCGGTCGCAAACCGGGTTCACCACGAGAGCAAAGCCACGGTGAGCCGTGTAACCACGGGTCCTCCGGATTCACTCGATCGTATGTTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78530.1
    ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAACATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGCGTTCTTCATCGATGCAAGAGTCAAGATATCCGTCTGCCGAGAGTCATCAAAATTCATTGGAATGAACATTAACACACCACACCTCCGACGTTGTGTGGGGCAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGGACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGATATCCCTAGCGAGCCAAATTACCACAAGTCCTCCAGATTCACTCAATCGTTTATTATGTTGTTTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78529.1
    TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAATATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGGGATTCTTCATCGATGCAAGAGCCAAGATATCCCTCGGCGAGAGTCATCAACAATTCATTGGCATGACCAATAACACACCACACCTCCGACTTTTTTGACAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGCACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGAAATCCCTAGCGAGCCAAATAACCACAAGTCCTCCAGATTCACTCAATCGTATATTCTGCTGTCTCAACAATGTCCTTCGGCAGCTCGCCGT
    Z78527.1
    GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACAATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACCCAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAAGAGCCAGATATCCGTTGGCCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAACACACTGCCACTCCAACTTTTGGTGACAGTAGTGCCATTTGCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGCCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGCGAGCCGAGTAACCACAGGTCATCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78526.1
    ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGGCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAGCACACTGCCACTCCAACTTTTGGTGACAATAATGCCATTGTTCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGGCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGGGAGCCGAGTAACCACAGGTCCTCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78525.1
    TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATCCAAGGGTCATTTTGCCTCCACTGGGGTTCAAACAAGTTATTATATAGGTCCGATAAATACGATAACATGGGGGATTTAAATGACAACCTTGTGGCAGCCAATGTCCTCCACCTCCTGCCAAGTTTCGGGCCATCAAAACATTCCTTAGACCCAACGCGCCAAGGCACAAGGGGCCAATCTTTCGCATCCGCACAATGCAAGCCTAGATATATGGACAGGCATTGGCAAGAGAGACACAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCCGATGCAAGAGCCAAGATATCCGCTTGTCGAGAGTCATCAAAATTATATTCGGCATGCACAGGAGGACACCGCCCCTCCAACTTTTTTGACAATAATGTCATTTGTCGGTGATGCTTCATATGACTTGGCACAAACTGTGCCGTGCTAGAGTTTGAAACCTGCCATGAAAGCTCCCGAGGAGGCCCAACGACAAATTTGGGTCACCACAAAAGCAAAGCCCTCGGCAAGCCGAATAACCACAGGTCCTCCGGATTCACTCGATGTATATTCTGCTATCTCAACA
    Z78524.1
    GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAACCGCCCAAGGGTCATTTTTCCTTCATTGGGGTTCAAACAAGTTCTTTTATAGGTTCAACAAATACGACAACATGGGGAATTCAAATGGCAGCCTTGCAGCTACCAACATTCTCCACCTCCTGCCGAGTTCCGAGCTATCAAAATACTAATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATATTTCATCCGCACAATGCTAACCTAGCTATATGCTGGGCAAGCATTGGNAGGAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGGATGTAAAGAGCCAAGTTCCGTTGCGAGAGTCATCAAAAAATTCATTAGCATGCACAACAGCATGCCGCCCCTCCGACTTTTTTAACAACAATGCCATTCATCAGGGATGCTTATATGACTTGGCACAAACTGCGCCGTGCTAAAGGTTTGAACCCGCCATGAAAGCTCCTAAGGAGGCCCAATTAAGGTCACCACAAAAGCAAAGCACTGACGAGCCAAGTAACCACATGTCCTCCATATTCACTCAATCATATATTCTACTATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78523.1
    CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGGTTCAAACGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCGCATCCACGCAATGCAAACCAATCTATATGATAGGAAAGCATTGATAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACATATGTCGAGAGTCAATAAAAAATTCATTATGCATGCATGCAAGAGCACGCTTCCACTCTAGCTTTTGGGACAATAATGTCATTCCTCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCACACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAGTTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTGCCGCAGGTTCACCTACGGAAACCTGGTTACG
    Z78522.1
    CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCCCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACGCAATGCAAACCTGTCTATATGACGTATAACCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTTCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCGATCACACCACGATATGTCGAGAGTCATAAAAAATTCATTTGCATGTATGCAATAGCACGCTTCCACTCCAACTTTTTGGACAATAATGTCATTCATCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGGGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78521.1
    GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAGTTGTTTGCCTCCATTGGGATTCAAACATGTTTTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGTACCTAGGCACAAGGGGCCATTCTTTCACATCCACGCAATGCAAACCTATCTATATGACATGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTTGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAATTCATTTGCATGCATGCAATAGCACGCTTCCACTCCAACTTTTTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGGCAAATTAGGGTCACCGCCAAAGCCAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78520.1
    AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCCAGGTCCGACAAACACGACAATATGGGGGATTCACATGACAACCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCGATATGACAGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCACATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTTCGCTGCGTTCTTCATCCGATGCAAGAGCCAAGATATCCCTTGTCGAGAGTCATAAAATATTCTATTGCATGCATGCAATAGCACGCTTCTACTCCAACTTTTTGGACAATAATGTCATTCATCGGTGATGATTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78519.1
    TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCAACAAACGCGACAATATGGGGGATTCAAATGATAGCCTTATATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACATTAAGTTCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCTATATGATGGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAAATTCATTTGCATGCATGCAGTAGCACGCTTCCACTCCAACTTTTGGACAATAATGTCATTCATCGGCAATGCTTCATATGACTTGGTGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCAGATTCACTCGATCATAT
    Z78518.1
    GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCAAGGGTTGTTTTGCCTCCATAGGGGTTCAAATAGGTTCTTCTGTGGGTCCGGCAAATACGACAACATGTGGGATTTGAACGACAGCCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTTGGGCCATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGAGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGTATGCGTAACACACACCGTCCATCCGAGTTTGGGGGACAATAATGTAATTCGTCGGCGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGCTCAAACCCGCCATAAAAGCTCTCGAGGATGCCCAACTACAAATCAGGGTCACCACAAAAGTTAAGCCTTCGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCGAATATTCTACTATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78517.1
    GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAATGACAGTCTTGCGACTATCAACATGATCCACCTCCTGCCGGGGTTCGGGCCACCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGTCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGGAAGAGCAAGATATCCGTTGGCCGAGAGTCATCAAAAATTTATTGGCATGCACAACAGCACACCGCCCATCCGACTTTTGTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCAAACTGCACCGTGATAGAGGTTCAAACCCGCCATAAAATCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCTGGATTCACTCGATCGTATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78516.1
    TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCAATCACTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAACGATAGTCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTCGGGACATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTAAGATATCCGTTGACGAGAGTCATCAAAAAATTATTGGTATGCACAATAGCACACCGCCCATCCGACTTTGGTGACAATAATGTCATTCGTCGGTGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGTTCAAACCCGCCATAAAAGCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCATATATTATACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78515.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAATCAACCGCCCATGGGTTGGGTTGATGTGGCACCAATGGGGTTCCGAAGGGTCCTACGAATTATGATGGCCCATCAAATCGGACAACATACGGCATTCGCACAACAGCTTTGGCTCAGCCATTGCCCATCCACCTCTTGTCGAATTTTGGACCATCAAAAGCCCGAGGTCTTAGACCCATTGAACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGACCGAACGGTATTCTTGGACAACCATTGGTGGGAGAGACACAGTACACAACGCCCAGGCAGGCGTGCCCTTAGCCTGACGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTAAAAATGTAACTTGGGGGGCACGAATCAAGTGCCACCCCCACCGGCTACAAGGGAACACCCCTATGCCGATTCTTGTCTCAAAAGACTTGGCGCGAAGCTGCGCCGGGCTAGAAGTTTGATCGCCTTCAAGAACACTCCCAAGGAGGCCCGATGGCAGATCAAAGGCCACCACAGTAGCGAAAGCCCCTGGGAGATCGAAGTACACCGGTCCTCCAGATTAACTCGATCGTTTATTCTACGGTCTCAGCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78514.1
    TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCCATGGGTTGGGTCGACGTGCCTCCGACGGGGTTCGAAAGGGATCATTCCGGCCCATCAAATACGACAACATAGGGCATTCGCACAACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGACGAGAGTTGTCAGACGAATCTTAAGGGGGCACAGCAGCACGCCGGCCCTCCCCGGTCCATCCATGGAGGACGAAGACCCTTGTCGGTTCTATCTCATTTGACTTGGAGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATAGAGAACTCTCCCAAGGAGGTTCGATGGGAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTATCTCGATCGTATATTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78513.1
    CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGGTTGGGTCGACGTGCCTCCGAAGGGGTTCGAAAGGGATCGTTCTGGCACATCAAATACGACAACATAGGGCATTCGCACAAAAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACAGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCGAAAATTCTTTCGGGCACAGCAGCATGCCGGCCTCCCCGGCTCCATGGAGGACGATGCCCTTGCCGGTTCTATCTCATTTGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATGGAGAAATCTCCCAATGAGGCTCGATGGCAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTAACTCGATCGTGTATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78512.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAGCTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCGAGAGTTGTCAAAAAATAATTTCATTGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGGGGACGATGCCCTCCGCCGGTTCCATCTCATTTGATTTGGGGCGAAACTGCGCCGGGCCAAGGGTTTCATTTGCCATCAAGAATTCCTCCCAGAGGGTTCGACGGAAATTCACGGTCACCACAGTAGCCAAGGCCCCTGGGAGACCAAACTACACCGGCCCTCCGGATTAATTCGATCGTTTTTTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78511.1
    TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGCTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTTCATCGGATGCAAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAAAAATTCATGGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAATTCTCCCAAGGAGGTTCGACGGAAATTCACGGCCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78510.1
    TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAAAACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAATCTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCCATCCGATGCAAAGAGGCGAGATATCCGTTGCGAGAGTTGTCAAAAAATCCGAGGGGGGAACGGAAGAATGCCGGCCCTCCCCGGTTCCATGGGGGGGGACGCCCCCTGCCGGTCCTATCTCATTTGATTGGGGGGGAAATTGCGCCGGGCCAAGGGTCCATTGGCCACCAAGAATTCTCCCAGGGGGGTTCGATGGAAATTCACGGCCACCAAGGGGGGAAAGCCCCTGGGGAGACCAAATTAAACCGGTCCTCCGGTTTAATTCGATCGTTTTTTCGGGGGCTTAAAAAGGAATCCTCCCGAAGGTCACCTCGGAACCCTGGTTAG
    Z78509.1
    TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACCCCCCATTGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAAAACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAACCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCTGCTTATCACATTTAAACCTTGTTACGTCCGTTGCCGAGAGTTGTCAGAAAAATTCATGGGGGGGCACGGAGCATGCCGGCCCTCCCCGCTCCATGGAGGACGACGCCCTCCGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78508.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAATTCATTGGGGGCACGGCAGCATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78507.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGTGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAAATTCATTGGGGGACGGAAGCACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGAAATTCACGGTCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78506.1
    TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTATGAGCTGAGGGCTCCTTGAACGAGAGTTGTCAGCCCGATTCAAAGGGGGTACGGCGGCATGCCGGTCCTCCCCGGTTCCATGGAGGACGAAGCCCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78505.1
    AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGACCGAACAGTATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTGACGAGAGTTGTCACAAAGATTCATTGGGGGTATGGTGGTATGCCGGGCCTCCCCGGTTCCATGGAGGACGAAGACCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCACGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAATATCATGGTCACCGCAGTAGCCACCGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78504.1
    TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGTTTTGTCGCAGCCACCGTCCGTCCACCTCTTGTCGGATTTGGTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGTCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAGATTCATTGGGGGCACGGCGGGCATGCCGGCCCTCCCCGGCTCCATGGAGGGACGATGCCCTCTGACGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78503.1
    TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAAAGAGCTGAGACTCCGTTGCCGAGAGTTGTCAAAATAATTCATTGGGGGTACGGTAGCATGCCGGCCCTCCCCGGTTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGTAAATCACGGACACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAAGGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78502.1
    GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTACGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTTCGCTGCGTCCTTCCATCGATGCAAGAGCCGAGATTCCGGCCGAGAGTTGTCAATATAAAATTCATTGGGGGGACGGTAGTATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGACCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGTCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGGGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78501.1
    TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCTCGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACAACACTTATCACAATTTCGCTGCGTCCTTCCATCCGATGCAAGGAGCCGAGATTTCCCGTTTGGCCGAGAGTTGCCCAAAAAAAAATTCATTGGGGGCACGGCAACACGCCGGCCCTCCCCGGCTCCATGGAGGACGATTCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTTCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78500.1
    CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCGAGATTCCCGTTGGCCGAGAAGTTTTCCCAAAGAAAAATTCATTGGGGGCACGGCAGGACGCCGGCCCTCCCCGGCTCCATGGAGGGCGATGCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCAATTATTTTGCGGTCTCAACAATGAGCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78499.1
    GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACACTTATCGCATTTCGCTGCGTCCTCCATCCGATGCAAAGAGCCGAGATATCCGTTGGCTGAGAGTTGTCAAAAAAATAATTGGGAGAGGTCAAGGCCACATGCCGTCCCCTCCTCGAATTCAACGAAGGGGGTATGTCCTTCACCATTTATGTGTCATATGACTTGGTGCAAAACTGCGACGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGAAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGCGATCTCAACAATGATCCCTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78498.1
    GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAACTGATCACGGGTTGATGTGCCTCCAGTGGGGTTCACAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAGTTTTGTTGTAATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGTACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCCAAAAATAATTTGGAGAGGGGGAGTGATTGTAGGGTCAAGGGCAAAATGCCGCCCCCTCCTTCGCCATTATTGTGTCATATGACTTGGGGCAAAACTGCCCCGGGCAAGGGTAAGATCGCCGGCAAGAAAACTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCACCCATGGGTGACCGGAGTAAACTGATCCTCTGGAGTAACTCGATCAATTATTATGTGATCTCAACAATGACCTTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78497.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTGTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCGCGTTGCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACACCCGCACCGGCTGAACAGTATGCCTGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAGAGCTGAGATCTCCGTTGCTGAGAGTTGTTTCACAAAAAATAATTTGGAGAGGGGGAGGGAGTGTAGGCTCAAGGCCACATGCCCTTCACCTCCTTGAATTCACAAAGGGCCGTGCCCTTCACCATTTATGTGTCATGTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCGCCCATGGGCGACCAAAGTAAACTGATCCTCTGGATTCACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78496.1
    GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGGTGCGTCCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTCAAAAAATTAATTGGAGAGGTCAAGGCCACATGCCGGCCCCTCCTCGAATTCACGAAGGGATATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGCAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78495.1
    CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTGGATTATGCCGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAACAGTGTTGTTGTAGCCTGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCTGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTGAGGGAGAGAGATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTCATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGGATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTCAGAAAAATACTTTGGAGAGGGGGAGAGCGAGTGCAGTGTCAAGGCCACATGCCGCCCCCCCTCCTTGAATACACGAAGGGCTATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACCGCGCCGGACAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCTCCATGGCAACTCTAGGTCACCGCAATAGCAAGCGCCCATGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78494.1
    CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTACGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGACAGGCGTTCCCTTGGCCTGTTGGCCTCGGGCCCAACTTGCGTTCAAAGACTCGATGTTAACGGGATTCTGCAATTCACACCTGAGATATCCGAAGTTGAGAGTGTTACCTCAATAATGGGGAGTTGGTCAAGGCAGTATGCCGTCCCCCTTCATTAATTATGGGTCATTTGATTTGGCGCATTACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGTTCCCAAGGAGGCCCGATGGTAAATCTAGGTCACTGCAACAGTAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGACCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78493.1
    GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCACACCGGATGACCAGTACGTTTGGACAGCCTCATTAAGGGAGAGATATGAGTCGCAATGTCCAGGCAGGCGTCCCCTTGGGCTGATGGCCTCGCGCGCAACTTGCGTTCAAAGATTCGATGGTTCACGGGATTCTGCAAGGCACACCATTTATCGAATCTCGCTGCGTTCCTTCATCGATGCATGAGCCGAGATATCCGTTGCTGAGAGTTGTTACCTAAATACTTGGGGGCTGGTCAAGGAAGAATGCCGCCCCCCTTCACCACTTATGTAAAATTTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGCTCCCAGGGAGGCCCGATGGAAATTCTAGGTCACTGCAACAGAAAATGCCCATGGGTGACCAAAGTAAACCGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78492.1
    TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCTGATTAAGCTGGCCCATGTAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAATTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCGTGCAGCTCTTAGACCCCCCGTCCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTCCACGGGATTCTGCAAATCACACCAGTATCGCATTTCGCTGGGTTCTACATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGAGTGGGTCAAGGCAGTACGCTCCCCCCCTTCACCAATTATGTGACATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78491.1
    GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTCCAAAGGGTCTTCTGATTATGCTGGTCCATGTAACACAACAACCCGGGGCATTCGCACAACATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAATATAATTTGGAGCGGGTCAAGGCAGCACGCCTCCCCCCAACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCTAAGGAGGCCCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78490.1
    TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGACGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAAATTATGCTGGCCCATCTAATACGACAACGCGGGGCATTCGCACAACACTTGTTGCAGAACAGTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCGTTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCGAGGCAGCATGCCGCCCCCTTCACCAATTATGTGCCATATGACTTGGCGCAAAACTGCGCCGGACTAGGGTTCAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78489.1
    GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCATGGGGTTCAAAAGGGTCTTCAGATATGCTGGGCCCATCTAATACGACAACTCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAATCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACCCACATCCCCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATCTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTCGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCTCGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78488.1
    AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGGGTTGATGGGCCTCCGATGGGGTTCAAAAGGGTCTTCAGATTATGATGGCCCATCTAATACGACAACCCGGGGGATCCGCACAACAGTTTGGTTGTGGCGTTGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGGCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTCCATCCCGTTGATGAGAGTTGTTAGAAAATAATTAGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCACCAATTATGTGTCGTATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTAAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAACAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTGCGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACAG
    Z78487.1
    TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTACGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGCATTTAGGACCATCGATAGCCCATGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCGAACTCACATCCGCACCGGCTGAACAATATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTAAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGAGTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78486.1
    TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGGCACTGCAATAGCAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATAATCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGATTCACCTACGGAAACCTCGTGACG
    Z78485.1
    TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAGGGGTCTTTAGATTGTGCTGGCCCGTCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGGCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGCCCAAACTCACATCCGCACCGGCTGCACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCCCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGATATCTCAACAATGATCCATCCGCAGATTCACCTTCGGACACCAGGTTCAG
    Z78484.1
    AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTGTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGGGAGGGTAGGAAACATGCCGCCCCTTCACAATTATGTGTCATATGACTTGGGCCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGGAGGCTCGATGGTAAACTCTCGGTCACTGCAATAGCCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCCCAGGTTCACCTACGGAAACCTTGTTACG
    Z78483.1
    TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGATATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATTGTTCACGGGATTCTGCAATTCACACCATTTATGATTTCGGTTTGGTCTTCATCCATGCAAGAGGCCAGATTTCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78482.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACGCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGNTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGAAGAGGCGAGATATCCGTTGCNGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTACGGTTTAGATCGCCAGCAAGAAAGCTCCCAGGAGGCTCGATGGCAAATCTCGGTCACTGCAGTAGA
    Z78481.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTCCGTTGCTGAGAGTTGTTAAAAAAATGATTTGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78480.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCGCACAACATTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78479.1
    ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTGATCTTGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAGACTCGATGGTTCACGGGATCTCTGTATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGGCAAGGCAAGGCAGCATCCCTCCCCTTCCAGTAATGTCTCATATAAATTGGCGCAAATCTGCGCCGGGCAAGGGTTTAGTTCGGCTCGATGGCAAATCTAGGTCACTGAAACAGCAAATGCCCATGGGTGACCAAAGTAACCTGATCCTCCTGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78478.1
    GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAACTGATAATGGGTTGATGGGCCTCCAATGGGGGTCAAAAGGGTCTTTAGATTATGCTGACCCATCTAAACGACAACCGGGGCATTCGCACAACAGTTCTGTTGTCACATAGCTCGTCCACCTCTTGCCGTATTTAGGACCATCCAACCCATACAGCTCTTAGACCCCCGTACCAAAGGACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGCAGCCTCATTAAGGGAGAGATATGTCTCACAATTTCCAGGCAGGCGTCCCCTTGACCTGATGGCCTCGGACGCAACATGCGTTCAAAACTCGATGTTCAGGGATCTGCTATCACGCCATTATCAATATTTGGGGGAGCCATGGCAGCATGCCGCCCCTTCCAGTTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAACTCTAGGCCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCCCCAGATTAACTCGATCAATTATTATGTGATCTCAACACTGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78477.1
    GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCCGCAACTTGCGTTCCAAAGACTCGATGGTTCAGGGATTCTTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAGAGCCGAGATACTCGTTGCTGAGAGTTGGTTAACAAATATTTGGGGAGGGCAAGGCAGCATGCCGCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGCAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78476.1
    GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGATTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCGCAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTGGATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCGGATGGCCTAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGCCAAGGCAGCATGCCGCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78475.1
    ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCATTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGACGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78474.1
    AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAGAAGGGTCTGTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGTGTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTCAAAAATATAATTTGGGGAGGGTCAAGGCAGGATGCCGCCCTTCCAATTATGTGTCATATGATTGGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGAAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTACGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78473.1
    CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGGCAAGAACAAGGGGCCAAACTCACATCCGCACCGACTGAACAGTATGTATGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAACAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCACCAAGAAAGCTCCCATGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78472.1
    GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATTGGCCTCTAAAGGGGTTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGAACTACAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAAGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATAGACAGCCTCATTAAAGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATTCCGTTGCTGAGAGTTATCAAAAAATAATTTGGGGAGGGTCTTAGACTGCATGCCGCCCCCTTCCAATTATGTGTGAAATGACTTGGCGCAAGACTGCGCCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78471.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTACTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGATGTTGAGAGTTGTCAAAAAATTAATTTGGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCATTATGTGTCGGGGTGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGACAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACCGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78470.1
    AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTAGTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTAAAAAAATAATTTGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78469.1
    AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGCACTATCATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCTAGACCGCATGCCGCCCCCTTCCAATTATGTGTGATATGACTTGGCGCAAGACTGCGCCGGGCAACGAATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78468.1
    AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCACACAATAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTCCATCAAAAGCCCAGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGAGTTCCGGAGCTGAGAGTCGTTAAAAAAATGATATGGGGAGGGTCAAGGCAGCATGCTGCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTCAGATCGCCATCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATTCCCATGAGTGACCAAAGTAAACTGATCCTCCAGATTAACTCAATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78467.1
    TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAGAGGGTCTTTAGATTATGCTGGCCCAGCTAATACGACAACCCGGGTCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCACCAAAGGCACATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTACCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78466.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACGCAAGAACAAGGGGCCAAACTCACATCCGCACAGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCAAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTACAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGAAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78465.1
    GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACTCGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTTCGTTCAAAGACTCGATGGTTTACGGGATTCTTCAATTCACACCAATTATCGCATTTCGCTGCGTTCCTCATCGATTCAAGAACCGAGAAATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78464.1
    GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATTACGAAAGCCCGGGGAATCCGAACAGAAATTTGGTTGAAGAATAGTCCGCCCACCTCTTGCCGTTTTTAGGACCACCAAAAGCCCATGAAGTTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGAACCGGTTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATTTGATTCGAAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGAACTTGGCGTTCAAAGACTCGATGGTTCACGGGATTCTGAAATTCACACCATTTATCGGATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCTCCGGGCAAGGGTTTAGATCTCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGATCAATTATTATGTGATCTCAACAATGACCCTTCCGCTCACCTACGGAAACCTTGTTACG
    Z78463.1
    GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGTATTCGCACAGCAATTTTGTTGCAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAACATGCCGCCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTCCCCCGGGCGAGGGTTTAGATCCCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGAACAATGATTATGTGATCTCAACAATGAACCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78462.1
    ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTGCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGGTCACATCCGGAGACCTCGTGACG
    Z78461.1
    TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCATGGGTTGATGGGCTCCAATGGGGTTCAAAGGGTCTTCAGATTACGGTGGCCCATCTTATACGACAACCCGGGGGCTTTGGACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGGATGTATGGACAGGCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGGAATTCACACCATTTATCGGATTTCGGTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTGAAAAAATTATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGTCATTTGACTGGGGGAAAAATTGCGCCGGGTAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTCCAGCAGCAACTGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78460.1
    TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTGGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78459.1
    AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAATTTATCGCAATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78458.1
    CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCGTACGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78457.1
    CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACGATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGGTCCTTCATCCGATGAAAGAGCCGAGATATCCGTTGTTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGCCATTTGACTGGGGGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTACAACAGAAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78456.1
    GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGACCGGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTCCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78455.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGAAGCTCTTAGACCCCCCGTGCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCAACCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGGATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGAAAGGGTTTAGTTCTCGCCAAAAAGAAAGTTCCCAAGGGGGTTCGATGGAAATTCTAGGTCACTCCAAAAGAAATTGTTCATGGGTGACCAAAGTAAAAAGTTCCTCCAGATTAACTCGATCAATTTTTATGTGATCTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTGGTTACG
    Z78454.1
    GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGGTTGATGGGCCTCCAATGGGTTCAAAAGGGGCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAGAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGTCTGATGGCCTCGGCCGCAACTTGCGTTCACAGACTCGATGGTTCACGGGATTCTGCAAATCACACCATTTATCGCATGAGTTCCGTTGATGAGAGTTGTTAACAATAGTGGGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCGATTATGTGTCAAATGACTTGGAGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGTAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78453.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAGTTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78452.1
    TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACAGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCATCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATTCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAAAAAGAAAGCTCCCAAGGAGGCTTGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAATAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78451.1
    GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCTATAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGGCTTCGCACAACAATTTTATTGTAGGGATAGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAATCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGGTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATTTGACTCGCGATGCCCAGGGAGGGGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGATTTCGCTGGTTCTTCATCGATGAAAGAGCGAGATTTCCGTTGTTGAGAGTTGCAAAAAATACTTTGGGGAGGGTCAGGGCAGGATGCCCGCCCCCTACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCACGGGTTTAGATCTCGCCAGCAAGAAAGCTCCCAGGGGGGCTCCATGGAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTACACCTACGGAAACCTTGTTACG
    Z78450.1
    CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGTCCTCCAATGGGGTTCAACAGGGTCTACAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGCACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCAGTTTTGAGGACCATCATATGCCCATGCAGTTCTTAGACCCCCCGGACCAAAGAACAAGGGGCCACACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAGGGGAGAGATATGACTCGGAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGGCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGTAAATGCTCATGGATGACCAAAGTAAACAGATCCTCCAGCTTAACTCGATCAATTATTATGTGATATCAGCAATGATCCTTCC
    Z78449.1
    GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTAGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGGAGGGTCAGGGCAGGATGCCACCCTTCACCAATTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78448.1
    CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTGATCATGGGCTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTNTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGCGTTCTTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAAGGGAAGGATGCCACCCCCCTTCACCAATTATGTGTCATACGACTTGGCGCAAAACTGCGCCGGGAAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78447.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATAGGATGCCACCCCCTTCACCACTTATGTGTCATTTGACTGGGGGCAAAATTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAAAGCTCCCAGGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGGAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78446.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGAGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGGCCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCTGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTATCGCAATTTCGCTGGTCCTCCATCGATGAAGAGCCGAGTATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGGAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78445.1
    ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTGCACGGGATTCTGCAATTCACACCATTTATCGCATTTGCCGAAGACTGAGTCTCCGTTGCTGAGAGTTGTTAAAAAAATAGTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCAAACGACTTGGAGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78444.1
    AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCTTCAGATATCGGTGGCCCATCTTATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCATCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTAAGGGAAGAATGCCACCCCCTCCACCAATTTTGTGTAATTTGATTGGGGGAAAAATTGAGGGTTTAGATATCGCCAACAAGGATGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGTTCACCCTACGGAAACCTTGTTACG
    Z78443.1
    CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTATGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGGTGTGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGGCCCAAAACTCTGCCGGGCAAGGGTTTAGATATCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78442.1
    ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTTGTTGTGGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCGGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTTTGGATAGCCTCGTTAGGGGAGAGATATGATTCCCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCAAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78441.1
    CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCACAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAACTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATATCGGCAATGACCTTCC
    Z78440.1
    TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGCTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCAACGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGATATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAATGCTCCCAAGGAGGGTCGATGGAAATTCTAGGTCACTACAACAGAAATTGCCCATGGGTGACCAAAATATGCAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGGAGGTCCACCTACGGAAACCTTGTTACG
    Z78439.1
    GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTTCCGTATTTAGGACCATCCAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATTACCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCCAACTCCGCCGGGCAGGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATG



```python
# Create records list with reverse complement
# Read sequences from FASTA file
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
]
```


```python
# Print records length
len(records)
```




    94




```python
# Create records list with reverse complement
# Read sequences from FASTA file
# Limit lis to those shorter than 700 bases
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]
```


```python
# Print records length 
len(records)
```




    18




```python
# Generate records list with reverse complement
# Read sequences from FASTA file
# Limit lis to those shorter than 700 bases
records = (
rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
)

SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    18


## Multiple Sequence Alignment 
### Sequence Alignment
In this analysis, we used Biopython for sequence alignment reading, manuipulations, and conversions while generating matricecs and a phylogenetic tree.

```python
# Website for data import
#https://raw.githubusercontent.com/bipython/biopython/master/Doc/examples/PF05371_seed.sth
```


```python
# import AlignIO
from Bio import AlignIO
```


```python
# Read alignment 
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
# Print alignment
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
# Read alignment
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
# Print alignment length
print("Alignment length %i" % alignment.get_alignment_length())
```

    Alignment length 52



```python
# Print record in alignment with ID for each sequence
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
# Look for database cross-references 
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
# Print record
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
# Import Seq
# Import SeqRecord
# Import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
# Create three sequence alignment
align1 = MultipleSeqAlignment(
    [
       SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
       SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
       SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma"),
     ]     
)

align2 = MultipleSeqAlignment(
    [
       SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
      SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
      SeqRecord(Seq("GTCAGCTAC"), id="Zeta"),
  ]
)

align3 = MultipleSeqAlignment(
    [
     SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
     SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
    SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
 ]
)
```


```python
# Create alignment list 
my_alignments = [align1, align2, align3]
```


```python
# Print alignment
my_alignments
```




    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f932cb61910>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f932cb61850>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f932cb61410>]




```python
# Import AlignIO
# Write and save alignments
from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
# Read PHYLIP file
alignments = AlignIO.parse("my_example.phy", "phylip")

```


```python
# Print alignments
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAC Zeta
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota
    



```python
# Read PHYLIP file
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
# Grab last item in alignments
last_align = alignments[-1]
```


```python
# Print last alignment
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
# Grab item in first alignment
first_align = alignments[0]
```


```python
# Print first alignment
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma



```python
# Import AlignIO
from Bio import AlignIO
```


```python
# Read File
# Converts alignments 
count = AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.aln", "clustal")
```


```python
# Print amount of converted alignments
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
# Read file
alignments = AlignIO.parse("PF05371_seed.sth.txt", "stockholm")
```


```python
# Convert to clustal format
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
```


```python
# Print amount converted
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
# Convert alignment to PHYLIP format 
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python
# Convert alignment to PYHLIP relaxed format
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip-relaxed")
```




    1




```python
# Read alignment
# Create ID dictonary for names
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" % i
```


```python
# Print names
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
# Write alignment
AlignIO.write([alignment], "PF05371_seed.phy", "phylip")
```




    1




```python
# Read alignment
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
# Print alignment length
print("Number of rows: %i" % len(alignment))
```

    Number of rows: 7


```python
# Print alignment record
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
# Print alignments from row 3 to 7
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
# Print alignment in row 2, column 6
print(alignment[2,6])
```

    T



```python
# Save as code as last line
print(alignment[2].seq[6])
```

    T



```python
# Print all rows, column 6
print(alignment[:,6])
```

    TTT---T



```python
# Print rows 3 to 5, column 6
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
# Print all rows, first 6 columns
print(alignment[:, :6])
```

    Alignment with 7 rows and 6 columns
    AEPNAA COATB_BPIKE/30-81
    AEPNAA Q9T0Q8_BPIKE/1-52
    DGTSTA COATB_BPI22/32-83
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49
    FAADDA COATB_BPIF1/22-73



```python
# Print al rows, columns 6 to 8
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
# Print all rows, cloumns 9 to end
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
# Remove columns 6 to 8 
edited = alignment[:, :6] + alignment[:, 9:]
```


```python
# Print new alignment
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
# Sort sequence
edited.sort()
```


```python
# Print sequence
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49



```python
# Import SeqRecord
from Bio.SeqRecord import SeqRecord
```


```python
# Import MultipleSeqAlignment
from Bio.Align import MultipleSeqAlignment
```


```python
# Import Seq
from Bio.Seq import Seq
```


```python
# Creat 4 alignments
alignment = MultipleSeqAlignment(
  [
      SeqRecord(Seq("ACTCCTA"), id="seq1"),
      SeqRecord(Seq("AAT-CTA"), id="seq2"),
      SeqRecord(Seq("CCTACT-"), id="seq3"),
      SeqRecord(Seq("TCTCCTC"), id="seq4"),
  ]
)
```


```python
# Print alignment
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
# Create substitutions matrix
substitutions = alignment.substitutions
```


```python
# Print substitutions
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
# Save M as subtitutions for "ATCG"
m = substitutions.select("ATCG")
```


```python
# Print m
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
# Save M as subtitutions for "ACTG"
m = substitutions.select("ACTG")
```


```python
# Print M
print(m)
```

        A    C    T   G
    A 2.0  4.5  1.0 0.0
    C 4.5 10.0  0.5 0.0
    T 1.0  0.5 12.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
# Import library
import Bio.Align.Applications
```


```python
# Print directory 
dir(Bio.Align.Applications)
```




    ['ClustalOmegaCommandline',
     'ClustalwCommandline',
     'DialignCommandline',
     'MSAProbsCommandline',
     'MafftCommandline',
     'MuscleCommandline',
     'PrankCommandline',
     'ProbconsCommandline',
     'TCoffeeCommandline',
     '_ClustalOmega',
     '_Clustalw',
     '_Dialign',
     '_MSAProbs',
     '_Mafft',
     '_Muscle',
     '_Prank',
     '_Probcons',
     '_TCoffee',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__']




```python
# Import ClustualwCommandline
from Bio.Align.Applications import ClustalwCommandline
```


```python
# Reference for data import
#https://github.com/biopython/biopython/blob/master/Doc/examples/opuntia.aln
```


```python
#Import AlignIO
from Bio import AlignIO
```


```python
# Retreieve alignment
align = AlignIO.read("opuntia.aln", "clustal")
```


```python
# Print alignment
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
# Import Phylo
from Bio import Phylo
```


```python
# Stores data as phylogenetic tree
tree = Phylo.read("opuntia.dnd", "newick")
```


```python
# Draw phylogenetic tree
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
    

### Pairwise Alignment
In this analysis, we used Biopython to compare DNA sequence alignments in the local and global sphere with creating score, adjusting the scores with evaluating forwards and backwards alignments.

```python
# Import Align
from Bio import Align
```


```python
# Create aligner
aligner = Align.PairwiseAligner()
```


```python
# Create aligner with match score 1
aligner = Align.PairwiseAligner(match_score = 1.0)
```


```python
# Define target sequence
target = "GAACT"
```


```python
# Define query sequence
query = "GAT"
```


```python
# Calculate score for target and query sequences
score = aligner.score(target,query)
```


```python
# Print score
score
```




    3.0




```python
# Create alignment
alignments = aligner.align(target, query)
```


```python
# Print alignments
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
# Change mode to local
aligner.mode = "local"
```


```python
# Create target sequence
target = "AGAACTC"
```


```python
# Create query sequence
query = "GAACT"
```


```python
# Score alignments
score = aligner.score(target,query)
```


```python
# Print score
score
```




    5.0




```python
# Create alignment
alignments = aligner.align(target,query)
```


```python
# Print alignment
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
# Print aligner
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
# Show aligner algorithm 
aligner.algorithm
```




    'Smith-Waterman'




```python
# Print epsilon alinger
aligner.epsilon
```




    1e-06




```python
# Import Align 
from Bio import Align
```


```python
# Create alignment 
aligner = Align.PairwiseAligner()
```


```python
# Create target sequence
target = "GAACT"
```


```python
# Create query sequence
query = "GAT"
```


```python
# Create alignment with sequences
alignments = aligner.align(target,query)
```


```python
# Grab first alignment
alignment = alignments[0]
```


```python
# Print alignment
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7f6008242710>




```python
# Print alignment score
alignment.score
```




    3.0




```python
# Print target sequence
alignment.target
```




    'GAACT'




```python
# Print query sequence
alignment.query
```




    'GAT'




```python
# Print alignment
print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    



```python
# Print alignment coordinates
alignment.coordinates
```




    array([[0, 2, 4, 5],
           [0, 2, 2, 3]])




```python
# Print alignment length
len(alignment)
```




    2




```python
# Print alignment shape
alignment.shape
```




    (2, 5)




```python
# Change mode to local 
aligner.mode = "local"
```


```python
# Create alignment 
local_alignments = aligner.align("TGAACT", "GAC")
```


```python
# Create best match for alignment 
local_alignment = local_alignments[0]
```


```python
# Print alignment 
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
# Print alignment shape
local_alignment.shape
```




    (2, 4)




```python
# Change mode to global
aligner.mode = "global"
```


```python
# Create aligner
aligner = Align.PairwiseAligner(match = 1.0, mismatch_score = -10)
```


```python
# Print aligner
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
# Create alignments
alignments = aligner.align("AAACAAA", "AAAGAAA")
```


```python
# Print alignment length
len(alignments)
```




    2




```python
# Print first alignment 
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
# Print second alignment 
print(alignments[1])
```

    target            0 AAA-CAAA 7
                      0 |||--||| 8
    query             0 AAAG-AAA 7
    



```python
# Print alignment 
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
# Sort alignment
local_alignment.sort()
```


```python
# Print alignment 
print(local_alignment)
```

    target            0 GA-C 3
                      0 ||-| 4
    query             1 GAAC 5
    



```python
# Import Align
from Bio import Align
```


```python
# Import reverse complement 
from Bio.Seq import reverse_complement
```


```python
# Create target sequence
target = "AAACCC"
```


```python
# Create query sequence
query = "AACC"
```


```python
# Create aligner
aligner = Align.PairwiseAligner(mismatch_score= -1, internal_gap_score= -1)
```


```python
# Score aligner
aligner.score(target, query)
```




    4.0




```python
# Score aligner againt reverse complement 
aligner.score(target, reverse_complement(query))
```




    0.0




```python
# Score aligner with reverse complement 
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
# Score aligner
aligner.score(target, query, strand = "-")
```




    0.0




```python
# Create alignment
alignments = aligner.align(target,query)
```


```python
# Print aligment length
len(alignments)
```




    1




```python
# Print  first alignment 
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# Print first alignment in BED format
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
# Create aligner with reverse complement 
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
# Print alignment in BED format
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	-	1	5	0	1	4,	0,
    



```python
# Create alignment 
alignments = aligner.align(target, query, strand = "-")
```


```python
# Print alignment length
len(alignments)
```




    2




```python
# Print first alignment 
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
# Print second alignment 
print(alignments[1])
```

    target            0 ----AAACCC  6
                      0 ---------- 10
    query             4 GGTT------  0
    



```python
# Change left gap score
aligner.left_gap_score = -0.5
```


```python
# Change right gap score
aligner.right_gap_score = -0.2
```


```python
# Score aligner
aligner.score(target,query)
```




    3.3




```python
# Create alignments
alignments = aligner.align(target, query)
```


```python
# Print alignment length
len(alignments)
```




    1




```python
# Print first alignment 
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
# Create alignment with reverse complement
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
# Print alignment 
print(alignments)
```

    <Bio.Align.PairwiseAlignments object at 0x7f600819e890>



```python
# Score aligner with reverse complement
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
# Print first alignment
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             4 -AACC- 0
    



```python
# Print score for forward strand
aligner.score(target, query, strand = "+")
```




    3.3




## Blast and Challenge
### Blast 
In this analysis, we used Biopython with BLAST to take DNA from a BLAST database and print it through, looking for sequences mathcing a specific E-value. 

```python
# Import NCBIWW
from Bio.Blast import NCBIWWW
```


```python
# Register email
NCBIWWW.email = "nml020@email.latech.edu"
```


```python
# Search database
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
```


```python
# Reference for data
#https://github.com/biopython/biopython/blob/master/Doc/examples/m_cold.fasta
```


```python
# Import Seq IO
from Bio import SeqIO
```


```python
# Create record from FASTA file
record = SeqIO.read("m_cold.fasta", format = "fasta")
```


```python
# Print record
print(record)
```

    ID: gi|8332116|gb|BE037100.1|BE037100
    Name: gi|8332116|gb|BE037100.1|BE037100
    Description: gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
    Number of features: 0
    Seq('CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...TTC')



```python
# Submit sequence 
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
# Write BLAST
with open("m_cold.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
# Import NCBIXML
from Bio.Blast import NCBIXML
```


```python
# Open BLAST file
result_handle = open("my_blast.xml")
```


```python
# Read file
blast_record = NCBIXML.read(result_handle)
```


```python
# Define E-value 
E_VALUE_THRESH = 0.04
```


```python
# Check alignments
# Print alignments
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gi|2618480339|ref|XM_048479995.2| PREDICTED: Ziziphus jujuba cold-regulated 413 plasma membrane protein 2 (LOC107424728), mRNA
    length: 1028
    e value: 3.52915e-72
    GAA-CAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGA...
    ||| || ||||||||||| |  ||| |||  ||||| |||| |||||||| |   |||  |||| |  ||||  |...
    GAAGCA-AAAATGGGGAG-G--ATGGAGTTTTTGGCTATGAGAACTGATCCA---GCCACGGCTGACTTGATAAA...
    ****ALIGNMENT****
    sequence: gi|1227938481|ref|XM_022049453.1| PREDICTED: Carica papaya cold-regulated 413 plasma membrane protein 2-like (LOC110820077), mRNA
    length: 1009
    e value: 1.27832e-66
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCG...
    ||||||||||||| | |  ||| || ||||| ||||| ||||||||   ||||   ||| || | |||  ||| |...
    AGAAAATGGGGAG-G-ATGGAA-TATTTGGCTATGAAGACTGATCA---GGCCACTGCTGATCTCATCACTTCTG...
    ****ALIGNMENT****
    sequence: gi|2395983796|ref|XM_025094967.2| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X1, mRNA
    length: 980
    e value: 3.57944e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|1350315641|ref|XM_024180293.1| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X4, mRNA
    length: 868
    e value: 3.57944e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|2395983797|ref|XM_006466624.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X2, mRNA
    length: 968
    e value: 3.57944e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|1350315638|ref|XM_006425719.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X3, mRNA
    length: 893
    e value: 3.57944e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|2395983800|ref|XM_006466626.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X5, mRNA
    length: 913
    e value: 3.57944e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|1350315636|ref|XM_006425716.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X2, mRNA
    length: 881
    e value: 3.57944e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|2395983798|ref|XM_006466623.4| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X3, mRNA
    length: 1052
    e value: 3.57944e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|2395983799|ref|XM_006466625.3| PREDICTED: Citrus sinensis cold-regulated 413 plasma membrane protein 2 (LOC102620025), transcript variant X4, mRNA
    length: 978
    e value: 3.57944e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|1350315634|ref|XM_006425717.2| PREDICTED: Citrus clementina cold-regulated 413 plasma membrane protein 2 (LOC18037141), transcript variant X1, mRNA
    length: 952
    e value: 3.57944e-62
    AAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATA...
    |||| ||| ||| | |||| || ||||| |||||||||||| |   ||  |  ||  |  |||||   || ||| ...
    AAAT-GGG-GAG-ATTGAATTATTTGGCTATGAAAACTGATGATCAGGTTGCAGCAGAGTTGATCAGCTCTGATT...
    ****ALIGNMENT****
    sequence: gi|2526866810|ref|XM_057645500.1| PREDICTED: Actinidia eriantha cold-regulated 413 plasma membrane protein 2-like (LOC130785340), mRNA
    length: 1152
    e value: 1.28739e-61
    AAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGAT...
    ||||| ||| ||| |||| | || ||||| ||||| || |||| |  |||  || | ||| ||||| |||| || ...
    AAAAT-GGG-GAG-AATGGATTATTTGGCGATGAAGACCGATCCAGCGGC--TGCCGAAT-TGATCAATTCGGAC...
    ****ALIGNMENT****
    sequence: gi|1768569081|ref|XM_031406607.1| PREDICTED: Pistacia vera cold-regulated 413 plasma membrane protein 2-like (LOC116120644), mRNA
    length: 982
    e value: 3.60485e-57
    TCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATTACGGGTTT...
    || |||||||||||||||||  | ||  ||| || ||| || ||  | || | ||  ||||||  |   || |||...
    TCTGATATCAATGAGCTTAAGCTTGCTGCAAAGAAGCTTATTAACCACGCAACTAAACTCGGTGGTCTTGGCTTT...
    ****ALIGNMENT****
    sequence: gi|1954740698|ref|XM_038867092.1| PREDICTED: Tripterygium wilfordii cold-regulated 413 plasma membrane protein 2 (LOC120014952), mRNA
    length: 999
    e value: 3.60485e-57
    GAACAGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGAT...
    ||| ||||||| ||| ||| || | | || ||||| ||||| |||||||    ||  ||||   || |||||   ...
    GAAAAGAAAAT-GGG-GAG-AACGGATTATTTGGCGATGAAGACTGATC---CGGTTGTGGACGATTTGATCAGC...
    ****ALIGNMENT****
    sequence: gi|2526831378|ref|XM_057642435.1| PREDICTED: Actinidia eriantha cold-regulated 413 plasma membrane protein 2-like (LOC130782944), mRNA
    length: 1214
    e value: 3.60485e-57
    GAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGA...
    |||||| ||| ||| |||| | || ||||| ||||| || |||| |   | | || | ||| ||||| |||||||...
    GAAAAT-GGG-GAG-AATGGATTATTTGGCGATGAAGACCGATCCAGCTG-C-TGCCGAAT-TGATCAATTCCGA...
    ****ALIGNMENT****
    sequence: gi|1350280614|ref|XM_024170292.1| PREDICTED: Morus notabilis cold-regulated 413 plasma membrane protein 2 (LOC21394987), mRNA
    length: 1020
    e value: 2.16956e-54
    TACTTGGCCATGAAAACTGATCAATTGGCCGTGGCTAATATGATCGATTCCGATATCAATGAGCTTAAAATGGCA...
    || |||||||||||||| || | |   |||  |||| || ||||  |||| |||||||| ||||| || || || ...
    TATTTGGCCATGAAAACGGACCCA---GCCACGGCTGATTTGATAAATTCTGATATCAACGAGCTGAAGATCGCT...
    ****ALIGNMENT****
    sequence: gi|2350562715|ref|XM_052350285.1| PREDICTED: Diospyros lotus cold-regulated 413 plasma membrane protein 2-like (LOC127810688), mRNA
    length: 1265
    e value: 1.0094e-52
    GGGTTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAAC...
    || || ||||||   |||||||||||| | |||| | |||| |||||||||||||||||||||||||||||||||...
    GGCTTCGGCACTTCCTTCCTCAAATGGATTGCCTCCTTTGCTGCTATTTACTTGTTGATATTGGATCGAACAAAC...
    ****ALIGNMENT****
    sequence: gi|1585724761|ref|XM_028202722.1| PREDICTED: Camellia sinensis cold-regulated 413 plasma membrane protein 2-like (LOC114262355), mRNA
    length: 910
    e value: 3.63044e-52
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAAT-TGGCCGTGGCTAATATGATCGATTCC...
    |||||||||||||  ||||| |||| ||||| ||||| || |||||    | ||    | |||  |||  |||||...
    AGAAAATGGGGAGGAAAATGGAGTATTTGGCAATGAAGACCGATCATCCAGCCCCAACCCAAT-CGATGAATTCC...
    ****ALIGNMENT****
    sequence: gi|2583747300|ref|XM_059787294.1| PREDICTED: Cornus florida cold-regulated 413 plasma membrane protein 2-like (LOC132285128), mRNA
    length: 1126
    e value: 1.30574e-51
    AGAAAATGGGGAGAGAAATGAAGTACTTGGCCATGAAAACTGATCAATTGGCC-GTGGCTAATATGATCGATTCC...
    ||||||| ||| |||||| | |||| ||||| |||||||||||||    ||||   | | ||| ||||| |||||...
    AGAAAAT-GGG-GAGAAA-GGAGTATTTGGCTATGAAAACTGATC---CGGCCACAGCCGAAT-TGATCAATTCC...
    ****ALIGNMENT****
    sequence: gi|2118882425|ref|XM_044613294.1| PREDICTED: Mangifera indica cold-regulated 413 plasma membrane protein 2-like (LOC123198583), mRNA
    length: 1083
    e value: 4.69627e-51
    CTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAACTGGAGAACCA...
    ||| ||| ||||||||| | |||| | |||| ||||||||||||||||||||||||||||||||||||||||| |...
    CTC-TTTTCTCAAATGGGTTGCCTCCTTTGCTGCTATTTACTTGTTGATATTGGATCGAACAAACTGGAGAACAA...
    ****ALIGNMENT****
    sequence: gi|2537875835|ref|XM_021822463.2| PREDICTED: Hevea brasiliensis cold-regulated 413 plasma membrane protein 2 (LOC110663212), transcript variant X1, mRNA
    length: 940
    e value: 2.82643e-48
    CGGGTTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAA...
    ||| ||||| |||  ||| ||||||||| | ||||   |||| || || || ||| |||||||||||||||||||...
    CGGCTTTGGTACTTCTTTTCTCAAATGGGTTGCCTCTTTTGCTGCCATATATTTGCTGATATTGGATCGAACAAA...
    ****ALIGNMENT****
    sequence: gi|1220047144|ref|XM_021953815.1| PREDICTED: Prunus avium cold-regulated 413 plasma membrane protein 2-like (LOC110753022), mRNA
    length: 894
    e value: 2.82643e-48
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    ||||  |||| || |||||||| || || || ||| | | ||  |||||||| ||||| | || ||| |||   |...
    TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAAAGATGCCACTAAGCTTGGTGGGT...
    ****ALIGNMENT****
    sequence: gi|1220094463|ref|XM_021978417.1| PREDICTED: Prunus avium cold-regulated 413 plasma membrane protein 2-like (LOC110773902), mRNA
    length: 896
    e value: 2.82643e-48
    TGATCGATTCCGATATCAATGAGCTTAAAATGGCAACAATGAGGCTCATCAATGATGCTAGTATGCTCGGTCATT...
    ||||  |||| || |||||||| || || || ||| | | ||  |||||||| ||||| | || ||| |||   |...
    TGATAAATTCAGACATCAATGATCTCAAGATTGCAGCCAAGAAACTCATCAAAGATGCCACTAAGCTTGGTGGGT...
    ****ALIGNMENT****
    sequence: gi|2118854713|ref|XM_044609109.1| PREDICTED: Mangifera indica cold-regulated 413 plasma membrane protein 2-like (LOC123195424), mRNA
    length: 742
    e value: 1.71315e-40
    GGGTTTGGCA-CTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAA...
    || ||||| | ||| ||| ||| ||||||| |||| | |||| |||||||| |||||||||||||||||||||||...
    GGCTTTGGAATCTC-TTTTCTCCAATGGCTTGCCTCCTTTGCTGCTATTTATTTGTTGATATTGGATCGAACAAA...
    ****ALIGNMENT****
    sequence: gi|2537875838|ref|XM_021822464.2| PREDICTED: Hevea brasiliensis cold-regulated 413 plasma membrane protein 2 (LOC110663212), transcript variant X2, mRNA
    length: 924
    e value: 4.83107e-31
    CGGGTTTGGCACTCATTTCCTCAAATGGCTCGCCTGCCTTGCGGCTATTTACTTGTTGATATTGGATCGAACAAA...
    ||| ||||| |||  ||| ||||||||| | ||||   |||| || || || ||| |||||||||||||||||||...
    CGGCTTTGGTACTTCTTTTCTCAAATGGGTTGCCTCTTTTGCTGCCATATATTTGCTGATATTGGATCGAACAAA...
    ****ALIGNMENT****
    sequence: gi|2739574736|gb|CP157410.1| Ziziphus jujuba isolate TW-2024b chromosome 08
    length: 27596186
    e value: 1.06064e-17
    TTACTTGTTGATATTGGATCGAACAAACTGGAGAACCAACATGCTCACGTCACTTTTAGTCCCTTACATATTCCT...
    ||| |||||||||||||||||||||||||||||||| || ||||| || ||||||||||| ||||| || ||| |...
    TTATTTGTTGATATTGGATCGAACAAACTGGAGAACAAATATGCTTACTTCACTTTTAGTTCCTTATATCTTCTT...
    ****ALIGNMENT****
    sequence: gi|3129608210|emb|OZ375381.1| Ziziphus lotus genome assembly, chromosome: 8
    length: 22133894
    e value: 4.9347e-16
    TTACTTGTTGATATTGGATCGAACAAACTGGAGAACCAACATGCTCACGTCACTTTTAGTCCCTTACATATTCCT...
    ||| |||||||||||||||||||||||||||||||| || ||||| || ||||| ||||| ||||| || ||| |...
    TTATTTGTTGATATTGGATCGAACAAACTGGAGAACAAATATGCTTACTTCACTCTTAGTTCCTTATATCTTCTT...


### Challenge 
In this analysis, I found a cancer gene of ESR1 and ran it through BLAST finding the E-value and comparing it to that of the chimp.
```python
from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = "nml020@email.latech.edu"
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "NM_005228.5")
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("ESR1.txt", format = "fasta")
```


```python
print(record)
```

    ID: ESR1
    Name: ESR1
    Description: ESR1
    Number of features: 0
    Seq('AGGAGCTGGCGGAGGGCGTTCGTCCTGGGACTGCACTTGCTCCCGTCGGGTCGC...AAA')



```python
with open("EGFR", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("EGFR")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
with open("EGFR") as result_handle:
    blast_record = NCBIXML.read(result_handle)
```


```python
E_VALUE_THRESH = 100
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gi|1519245592|ref|NM_005228.5| Homo sapiens epidermal growth factor receptor (EGFR), transcript variant 1, mRNA
    length: 9905
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ****ALIGNMENT****
    sequence: gi|2694467815|ref|XM_034964127.3| PREDICTED: Pan paniscus epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 10112
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCACACGGTGTGAGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|2697829648|ref|XM_519102.9| PREDICTED: Pan troglodytes epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 10049
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||| |||||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||| ||||...
    AGACGACCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCACACGGTGTGAGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|2217365841|ref|XM_047419952.1| PREDICTED: Homo sapiens epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 48454
    e value: 0.0
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    ****ALIGNMENT****
    sequence: gi|2217365843|ref|XM_047419953.1| PREDICTED: Homo sapiens epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 10631
    e value: 0.0
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    ****ALIGNMENT****
    sequence: gi|1676318675|ref|NM_001346900.2| Homo sapiens epidermal growth factor receptor (EGFR), transcript variant 8, mRNA
    length: 9676
    e value: 0.0
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    ****ALIGNMENT****
    sequence: gi|2462612701|ref|XM_054357417.1| PREDICTED: Homo sapiens epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 9577
    e value: 0.0
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    ****ALIGNMENT****
    sequence: gi|2462612703|ref|XM_054357418.1| PREDICTED: Homo sapiens epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 10066
    e value: 0.0
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    ****ALIGNMENT****
    sequence: gi|2695076151|ref|XM_031012805.3| PREDICTED: Gorilla gorilla gorilla epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 10067
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|2695861090|ref|XM_054494705.2| PREDICTED: Pongo pygmaeus epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 10254
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|2697829649|ref|XM_054687551.2| PREDICTED: Pan troglodytes epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9642
    e value: 0.0
    AAGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AAGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGA...
    ****ALIGNMENT****
    sequence: gi|2694660314|ref|XM_024249444.3| PREDICTED: Pongo abelii epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 9886
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|1890245545|ref|NM_001346899.2| Homo sapiens epidermal growth factor receptor (EGFR), transcript variant 7, mRNA
    length: 9770
    e value: 0.0
    GCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGGGGTGCAGGAGAGGAGAACTGCCAGAAACTGACCA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGGGGTGCAGGAGAGGAGAACTGCCAGAAACTGACCA...
    ****ALIGNMENT****
    sequence: gi|1890245545|ref|NM_001346899.2| Homo sapiens epidermal growth factor receptor (EGFR), transcript variant 7, mRNA
    length: 9770
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ****ALIGNMENT****
    sequence: gi|2694557929|ref|XM_055293115.2| PREDICTED: Symphalangus syndactylus epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 10025
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||| |||||||||||||||||||||||||||||||| |||||||||||||||||| | ||...
    AGACGTCCGGGCAGCCCCTGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGAGCGCCCGCCTCG...
    ****ALIGNMENT****
    sequence: gi|2695861091|ref|XM_063667305.1| PREDICTED: Pongo pygmaeus epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 12014
    e value: 0.0
    AAAGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AAAGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGG...
    ****ALIGNMENT****
    sequence: gi|2695861091|ref|XM_063667305.1| PREDICTED: Pongo pygmaeus epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 12014
    e value: 8.78936e-165
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|2695861095|ref|XM_063667307.1| PREDICTED: Pongo pygmaeus epidermal growth factor receptor (EGFR), transcript variant X4, mRNA
    length: 9598
    e value: 0.0
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    ****ALIGNMENT****
    sequence: gi|2695861093|ref|XM_063667306.1| PREDICTED: Pongo pygmaeus epidermal growth factor receptor (EGFR), transcript variant X3, mRNA
    length: 10626
    e value: 0.0
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    ****ALIGNMENT****
    sequence: gi|1788696280|ref|XM_023232039.3| PREDICTED: Piliocolobus tephrosceles epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 9946
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||| ||||||||||| ||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCTGCGCAGCGCGGACGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|1675033261|ref|NM_001346941.2| Homo sapiens epidermal growth factor receptor (EGFR), transcript variant EGFRvIII, mRNA
    length: 9104
    e value: 0.0
    GTAATTATGTGGTGACAGATCACGGCTCGTGCGTCCGAGCCTGTGGGGCCGACAGCTATGAGATGGAGGAAGACG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GTAATTATGTGGTGACAGATCACGGCTCGTGCGTCCGAGCCTGTGGGGCCGACAGCTATGAGATGGAGGAAGACG...
    ****ALIGNMENT****
    sequence: gi|1675033261|ref|NM_001346941.2| Homo sapiens epidermal growth factor receptor (EGFR), transcript variant EGFRvIII, mRNA
    length: 9104
    e value: 1.48707e-174
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ****ALIGNMENT****
    sequence: gi|1777319992|ref|XM_003895963.5| PREDICTED: Papio anubis epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 10029
    e value: 0.0
    TCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCGGCCGA...
    |||||||||||||||||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||||| |...
    TCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCGGCCAA...
    ****ALIGNMENT****
    sequence: gi|1751284283|ref|XM_010380548.2| PREDICTED: Rhinopithecus roxellana epidermal growth factor receptor (EGFR), mRNA
    length: 10612
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||| ||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCTGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|1411143690|ref|XM_025378159.1| PREDICTED: Theropithecus gelada epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 9927
    e value: 0.0
    GTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCGGCCG...
    ||||||||||||||||||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||||| ...
    GTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCGGCCA...
    ****ALIGNMENT****
    sequence: gi|1059138040|ref|XM_017862139.1| PREDICTED: Rhinopithecus bieti epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 10646
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||| ||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCTGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|3129478328|ref|XM_015133436.3| PREDICTED: Macaca mulatta epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 9870
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|1825932383|ref|XM_033194501.1| PREDICTED: Trachypithecus francoisi epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 10680
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||| ||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCTGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|2898751211|ref|XM_011760496.2| PREDICTED: Macaca nemestrina epidermal growth factor receptor (LOC105493061), transcript variant X1, mRNA
    length: 9867
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|2309461522|ref|XM_050783477.1| PREDICTED: Macaca thibetana thibetana epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 9864
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|2978663125|ref|XM_005549559.5| PREDICTED: Macaca fascicularis epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 9869
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|795268648|ref|XM_012065405.1| PREDICTED: Cercocebus atys epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9858
    e value: 0.0
    CGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCGGCCGAGG...
    |||||||||||||||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||||| | |...
    CGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCGGCCAAAG...
    ****ALIGNMENT****
    sequence: gi|2954769443|ref|XM_007981280.3| PREDICTED: Chlorocebus sabaeus epidermal growth factor receptor (EGFR), mRNA
    length: 9871
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|1788696281|ref|XM_023232041.2| PREDICTED: Piliocolobus tephrosceles epidermal growth factor receptor (EGFR), transcript variant X3, mRNA
    length: 9822
    e value: 0.0
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    |||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    ****ALIGNMENT****
    sequence: gi|1777319994|ref|XM_009202927.4| PREDICTED: Papio anubis epidermal growth factor receptor (EGFR), transcript variant X3, mRNA
    length: 9853
    e value: 0.0
    AAGAA--AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    |||||  |||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||...
    AAGAATCAGTTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    ****ALIGNMENT****
    sequence: gi|1411143694|ref|XM_025378161.1| PREDICTED: Theropithecus gelada epidermal growth factor receptor (EGFR), transcript variant X3, mRNA
    length: 9814
    e value: 0.0
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    |||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    ****ALIGNMENT****
    sequence: gi|1059138042|ref|XM_017862140.1| PREDICTED: Rhinopithecus bieti epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9608
    e value: 0.0
    AAGAA--AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    |||||  |||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||...
    AAGAATCAGTTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    ****ALIGNMENT****
    sequence: gi|795296714|ref|XM_011966774.1| PREDICTED: Mandrillus leucophaeus epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9572
    e value: 0.0
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    |||||||||||||||||||||||| ||||| ||||||||||||||||||||||||||||||||||||||||||||...
    AGTTTGCCAAGGCACGAGTAACAAACTCACACAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    ****ALIGNMENT****
    sequence: gi|795332916|ref|XM_011925699.1| PREDICTED: Colobus angolensis palliatus epidermal growth factor receptor (EGFR), mRNA
    length: 9805
    e value: 0.0
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    |||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    ****ALIGNMENT****
    sequence: gi|795296705|ref|XM_011966773.1| PREDICTED: Mandrillus leucophaeus epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 9801
    e value: 0.0
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    |||||||||||||||||||||| ||||| ||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGCCAAGGCACGAGTAACAAACTCACACAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    ****ALIGNMENT****
    sequence: gi|2898751282|ref|XM_011760498.3| PREDICTED: Macaca nemestrina epidermal growth factor receptor (LOC105493061), transcript variant X3, mRNA
    length: 11725
    e value: 0.0
    AAGAA--AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    |||||  |||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||...
    AAGAATCAGTTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    ****ALIGNMENT****
    sequence: gi|2898751212|ref|XM_071095193.1| PREDICTED: Macaca nemestrina epidermal growth factor receptor (LOC105493061), transcript variant X2, mRNA
    length: 9541
    e value: 0.0
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    |||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    ****ALIGNMENT****
    sequence: gi|2978663126|ref|XM_005549560.5| PREDICTED: Macaca fascicularis epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9592
    e value: 0.0
    AAGAA--AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    |||||  |||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||...
    AAGAATCAGTTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    ****ALIGNMENT****
    sequence: gi|2978663127|ref|XM_045389159.3| PREDICTED: Macaca fascicularis epidermal growth factor receptor (EGFR), transcript variant X3, mRNA
    length: 10129
    e value: 0.0
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    |||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGTTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    ****ALIGNMENT****
    sequence: gi|795268686|ref|XM_012065409.1| PREDICTED: Cercocebus atys epidermal growth factor receptor (EGFR), transcript variant X6, mRNA
    length: 9670
    e value: 0.0
    AAGAA--AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    |||||  |||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||...
    AAGAATCAGTTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCC...
    ****ALIGNMENT****
    sequence: gi|795268677|ref|XM_012065408.1| PREDICTED: Cercocebus atys epidermal growth factor receptor (EGFR), transcript variant X5, mRNA
    length: 9813
    e value: 0.0
    AGTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    |||||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||...
    AGTTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGAT...
    ****ALIGNMENT****
    sequence: gi|795268638|ref|XM_012065404.1| PREDICTED: Cercocebus atys epidermal growth factor receptor (EGFR), transcript variant X1, mRNA
    length: 14852
    e value: 0.0
    GTTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATG...
    ||||||||||||||||||||||| |||||||||||||||||||||||||||||||||||||||||||||||||||...
    GTTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATG...
    ****ALIGNMENT****
    sequence: gi|795268668|ref|XM_012065407.1| PREDICTED: Cercocebus atys epidermal growth factor receptor (EGFR), transcript variant X4, mRNA
    length: 9540
    e value: 0.0
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    |||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    ****ALIGNMENT****
    sequence: gi|795268659|ref|XM_012065406.1| PREDICTED: Cercocebus atys epidermal growth factor receptor (EGFR), transcript variant X3, mRNA
    length: 17529
    e value: 0.0
    TTTGCCAAGGCACGAGTAACAAGCTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    |||||||||||||||||||||| ||||||||||||||||||||||||||||||||||||||||||||||||||||...
    TTTGCCAAGGCACGAGTAACAAACTCACGCAGTTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGT...
    ****ALIGNMENT****
    sequence: gi|1777319993|ref|XM_009202926.4| PREDICTED: Papio anubis epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 10180
    e value: 0.0
    AGACCATCCAGGAGGTGGCTGGTTATGTCCTCATTGCCCTCAACACAGTGGAGCGAATTCCTTTGGAAAACCTGC...
    |||||||||||||||||||||||||||||||||| |||||||||||||||||||| |||||||||||||||||||...
    AGACCATCCAGGAGGTGGCTGGTTATGTCCTCATCGCCCTCAACACAGTGGAGCGGATTCCTTTGGAAAACCTGC...
    ****ALIGNMENT****
    sequence: gi|2694557930|ref|XM_055293116.2| PREDICTED: Symphalangus syndactylus epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9890
    e value: 0.0
    GCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGGGGTGCAGGAGAGGAGAACTGCCAGAAACTGACCA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGGGGTGCAGGAGAGGAGAACTGCCAGAAACTGACCA...
    ****ALIGNMENT****
    sequence: gi|2694557930|ref|XM_055293116.2| PREDICTED: Symphalangus syndactylus epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9890
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    |||||||||||||||||| |||||||||||||||||||||||||||||||| |||||||||||||||||| | ||...
    AGACGTCCGGGCAGCCCCTGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGAGCGCCCGCCTCG...
    ****ALIGNMENT****
    sequence: gi|1297699725|ref|XM_023232040.1| PREDICTED: Piliocolobus tephrosceles epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9810
    e value: 0.0
    GCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGGGGTGCAGGAGAGGAGAACTGCCAGAAACTGACCA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGGGGTGCAGGAGAGGAGAACTGCCAGAAACTGACCA...
    ****ALIGNMENT****
    sequence: gi|1297699725|ref|XM_023232040.1| PREDICTED: Piliocolobus tephrosceles epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9810
    e value: 0.0
    AGACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCG...
    ||||||||||||||||||| ||||||||||| ||||||||||||||||||| |||||||||| ||||||| ||||...
    AGACGTCCGGGCAGCCCCCTGCGCAGCGCGGACGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCG...
    ****ALIGNMENT****
    sequence: gi|1411143692|ref|XM_025378160.1| PREDICTED: Theropithecus gelada epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9794
    e value: 0.0
    GCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGGGGTGCAGGAGAGGAGAACTGCCAGAAACTGACCA...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    GCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGGGGTGCAGGAGAGGAGAACTGCCAGAAACTGACCA...
    ****ALIGNMENT****
    sequence: gi|1411143692|ref|XM_025378160.1| PREDICTED: Theropithecus gelada epidermal growth factor receptor (EGFR), transcript variant X2, mRNA
    length: 9794
    e value: 0.0
    ACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCCGCACGGTGTGAGCGCCCGACGCGGC...
    ||||||||||||||||||||||||||||||||||||||||||||||||| |||||||||| ||||||| ||||||...
    ACGTCCGGGCAGCCCCCGGCGCAGCGCGGCCGCAGCAGCCTCCGCCCCCGGCACGGTGTGCGCGCCCGCCGCGGC...



## Open CV
### Open CV 1: "OpenCVBasics"
In this analysis, we used an image, particularly of a sloth, to learn the basics of image processing with modifications including color conservion and image rotation. 

```python
# Import Libraries
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
import cv2
```


```python
# Import image
img = cv2.imread("sloth.jpeg")
```


```python
# Check data type is loaded
type(img)
```




    numpy.ndarray




```python
# Test and try to run image from wrong path
img_wrong = cv2.imread('wrong/path/doesnt/abcdegh.jpg')
```


```python
type(img_wrong)
```




    NoneType




```python
# Display image
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7efcdf41fd50>

<img width="555" height="406" alt="image" src="https://github.com/user-attachments/assets/1b53c8bf-e74c-441b-93ab-a6d321decf23" />





```python
# Convert image to correct color display 
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
# Display image
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7efcdf3d6dd0>

<img width="588" height="414" alt="image" src="https://github.com/user-attachments/assets/82be48ec-93f4-4681-aa6b-d95b4c360425" />




```python
# Give image dimensions
img_gray = cv2.imread("sloth.jpeg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (190, 266)




```python
# Display image
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7efcddb49b10>


<img width="582" height="397" alt="image" src="https://github.com/user-attachments/assets/290ceeae-fb1e-4d9d-9494-f2b3c8397d69" />





```python
# Convert image to gray and display
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7efcddabc710>


<img width="580" height="399" alt="image" src="https://github.com/user-attachments/assets/364537b2-7422-4603-a4c5-e3fd07f097c3" />





```python
# Give image dimensions
fix_img.shape
```




    (190, 266, 3)




```python
# Resize image with new dimensions
new_img = cv2.resize(fix_img,(1000,400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7efcdc224f10>

<img width="636" height="279" alt="image" src="https://github.com/user-attachments/assets/61590f54-f4bb-4b79-8fa6-3a3061be5943" />





```python
# Give image dimensions
new_img.shape
```




    (400, 1000, 3)




```python
# Resize image with new height and width ratios
w_ratio = 0.5
h_ratio = 0.5

new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
# Display image
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7efcdc207a10>


<img width="573" height="412" alt="image" src="https://github.com/user-attachments/assets/bfa94e2e-d0fb-45e6-b951-812bc013e8a5" />





```python
# Give dimensions
new_img.shape
```




    (95, 133, 3)




```python
# Vertical Flip image
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7efcdc174690>


<img width="576" height="399" alt="image" src="https://github.com/user-attachments/assets/3b093fce-d776-4009-b07c-175e6d6ee037" />





```python
# Horizontal image flip
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7efcdc0e3790>

<img width="615" height="414" alt="image" src="https://github.com/user-attachments/assets/4886179e-2463-475e-b919-ee6086043cc0" />





```python
# Check image is saved in numpy
type(fix_img)
```




    numpy.ndarray




```python
# Save modified image
cv2.imwrite('sloth_fixed_image.jpeg', flip_img)
```




    True

## Open CV Pt. 2
In this analysis, we continued with color conversions with the sloth image but added a new concept of image mixing and overlay. 
```python
# Import libraries
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# Import image
img = cv2.imread("sloth.jpeg")
```


```python
# Display image 
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f371c22f090>


<img width="623" height="410" alt="image" src="https://github.com/user-attachments/assets/cc022b08-9cac-4b73-8310-ef6de6c76dbe" />




```python
# Convert image to correct coloring 
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
# Display image
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f3717ab0bd0>


<img width="578" height="413" alt="image" src="https://github.com/user-attachments/assets/6fa3d5d7-31cb-4587-a564-7c146a8d3baa" />



```python
# Convert image colors using HSV coloring
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
# Display image
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f3717a30810>


<img width="586" height="408" alt="image" src="https://github.com/user-attachments/assets/af995dc5-2cd2-4e6a-8b2e-547d418a8d7b" />




```python
# Convert image colors using HLS coloring
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
# Display image
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7f3717a15fd0>


<img width="664" height="413" alt="image" src="https://github.com/user-attachments/assets/da36ccdf-dc9e-45b1-aaf9-29aa420df219" />




```python
# Load other image
img1 = cv2.imread('do-not-copy.jpeg')
img2 = cv2.imread("sloth.jpeg")
```


```python
# Display image
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f3717993810>


<img width="657" height="413" alt="image" src="https://github.com/user-attachments/assets/fea3460c-805e-40a9-9a4c-92203cb527dd" />



```python
# Covert image to correct coloring 
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
# Display image
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7f37178ffc90>

<img width="605" height="401" alt="image" src="https://github.com/user-attachments/assets/a97d9c58-ada3-40a5-8053-66f70aef9005" />




```python
# Display first image
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7f3717878090>


<img width="590" height="423" alt="image" src="https://github.com/user-attachments/assets/6535a8a8-dc58-4af4-9246-7a0fe158cff2" />




```python
# Resize images to the same size
img1 = cv2.resize(img1, (1200,1200))
img2 = cv2.resize(img2, (1200, 1200))
```


```python
# Set blending weights of images
alpha = 0.5
beta = 0.5
```


```python
# Blend images together
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


```python
# Display blended image
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7f3716381c90>


<img width="540" height="390" alt="image" src="https://github.com/user-attachments/assets/63d0f2b0-c4e5-400f-9913-7e75d4cdb7d8" />




```python
# Change blending weights
# Blend image
# Display blended image
alpha = 0.2
beta = 0.8

blended1 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended1)
```




    <matplotlib.image.AxesImage at 0x7f371574e5d0>


<img width="471" height="400" alt="image" src="https://github.com/user-attachments/assets/6cf4fd48-96ef-4af6-873b-8029fcde3134" />




```python
# Reload original images
# Convert images to RGB coloring
# Resize sloth image
img1 = cv2.imread('do-not-copy.jpeg')
img2 = cv2.imread('sloth.jpeg')

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)

img1 = cv2.resize(img1, (100,100))
```


```python
# Give images different sizes
# Set image placements
# Create overaly boundary
# Paste smaller image onto larger image
# Display image combination
large_img = img2
small_img = img1

x_offset = 0
y_offset = 0

x_end = x_offset + small_img.shape[0]
y_end = y_offset + small_img.shape[1]

large_img[y_offset:y_end, x_offset:x_end] = small_img

plt.imshow(large_img)
```




    <matplotlib.image.AxesImage at 0x7f37177f6710>

<img width="556" height="419" alt="image" src="https://github.com/user-attachments/assets/b45a6528-5db5-473c-bfc0-9018f621d1ca" />

### Open CV Pt. 3
In this analysis, we finished off images revisions with grayscaling and applying various thresholds to images in comparing pixelation intensity. 
```python
# Import libraries
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# Load image
img = cv2.imread('rainbow.jpg')
```


```python
# Display image
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7f35d63e1150>


<img width="271" height="263" alt="image" src="https://github.com/user-attachments/assets/09070584-a163-43e2-9361-da3160f1624c" />




```python
# Load image in grayscale
img = cv2.imread('rainbow.jpg', 0)
```


```python
# Display image
plt.imshow(img, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f35d6382810>


<img width="286" height="254" alt="image" src="https://github.com/user-attachments/assets/b58aa5da-78ad-4dee-bd90-828471bcd0b6" />




```python
# Apply binary threshold
ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
```


```python
# Show threshold value used
ret1
```




    127.0




```python
# Display image
plt.imshow(thresh1, cmap= "gray")
```




    <matplotlib.image.AxesImage at 0x7f35d42ebd10>


<img width="306" height="257" alt="image" src="https://github.com/user-attachments/assets/bc50dbcb-8a70-487e-987f-b9ad43be6e55" />




```python
# Apply truncation threshold
# Display image
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f35d4257650>

<img width="281" height="254" alt="image" src="https://github.com/user-attachments/assets/44c96e32-eb2e-48ec-9fd1-5a0f42638a89" />




```python
# Apply to-zero threshold
# Display image
img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f35d42440d0>


<img width="262" height="256" alt="image" src="https://github.com/user-attachments/assets/f8a15103-0863-44e6-9bd9-896de92d2d77" />




```python
# Import new image
# Load image in grayscale
# Display image
img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f35d419fcd0>


<img width="238" height="257" alt="image" src="https://github.com/user-attachments/assets/4bdb3080-a2b7-4b34-8d63-b902829947e5" />




```python
# Define functions for image display
# Make image larger
# Add subplot
# Load image in grayscale
def show_pic(img):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
# Display image
show_pic(img_r)
```

<img width="636" height="858" alt="image" src="https://github.com/user-attachments/assets/5b8c235c-b2f9-4623-a2ff-b227413333de" />




```python
# Apply binary threshold to new image
# Display image
ret, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```
<img width="583" height="861" alt="image" src="https://github.com/user-attachments/assets/633e5096-98fb-422b-8bfd-84e03136707c" />




```python
# Apply higher binary threshold
# Display image
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```

<img width="646" height="862" alt="image" src="https://github.com/user-attachments/assets/f4c31f00-d3e2-4e45-96cc-973baa81e28c" />




```python
# Apply adaptive mean threshold
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
```


```python
# Display image
show_pic(th2)
```
<img width="595" height="868" alt="image" src="https://github.com/user-attachments/assets/f49b8790-5c58-4abf-9035-a714f695eebb" />




```python
# Blend threshold images
# Display image
blended = cv2.addWeighted(src1 = th1, alpha =0.6,
                          src2 = th2, beta =0.4, gamma = 0)
show_pic(blended)
```
<img width="594" height="875" alt="image" src="https://github.com/user-attachments/assets/3b20338a-a52f-4575-84cc-584927ed3044" />




```python
# Apply adaptive mean threshold again
# Blend images
# Display image
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th3, beta = 0.4, gamma = 0)

show_pic(blended)
```
<img width="613" height="868" alt="image" src="https://github.com/user-attachments/assets/aa955b24-62a7-4573-976d-62e161d1a886" />




## Aspect Detection
### Corner Detection
In this analysis, we used to chessboard images to detect and visualize their corner features using Harris and Shi-Tomasi corner detection methods. 
```python
# Import libraries
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# Load flat chessboard image
# Display image in color
flat_chess = cv2.imread('chessboard_green.png')
flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2RGB)
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7f071da3cf90>

<img width="450" height="403" alt="image" src="https://github.com/user-attachments/assets/f64e5d4b-1975-404f-9ef1-d222fca6965c" />




```python
# Display image in gray
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7f071d9aaa10>


<img width="457" height="411" alt="image" src="https://github.com/user-attachments/assets/e25f604c-e734-435b-bdb1-a1e7eb4cc2a2" />




```python
# Load real chessboard image
real_chess = cv2.imread("chessboard.jpg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
# Display image
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7f071c101c90>

<img width="643" height="413" alt="image" src="https://github.com/user-attachments/assets/6e1efa85-5e6d-4ac0-9d0d-4d25a91c80d3" />




```python
# Display image in gray
gray_real_chess = cv2. cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7f071c0950d0>


<img width="662" height="427" alt="image" src="https://github.com/user-attachments/assets/72681734-abe9-4e3c-8110-fa6cbf4e38ab" />




```python
# Convert flat chessboard to grayscale
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)

dst = cv2.dilate(dst, None)
```


```python
# Display image
flat_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7f071c0df9d0>

<img width="466" height="421" alt="image" src="https://github.com/user-attachments/assets/8babc0e1-6b37-40c9-9165-d10a25d06e19" />




```python
# Use of Harris Detection
# Convert real chessboard to grayscale
# Display image
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k = 0.04)
dst = cv2.dilate(dst, None)

real_chess[dst>0.01*dst.max()] = [255, 0, 0]

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7f0715f3d650>


<img width="642" height="407" alt="image" src="https://github.com/user-attachments/assets/1cdb3322-131c-4b88-9237-cad8a9621a06" />




```python
#Shi-Tomasi Corner Detection

corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
```


```python
# Mark detected corners
# Display image
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y,), 3, (255,0,0), -1)
    
plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7f0715ec0b90>


<img width="446" height="405" alt="image" src="https://github.com/user-attachments/assets/2bbc80e4-33ed-4fc6-8e2b-6e3a3c0f616e" />



```python
# Detect strongest 100 corners
# Display image
corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)
    
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7f0715e97b50>

<img width="712" height="407" alt="image" src="https://github.com/user-attachments/assets/7ae4a27b-be2b-43d7-9ccd-c35ca37a7466" />




### Edge Detection
In this analysis, we applied Canny edge detection to an image, sloth image, and applied varying parameters in vieweing the edges.

```python
# Import libary 
import cv2
```


```python
# Import numpy
import numpy as np
```


```python
# Import matplotlib
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# Read imported image
# Print image
img = cv2.imread("sloth.jpeg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fb021e6ff50>




<img width="546" height="385" alt="image" src="https://github.com/user-attachments/assets/ead241b9-a1fc-4306-9224-05226cfc2856" />




```python
# Apply Canny edge detection
# Fix threshold
# Print image
edges = cv2.Canny(image = img, threshold1= 127, threshold2 = 127)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb020dad6d0>




<img width="526" height="379" alt="image" src="https://github.com/user-attachments/assets/29f2a134-a823-4301-b8ba-a1da5366aa71" />




```python
# Calculate median pixel intensity of image
med_value = np.median(img)
med_value
```




    114.0




```python
# Creates upper and lower threshold
# Apply Canny detection with thresholds
# Print image
lower = int(max(0, 0.7*med_value))
upper = int(min(255, 1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb020ca9d10>



<img width="567" height="403" alt="image" src="https://github.com/user-attachments/assets/e1eda57a-13a0-4d15-bde5-6829d08e452e" />




```python
# Increase upper threshold
# Print image
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper +100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb020c145d0>




<img width="573" height="385" alt="image" src="https://github.com/user-attachments/assets/4a30a519-d5f9-4aa7-920d-d7517688971e" />




```python
# Apply blur to image
# Print image
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb020bff590>




<img width="583" height="388" alt="image" src="https://github.com/user-attachments/assets/e87e4fb9-8d02-4a03-8074-2216a0ad11b8" />




```python
# Apply stronger blur
# Print image
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb020b6d350>



<img width="545" height="395" alt="image" src="https://github.com/user-attachments/assets/5508166c-9216-4a6b-9436-30ff4e7fb3ed" />



```python
# Apply blur
# Increase upper threshold
# Print image
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper +50)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb020add0d0>



<img width="532" height="382" alt="image" src="https://github.com/user-attachments/assets/4c4ebe38-6326-47d7-b614-d50531eba218" />




```python
# Apply blur
# Increase upper threshold again
# Print image
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper +100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb020ac2dd0>




<img width="544" height="391" alt="image" src="https://github.com/user-attachments/assets/36b0cdeb-1e56-4d37-88f4-a466797e373d" />




```python
# Apply stronger blur
# Increase upper threshold
# Print image
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 60)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb0209a2910>



<img width="517" height="399" alt="image" src="https://github.com/user-attachments/assets/9b744542-e1b5-4e8c-97c5-5466cfee96a2" />






```python
# Apply stronger blur
# Print image
blurred_img = cv2.blur(img, ksize = (8,8))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fb02098f4d0>



<img width="561" height="400" alt="image" src="https://github.com/user-attachments/assets/f4ed68d3-45fc-42cb-bddd-d4522941907c" />



## Feature Detection
### Feature Matches
In this analysis, we took images of Apple Jacks and a variety of cereal and used ORB and SIFT to detect and match, using BF and FLANN matching, within one another with visualizations through line connections.

```python
# Import libaries
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
# Define functions for image display
def display(img,cmap = 'gray'):
    fig = plt.figure(figsize =(12, 10))
    ax = fig.add_subplot(111)
    ax.imshow(img,cmap = 'gray')
```


```python
# Image display in gray
apple_jacks = cv2.imread("apple_jacks.jpg", 0)
display(apple_jacks)
```
<img width="609" height="930" alt="image" src="https://github.com/user-attachments/assets/dd5f541a-3d30-4f4d-8249-94a392a86567" />




```python
# Image display in gray
cereals = cv2.imread('all_cereal.jpg', 0)
display(cereals)
```
<img width="1168" height="857" alt="image" src="https://github.com/user-attachments/assets/28e6a75d-d889-4ac8-8349-e431ce530640" />




```python
# Creation of ORB detection
orb = cv2.ORB_create()

kp1,des1 = orb.detectAndCompute(apple_jacks, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
```


```python
# Creation of BF matching
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches = bf.match(des1, des2)
```


```python
# Sort matches
matches = sorted(matches, key = lambda x:x.distance)
```


```python
# Top 25 matches made
apple_jacks_matches = cv2.drawMatches(apple_jacks, kp1, cereals, kp2, matches[:25], None, flags = 2)
```


```python
# Display image with matches
display(apple_jacks_matches)
```

<img width="1160" height="714" alt="image" src="https://github.com/user-attachments/assets/2404ff72-3e52-41d4-8eec-b7001ed0cce8" />



```python
# Create SIFT detection
sift = cv2.SIFT_create()
```


```python
# Give keypoints and descriptors for matching
kp1 = sift.detect(apple_jacks, None)
kp1, des1 = sift.compute(apple_jacks, kp1)

kp2 = sift.detect(cereals, None)
kp2, des2 = sift.compute(cereals, kp2)
```


```python
# Apply BF matching
bf = cv2.BFMatcher()
macthes = bf.knnMatch(des1, des2, k=2)
```


```python
# Filter good images
good = []

for match1, match2 in macthes:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
# Print statistics from matching
print('Length of total matches:', len(macthes))
print('Length of good matches:', len(good))
```

    Length of total matches: 4316
    Length of good matches: 120



```python
# Draw matches
sift_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags = 2)
display(sift_matches)
```

<img width="1157" height="708" alt="image" src="https://github.com/user-attachments/assets/e2b929ff-a0ba-47ef-a11c-539350c35a81" />



```python
# Create FLANN matching
sift = cv2.SIFT_create()

kp1 = sift.detect(apple_jacks, None)
kp1, des1 = sift.compute(apple_jacks, kp1)

kp2 = sift.detect(cereals, None)
kp2, des2 = sift.compute(cereals, kp2)
```


```python
# Indiccate parameters
flann_index_KDtree = 0
index_params = dict(algorithm=flann_index_KDtree, trees = 5)
search_params = dict(checks=50)
```


```python
# Load and Perform FLANN matching
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k=2)

good = []

for match1, match2, in matches:
    if match1.distance < 0.75*match2.distance:
        good.append([match1])
```


```python
# Draw matches
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, good, None, flags = 0)
display(flann_matches)
```

<img width="1271" height="721" alt="image" src="https://github.com/user-attachments/assets/6263de34-2246-4797-b8c0-96654915177e" />



```python
# Create SIFT detection again 
sift = cv2.SIFT_create()

kp1 = sift.detect(apple_jacks, None)
kp1, des1 = sift.compute(apple_jacks, kp1)

kp2 = sift.detect(cereals, None)
kp2, des2 = sift.compute(cereals, kp2)
```


```python
# Define FLANN parameters
flann_index_KDtree = 0
index_params = dict(algorithm = flann_index_KDtree, trees = 5)
search_param = dict(checks = 50)
```


```python
# Create FLANN matching
flann = cv2.FlannBasedMatcher(index_params, search_params)

matches = flann.knnMatch(des1, des2, k = 2)
```


```python
matchesMask = [[0,0] for i in range(len(matches))]
```


```python
# Define parameters
for i, (match1, match2) in enumerate(matches):
    if match1.distance <0.75*match2.distance:
        matchesMask[i] = [1,0]

draw_params = dict(matchColor = (0,255,0),
                  singlePointColor = (255,0,0),
                  matchesMask = matchesMask,
                  flags = 0)
```


```python
# Display matches
flann_matches = cv2.drawMatchesKnn(apple_jacks, kp1, cereals, kp2, matches, None, **draw_params)

display(flann_matches)
```
<img width="1152" height="716" alt="image" src="https://github.com/user-attachments/assets/7d4fe00e-5178-424e-99e2-5f41adc29208" />



### Object Detection
In this analysis, we took images of sunflowers for finding smaller train images in larger test images using using matching methods with heat map visualization and rectangles for specific location detection.

```python
# Import libraries
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
```


```python
%matplotlib inline
```


```python
# Load training image 
full = cv2.imread('training_sunflower.jpg')
```


```python
# Convert image to color
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
# Display image
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7ff6c62b0e10>

<img width="441" height="474" alt="image" src="https://github.com/user-attachments/assets/580f3481-9e16-4bc8-9e16-6f58a169f79f" />



```python
# Load testing image
test = cv2.imread('testing_sunflower.jpg')
```


```python
# Convert image to color
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
```


```python
# Display image
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7ff6c621dc10>
<img width="466" height="351" alt="image" src="https://github.com/user-attachments/assets/4c6ddeb7-8e71-4159-88a5-9d4211d713c8" />




```python
# Print image dimensions
print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Test image shape: (1555, 1404, 3)
    Training image shape: (310, 510, 3)



```python
# Define matching methods
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_CCORR_NORMED', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']

```


```python
# Instruction loop for each macthing method: copy, method, template, location, box detection
# Heatmap display with detection results
for m in methods:
    plt.figure()
    
    test_copy = test.copy()
    method = eval(m)
    
    res= cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255,0,0), 10)
    
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title("Detection of template")
    
    plt.suptitle(m)
    
    plt.show
    print('\n')
    print('\n')
```
<img width="659" height="895" alt="image" src="https://github.com/user-attachments/assets/ac9f6ba4-a502-4045-866d-8476ae7d6c5e" />

<img width="662" height="893" alt="image" src="https://github.com/user-attachments/assets/019cdd06-bbf9-49e8-a24f-5ed5fd717b23" />

<img width="635" height="932" alt="image" src="https://github.com/user-attachments/assets/e0d2e34b-6392-4602-a9ee-d80e3359b8f7" />

