# NatalieLeBlanc_Python2_Portfolio
This is the portfolio of the Python 2 codes that I learned through Winter 2025-2026 in the Python 2 class (BISC 4503).


## Sequence Objects
In this analysis, 

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
# Pint sequence
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
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
translate(my_string)
```




    'AVMGRWKGGRAAG*'





## Sequence Annotations


## Sequence I/O


## Multiple Sequence Alignment 



## Blast 


## Challenge 1


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
# Detect starongest 100 corners
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

