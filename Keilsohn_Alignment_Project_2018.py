# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 21:07:11 2018

@author: William Keilsohn
"""

# Import Packages
#import Bio #http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc4
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


# Load in data
## European Earwig CO1 data
dataDic = {}


for seq_record in SeqIO.parse('sequence.gb', 'genbank'): #Loops through the file and pulls out all the sequence data
    dataDic[seq_record.id] = seq_record.seq._data
    
# Account for unequal sequence length
## https://stackoverflow.com/questions/31152011/multiple-sequence-alignment-with-unequal-string-length
lenLis = []
for k, v in dataDic.items():
    lenLis.append(len(v))
maxLength = max(lenLis)

extendedLis = []
for k, v in dataDic.items():
    extendedLis.append(v.ljust(maxLength, '-'))

# Apply the aplhabet conotations
alphaLis = []
for i in extendedLis:
    alphaLis.append(Seq(i, IUPAC.unambiguous_dna))

# Make a new dictionary with the adjusted sequences
seqDic = {}
idLis = []

for k, v in dataDic.items():
    idLis.append(k)

dataLen = len(idLis)
    
for x in range(0, dataLen):
    seqDic[idLis[x]] = alphaLis[x]

# Alighn the Data
dataLis = []
for key, value in seqDic.items():
    dataLis.append(SeqRecord(value, id = key))

# "Standard" Multiple Alignment according to Biopython   
dataAlign = MultipleSeqAlignment(dataLis) #Align the data
AlignIO.write(dataAlign, 'Alignment.phy', 'phylip') #For now this seems to work... I think.

### For the following alignments, the data must be in a FASTA file format
SeqIO.write(dataLis, "dataFile.fasta", 'fasta')
