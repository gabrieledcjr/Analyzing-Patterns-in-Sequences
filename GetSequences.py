#!/usr/bin/python -tt
'''
Create on January 2, 2014

@author: Gabriel de la Cruz
'''

# [BioPython] third party modules use to manipulate sequence files
from Bio import SeqIO

import GlobalVars
from AssembleFASTAFiles import WriteToFileStopCodons


def GetSequences(fileName="", fileType="fasta", path="", isPatterns=False, name=""):

  ctr = 0          # counter
  seqIds = []      # allele pattern's id names
  sequences = []   # allele pattern's amino acid sequences

  ctr2 = 0         # counter for sequences with stop codons
  seqsWithCodons = []

  # Parse through the file, normally a fasta file
  # and retrive all sequences
  for seq_record in SeqIO.parse(path + fileName, fileType):

    # Allele Patterns
    if isPatterns:
      seqIds.insert(ctr, seq_record.id)
      sequences.insert(ctr, str(seq_record.seq))

    # Sample Sequences
    # Don't allow duplicates from input sample and
    # write to separate file sequences with stop codons
    else:

      if '*' in str(seq_record.seq):
        seqsWithCodons.insert(ctr2, [seq_record.id, str(seq_record.seq)])
        ctr += 1
        ctr2 += 1
        continue

      isAdded = False  # flag variable to check whether sequence has been added

      # loop through the parsed sequences to check if current sequence exist
      for i in range(0, len(seqIds)):
        if str(seq_record.seq) == sequences[i]:
          seqIds[i] = seqIds[i] + ',' + seq_record.id
          isAdded = True
          break

      if not isAdded:
        seqIds.insert(ctr, seq_record.id)
        sequences.insert(ctr, str(seq_record.seq))

    ctr += 1

  if GlobalVars.DEBUG:
    print

  # Write to file sequences with stop codons
  if not isPatterns and ctr2 > 0:
    WriteToFileStopCodons(fileName.split('/')[-1].split('.')[0], seqsWithCodons)

  return seqIds, sequences, ctr, ctr2
# end GetSequences()
