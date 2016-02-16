

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

transcripts = []
with open("TR2012b.fa.txt", 'r') as fasta:
    for seq_record in SeqIO.parse(fasta, 'fasta', alphabet=IUPAC.ambiguous_dna):
        transcripts.append(seq_record)
        
open_frames = []

for i, transcript in enumerate(transcripts):
    if transcript.seq.tostring().find("AAAAAAAAAAAA") == -1: #if no poly-A, do all six frames
        for_1 = transcript.seq.translate(to_stop=True)
        for_2 = transcript.seq[1:].translate(to_stop=True)
        for_3 = transcript.seq[2:].translate(to_stop=True)
        rev_1 = transcript.seq.reverse_complement().translate(to_stop=True)
        rev_2 = transcript.seq[:-1].reverse_complement().translate(to_stop=True)
        rev_3 = transcript.seq[:-2].reverse_complement().translate(to_stop=True)
    
        frames = [for_1, for_2, for_3, rev_1, rev_2, rev_3]
        longestORF = max(frames, key=len)
        seqrecord = SeqRecord(longestORF, id=transcript.id)
        open_frames.append(seqrecord)
        
    else:  #if poly-A tail is found, only do forward translation
        for_1 = transcript.seq.translate(to_stop=True)
        for_2 = transcript.seq[1:].translate(to_stop=True)
        for_3 = transcript.seq[2:].translate(to_stop=True)
        frames = [for_1, for_2, for_3]
        longestORF = max(frames, key=len)
        seqrecord = SeqRecord(longestORF, id=transcript.id)
        open_frames.append(seqrecord)
        
with open("spodoptera_proteome_unambig_long.fasta", "w") as proteome:
    SeqIO.write([record for record in open_frames if record.seq.tostring().find("X") == -1
                and record.seq.tostring().find("B") == -1
                and record.seq.tostring().find("Z") == -1
                and len(record) > 100
                ], proteome, "fasta")
