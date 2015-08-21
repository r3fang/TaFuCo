import sys
import collections
import random 
from itertools import groupby

genes_names = [];
transcripts_names = collections.defaultdict(list);
transcripts_seq = {};
gene_pair = [];

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

def main():
    for (head, seq) in fasta_iter("transcript.fa"):
        gene_name = head.replace(">", "").split('|')[0] 
        transcript_id = head.strip().split('|')[1]
        transcripts_names[gene_name].append(transcript_id)
        transcripts_seq[head.replace(">", "")] = seq
        if gene_name not in genes_names:
            genes_names.append(gene_name)

    for i in xrange(1, 50):
        # randomly choose 2 genes
        gene1 = genes_names[random.randint(0, len(genes_names)-1)]
        gene2 = genes_names[random.randint(0, len(genes_names)-1)]
        # randomly choose 2 genes and the gene pair has to unique        
        if gene1 != gene2 and tuple(sorted([gene1, gene2])) not in gene_pair:
            (gene1, gene2) = tuple(sorted([gene1, gene2])) # change order
            name1 = gene1 + "|" + random.choice(transcripts_names[gene1])
            name2 = gene2 + "|" + random.choice(transcripts_names[gene2])
            seq1 = transcripts_seq[name1]
            seq2 = transcripts_seq[name2]
            strlen1 = random.randint(50, len(seq1)-1)
            strlen2 = random.randint(0, len(seq2)-50)
            seq = seq1[:strlen1]+seq2[strlen2:]
            print ">"+gene1+"_"+gene2+'\t'+str(strlen1)+'\t'+str(len(seq))
            print seq

if __name__=="__main__":
    main()
    