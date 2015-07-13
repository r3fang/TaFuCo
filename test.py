import foo, itertools
#seqs = dict(foo.FastaReader('sample_data/test.seq'))
seq = [1, 3, 4, 5, 6, 8]
#print foo.ReverseComplement('aaaattGGc')
# foo.index('sample_data/exons.fa.gz', 15)
foo.predict('sample_data/exons.fa.gz', 'sample_data/U2OS-1-ABGHI_S29_L001_R1_001.fastq.gz', 'sample_data/U2OS-1-ABGHI_S29_L001_R2_001.fastq.gz')
#foo.TRY()
# {BAG}-[sample_data]:(add large exons.fa)/[src/index.c]:(fix bugs if seq shorted than k)
