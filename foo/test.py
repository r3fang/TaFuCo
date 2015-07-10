import foo, itertools
seqs = dict(foo.FastaReader('sample_data/test.seq'))
print seqs
seqs = dict(foo.FastaReader('sample_data/test.fa'))
print seqs
seq = [1, 3, 4, 5, 6, 8]
print foo.totalIter(seq)
print foo.totalIter(xrange(10))
print foo.ReverseComplement('aaaattGGc')
foo.index('sample_data/test.fa', 30)
foo.predict('sample_data/test.fa', 'sample_data/U2OS-1-ABGHI_S29_L001_R1_001.fastq.gz')
