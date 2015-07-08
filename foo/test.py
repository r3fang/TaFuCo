import foo, itertools
seqs = dict(foo.FastaReader('sample_data/test.seq'))
print seqs
seq = [1, 3, 4, 5, 6, 8]
print foo.totalIter(seq)
print foo.totalIter(xrange(10))
print foo.ReverseComplement('aaaattGGc')
