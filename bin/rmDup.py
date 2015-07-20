import sys
import collections

BAG = collections.defaultdict(list)
fin = sys.stdin

line = ""
while True:
    line = fin.readline() 
    if len(line) == 0:
        break
    if '>' in line:
        [name, num] = line.strip().split()
        num = int(num)
        for i in xrange(num):
            read = fin.readline().strip() 
            BAG[name].append(read)

for name in BAG:
    tmp = set(BAG[name])
    print name, len(tmp)
    for read in tmp:
        print read
    
    