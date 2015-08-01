import sys
import random 

def rev_com(s):
    return s.upper().replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]

def main():
    k = int(sys.argv[1])
    seq = open(sys.argv[2]).readlines()[1]
    seq_rev = rev_com(seq)    
    name = sys.argv[3]
    fout1 = open(name+"_R1.fq", 'w')
    fout2 = open(name+"_R2.fq", 'w')
    for i in xrange(len(seq) - k):
        if(random.randint(1, 10)<=4):
            fout2.write("@"+str(i)+'\n')
            fout2.write(seq[i:i+k]+'\n')
            fout2.write('+\n')
            fout2.write('+'*k+'\n')    
            fout1.write("@"+str(i)+'\n')
            fout1.write(seq[i:i+k]+'\n')
            fout1.write('+\n')
            fout1.write('+'*k+'\n')            
        
    fout1.close()
    fout2.close()
if __name__ == '__main__':
    main()
    