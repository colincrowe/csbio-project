import random
import numpy as np
import matplotlib.pyplot as plt

# helper functions

def get_rc(sequence):
    '''Takes a k-mer and returns the reverse compliment of that k-mer'''
    
    rc = ""
    length = len(sequence)
    for i in range(length)[::-1]:
        nc = sequence[i]
        if nc == 'A' or nc == 'a':
            rc += 'T'
        if nc == 'T' or nc == 't':
            rc += 'A'
        if nc == 'C' or nc == 'c':
            rc += 'G'
        if nc == 'G' or nc == 'g':
            rc += 'C'
            
    return rc

def get_kmer_rc_count(sequence, kmer):
    '''Takes a sequence and a kmer and returns the frequency of the kmer and its reverse complement in that sequence'''
    
    rc = get_rc(kmer)
    length = len(kmer)
    kmer_count = sequence.count(kmer)
    rc_count = sequence.count(rc)
        
    return kmer_count, rc_count

def at_cg_count(sequence):
    '''Takes a sequence and returns the frequency of A + T nucelotides and C + G nucleotides'''
    at = 0
    cg = 0
    for i in sequence:
        if i == 'C' or i == 'G':
            cg += 1
        else:
            at += 1
    return at, cg

def random_sequence(length):
    '''Returns a random DNA sequence of length specified by input'''

    random.seed()
    
    sequence = ""
    ncs = ['A','T','C','G']

    for i in range(length):
        select = random.randint(0, 3)
        sequence += ncs[select]
        
    return sequence
    
def random_sequence_cg(length, CG):
    '''Returns a random DNA sequence where CG nucleotides show up with frequency specified by input'''
    
    sequence = ""
    cg = ['C','G']
    at = ['A','T']
    
    for i in range(length):
        select = random.random()
        select2 = random.randint(0,1)
        if select <= CG:
            sequence += cg[select2]
        else:
            sequence += at[select2]
            
    return sequence

def random_inversion(sequence, num, length):
    '''Replaces NUM random k-mers of length LENGTH with their RCs'''
    for i in range(num):
        select = random.randint(0, len(sequence) - length)
        kmer = sequence[select : select + length]
        first_half = sequence[0 : select]
        second_half = sequence[select + length -1 : -1]
        rc = get_rc(kmer)
        sequence = first_half + rc + second_half
        
    return sequence

def get_3mers():
    '''Returns a list containing half the 3-mers; 
       (if a k-mer is in the list, its RC will not be in the list,
       but other than that it is comprehensive)'''
    
    nctides = ['A', 'T', 'C', 'G']
    three_mers = []
    rcs = []

    for i in range(4):
        for j in range(4):
            for k in range(4):
                three_mer = nctides[i] + nctides[j] + nctides[k]
                if three_mer not in rcs:
                    three_mers.append(three_mer)
                    rcs.append(three_mer)
                    rcs.append(get_rc(three_mer))
                    
    return three_mers

def read_BRCA1():
    '''Opens the file containing the BRAC1 genome
       and returns a string containing the genome'''

    BRCA1= ""

    with open('gene.fna') as f:
        meta_data = f.readline()
        while True:
            seq = f.readline().rstrip()
            if len(seq) == 0:
                break
            if seq[0] == '>':
                continue
            
            BRCA1+= seq

    return BRCA1

def count_kmers(sequence, kmers):
    '''Returns how many times each KMERS occurs in SEQUENCE'''
    tmer_count = []
    rc_count = []
    for i in kmers:
        tmer, rc = get_kmer_rc_count(sequence, i)
        tmer_count.append(tmer)
        rc_count.append(rc)

    return tmer_count, rc_count

def dump_graph():
    '''Creates the graph showing 3-mer and RC counts in BRCA1'''

    BRCA1 = read_BRCA1()

    three_mers = get_3mers()
    tmer_count, rc_count = count_kmers(BRCA1, three_mers)
    
    empty = []
    for i in range(32):
        empty.append(0)

    data = [tmer_count, rc_count]
    X = np.arange(32)
    fig = plt.figure()
    ax = fig.add_axes([0,0,2,2])
    ax.bar(X - 0.2, data[0], color = '#7eb54e', width = 0.4)
    ax.bar(X + 0.2, data[1], color = 'g', width = 0.4)
    ax.bar(three_mers, empty)
    ax.set_ylabel('Frequency')
    plt.legend(['3-mer Frequency', 'RC Frequency'],loc=1,prop={'size': 14})
    plt.ylim([0, 4200])

    return BRCA1