#
# UTILITY CODE FOR GENERATING CODEBOOK AND ENCODING/DECODING READS
# WEA (wallen@stanford.edu)
# Deisseroth Lab
# (C) Stanford University 2017
#

import numpy as np
import itertools


colors = [0, 1, 2, 3]  # red, orange, green, blue
bases = ['A', 'T', 'C', 'G']

colormap = {
    'AT': 3, 'CT': 2, 'GT': 1, 'TT': 0,
    'AG': 2, 'CG': 3, 'GG': 0, 'TG': 1,
    'AC': 1, 'CC': 0, 'GC': 3, 'TC': 2,
    'AA': 0, 'CA': 1, 'GA': 2, 'TA': 3
}

reverse_map = {
    0: ['AA', 'CC', 'GG', 'TT'],
    1: ['AC', 'CA', 'GT', 'TG'],
    2: ['AG', 'CT', 'GA', 'TC'],
    3: ['AT', 'CG', 'GC', 'TA']
}

start_pos = [1, 0, 4, 3, 2]

alexa_map = {0 : 'G', 1 : 'T', 2 : 'A', 3 : 'C'}

def hamming_distance(s1, s2):
    # Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def SeqToColor(seq):
    code = {'G' : 1,
            'T' : 2,
            'A' : 3,
            'C' : 4}
    return [code[v] for v in seq]

def DecodeAlexa(color_seq):
    return ''.join([alexa_map[int(c)] for c in color_seq])

def ColorspaceToString(x):
    return "".join([str(i) for i in SeqToColor(x)])

def DecodeSOLID(color_seq, start_base):
    color_seq = [int(i) for i in color_seq]
    seq = ""
    possible_bases = reverse_map[color_seq[0]]
    for p in possible_bases:
        if p[0] == start_base:
            seq = p
    for i in range(1, len(color_seq)):
        possible_bases = reverse_map[color_seq[i]]
        for p in possible_bases:
            if p[0] == seq[i]:
                seq += p[1]
    return seq


def ReactionsToSOLID(rxns):
    """ Given list of primers x ligations, re-order into proper solid colorspace """
    primers = len(rxns)
    ligations = len(rxns[0])
    colors = np.zeros((primers * ligations, 1))
    for r in range(len(rxns)):  # iterate over primers
        for l in range(len(rxns[r])):
            curr_start = start_pos[r] + 5 * l
            colors[curr_start] = rxns[r][l]
    return colors.flatten()


def EncodeSOLID(seq):
    color_seq = [colormap[seq[0:2]]]
    for i in range(2, len(seq)):
        color_seq.append(colormap[seq[(i - 1):(i + 1)]])
    return color_seq


def SimulateEncodeSOLID(seq):
    colors = []
    seq_rounds = int(np.floor(len(seq) / 5.))  # number of ligations/primer
    nprimers = 5
    ordered_colors = np.zeros((seq_rounds * nprimers, 1))
    for i in range(nprimers):  # primer rounds
        curr_colors = []
        for j in range(seq_rounds):
            curr_start = start_pos[i] + 5 * j
            curr_pair = colormap[seq[curr_start:(curr_start + 2)]]
            ordered_colors[curr_start] = curr_pair
            curr_colors.append(curr_pair)
        colors.append(curr_colors)
    return (ordered_colors.flatten(), colors)


def test_solid():
    x = 'ATCCGGATCCGTACTCGTAATGCTAT'
    (y, _) = SimulateEncodeSOLID(x)
    y = EncodeSOLID(x)
    z = DecodeSOLID(y, 'A')
    print(x == z)
    x = "ATCCGGATCCGT"
    (o, c) = SimulateEncodeSOLID(x)
    print(o == ReactionsToSOLID(c))


def SampleCode(nlen, alphabet=['A', 'T', 'C', 'G']):
    """ Naive way to sample random sequences... """
    return "".join([alphabet[i] for i in np.random.random_integers(0, len(alphabet) - 1, nlen)])


def SampleMultiCodes(ncodes, nlen, min_dist=0):
    """ Uniformly sample codes with a Hamming distance constraint. Naively will just continue sampling codes until fills up codebook. """
    all_codes = [SampleCode(nlen)]
    k = 1
    while k < ncodes:
        curr_code = SampleCode(nlen)
        do_add = True
        for c in all_codes:
            if hamming_distance(c, curr_code) <= min_dist:
                do_add = False
                break
        if do_add:
            all_codes.append(curr_code)
            k += 1
    return all_codes

def SaveCodebook(codes, fname):
    solid_colors = [np.array(EncodeSOLID(x)) for x in codes]
    good_codes = []
    for i in range(len(solid_colors)):
        if not np.all(solid_colors[i] == 0) and \
            not np.all(solid_colors[i] == 1) and \
            not np.all(solid_colors[i] == 2) and \
            not np.all(solid_colors[i] == 3):
            good_codes.append(codes[i])
        else:
            print(solid_colors[i])
            print(DecodeSOLID(solid_colors[i], 'C'))
    f = open(fname,'w')
    i = 0
    for code in good_codes:
        i += 1
        #print i
        f.write("%s\n" % code)
    f.close()

def LoadCodebook(fname):
    codes = []
    with open(fname, 'r') as f:
        for line in f:
            codes.append(line.strip())
    return codes

def GenerateAllCodes(nlen,bases=['A','T','G','C']):
    return [''.join(p) for p in itertools.product(bases, repeat=nlen)]


def EncodeSOLIDToString(s):
    return "".join([str(x) for x in EncodeSOLID(s)])

def FindCodebookInColorspace(ncodes, nlen, min_dist=0):
    all_codes = np.random.permutation(GenerateAllCodes(nlen))
    k = 1
    i = 1
    good_codes = [all_codes[0]]
    for niter in range(3): # give it 5 attempts
        print(len(good_codes))
        for curr_code in all_codes:
            gc_fraction = float(curr_code.count('G')+curr_code.count('C'))/len(curr_code)
            do_add = True
            for c in good_codes:
                if hamming_distance(EncodeSOLIDToString(c),EncodeSOLIDToString(curr_code)) <= min_dist or \
                        gc_fraction < 0.2 or gc_fraction > 0.8:
                    do_add = False
            if do_add:
                good_codes.append(curr_code)
                k += 1
            i += 1
            if k == ncodes:
                break
    return good_codes

def FindCodebook(ncodes, nlen, min_dist=0):
    all_codes = np.random.permutation(GenerateAllCodes(nlen))
    k = 1
    i = 1
    good_codes = [all_codes[0]]
    for niter in range(10): # give it 5 attempts
        for curr_code in all_codes:
            gc_fraction = float(curr_code.count('G')+curr_code.count('C'))/len(curr_code)
            do_add = True
            for c in good_codes:
                if hamming_distance(c,curr_code) <= min_dist or \
                        gc_fraction < 0.2 or gc_fraction > 0.8:
                    do_add = False
            if do_add:
                good_codes.append(curr_code)
                k += 1
            i += 1
            if k == ncodes:
                break
    return good_codes

def CodebookToSOLID(codebook):
    return [EncodeSOLID(c) for c in codebook]
