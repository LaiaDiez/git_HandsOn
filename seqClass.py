#!/usr/bin/env python

import sys, re
from argparse import ArgumentParser

parser = ArgumentParser(description = 'Classify a sequence as DNA or RNA')
parser.add_argument("-s", "--seq", type = str, required = True, help = "Input sequence")
parser.add_argument("-m", "--motif", type = str, required = False, help = "Motif")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()
args.seq = args.seq.upper()                 
if re.search('^[ACGTU]+$', args.seq): 
    if re.search('^[ACGT]+$', args.seq): # If the sequence have T and not U is DNA
        print ('The sequence is DNA') 
    elif re.search('^[ACGU]+$', args.seq): # If the sequence have U and not T is RNA
        print ('The sequence is RNA') 
    else:
        print ('The sequence is not DNA nor RNA') # If we see both, then it can't be RNA or DNA
else:
    print ('The sequence is not DNA nor RNA') 
if args.motif:
    args.motif = args.motif.upper()
    print(f'Motif search enabled: looking for motif "{args.motif}" in sequence "{args.seq}"... ', end = '')
    if re.search(args.motif, args.seq):
        print("FOUND") #if it find a motif it prints FOUND
    else:
        print("NOT FOUND") #if it doesn't find a motif prints NOT FOUND

