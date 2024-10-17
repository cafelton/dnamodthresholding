import pysam, argparse, math
import numpy as np
from sklearn import metrics

def getcomp(seq):
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq)#newseq[::-1]

typesOfMods = {'5mC':[('C', 0, 'm')], '5hmC': [('C', 0, 'h')], '5fC': [('C', 0, 'f')], '5caC': [('C', 0, 'c')],
               '5hmU': [('T', 0, 'g')], '5fU': [('T', 0, 'e')], '5caU': [('T', 0, 'b')],
               '6mA': [('A', 0, 'a'), ('A', 0, 'Y')], '8oxoG': [('G', 0, 'o')], 'Xao': [('N', 0, 'n')]}

compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

def getMods(posstag, s):
    ##this section is just retrieving the modification information in a robust way
    if s.is_reverse: posstag = [(x[0], 1, x[2]) for x in posstag]
    ml = None
    if not s.modified_bases:  # isinstance(s.modified_bases, dict):
        return None, None
    for t in posstag:
        if t in s.modified_bases:
            ml = s.modified_bases[t]
            break
    if not ml:
        print(readname, 'does not have modification information', s.modified_bases.keys())
        return None, None

    ###this determines whether As with no score are registered as unknown or unmodified
    if s.has_tag('MM'):
        skippedBase = -1 if s.get_tag('MM').split(',', 2)[0][-1] == '?' else 0
    elif s.has_tag('Mm'):
        skippedBase = -1 if s.get_tag('Mm').split(',', 2)[0][-1] == '?' else 0
    else:
        return None, None

    seq = s.query_sequence  # s.get_reference_sequence().upper() #s.query_sequence
    if s.is_reverse:  ###need to get compliment of sequence, but not reverse!!
        seq = getcomp(seq)

    ##get the position of every A in the read sequence
    seqApos = set()
    c = 0
    for b in seq:
        if b == 'A':
            seqApos.add(c)
        c += 1

    ###if As are not already accounted for in the modification calls, add the determined skippedBase code at that position
    ml = dict(ml)
    for i in seqApos:
        if i not in ml:
            ml[i] = skippedBase

    return ml, seq

def getScores(filename, expVal):
    samfile = pysam.AlignmentFile(filename, "rb")
    d = 0

    pred_list = []
    exp_list = []
    ROCdict = {}
    for s in samfile:
        base_qualities = s.query_qualities
        if base_qualities is not None and sum(base_qualities)/len(base_qualities) > 10:
            if s.is_mapped and not s.is_secondary and not s.is_supplementary:
                modinfo, readseq = getMods(posstag, s)
                if not modinfo: continue

                seqlen = len(readseq)

                for i in range(halfkmersize-1, seqlen - halfkmersize):
                    if i in modinfo:
                        thiskmer = readseq[i - (halfkmersize - 1):i + halfkmersize]
                        modpred = modinfo[i]/255

                        if thiskmer not in ROCdict:
                            ROCdict[thiskmer] = ([], [])

                        ROCdict[thiskmer][0].append(modpred)
                        ROCdict[thiskmer][1].append(expVal)
                d += 1
                if d%1000 == 0: print(d, 'reads processed')
    for kmer in ROCdict:
        modpreds = ROCdict[kmer][0]
        modexps = [expVal] * len(modpreds)
        ROCdict[kmer] = (modpreds, modexps)
    samfile.close()
    print(filename + ' processed')
    print()
    return ROCdict

def calcROC(outFile, posDict, negDict):
    with open(str(outFile), 'w') as outfile: 
        for kmer in posDict.keys():
            if kmer in negDict.keys():
                combined_predictedScore = negDict[kmer][0] + posDict[kmer][0]
                combined_expectedScore = negDict[kmer][1] + posDict[kmer][1]
                fpr, tpr, thresholds = metrics.roc_curve(combined_expectedScore, combined_predictedScore, pos_label=1)
                auc = metrics.auc(fpr, tpr)
                
                diff = tpr - fpr
                maxDiff = 0
                for d, fp, tp, th in zip(diff, fpr, tpr, thresholds):
                    if d > maxDiff:
                        maxDiff = d
                        maxFPR = fp
                        maxTPR = tp
                        maxThreshold = th
                outfile.write(f"{kmer}\t{maxThreshold}\t{maxTPR}\t{maxFPR}\t{auc}\n")
    print('done')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--negativeControl','-n', type=str,action='store',help='path to negative_control.bam')
    parser.add_argument('--positiveControl','-p', type=str,action='store',help='path to postive_control.bam')
    parser.add_argument('--kmerSize', '-k', default = 9, type=int, action='store',help='size of kmer (default: 9')
    parser.add_argument('--modType','-m', default = '6mA', type=str, action='store',help='modification type (default: 6mA)')
    parser.add_argument('--outFile', '-o', type= str, action ='store', help = 'name of output file (.txt)')

    args = parser.parse_args()

    if not args.outFile or '.txt' not in args.outFile:
        print('Error: Please enter a valid output file name ending with .txt')
        exit(1) 
    elif not args.negativeControl or not args.positiveControl or '.bam' not in args.negativeControl or '.bam' not in args.positiveControl:
        print('Error: Please enter valid positive and negative control files ending with .bam')
        exit(1)

    global kmersize, modtype, posstag, halfkmersize
    modtype = args.modType
    posstag = typesOfMods[modtype]
    kmersize = args.kmerSize 
    halfkmersize = math.ceil(kmersize/2)

    negDict = getScores(args.negativeControl, 0)
    posDict = getScores(args.positiveControl, 1)
    calcROC(args.outFile, posDict, negDict)

if __name__ == "__main__":
    main()