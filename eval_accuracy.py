import pysam, argparse, math
from sklearn import metrics
import matplotlib.pyplot as plt 


compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
def getcomp(seq):
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq)#newseq[::-1]

typesOfMods = {'5mC':[('C', 0, 'm')], '5hmC': [('C', 0, 'h')], '5fC': [('C', 0, 'f')], '5caC': [('C', 0, 'c')],
               '5hmU': [('T', 0, 'g')], '5fU': [('T', 0, 'e')], '5caU': [('T', 0, 'b')],
               '6mA': [('A', 0, 'a'), ('A', 0, 'Y')], '8oxoG': [('G', 0, 'o')], 'Xao': [('N', 0, 'n')]}

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

def nothresh(score): #threshold of 0
    return score

def simple_threshold(score): #threshold of 0.5
    if score >= threshold:
        return 1
    else:
        return 0

def dynamic_threshold(thiskmer, score, trainingDict): #excluding kmers with a TPR-FPR in bottom 10%
    if thiskmer in trainingDict: 
        if score >= trainingDict[thiskmer]:
            return 1
        else:
            return 0

def qual_thresh(qualscore):
    if qualscore < 25: 
        return 1
    else: 
        return 0

def getScores(filename, trainingFile):
    if thresholdtype == 'unfiltered-dynamic' or thresholdtype == 'filtered-dynamic':
        training_dict = {}
        filtered_training_dict = {}
        to_filter = []

        with open(trainingFile) as training_data:
            for line in training_data:
                split_line = line.strip().split()
                kmer, threshold, tpr, fpr = split_line[0], split_line[1], split_line[2], split_line[3]
                training_dict[kmer] = float(threshold)
                score = float(tpr)-float(fpr)
                to_filter.append((kmer, score, threshold))
        to_filter.sort(key=lambda x: x[1])
        cutoff = int(0.1 * len(to_filter))

        for i, (kmer, score, threshold) in enumerate(to_filter):
            if i >= cutoff:
                filtered_training_dict[kmer] = float(threshold)

    samfile = pysam.AlignmentFile(filename, "rb")
    d = 0

    pred_list = []
    expected_list = []
    best_list = []
    for s in samfile:
        base_qualities = s.query_qualities
        if base_qualities is not None and sum(base_qualities)/len(base_qualities) > 10:
            if s.is_mapped and not s.is_secondary and not s.is_supplementary:
                d += 1
                if d%100 == 0: print(d, 'reads processed')

                modinfo, readseq = getMods(posstag, s)
                if not modinfo: continue
                seqlen = len(readseq)

                if thresholdtype == 'none':
                    for seqpos in modinfo:
                        modpred = modinfo[seqpos]
                        predscore = modpred/255 
                        pred_list.append(nothresh(predscore))                    

                elif thresholdtype == 'simple':
                    for seqpos in modinfo:
                        modpred = modinfo[seqpos]
                        predscore = modpred/255 
                        pred_list.append(simple_threshold(predscore)) 

                elif thresholdtype == 'quality-all':
                    for i in range(seqlen):
                        qual = base_qualities[i]
                        pred_list.append(qual_thresh(qual))

                elif thresholdtype == 'quality-a':
                    for i in range(seqlen):
                        if readseq[i] == 'A':
                            qual = base_qualities[i]
                            pred_list.append(qual_thresh(qual))

                elif thresholdtype == 'unfiltered-dynamic':
                    for i in range(halfkmersize-1, seqlen - halfkmersize):
                        if i in modinfo:
                            thiskmer = readseq[i - (halfkmersize - 1):i + halfkmersize]
                            modpred = modinfo[i]
                            predscore = modpred/255
                            pred_list.append(dynamic_threshold(thiskmer, predscore, training_dict))

                elif thresholdtype == 'filtered-dynamic':
                    for i in range(halfkmersize-1, seqlen - halfkmersize):
                        if i in modinfo:
                            thiskmer = readseq[i - (halfkmersize - 1):i + halfkmersize]
                            modpred = modinfo[i]
                            predscore = modpred/255
                            pred_list.append(dynamic_threshold(thiskmer, predscore, filtered_training_dict))

                else:
                    print('Error: Enter valid thresholding method: none, simple, quality-all, quality-a, unfiltered-dynamic, filtered-dynamic')
                    exit(1)

    samfile.close()
    return pred_list


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--negativeControl','-n', type=str,action='store',help='path to negative_control_testset.bam')
    parser.add_argument('--positiveControl','-p', type=str,action='store',help='path to postive_control_testset.bam')
    parser.add_argument('--trainingSet','-s', type=str,action='store',help='path to training file (.txt) generated in calcROCmetrics.py (dynamic thresholding only)')
    parser.add_argument('--kmerSize', '-k', default = 9, type=int, action='store',help='size of kmer (default: 9')
    parser.add_argument('--modType','-m', default = '6mA', type=str, action='store',help='modification type (default: 6mA)')
    parser.add_argument('--thresholdingMethod', '-t', type = str, action='store', help = 'choose thresholding method: none, simple, quality-a, quality-all, unfiltered-dynamic, filtered-dynamic')

    args = parser.parse_args()

    if not args.negativeControl or not args.positiveControl or '.bam' not in args.negativeControl or '.bam' not in args.positiveControl:
        print('Error: Please enter valid positive and negative control files ending with .bam (-n neg_test.bam -p pos_test.bam)')
        exit(1)

    global kmersize, modtype, posstag, halfkmersize, thresholdtype, threshold
    modtype = args.modType
    posstag = typesOfMods[modtype]
    kmersize = args.kmerSize 
    halfkmersize = math.ceil(kmersize/2)
    thresholdtype = args.thresholdingMethod

    if thresholdtype == 'simple':
        threshold = float(input('Enter threshold: '))
    elif thresholdtype == 'unfiltered-dynamic' or thresholdtype =='filtered-dynamic':
        if not args.trainingSet or '.txt' not in args.trainingSet:
            print('Error: Please enter valid training set ending with .txt (-s train.txt)')
            exit(1)


    predListNeg = getScores(args.negativeControl, args.trainingSet) 
    predListPos = getScores(args.positiveControl, args.trainingSet) 

    predListNeg = [x for x in predListNeg if x is not None]
    predListPos = [x for x in predListPos if x is not None]

    expListNeg = sum(predListNeg)/len(predListNeg)
    expListPos = sum(predListPos)/len(predListPos)

    print('fpr', expListNeg)
    print('tpr', expListPos)

if __name__ == "__main__":
    main()
