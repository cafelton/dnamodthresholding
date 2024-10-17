import pysam, math, random, sys, argparse
from hmmlearn import hmm
import numpy as np

####Usage: python3 test.py file.bam (must be sorted and indexed)

# sys.argv = ['text.py', 'file.bam']

compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

typesOfMods = {'5mC':[('C', 0, 'm')], '5hmC': [('C', 0, 'h')], '5fC': [('C', 0, 'f')], '5caC': [('C', 0, 'c')],
               '5hmU': [('T', 0, 'g')], '5fU': [('T', 0, 'e')], '5caU': [('T', 0, 'b')],
               '6mA': [('A', 0, 'a'), ('A', 0, 'Y')], '8oxoG': [('G', 0, 'o')], 'Xao': [('N', 0, 'n')]}

def getcomp(seq):
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq)#newseq[::-1]

def processBedFile(tssFile):
    tssscores = []
    for i in range(2001):
        tssscores.append([])

    tsspos = {}
    file = open(tssFile)
    for line in file:
        if line[:5] != 'track':
            line = line.strip().split()
            chr, dir = line[0], line[4] #line[5] - changed for +1 file
            pos = int(line[1]) if dir == '+' else int(line[2])
            if chr not in tsspos:
                tsspos[chr] = []
            tsspos[chr].append((pos-1000, pos+1000, dir))
    return tsspos, tssscores

def createTrainingDict(trainingFile):
    training_dict = {}
    filtered_training_dict = {}
    to_filter = []

    with open(trainingFile) as training_data: # set up to filter by auc
        for line in training_data:
            kmer, thresh, tpr, fpr, auc = line.rstrip().split()
            score = float(tpr)-float(fpr)
            training_dict[kmer] = (float(thresh))
            if float(auc) >= 0.8:
                filtered_training_dict[kmer] = (float(thresh))

        ### Filter by bottom % of TPR-FPR:
        #     to_filter.append((kmer, score, thresh))
        # to_filter.sort(key=lambda x: x[1])
        # cutoff = int(0.05 * len(to_filter)) # 0.1 for 10%s

        # for i, (kmer, score, threshold) in enumerate(to_filter):
        #     if i >= cutoff:
        #         filtered_training_dict[kmer] = (float(thresh))

    if thresholdtype == 'unfiltered-dynamic':
        return training_dict

    elif thresholdtype == 'filtered-dynamic':
        return filtered_training_dict

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

def qual_all_threshold(qualscore):
    thisscore = 0 if (50-qualscore)/50. < 0.58 else 1
    return thisscore

def qual_nothresh(qualscore):
    thisscore = 50-qualscore
    return thisscore

def qual_a_threshold(qualscore):
    thisscore = 0 if (50-qualscore)/50. < 0.6 else 1
    return thisscore

def dynamicThreshold(scores, seq, halfkmersize, kmerToThresh):
    outscores = {}
    for i in range(halfkmersize-1, len(seq) - halfkmersize):
        if i in scores:
            thiskmer = seq[i - (halfkmersize - 1):i + halfkmersize]
            if thiskmer in kmerToThresh:
                thisthresh = kmerToThresh[thiskmer]
                thisscore = 0 if scores[i] < thisthresh else 1
                outscores[i] = thisscore
    return outscores

def simplethreshold(scores, threshold):
    outscores = {}
    for i in scores:
        thisscore = 0 if scores[i] / 255. < threshold else 1
        outscores[i] = thisscore
    return outscores

def no_threshold(scores):
    return(scores)

def loadHMM(hmmmatrix):
    startmat, emmat, transmat = [], [], []
    for line in open(hmmmatrix):
        line = line.rstrip().split('\t')
        vals = [[float(y) for y in x.split(',')] for x in line[1:]]
        if line[0] == 'Starting Matrix': startmat = np.array(vals[0])
        elif line[0] == 'Transmission Matrix': transmat = np.array(vals)
        elif line[0] == 'Emission Matrix': emmat = np.array(vals)

    print(startmat, transmat, emmat)
    model = hmm.CategoricalHMM(n_components=2)
    model.startprob_ = startmat
    model.transmat_ = transmat
    model.emissionprob_ = emmat

    return model

def processSamfile(chromatinFile, model, tsspos, tssscores, training_dict):
    samfile = pysam.AlignmentFile(chromatinFile, "rb")
    c = 0
    posMod = {}
    qualscores = []
    for s in samfile:#.fetch('chrIII'): ####If wanting to look at whole genome, just: for s in samfile:
        chr = s.reference_name
        if not s.is_secondary and chr in tsspos: #and not s.reference_name == 'chrM'
            c += 1
            if c % 1000 == 0: print(c)
            alignstart, alignend = s.reference_start, s.reference_end
            hassite = False
            readname = s.query_name
            base_qualities = s.query_qualities
            if sum(base_qualities) / len(base_qualities) <= 10: continue
            cigar = s.cigartuples

            if thresholdtype == 'quality-all':
                ml = {x:qual_all_threshold(base_qualities[x]) for x in range(len(base_qualities))}

            elif thresholdtype == 'quality-a':
                ml = {}
                ogscores, readseq = getMods(posstag, s)
                for x in range(len(readseq)):
                    if readseq[x] == 'A':
                        ml[x] = qual_a_threshold(base_qualities[x])

            elif thresholdtype == 'quality-none':
                ml = {x:qual_nothresh(base_qualities[x]) for x in range(len(base_qualities))}

            else:
                ogscores, readseq = getMods(posstag, s)
                if thresholdtype in ['unfiltered-dynamic', 'filtered-dynamic']:
                    ml = dynamicThreshold(ogscores, readseq, halfkmersize, training_dict)
                elif thresholdtype == 'none':
                    ml = no_threshold(ogscores)
                else:
                    ml = simplethreshold(ogscores, threshold)

            knownstates, knownpos = [], []
            sorted_ml = dict(sorted(ml.items()))
            for position,score in sorted_ml.items():
                knownstates.append(score)
                knownpos.append(position)

            knownstates = np.array(knownstates).reshape(-1, 1)
            if len(knownstates) != 0:
                hiddenstates = model.predict(knownstates)
                newmodinfo = {knownpos[i]:hiddenstates[i] for i in range(len(knownpos))}

            ref, quer = 0, 0
            posOnGenome = []
            for block in cigar:
                if block[0] in {0, 7, 8}:  # match, consumes both
                    for i in range(block[1]):
                        if quer in newmodinfo: posOnGenome.append([ref + alignstart, newmodinfo[quer]])
                        ref += 1
                        quer += 1
                elif block[0] in {1, 4}:  # consumes query
                    quer += block[1]
                elif block[0] in {2, 3}:  # consumes reference
                    ref += block[1]
            dirtowrite = '-' if s.is_reverse else '+'
            posdict = dict(posOnGenome)

            for coord in tsspos[chr]:
                start, stop, dir = coord[0], coord[1], coord[2]
                if start > alignstart and stop < alignend:
                    for pos in range(start, stop+1):
                        if pos in posdict:
                            thisscore = posdict[pos]

                            if dir == '+':
                                tssscorepos = pos-start
                            else:
                                tssscorepos = 2000 - (pos-start)
                            tssscores[tssscorepos].append(thisscore)

    samfile.close()

    tssscores = [sum(x)/len(x) if len(x) > 0 else 0 for x in tssscores]
    return tssscores

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--annotationFile','-a', type=str, default = '+1Nucs.bed', action='store',help='path to annotated BED file')
    parser.add_argument('--chromatinFile','-c', type=str,action='store',help='path to chromatin.bam alignment file')
    parser.add_argument('--hmmMatrix', '-hmm', type=str,action='store',help='path to hmm starting matrix (output from trainHMM.py)')
    parser.add_argument('--trainingFile', '-s', type=str,action='store',help='path to training file (.txt) generated in calcROCmetrics.py (dynamic thresholding only)')
    parser.add_argument('--kmerSize', '-k', default=9, type=int, action='store',help='size of kmer (default: 9)')
    parser.add_argument('--modType','-m', default = '6mA', type=str, action='store',help='modification type (default: 6mA)')
    parser.add_argument('--thresholdingMethod', '-t', type = str, action='store', help = 'choose thresholding method: none, simple, unfiltered-dynamic, filtered-dynamic, quality-none, quality-a, quality-all')
    parser.add_argument('--outFile', '-o', type = str, action='store', help = 'name of output file')

    args = parser.parse_args()

    if args.outFile is None:
        args.outFile = f"{args.annotationFile.split('.bed')[0]}_{args.thresholdingMethod}_tssScores.txt"

    global modtype, posstag, kmersize, halfkmersize, thresholdtype, threshold
    modtype = args.modType
    posstag = typesOfMods[modtype]
    kmersize = args.kmerSize
    halfkmersize = math.ceil(kmersize/2)
    thresholdtype = args.thresholdingMethod

    model = loadHMM(args.hmmMatrix)
    tsspos, tssscores = processBedFile(args.annotationFile)
 
    if thresholdtype == 'simple':
        threshold = float(input('Enter threshold: '))

    training_dict = None
    if thresholdtype == 'unfiltered-dynamic' or thresholdtype =='filtered-dynamic':
        if not args.trainingFile or '.txt' not in args.trainingFile:
            print('Error: Please enter valid training set ending with .txt (-s train.txt)')
            exit(1)
        else:
            training_dict = createTrainingDict(args.trainingFile)
        
    tssfinalscores = processSamfile(args.chromatinFile, model, tsspos, tssscores, training_dict)

    out = open(args.outFile + '.tsv', 'w')
    for l in tssfinalscores:
        out.write(str(l) + '\n')

if __name__ == '__main__':
    main()


