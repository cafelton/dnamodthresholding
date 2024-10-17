# example usage [python trainHMM.py -c ys18_chrom_150U_subsample.bam -s 9merROC_metrics_newtrain.txt -o test -t filtered-dynamic]
# don't need -s with simple thresholding


import pysam, math, random, argparse
import matplotlib.pyplot as plt
from hmmlearn import hmm
import numpy as np

compbase = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}

typesOfMods = {'5mC':[('C', 0, 'm')], '5hmC': [('C', 0, 'h')], '5fC': [('C', 0, 'f')], '5caC': [('C', 0, 'c')],
               '5hmU': [('T', 0, 'g')], '5fU': [('T', 0, 'e')], '5caU': [('T', 0, 'b')],
               '6mA': [('A', 0, 'a'), ('A', 0, 'Y')], '8oxoG': [('G', 0, 'o')], 'Xao': [('N', 0, 'n')]}

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


def getcomp(seq):
    newseq = []
    for base in seq: newseq.append(compbase[base])
    return ''.join(newseq)#newseq[::-1]

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

def processSamfile(filename, thresholdtype, training_dict = None):
	samfile = pysam.AlignmentFile(filename, "rb")
	d = 0
	allknownstates = []
	for s in samfile:#.fetch('chrII'):
	    if s.is_mapped and not s.is_secondary and not s.is_supplementary:
	        d += 1
	        if d%10000 == 0: print(d, 'reads processed')
	        if d > 10: break
	        base_qualities = s.query_qualities


	        ##If scoring using modification
	        modinfo, readseq = getMods(posstag, s)
	        if not modinfo: continue
	        seqlen = len(readseq)

	        ###THIS WHOLE NEXT SECTION IS VARIOUS THRESHOLDING, REPLACE WITH YOUR METHOD

	        knownstates, knownpos = [], []
	        ##If do need to know kmer to threshold
	        if sum(base_qualities) / len(base_qualities) > 10:

	        	if thresholdtype == 'none':
	        		for pos in range(seqlen):
	        			if pos in modinfo:
	        				knownstates.append(modinfo[pos])

	        	elif thresholdtype == 'simple':
	        		for pos in range(seqlen):
	        			if pos in modinfo:
	        				modval = 1 if modinfo[pos]/255. > threshold else 0
	        				knownstates.append(modval)
	        				knownpos.append(pos)

	        	elif thresholdtype == 'quality-none':
	        		 for i in base_qualities:
	        		 	knownstates.append(50-i)

	        	elif thresholdtype == 'quality-all': #change threshold (0.6) to best qual score thresh
	        		 for i in range(seqlen):
	        		 	qualscore = base_qualities[i]
	        		 	normalized_score = (50 - qualscore) / 50
	        		 	if normalized_score > 0.58:
	        		 		knownstates.append(1)
	        		 	else: 
	        		 		knownstates.append(0)

	        	elif thresholdtype == 'quality-a': #change threshold (0.6) to best qual score thresh
	        		for seqpos in modinfo:
	        			if readseq[seqpos] == 'A':
	        				qualscore = base_qualities[seqpos]
	        				qualval = 1 if ((50-qualscore)/50) > 0.6 else 0
	        				knownstates.append(qualval)

	        	elif thresholdtype == 'unfiltered-dynamic' or thresholdtype == 'filtered-dynamic':
	        		for i in range(halfkmersize-1, seqlen - halfkmersize):
	        			if i in modinfo:
	        				thiskmer = readseq[i - (halfkmersize - 1):i + halfkmersize]
	        				modpred = modinfo[i]
	        				if thiskmer in training_dict: 
	        					modval = 1 if modpred/255. > training_dict[thiskmer] else 0
	        					knownstates.append(modval)

	        ##for training model
	        allknownstates.extend(knownstates)
	samfile.close()
	return allknownstates

def trainModel(allknownstates, out):
	##this is for training model
	nreads = len(allknownstates)//10000
	allknownstates2 = [allknownstates[(x*10000):(x*10000)+10000] for x in range(nreads)]
	print(nreads, len(allknownstates))

	random.shuffle(allknownstates2)
	X_train = np.array(allknownstates2[:int(len(allknownstates2)*.7)])
	X_validate = np.array(allknownstates2[int(len(allknownstates2)*.7):])

	best_score = best_model = None
	n_fits = 10
	for idx in range(n_fits):
	    model = hmm.CategoricalHMM(n_components=2, random_state=idx,init_params='')
	    t1, t2 = random.uniform(0.9, 1), random.uniform(0.9, 1)
	    model.transmat_ = np.array([[t1,1-t1],
	                                [1-t2,t2]])
	    model.startprob_ = np.array([.8, 0.2])
	    if thresholdtype in ['simple', 'unfiltered-dynamic', 'filtered-dynamic', 'quality-none', 'quality-a', 'quality-all']:
	    	model.emissionprob_ = np.array([[.7,.3],[.3, .7]])
	    elif thresholdtype == 'none':
	    	model.emissionprob_ = np.array([[.5] + [.5/255]*255, [.5/255]*255 + [.5]])
	    else:
	    	model.emissionprob_ = np.array([[.5] + [.5/49]*49, [.5/49]*49 + [.5]])
	    model.fit(X_train)
	    score = model.score(X_validate)
	    print(f'Model #{idx}\tScore: {score}')
	    if best_score is None or score > best_score:
	        best_model = model
	        best_score = score

	out = open(out + '.tsv', 'w') ###CHANGE THIS FILENAME
	out.write('Starting Matrix\t' + ','.join(str(x) for x in best_model.startprob_) + '\n')
	out.write('Transmission Matrix\t' + '\t'.join(','.join([str(y) for y in x]) for x in best_model.transmat_) + '\n')
	out.write('Emission Matrix\t' + '\t'.join(','.join([str(y) for y in x]) for x in best_model.emissionprob_) + '\n')
	out.close()

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--chromatinFile','-c', type=str,action='store',help='path to chromatin.bam')
	parser.add_argument('--trainingFile', '-s', type=str,action='store',help='path to training file (.txt) generated in calcROCmetrics.py (dynamic thresholding only)')
	parser.add_argument('--kmerSize', '-k', default=9, type=int, action='store',help='size of kmer (default: 9)')
	parser.add_argument('--modType','-m', default = '6mA', type=str, action='store',help='modification type (default: 6mA)')
	parser.add_argument('--thresholdingMethod', '-t', type = str, action='store', help = 'choose thresholding method: none, simple, unfiltered-dynamic, filtered-dynamic, quality-none, quality-a, quality-all')
	parser.add_argument('--outFile', '-o', type = str, action='store', help = 'name of output file')

	args = parser.parse_args()

	global modtype, posstag, kmersize, halfkmersize, thresholdtype, threshold
	modtype = args.modType
	posstag = typesOfMods[modtype]
	kmersize = args.kmerSize
	halfkmersize = math.ceil(kmersize/2)
	thresholdtype = args.thresholdingMethod
	outFile = args.outFile

	filename = args.chromatinFile
	trainingFile = args.trainingFile


	if thresholdtype == 'simple':
		threshold = float(input('Enter threshold: '))

	if thresholdtype == 'unfiltered-dynamic' or thresholdtype =='filtered-dynamic':
		if not args.trainingFile or '.txt' not in args.trainingFile:
			print('Error: Please enter valid training set ending with .txt (-s train.txt)')
			exit(1)
		else:
			training_dict = createTrainingDict(trainingFile)
			allknownstates = processSamfile(filename, thresholdtype, training_dict)
	else:
		allknownstates = processSamfile(filename, thresholdtype)

	trainModel(allknownstates, outFile)


if __name__ == "__main__":
	main()