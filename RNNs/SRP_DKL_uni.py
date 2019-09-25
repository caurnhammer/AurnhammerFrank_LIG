
###############################################################################
# Code by Christoph Aurnhammer, based on                                      #
#    https://github.com/pytorch/examples/tree/master/word_language_model      #
# Citation: Aurnhammer & Frank (2019), Neuropsychologia.                      #
# LIG_1:                                                                      #
# This code averages all instances of a model snapshot and                    #
# then computes surprisal and the Kullbach-Leibler Divergence                 #
# on the experimental stimuli from Frank (2013)                               #
###############################################################################

import argparse
import torch
from torch.autograd import Variable
import torch.nn as nn
import data
import pandas
from torch.nn.functional import log_softmax, softmax
from math import exp, log
import glob
import re

#########################################
# Parse user input (command line flags) #
#########################################
parser = argparse.ArgumentParser(description='PyTorch ENCOW Language Model')
parser.add_argument('--data', type=str, default='./corpus',
                    help='location of the data corpus')
parser.add_argument('--bptt', type=int, default=41,
                    help='sequence length')
parser.add_argument('--cuda', action='store_true',
                    help='use CUDA')
args = parser.parse_args()
# Notify user if a cuda decive could be used
if torch.cuda.is_available():
    if not args.cuda:
        print("WARNING: You have a CUDA device, so you should probably run with --cuda")

#####################################
# Define all project specific input #
#####################################
# Define the names of the snapshots and the directories where they are found
snapshot_names = ['_1000_', '_3000_', '_10000_', '_30000_', '_100000_', '_300000_', '_1000000_',
                  '_3000000_', '_6470000_']
# We are only interested
snapshot_names = ['LSTM'+name for name in snapshot_names]
# Shortened names for output
snapshot_output_names = ['1k', '3k', '10k', '30k', '100k', '300k', '1M', '3M', '6-47M']

# Directories where the models are stored (normal LM = nextword; two steps ahead = plusone)
nextword_dir = './output/'

# Get snapshots absolute paths from dir
def get_paths(dir):
    snapshots = []
    for snap_name in snapshot_names:
        for repetition in range(0, 6):
            found_checkpoints = glob.glob(dir + snap_name + str(repetition))
            def numericalSort(value):
                # Numerical sort from here on
                numbers = re.compile(r'(\d+)')
                parts = numbers.split(value)
                parts[1::2] = map(int, parts[1::2])
                return parts

            if found_checkpoints != []:
                for cp in sorted(found_checkpoints, key=numericalSort):
                    snapshots.append(cp)
    return snapshots

nextword_paths = get_paths(nextword_dir)
#plusone_paths = get_paths(plusone_dir)

# Define path of output file
outfile = './srp_entr/SRP_DKL_unigram.txt'
# Load items from the experimental stimuli
with open('./corpus/test.txt') as inputfile:
    inputfile = inputfile.read()
    inputfile = inputfile.replace(' ,', ',')
    inputfile = inputfile.replace(' n\'t', 'n\'t')
    inputfile = inputfile.replace(' \'', '\'')
    inputfile = inputfile.split('\n')
    del inputfile[-1]
    inputfile = [sentence.split(' ') for sentence in inputfile]

sent_nr = []
word_pos = []
words = []
for sent_ind, sentence in enumerate(inputfile):
    for word_ind, word in enumerate(sentence):
        sent_nr.append(sent_ind + 1)
        word_pos.append(word_ind + 1)
        words.append(word)

# Prepare output file
df = pandas.DataFrame()
df['sent_nr'] = sent_nr
df['word_pos'] = word_pos
df['word'] = words
df['item'] = pandas.read_csv('./input/itemnumbers_frank2013.csv', delimiter='\t', header = None)
df['ENCOW_log_freq'] = pandas.read_csv('./input/ENCOWfreqs_frank2013.csv', delimiter = '\t')
df.to_csv(outfile, sep='\t', index=False)

##########################
# Prepare classification #
##########################
# Load corpus into memory (calls data.py)
print('Loading corpus from {}'.format(args.data))
corpus = data.Corpus(args.data, args.bptt)
ntokens = len(corpus.dictionary)
# Define the criterion for evaluation
criterion = nn.CrossEntropyLoss()

# chr: define end of sentence index
if args.cuda:
    a = torch.cuda.LongTensor([corpus.dictionary.word2idx['<eos>']])
else:
    a = torch.LongTensor([corpus.dictionary.word2idx['<eos>']])

# Define helper functions
def batchify(data, bsz):
    # Work out how cleanly we can divide the dataset into bsz parts.
    nbatch = data.size(0) // bsz
    # Trim off any extra elements that wouldn't cleanly fit (remainders).
    data = data.narrow(0, 0, nbatch * bsz)
    # Evenly divide the data across the bsz batches.
    data = data.view(bsz, -1).t().contiguous()
    if args.cuda:
        data = data.cuda()
    return data

def repackage_hidden(h):
    """Wraps hidden states in new Variables, to detach them from their history."""
    if type(h) == Variable:
        return Variable(h.data)
    else:
        return tuple(repackage_hidden(v) for v in h)

def get_batch(source, i, evaluation=False):
    seq_len = min(args.bptt, len(source) - 1 - i)
    data = Variable(source[i:i + seq_len])
    target = Variable(source[i + 1:i + 1 + seq_len].view(-1))
    return data, target

######################
# Classify test data #
######################
# chr: Define list for surprisal values per sentence
def classify_sent(model, data):
    hidden = model.init_hidden(bsz=1)
    output, hidden = model(data, hidden)
    output_flat = output.view(-1, ntokens)
    output_sftmx = softmax(output_flat, dim=1)
    return output_sftmx

def average_sent(sentence_output):
    average_dist = []
    num_words = len(sentence_output[0])
    for word_outer in range(0, num_words):
        word_dists = []
        for model in sentence_output:
            word_dists.append(model[word_outer])
        word_avg = sum(word_dists)/len(word_dists)
        average_dist.append(word_avg)
    return average_dist

def compute_surprisal(dists, targets):
    # accept a list of probability distributions and the indices of the correct word and compute surprisal for each word
    sent_surprisal = []
    for target, prob in zip(targets.data, dists):
        sent_surprisal.append([corpus.dictionary.idx2word[target], round(log(float((prob[target]))), 4) * -1])
    return sent_surprisal

def compute_KLD(dist_nextword, dist_plusone, targets):
    # accept two lists of probability distributions and compute the Kullback-Leibler Divergence for each pair
    del dist_nextword[0]
    del dist_plusone[-1]
    # We can't compute the KLD for the first words of sentence
    sent_KLD = [[corpus.dictionary.idx2word[targets[0]], None]]
    targets = targets[1:]
    for dist_nw, dist_po, target in zip(dist_nextword, dist_plusone, targets):
        cross_entropy = -sum(dist_nw * dist_po.log())
        plusone_entropy = -sum(dist_nw * dist_nw.log())
        KLD = cross_entropy - plusone_entropy
        sent_KLD.append([corpus.dictionary.idx2word[target], round(KLD.item(), 4)])
    return sent_KLD

def compute_unigram_KLD(dist_nextword, freq_dist, targets):
    # Accept one lists of probability distributions and compute the Kullback-Leibler Divergence
    # for the prob dist for the next word (RNN) and the word afterwards (unigram based)
    # We can't compute the KLD for the first words of sentence
    del dist_nextword[0]
    sent_KLD = [[corpus.dictionary.idx2word[targets[0]], None]]
    targets = targets[1:]

    dist_po = freq_dist
    #dist_po = torch.cuda.FloatTensor([6.6908] * len(dist_nextword[0]))
    for dist_nw, target in zip(dist_nextword, targets):
        cross_entropy = -sum(dist_nw * dist_po.log())
        plusone_entropy = -sum(dist_nw * dist_nw.log())
        KLD = cross_entropy - plusone_entropy
        sent_KLD.append([corpus.dictionary.idx2word[target], round(KLD.item(), 8)])
    return sent_KLD

def add_to_df(input_lol, df, snap_name, metric_name):
    # Clean up sentences: remove commas and words with clitics (ensures equal lengths of sentences)
    for i, sentence in enumerate(input_lol):
        for j, word in enumerate(sentence):
            if word[0] == ',':
                del (input_lol[i][j])
            elif '\'' in word[0]:
                del (input_lol[i][j])
    # Add metrics to new column
    new_col = []
    for row in df.iterrows():
        word_value = input_lol[row[1][0] - 1][row[1][1] - 1]
        word = row[1][2].lower()
        word = word.strip(',')
        if word_value[0] == word and word_value[1] != None:
            new_col.append(float(word_value[1]))
        else:
            new_col.append(None)
    assert len(df) == len(new_col)
    df[metric_name + '_' + snap_name] = new_col
    return df

### Load vocabulary dictionary
def load_vocab(vocab_path):
    with open(vocab_path, encoding="utf-8") as vocab:
        x = eval(vocab.read())
    print('> chr: read vocabulary', type(x))
    return x

######################################
# Loop through snapshots, sentences  #
# Compute surprisal, KLD             #
# Add data to file for each snapshot #
######################################
test_data = batchify(corpus.test, bsz = 1)
for snap_name, snap_out_name in zip(snapshot_names, snapshot_output_names):
    snap_paths_nextword = [path for path in nextword_paths if snap_name in path]
    
    sent_counter = 0
    surprisal_list = []
    KLD_list = []

    # unigram dist: load vocabulary from py dicrtionary file
    frequencies = load_vocab('./input/ENCOW_train_vocab.txt')
    N = sum(frequencies.values())
    # new list to store a freq dist, that matches the order of vocab entries of the NN
    freqs_corpusorder = []
    # Loop through entries in NN vocab
    for i in range(0, len(corpus.dictionary)):
        # get the word type of the dict entry
        word = corpus.dictionary.idx2word[i]
        # get the freq of the dict entry
        try:
            freqs_corpusorder.append(frequencies[word] / N)
        except KeyError:
            freqs_corpusorder.append(0.00000001)

    # convert the ordered freq dist to torch tensor
    if args.cuda:
        freqs_corpusorder = torch.cuda.FloatTensor(freqs_corpusorder)
    else:
        freqs_corpusorder = torch.FloatTensor(freqs_corpusorder)

    for x in range(0, test_data.size(0) - 1, args.bptt):
        data, targets = get_batch(test_data, x, evaluation=True)
        # chr: cut off data at end of sentence
        for j in range(len(data)):
            if (data[j].data == a)[0] == True:
                sent_counter += 1
                data = data[:j, :1]
                targets = targets[:j]
                break

        # For next-word prediction models
        output_sftmx_nextword = []
        for model_path in snap_paths_nextword:
            with open(model_path, 'rb') as f:
                model = torch.load(f)
            model.eval()
            output_sftmx_nextword.append(classify_sent(model, data))
        # list of 6 models with each n_targets x vocabulary size
        avg_nextword = average_sent(output_sftmx_nextword)

        # For each word in the current sentence get surprisal and KLD
        # Compute surprisal for each word
        surprisal_list.append(compute_surprisal(avg_nextword, targets))
        # Compute KLD for each word
        KLD_list.append(compute_unigram_KLD(avg_nextword, freqs_corpusorder, targets))
        print("Classified {} sentences at snapshot {}".format(sent_counter, snap_out_name), end='\r')

    df = add_to_df(surprisal_list, df, snap_out_name, metric_name='srp')
    df = add_to_df(KLD_list, df, snap_out_name, metric_name='KLD')
    df.to_csv(outfile, sep='\t', index=False)

