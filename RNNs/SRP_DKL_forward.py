###############################################################################
# Code by Christoph Aurnhammer, based on                                      #
#    https://github.com/pytorch/examples/tree/master/word_language_model      #
# Citation: Aurnhammer & Frank (2019), Neuropsychologia.                      #
# LIG_2:                                                                      #
# This code averages all instances of a model snapshot and                    #
# then computes surprisal and the Kullbach-Leibler Divergence                 #
# on the experimental stimuli from Frank (2013)                               #
# (Christoph Aurnhammer, 05.04.2019                                           #
# for Aurnhammer, Frank (upcoming)                                            #
###############################################################################

import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.nn.functional import softmax
import argparse
import pandas
import glob
import re
import numpy as np
from math import log, exp, isnan
# script data.py
import data


def parse_args():
    parser = argparse.ArgumentParser(description='PyTorch ENCOW Language Model')
    parser.add_argument('--data', type=str, default='./corpus/',
                        help='location of the data corpus')
    parser.add_argument('--bptt', type=int, default=42,
                        help='sequence length')
    parser.add_argument('--cuda', action='store_true',
                        help='use CUDA')
    arguments = parser.parse_args()
    return arguments


# Get snapshots absolute paths from dir
def get_paths(directory, names_snapshot):
    snapshots = []
    for snap in names_snapshot:
        for repetition in range(0, 6):
            found_checkpoints = glob.glob(directory + snap + str(repetition))
            # If list is not empty
            if found_checkpoints:
                for cp in sorted(found_checkpoints, key=numerical_sort):
                    snapshots.append(cp)
    return snapshots


def numerical_sort(value):
    # Numerical sort from here on
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


def prepare_outfile(out_path):
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
    dataframe = pandas.DataFrame()
    dataframe['sent_nr'] = sent_nr
    dataframe['word_pos'] = word_pos
    dataframe['word'] = words
    dataframe['item'] = pandas.read_csv('./input/itemnumbers_frank2013.csv', delimiter='\t', header=None)
    dataframe['ENCOW_log_freq'] = pandas.read_csv('./input/ENCOWfreqs_frank2013.csv', delimiter='\t')
    dataframe.to_csv(out_path, sep='\t', index=False)
    return dataframe


def batchify(dt, bsz, arguments):
    # Work out how cleanly we can divide the dataset into bsz parts.
    nbatch = dt.size(0) // bsz
    # Trim off any extra elements that wouldn't cleanly fit (remainders).
    dt = dt.narrow(0, 0, nbatch * bsz)
    # Evenly divide the data across the bsz batches.
    dt = dt.view(bsz, -1).t().contiguous()
    if arguments.cuda:
        dt = dt.cuda()
    return dt


def repackage_hidden(h):
    """Wraps hidden states in new Variables, to detach them from their history."""
    if type(h) == Variable:
        return Variable(h.data)
    else:
        return tuple(repackage_hidden(v) for v in h)


def get_batch(source, i):
    seq_len = min(args.bptt, len(source) - 1 - i)
    dt = Variable(source[i:i + seq_len])
    target = Variable(source[i + 1:i + 1 + seq_len].view(-1))
    return dt, target


def get_eos(arguments):
    # chr: define end of sentence index
    if arguments.cuda:
        eos_tensor = torch.cuda.LongTensor([corpus.dictionary.word2idx['<eos>']])
    else:
        eos_tensor = torch.LongTensor([corpus.dictionary.word2idx['<eos>']])
    return eos_tensor


def forward_model(rnn, sequence, types, args):
    # Initialise hidden state of PLM to zeros for new sequence
    hidden_true = rnn.init_hidden(1)

    # List of probability distributions over w_t and w_t1
    w_t_list = []
    w_t1_list = []

    # For each word in the sentence (starting from <sos>)
    for item in sequence:
        # Reformat item (technicality)
        if args.cuda:
            item = torch.cuda.LongTensor([[int(item)]])
        else:
            item = torch.LongTensor([[int(item)]])

        # Model current item.
        # This is returns the "true" output / hidden states corresponding
        # to the actually occuring items in the stimuli
        output_true, hidden_true = rnn(item, hidden_true)

        # Collect current P(w_t|w_1..._t-1) probability distribution
        p_wt_dist = softmax(output_true, dim=-1).data[0][0]

        # For P(w_t+1|w_1...t):
        # Allocate array with vocabulary size as rows and columns
        if args.cuda:
            probs = torch.cuda.FloatTensor(np.empty([types, types]))
        else:
            probs = torch.FloatTensor(np.empty([types, types]))

        # For each possible possible w_t
        for wt in range(0, types):
            # Select probability of current w_t
            p_wt = p_wt_dist[wt]

            # Run using current wt and rnn hidden state produced after last true item
            if args.cuda:
                output_wt1, hidden_wt1 = rnn(torch.cuda.LongTensor([[wt]]), hidden_true)
            else:
                output_wt1, hidden_wt1 = rnn(torch.LongTensor([[wt]]), hidden_true)

            # Collect current P(w_t+1|w_1...t) distribution
            p_wt1_dist = softmax(output_wt1, dim=-1).data[0][0]

            # Enter as column into matrix
            # Each cell is the probability of the j-th w_t1
            # multiplied by the prob of the current w_t
            probs[:, wt] = p_wt1_dist * p_wt

        # Compute sum per row, leaving a single vector with
        # one probability per possible w_t1
        p_wt1_dist = probs.sum(dim=1)

        # Append to output lists
        w_t_list.append(p_wt_dist)
        w_t1_list.append(p_wt1_dist)
    return w_t_list, w_t1_list


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


def compute_surprisal(dists, seq_targets):
    # accept a list of probability distributions and the indices of the correct word and compute surprisal for each word
    sent_surprisal = []
    for target, prob in zip(seq_targets.data, dists):
        sent_surprisal.append([corpus.dictionary.idx2word[target], round(log(float((prob[target]))), 4) * -1])
    return sent_surprisal


def compute_kld(dist_nextword, dist_plusone, seq_targets):
    # accept two lists of probability distributions and compute the Kullback-Leibler Divergence for each pair
    del dist_nextword[0]
    del dist_plusone[-1]
    # We can't compute the KLD for the first words of sentence
    sent_kld = [[corpus.dictionary.idx2word[seq_targets[0]], None]]
    seq_targets = seq_targets[1:]
    for dist_nw, dist_po, target in zip(dist_nextword, dist_plusone, seq_targets):
        cross_entropy = -sum(dist_nw * dist_po.log())
        plusone_entropy = -sum(dist_nw * dist_nw.log())
        kld = cross_entropy - plusone_entropy
        sent_kld.append([corpus.dictionary.idx2word[target], round(kld.item(), 4)])
    return sent_kld


def add_to_df(input_lol, dataframe, snap, metric_name):
    # Clean up sentences: remove commas and words with clitics (ensures equal lengths of sentences)
    for i_index, sentence in enumerate(input_lol):
        for j_index, word in enumerate(sentence):
            if word[0] == ',':
                del (input_lol[i_index][j_index])
            elif '\'' in word[0]:
                del (input_lol[i_index][j_index])
    # Add metrics to new column
    new_col = []
    for row in dataframe.iterrows():
        word_value = input_lol[row[1][0] - 1][row[1][1] - 1]
        word = row[1][2].lower()
        word = word.strip(',')
        if word_value[0] == word and word_value[1] is not None:
            new_col.append(float(word_value[1]))
        else:
            new_col.append(None)
    assert len(df) == len(new_col)
    dataframe[metric_name + '_' + snap] = new_col
    return dataframe


def evaluate(surprisal_values):
    N = 0
    Psum = 0

    for surp in surprisal_values:
        if isnan(surp):
            pass
        else:
            N += 1
            Psum += -surp
    print("Evaluated: Perplexity {}".format(exp(-1 / N * Psum)))
    return exp(-1 / N * Psum)


def store_eval(wt_perf, wt1_perf):
    output = pandas.DataFrame()
    output['snapshots'] = ['1K', '3K', '10K', '30K', '100K', '300K', '1M', '3M', '6.47M']
    output['wt'] = wt_perf
    output['wt1'] = wt1_perf
    output.to_csv('srp_entr/PLM_ppl.csv')
    print(output)

if __name__ == "__main__":
    # Parse command line arguments
    args = parse_args()

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
    # We are only select LSTM models (in case other matching models are in dir)
    snapshot_names = ['LSTM' + name for name in snapshot_names]
    # Shortened names used for output
    snapshot_output_names = ['1k', '3k', '10k', '30k', '100k', '300k', '1M', '3M', '6-47M']

    # Directories where the models are stored (normal LM = nextword; two steps ahead = plusone)
    nextword_dir = './output/'
    # Get full paths
    nextword_paths = get_paths(nextword_dir, snapshot_names)

    # Prepare output file
    outfile_path = './srp_entr/SRP_DKL_snapshots.txt'
    df = prepare_outfile(outfile_path)

    # Load corpus into memory (requires data.py)
    print('Loading corpus from {}'.format(args.data))
    corpus = data.Corpus(args.data, args.bptt)
    ntypes = len(corpus.dictionary)
    eos = get_eos(args)
    test_data = batchify(corpus.test, bsz=1, arguments=args)

    # Evaluation output
    wt_ppl = []
    wt1_ppl = []
    ######################################
    # Loop through snapshots, sentences  #
    # Compute surprisal, KLD             #
    # Add data to file for each snapshot #
    ######################################
    for snap_name, snap_out_name in zip(snapshot_names, snapshot_output_names):
        snap_paths_nextword = [path for path in nextword_paths if snap_name in path]

        sent_counter = 0
        wt_surprisal_list = []
        KLD_list = []
        wt1_surprisal_list = []
        for x in range(0, test_data.size(0) - 1, args.bptt):
            data, targets = get_batch(test_data, x)
            # chr: cut off data at end of sentence
            for j in range(len(data)):
                if (int(data[j].data) == int(eos)) is True:
                    sent_counter += 1
                    data = data[:j, :1]

                    # normal targets
                    wt_targets = targets[:j]
                    # targets for wt1
                    wt1_targets = targets[1:j]
                    break

            # Wt list
            wt_list = []
            wt1_list = []

            # Forward modeling of two steps ahead probability
            for model_path in snap_paths_nextword:
                with open(model_path, 'rb') as f:
                    rnn_model = torch.load(f)
                rnn_model.eval()
                wt_seq, wt1_seq = forward_model(rnn_model, data, ntypes, args)  # return two lists of dists
                wt_list.append(wt_seq)
                wt1_list.append(wt1_seq)

            wt_avg = average_sent(wt_list)
            wt1_avg = average_sent(wt1_list)

            # For each word in the current sentence get surprisal and KLD
            # Compute surprisal for each word wt
            wt_surprisal_list.append(compute_surprisal(wt_avg, wt_targets))
            # Compute KLD for each word (targets used to make word step identifiable)
            KLD_list.append(compute_kld(wt_avg, wt1_avg, wt_targets))
            # Compute surprisal for each word wt1 (only for evaluation purposes)
            wt1_surprisal_list.append(compute_surprisal(wt1_avg, wt1_targets))
            print("Classified {} sentences at snapshot {}".format(sent_counter, snap_out_name), end='\r')

        # Add surprisal to output
        df = add_to_df(wt_surprisal_list, df, snap_out_name, metric_name='srp')
        df = add_to_df(KLD_list, df, snap_out_name, metric_name='KLD')
        df.to_csv(outfile_path, sep='\t', index=False)

        # List surprisal for evaluation
        wt_ppl.append(evaluate([y for x in wt_surprisal_list for y in x]))
        wt1_ppl.append(evaluate([y for x in wt1_surprisal_list for y in x]))

    # Store evaluation output
    store_eval(wt_ppl, wt1_ppl)
