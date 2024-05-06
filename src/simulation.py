import random
from Bio import SeqIO
from sys import argv
from textwrap import wrap
from pathlib import Path


def mutate_subseq(subseq, freq=(0.4, 0.8)):
    mut_frq = 1 - random.uniform(freq[0], freq[1])
    nucls = list(subseq)
    num_changes = int(len(nucls)*mut_frq)
    positions = [number for number in range(0, len(subseq))]
    for index in random.sample(positions, num_changes):
        nucls[index] = random.choice([x for x in "ACTG" if x != nucls[index].upper()])      
    return "".join(nucls)
            

def write_sequences(sequences, output_fhand):
    for name, seq in sequences.items():
        output_fhand.write(name+"\n")
        output_fhand.write(str(seq)+"\n")
        output_fhand.flush()    


def insert_into_nucleus(nuclear_fhand, sequences, nuclear_positions,out_fpath):
    insertions = []
    nucls = "".join([line.rstrip() for line in nuclear_fhand if not line.startswith(">")])
    index = 0
    for name, sequence in sequences.items():
        start = nuclear_positions[index]
        nucls = nucls[:start] + sequence + nucls[start:]
        insertions.append((start, name, len(sequence)))
        index += 1
    nucls = wrap("".join(nucls), 80)
    with open(out_fpath / "Nuclear_with_insertions.fasta", "w") as out_fhand:
        out_fhand.write(">Nuclear_with_insertions\n")
        for line in nucls:
            out_fhand.write(line+"\n")
            out_fhand.flush()
    return insertions


def get_sequences(input_fhand, insertion_positions, freq1, freq2, length, out_path):
    subseqs = {}
    record = SeqIO.read(input_fhand, "fasta")
    out_fhand = open(out_path / "subseqs_not_mutated.fasta", "w")
    out_fhand_mutated = open(out_path / "subseqs_mutated.fasta", "w")
    for position in insertion_positions:
        start = position
        end = start+length
        subseq = record.seq[start:end]
        mutated_subseq = mutate_subseq(subseq, freq=(freq1, freq2))
        subseq_name = ">{}_{}:{}".format(record.id, str(start), str(end))
        subseqs[subseq_name] = mutated_subseq
        out_fhand.write(subseq_name+"\n")
        out_fhand_mutated.write(subseq_name+"\n")
        out_fhand.write(str(subseq)+"\n")
        out_fhand_mutated.write(str(mutated_subseq+"\n"))
        out_fhand.flush()
        out_fhand_mutated.flush()
    out_fhand.close()
    out_fhand_mutated.close()
    return subseqs


def get_random_sequence_positions(input_sequence, num_insertions=100):
    record = SeqIO.read(input_sequence, "fasta")
    sequence_positions = list(range(0, len(record.seq)))
    insertion_positions = random.sample(sequence_positions, num_insertions)
    return insertion_positions
