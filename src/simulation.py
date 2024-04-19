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


def insert_into_nucleus(nuclear_fhand, sequences, out_fpath):
    insertions = []
    nucls = "".join([line.rstrip() for line in nuclear_fhand if not line.startswith(">")])
    positions_already_taken = []
    for name, sequence in sequences.items():
        invalid = True
        while invalid:
            start = random.randrange(len(nucls))
            end = start+(len(sequence))
            if end-start < len(nucls) and start not in positions_already_taken:
                invalid = False
                positions_already_taken.append(start)
        nucls = nucls[:start] + sequence + nucls[start:]
        insertions.append((start, name, len(sequence)))
    nucls = wrap("".join(nucls), 80)
    with open(out_fpath / "Nuclear_with_insertions.fasta", "w") as out_fhand:
        out_fhand.write(">Nuclear_with_insertions\n")
        for line in nucls:
            out_fhand.write(line+"\n")
            out_fhand.flush()
    return insertions


def get_sequences(input_fhand, starting_fhand, freq1, freq2, length, out_path):
    subseqs = {}
    record = SeqIO.read(input_fhand, "fasta")
    out_fhand = open(out_path / "subseqs_not_mutated.fasta", "w")
    out_fhand_mutated = open(out_path / "subseqs_mutated.fasta", "w")
    for line in starting_fhand:
        start = int(line.rstrip().split("_")[-1])
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


def main():
    freq1 = float(argv[4])
    freq2 = float(argv[5])
    length = int(argv[6])
    out_path = Path(argv[7])
    if not out_path.exists():
        out_path.mkdir(parents=True, exist_ok=True)
    with open(argv[1]) as input_fhand:
        with open(argv[2]) as starting_fhand:
            sequences = get_sequences(input_fhand, starting_fhand, freq1, freq2, length, out_path)
    with open(argv[3]) as nuclear_fhand:
        insertions = insert_into_nucleus(nuclear_fhand, sequences, out_path)
    with open(out_path / "summary.txt", "w") as out_fhand:
        for line in insertions:
            write_ = "{}\t{}\t{}\n".format(str(line[0]), str(line[1].replace(">", "")), str(line[2]))
            out_fhand.write(write_)


if __name__ == "__main__":
    main()
