from pathlib import Path
from subprocess import run
from sys import argv

def run_minimap2(genome, sequences, paf_fpath):
    cmd = ["minimap2", "-x map-pb", "--secondary no", genome,  
           sequences, ">", paf_fpath]
    run(" ".join(cmd), shell=True, capture_output=True)


def get_minimap_hits_on_insertions(insertions_source, minimap_hits):
    mapped_insertions = {insertion_name: [] for insertion_name in insertions_source}
    with open(minimap_hits) as minimap_fhand:
        for line in minimap_fhand:
            line = line.split()
            readName = line[0]
            readLength = line[1]
            readStart = line[2]
            readEnd = line[3]
            organelle_start = line[7]
            organelle_end = line[8]
            nuclMatches = line[9]
            alnLength = line[10]
            for insertion, positions in insertions_source.items():
                x = range(organelle_start, organelle_end)
                y = range(positions[0], positions[1])
                num_overlapping = len(set(x).intersection(y))
                if num_overlapping > 0:
                    aln_info = {"readName": readName, "readLength": readLength,
                                "readStart": readStart, "readEnd": readEnd,
                                "organelleStart": organelle_start, "organelleEnd": organelle_end,
                                "numMatches": nuclMatches, "alnLength": alnLength,
                                "overlappingNucls": num_overlapping}
                    mapped_insertions[insertion].append(aln_info)
    return mapped_insertions


def get_insertions_source(summary):
    insertions_source = {}
    with open(summary) as summary_fhand:
        for line in summary_fhand:
            insertion_name = line.rstrip().split()[1]
            start, end = insertion_name.split("_")[-1].split(":")
            insertions_source[insertion_name] = (start, end)
    return insertions_source

def main():
    out_path = Path(argv[1])
    summary = out_path / "summary.txt"
    genome_fpath = argv[2]
    sequences_fpath = out_path / "sd_0001.fastq"
    mapping_output = str(out_path / "reads_mapped_against_organelle.paf")
    insertions_source = get_insertions_source(summary)
    run_minimap2(genome_fpath, sequences_fpath, mapping_output)
    overlaps = get_minimap_hits_on_insertions(insertions_source, mapping_output)
    print(overlaps)


if __name__ == "__main__":
    main()