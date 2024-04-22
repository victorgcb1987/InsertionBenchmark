from subprocess import run

def run_minimap2(genome, sequences, paf_fpath):
    cmd = ["minimap2", "-cx map-pb", "-k 11", "--secondary no", "-t", "40", str(genome),  
           str(sequences), ">", str(paf_fpath)]
    print(" ".join(cmd))
    run(" ".join(cmd), shell=True)

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