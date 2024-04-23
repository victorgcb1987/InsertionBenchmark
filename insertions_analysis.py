#!/usr/bin/env python

import argparse
import pandas as pd
import sys

from pathlib import Path


from src.dataframes_load import (load_insertions_source_as_dataframe, 
                                 load_minimap2_hits_as_dataframe, 
                                 load_read_positions_from_maf_into_dataframe)
from src.minimap2 import run_minimap2


IDENTITIES = [(0.41, 0.5), (0.51, 0.60), 
              (0.61, 0.7), (0.71, 0.8),
              (0.81, 0.9), (0.91, 1)]


#Generating program options
def parse_arguments():
    desc = "Analyze mapping results"
    parser = argparse.ArgumentParser(description=desc)
    
    help_input_dir = "(Required) input dir"
    parser.add_argument("--input_dir", "-i", type=str,
                        help=help_input_dir)
    
    help_ref_genome = "(Required) Fasta file from the organelle"
    parser.add_argument("--organelle_ref", "-r", type=str,
                        help=help_ref_genome)
    
    help_skip_mapping = "(Optional) Skip mapping step"
    parser.add_argument("--skip_mapping", "-s", action="store_true",
                        help=help_skip_mapping)
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


#Parse and return values given to options when running this program
def get_arguments():
    parser = parse_arguments()
    input_dir = Path(parser.input_dir)
    ref_genome = Path(parser.organelle_ref)
    skip_mapping = parser.skip_mapping
    return {"input_dir": input_dir, 
            "ref_genome": ref_genome,
            "skip_mapping": skip_mapping}


def classify_minimap_hits(insertions_source, minimap_hits):
    mapped_insertions = {insertion_name: [] for insertion_name in insertions_source}
    with open(minimap_hits) as minimap_fhand:
        for line in minimap_fhand:
            line = line.split()
            readName = line[0]
            readLength = int(line[1])
            readStart = int(line[2])
            readEnd = int(line[3])
            strand = line[4]
            sequenceName = line[5]
            organelle_start = int(line[7])
            organelle_end = int(line[8])
            nuclMatches = int(line[9])
            alnLength = int(line[10])
            for insertion_name, positions in insertions_source.items():
                source_range = range(positions[0], positions[1])
                aln_range = range(organelle_start, organelle_end)
                num_overlaps = len(set(source_range).intersection(aln_range))
                if num_overlaps > 0:
                    aln_info = {"readName": readName, "readLength": readLength,
                                 "readStart": readStart, "readEnd": readEnd,
                                 "organelleStart": organelle_start, "organelleEnd": organelle_end,
                                 "numMatches": nuclMatches, "alnLength": alnLength, 
                                 "numOverlaps": num_overlaps, "strand": strand}
                    mapped_insertions[insertion_name].append(aln_info)
    return mapped_insertions


def get_insertions_source(summary):
    insertions_source = {}
    with open(summary) as summary_fhand:
        for line in summary_fhand:
            insertion_name = line.rstrip().split()[1]
            start, end = insertion_name.split("_")[-1].split(":")
            insertions_source[insertion_name] = (int(start), int(end))
    return insertions_source


def get_reads_in_insertions(minimap2_df, insertions_df):
    dict_to_dataframe = {"insertion": []}
    for row in insertions_df.itertuples():
        reads = minimap2_df.loc[~((row.organelleEnd <= minimap2_df["organelleStart"]) | (minimap2_df["organelleEnd"] <= row.organelleStart))]
        #reads = minimap2_df.loc[(minimap2_df["organelleStart"] >= row.organelleStart) & (minimap2_df["organelleEnd"] <= row.organelleEnd)]

def get_reads_from_insertions(insertions_df, sequences_in_nucleus_df):
    dict_to_dataframe = {"readName": [], "insertionStart": [],
                         "insertionEnd": [], "organellerStart": [],
                         "organelleEnd": [], "strand": []}
    for row in insertions_df.itertuples():
        reads = sequences_in_nucleus_df.loc[~((row.nuclearEnd <= sequences_in_nucleus_df["nuclearStart"]) | (sequences_in_nucleus_df["nuclearEnd"] <= row.nuclearStart))]
        if not reads.empty:
            for read in reads.itertuples():
                dict_to_dataframe["readName"].append(read.readName)
                dict_to_dataframe["organelleStart"].append(row.organelleStart)
                dict_to_dataframe["organelleEnd"].append(row.organelleEnd)
                dict_to_dataframe["nuclearStart"].append(row.nuclearStart)
                dict_to_dataframe["nuclearEnd"].append(row.nuclearEnd)
                dict_to_dataframe["mappingStart"].append(read.nuclearStart)
                dict_to_dataframe["mappingEnd"].append(read.nuclearEnd)
                dict_to_dataframe["strand"].append(read.strand)
    return pd.DataFrame.from_dict(dict_to_dataframe).sort_values("nuclearStart")



    #mirar cada read in sequences in nucles
    #si la posicion solapa con insertions:
    #####Tabla
    #readName NuclearInsertionStart, NuclearInsertionEnd, OrganelleStart, OrganelleEnd
    pass

    

def main():
    arguments = get_arguments()
    genome_fpath = arguments["ref_genome"]
    root_dir = arguments["input_dir"]
    for identity in IDENTITIES:
        in_fname = "{}_{}".format(identity[0], identity[1])
        in_fpath = root_dir / in_fname
        summary = in_fpath / "summary.txt"
        sequences_fpath = str(in_fpath / "sd_0001.fastq")
        mapping_output = str(in_fpath / "reads_mapped_against_organelle.paf")
        if not arguments["skip_mapping"]:
            run_minimap2(genome_fpath, sequences_fpath, mapping_output)
        insertions_df = load_insertions_source_as_dataframe(summary)
        minimap2_df = load_minimap2_hits_as_dataframe(mapping_output)
        sequences_in_nucleus_df = load_read_positions_from_maf_into_dataframe(in_fpath / "sd_0001.maf")
        df = get_reads_from_insertions(insertions_df, sequences_in_nucleus_df)
        print(df)

    # overlaps = classify_minimap_hits(insertions_source, mapping_output)
    # with open(out_path / "results.txt", "w") as results_fhand:
    #     results_fhand.write("Insertion\tNumOverlaps\n")
    #     for insertion, reads in overlaps.items():
    #         results_fhand.write(insertion+"\t"+str(len(reads))+"\n")

    
if __name__ == "__main__":
    main()
