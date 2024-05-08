#!/usr/bin/env python

import argparse
import pandas as pd
import sys

from pathlib import Path


from src.dataframes import (load_insertions_source_as_dataframe, 
                                 load_minimap2_hits_as_dataframe, 
                                 load_read_positions_from_maf_into_dataframe,
                                 get_reads_from_insertions,
                                 merge_minimap2_and_reference_nuclear,
                                 filter_by_readname)
from src.minimap2 import run_minimap2

from src.stats import get_mapping_stats


IDENTITIES = [(0.41, 0.5), (0.51, 0.60), 
              (0.61, 0.7), (0.71, 0.8),
              (0.81, 0.9), (0.91, 1)]

IDENTITIES = [(0.91, 1)]


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


def main():
    test_ir = [(84171, 110434), (128215, 154478)]
    arguments = get_arguments()
    genome_fpath = arguments["ref_genome"]
    root_dir = arguments["input_dir"]
    results_dict = {"identity": [], "total_reads": [],
                    "reads_from_insertions": [],
                    "reads_mapped": [], "reads_mapped_from_insertions": [],
                    "reads_mapped_not_from_insertion": []}
    for identity in IDENTITIES:
        in_fname = "{}_{}".format(identity[0], identity[1])
        in_fpath = root_dir / in_fname
        summary = in_fpath / "summary.txt"
        sequences_fpath = str(in_fpath / "sd_0001.fastq.gz")
        mapping_output = str(in_fpath / "reads_mapped_against_organelle.paf")
        if not arguments["skip_mapping"]:
           run_minimap2(genome_fpath, sequences_fpath, mapping_output)
        insertions_df = load_insertions_source_as_dataframe(summary, repetitive_regions=test_ir)
        insertions_df.to_csv("check.tsv", sep="\t")
        minimap2_df = load_minimap2_hits_as_dataframe(mapping_output, repetitive_regions=test_ir)
        minimap2_df.to_csv("minimap2.tsv", sep="\t", index=False)
        sequences_in_nucleus_df = load_read_positions_from_maf_into_dataframe(in_fpath / "sd_0001.maf.gz")
        total_sequences = set(sequences_in_nucleus_df['readName'].tolist())
        reads_from_insertions_df = get_reads_from_insertions(insertions_df, sequences_in_nucleus_df)
        from_insertions_readnames = set(reads_from_insertions_df['readName'].tolist())
        #reads_not_from_insertions_df = filter_by_readname(reads_from_insertions_df, from_insertions_readnames, mode="exclude")
        not_from_insertions_readnames = total_sequences - from_insertions_readnames
        mapped_reads = set(minimap2_df['readName'].tolist())
        mapped_reads_from_insertions = mapped_reads - not_from_insertions_readnames
        mapped_reads_not_from_insertions = mapped_reads - from_insertions_readnames
        results_dict["identity"].append(identity)
        results_dict["total_reads"].append(len(total_sequences))
        results_dict["reads_from_insertions"].append(len(from_insertions_readnames))
        results_dict["reads_mapped"].append(len(mapped_reads))
        results_dict["reads_mapped_from_insertions"].append(len(mapped_reads_from_insertions))
        results_dict["reads_mapped_not_from_insertion"].append(len(mapped_reads_not_from_insertions))
    results_df = pd.DataFrame.from_dict(results_dict)
    results_df.to_csv(root_dir/"results.tsv", sep="\t", index=False)
if __name__ == "__main__":
    main()
