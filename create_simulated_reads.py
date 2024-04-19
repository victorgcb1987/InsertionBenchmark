#!/usr/bin/env python

import argparse
import os
import sys

from pathlib import Path

from src.pbsim import run_pbsim


IDENTITIES = [(0.41, 0.5), (0.51, 0.60), 
              (0.61, 0.7), (0.71, 0.8),
              (0.81, 0.9), (0.91, 1)]

#Generating program options
def parse_arguments():
    desc = "Create in silico Long Reads from a fasta files"
    parser = argparse.ArgumentParser(description=desc)
    
    help_fasta_source = "(Required) input dir"
    parser.add_argument("--input_dir", "-i", type=str,
                        help=help_fasta_source)
    
    help_strategy = "(Optional) strategy used. By default is wgs"
    parser.add_argument("--strategy", "-s", type=str,
                        help=help_strategy, default="wgs")
    
    help_sequencing_depth = "(Optional) Sequencing Depth. Default is 100x"
    parser.add_argument("--sequencing_depth", "-d", type=int,
                        help=help_sequencing_depth, default=100)
    
    help_min_length = "(Optional) Min read length. Deafult is 1kb"
    parser.add_argument("--min_length", "-x", type=int,
                        help=help_min_length, default=1000)

    help_max_length = "(Optional) Max read length. Default is 30kb"
    parser.add_argument("--max_length", "-X", type=int,
                        help=help_max_length, default=30000)
    
    help_method = "(Optional) method used for read creation. Default is errhmm"
    parser.add_argument("--method", "-m", type=str,
                        help=help_method, default="errhmm")
    
    help_method_model = "(Required) method model file for read creation"
    parser.add_argument("--model_file", "-f", type=str,
                        help=help_method_model)
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


#Parse and return values given to options when running this program
def get_arguments():
    parser = parse_arguments()
    input_dir = Path(parser.input_dir)
    strategy = parser.strategy
    sequencing_depth = parser.sequencing_depth
    min_length = parser.min_length
    max_length = parser.max_length
    method = parser.method
    model_file = Path(parser.model_file)   
    return {"input_dir": input_dir, "strategy": strategy,
            "sequencing_depth": sequencing_depth, "min_length": min_length, 
            "max_length": max_length, "method": method, "model_file": model_file}


def main():
    arguments = get_arguments()
    for identity_range in IDENTITIES:
        freq1 = identity_range[0]
        freq2 = identity_range[1]
        out_dirname = "{}_{}".format(str(freq1), str(freq2))
        root_dir = arguments["input_dir"] / out_dirname        
        os.chdir(str(root_dir))
        run_pbsim(strategy=arguments["strategy"], depth=arguments["sequencing_depth"], 
                   min_length=arguments["min_length"], max_length=arguments["max_length"],
                   method=arguments["method"], method_model=arguments["model_file"], reference=root_dir / "Nuclear_with_insertions.fasta")
        

if __name__ == "__main__":
    main()