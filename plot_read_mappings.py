#!/usr/bin/env python

from argparse import ArgumentParser
from csv import DictReader
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def parse_arguments():
    desc = "Plot read mapping results from insertion_analysis output RESULTS"
    parser = ArgumentParser(description=desc)

    help_input_fpath = "(Required) Path to a TAB-file with insertion length and files paths"
    parser.add_argument("--input", 
                        "-i", type=str,
                        help=help_input_fpath,
                        required=True)
    help_output_fpath = "(Required) Output Dir"
    parser.add_argument("--output", 
                        "-o", type=str,
                        help=help_output_fpath,
                        required=True) 
    
    return parser


def get_arguments():
    options = parse_arguments().parse_args()
    results = {}
    with open(options.input) as input_fhand:
        for line in input_fhand:
            if line:
                line = line.rstrip().split()
                results[int(line[0])] = Path(line[1])
        out_path = Path(options.output)
        if not out_path.exists():
            out_path.mkdir(parents=True)
    return {"results": results,
            "output": out_path}



def get_values_from_results(length, results_fhand, results_dict):
    for line in DictReader(results_fhand, delimiter="\t"):
        results_dict["identity"].append(line["identity"])
        results_dict["reads_mapped"].append(float(line["reads_mapped_from_insertions"]) / float(line["reads_from_insertions"]))
        results_dict["length"].append(length)
    return results_dict


def plot_reads_mapped(df):
    #sns.color_palette("bright")
    sns.set_theme(style="ticks")
    sns.set_theme(style="darkgrid")
    sns.lineplot(x="identity", y="reads_mapped",
             hue="length", legend="full", palette=["#039BBA", "#D60270", "#FFC72C", 
                                                   "#A1E6EF", "#99B0EA", "#000000",
                                                   "#002497", "#E57200", "#009A44"],
             data=df)
    plt.show()

def main():
    sns.set_theme(style="ticks")
    results_dict = {"identity": [], "reads_mapped": [],
                    "length": []}
    arguments = get_arguments()
    for length, results_fpath in arguments["results"].items():
        with open(results_fpath) as results_fhand:
            results_dict = get_values_from_results(length, results_fhand, results_dict)   
    df = pd.DataFrame.from_dict(results_dict)
    print(df)
    plot_reads_mapped(df)
    

if __name__ == "__main__":
    main()