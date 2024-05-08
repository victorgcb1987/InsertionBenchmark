import gzip
import pandas as pd



def load_minimap2_hits_as_dataframe(minimap_hits, repetitive_regions=()):
    dict_to_dataframe = {"readName": [], "readLength": [],
                         "readStart": [], "readEnd": [],
                         "strand": [], 
                         "organelleStart": [], "organelleEnd": [],
                         "numMatches": [], "alnLength": [],
                         "IR": []}
    strand_conversion = {"+": "-", "-": "+"}
    with open(minimap_hits) as minimap_fhand:
        for line in minimap_fhand:
            repetitive = False
            line = line.split()
            read_name = line[0]
            read_length = int(line[1])
            read_start = int(line[2])
            read_end = int(line[3])
            strand = line[4]
            organelle_start = int(line[7])
            organelle_end = int(line[8])
            num_matches = int(line[9])
            aln_length = int(line[10])
            
            if organelle_start >= repetitive_regions[0][0] and organelle_end <= repetitive_regions[0][1]:
                rev_organelle_start = repetitive_regions[1][0] + (repetitive_regions[0][1] - organelle_end)
                rev_organelle_end = repetitive_regions[1][1] - (organelle_start - repetitive_regions[0][0])
                repetitive = True
            elif organelle_start >= repetitive_regions[1][0] and organelle_end <= repetitive_regions[1][1]:
                rev_organelle_start = repetitive_regions[0][0] + (repetitive_regions[1][1] - organelle_end)
                rev_organelle_end = repetitive_regions[0][1] - (organelle_start -  repetitive_regions[1][0])
                repetitive = True
            if repetitive:
                dict_to_dataframe["readName"].append(read_name)
                dict_to_dataframe["readLength"].append(read_length)
                dict_to_dataframe["readStart"].append(read_start)
                dict_to_dataframe["readEnd"].append(read_end)
                dict_to_dataframe["strand"].append(strand_conversion[strand])
                dict_to_dataframe["organelleStart"].append(rev_organelle_start)
                dict_to_dataframe["organelleEnd"].append(rev_organelle_end)
                dict_to_dataframe["numMatches"].append(num_matches)
                dict_to_dataframe["alnLength"].append(aln_length)
                dict_to_dataframe["IR"].append(repetitive)
            
            dict_to_dataframe["readName"].append(read_name)
            dict_to_dataframe["readLength"].append(read_length)
            dict_to_dataframe["readStart"].append(read_start)
            dict_to_dataframe["readEnd"].append(read_end)
            dict_to_dataframe["strand"].append(strand)
            dict_to_dataframe["organelleStart"].append(organelle_start)
            dict_to_dataframe["organelleEnd"].append(organelle_end)
            dict_to_dataframe["numMatches"].append(num_matches)
            dict_to_dataframe["alnLength"].append(aln_length)
            dict_to_dataframe["IR"].append(repetitive)
    return pd.DataFrame.from_dict(dict_to_dataframe)


def load_read_positions_from_maf_into_dataframe(maf_file):
    dict_to_dataframe = {"readName": [], "nuclearStart": [],
                         "nuclearEnd": [], "strand": []}
    with gzip.open(maf_file, "rt") as maf_fhand:
        for line in maf_fhand:
            if "ref" in line:
                nuclear_positions = (int(line.split()[2]), int(line.split()[2])+int(line.split()[3]))
            elif line.startswith("s"):
                readname = line.split()[1]
                strand = line.split()[4]
                dict_to_dataframe["readName"].append(readname)
                dict_to_dataframe["nuclearStart"].append(nuclear_positions[0])
                dict_to_dataframe["nuclearEnd"].append(nuclear_positions[1])
                dict_to_dataframe["strand"].append(strand)
    return pd.DataFrame.from_dict(dict_to_dataframe)


def load_insertions_source_as_dataframe(summary, repetitive_regions=()):
    dict_to_dataframe = {"nuclearStart": [], "nuclearEnd": [], 
                         "organelleStart":[], "organelleEnd": [],
                         "IR": []}
    with open(summary) as summary_fhand:
        for line in summary_fhand:
            repetitive = False
            line = line.rstrip().split()
            length = int(line[-1])
            nuclear_start = int(line[0])
            nuclear_end = nuclear_start+length
            organelle_start = int(line[1].split("_")[-1].split(":")[0])
            organelle_end = int(line[1].split("_")[-1].split(":")[1])

            if organelle_start >= repetitive_regions[0][0] and organelle_end <= repetitive_regions[0][1]:
                rev_organelle_start = repetitive_regions[1][0] + (repetitive_regions[0][1] - organelle_end)
                rev_organelle_end = repetitive_regions[1][1] - (organelle_start - repetitive_regions[0][0])
                repetitive = True
            elif organelle_start >= repetitive_regions[1][0] and organelle_end <= repetitive_regions[1][1]:
                rev_organelle_start = repetitive_regions[0][0] + (repetitive_regions[1][1] - organelle_end)
                rev_organelle_end = repetitive_regions[0][1] - (organelle_start -  repetitive_regions[1][0])
                repetitive = True
            if repetitive:
                dict_to_dataframe["nuclearStart"].append(nuclear_start)
                dict_to_dataframe["nuclearEnd"].append(nuclear_end)
                dict_to_dataframe["organelleStart"].append(rev_organelle_start)
                dict_to_dataframe["organelleEnd"].append(rev_organelle_end)
                dict_to_dataframe["IR"].append(repetitive)

            dict_to_dataframe["nuclearStart"].append(nuclear_start)
            dict_to_dataframe["nuclearEnd"].append(nuclear_end)
            dict_to_dataframe["organelleStart"].append(organelle_start)
            dict_to_dataframe["organelleEnd"].append(organelle_end)
            dict_to_dataframe["IR"].append(repetitive)

    return pd.DataFrame.from_dict(dict_to_dataframe).sort_values("nuclearStart")

def merge_minimap2_and_reference_nuclear(ref_df, minimap2_df):
    merge = ref_df.merge(minimap2_df, how="outer", on="readName")
    filtered_merge = merge.loc[~((merge["organelleEnd_y"] <= merge["organelleStart_x"]) | (merge["organelleEnd_x"] <= merge["organelleStart_y"]))]
    filtered_merge.to_csv("check2.tsv", sep="\t", index=False, na_rep='NULL')
    return filtered_merge


def get_reads_from_insertions(insertions_df, sequences_in_nucleus_df):
    dict_to_dataframe = {"readName": [],  "organelleStart": [],
                         "organelleEnd": [], "strand": [],
                         "refMappingStart": [], "refMappingEnd": [],
                         "pointOfInsertionStart": [], "pointOfInsertionEnd": []}
    for row in insertions_df.itertuples():
        reads = sequences_in_nucleus_df.loc[~((row.nuclearEnd <= sequences_in_nucleus_df["nuclearStart"]) | (sequences_in_nucleus_df["nuclearEnd"] <= row.nuclearStart))]
        if not reads.empty:
            for read in reads.itertuples():
                dict_to_dataframe["readName"].append(read.readName)
                dict_to_dataframe["organelleStart"].append(row.organelleStart)
                dict_to_dataframe["organelleEnd"].append(row.organelleEnd)
                dict_to_dataframe["pointOfInsertionStart"].append(row.nuclearStart)
                dict_to_dataframe["pointOfInsertionEnd"].append(row.nuclearEnd)
                dict_to_dataframe["refMappingStart"].append(read.nuclearStart)
                dict_to_dataframe["refMappingEnd"].append(read.nuclearEnd)
                dict_to_dataframe["strand"].append(read.strand)
    return pd.DataFrame.from_dict(dict_to_dataframe).sort_values("pointOfInsertionStart")

def filter_by_readname(dataframe_to_filter,  readnames, mode="exclude"):
    if mode == "exclude":
        filtered_dataframe = dataframe_to_filter.loc[~dataframe_to_filter["readName"].isin(readnames)]
    elif mode == "include":
        filtered_dataframe = dataframe_to_filter.loc[dataframe_to_filter["readName"].isin(readnames)]
    return filtered_dataframe

