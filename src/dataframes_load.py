import pandas as pd


def load_minimap2_hits_as_dataframe(minimap_hits):
    dict_to_dataframe = {"readName": [], "readLength": [],
                         "readStart": [], "readEnd": [],
                         "strand": [], 
                         "organelleStart": [], "organelleEnd": [],
                         "numMatches": [], "alnLength": []}
    with open(minimap_hits) as minimap_fhand:
        for line in minimap_fhand:
            line = line.split()
            dict_to_dataframe["readName"].append(line[0])
            dict_to_dataframe["readLength"].append(int(line[1]))
            dict_to_dataframe["readStart"].append(int(line[2]))
            dict_to_dataframe["readEnd"].append(int(line[3]))
            dict_to_dataframe["strand"].append(line[4])
            dict_to_dataframe["organelleStart"].append(int(line[7]))
            dict_to_dataframe["organelleEnd"].append(int(line[8]))
            dict_to_dataframe["numMatches"].append(int(line[9]))
            dict_to_dataframe["alnLength"].append(int(line[10]))
    return pd.DataFrame.from_dict(dict_to_dataframe)


def load_read_positions_from_maf_into_dataframe(maf_file):
    dict_to_dataframe = {"readName": [], "nuclearStart": [],
                         "nuclearEnd": []}
    with open(maf_file) as maf_fhand:
        for line in maf_fhand:
            if "ref" in line:
                nuclear_positions = (int(line.split()[2]), int(line.split()[2])+int(line.split()[3]))
            elif line.startswith("s"):
                readname = line.split()[1]
                dict_to_dataframe["readName"].append(readname)
                dict_to_dataframe["nuclearStart"].append(nuclear_positions[0])
                dict_to_dataframe["nuclearEnd"].append(nuclear_positions[1])
    return pd.DataFrame.from_dict(dict_to_dataframe)

def load_insertions_source_as_dataframe(summary):
    dict_to_dataframe = {"nuclearStart": [], "nuclearEnd": [], 
                         "organelleStart":[], "organelleEnd": []}
    with open(summary) as summary_fhand:
        for line in summary_fhand:
            line = line.rstrip().split()
            length = int(line[-1])
            dict_to_dataframe["nuclearStart"].append(int(line[0]))
            dict_to_dataframe["nuclearEnd"].append(int(line[0])+length)
            dict_to_dataframe["organelleStart"].append(int(line[1].split("_")[-1].split(":")[0]))
            dict_to_dataframe["organelleEnd"].append(int(line[1].split("_")[-1].split(":")[1]))
    return pd.DataFrame.from_dict(dict_to_dataframe).sort_values("nuclearStart")