
def get_mapping_stats(merged_df):
    reads_from_simulation_and_insertion = merged_df.loc[~(merged_df["organelleStart_x"].isnull())]
    reads_mapped = merged_df.loc[~(merged_df["organelleStart_x"].isnull()) & ~(merged_df["organelleStart_y"].isnull())]
    reads_from_insertion_not_mapped =  merged_df.loc[~(merged_df["organelleStart_x"].isnull()) & merged_df["organelleStart_y"].isnull()]
    reads_mapped_not_from_insertion = merged_df.loc[(merged_df["organelleStart_x"].isnull() & ~(merged_df["organelleStart_y"].isnull()))]
    reads_from_simulation_and_insertion.to_csv("reads_from_simulation_and_insertion.tsv", sep="\t", index=False, na_rep='NULL')
    reads_mapped.to_csv("reads_mapped.tsv", sep="\t", index=False, na_rep="NULL")
    reads_from_insertion_not_mapped.to_csv("reads_from_insertion_not_mapped.tsv", sep="\t", index=False, na_rep='NULL')
    reads_mapped_not_from_insertion.to_csv("reads_mapped_not_from_insertion.tsv", sep="\t", index=False, na_rep='NULL')
    return ({"reads_from_insertions": reads_from_simulation_and_insertion,
            "reads_mapped": reads_mapped,
            "reads_from_insertions_not_mapped": reads_from_insertion_not_mapped,
            "reads_not_from_insertions_mapped_in_organelle": reads_mapped_not_from_insertion})

def get_reads_in_insertions(minimap2_df, insertions_df):
    dict_to_dataframe = {"insertion": []}
    for row in insertions_df.itertuples():
        reads = minimap2_df.loc[~((row.organelleEnd <= minimap2_df["organelleStart"]) | (minimap2_df["organelleEnd"] <= row.organelleStart))]
        #reads = minimap2_df.loc[(minimap2_df["organelleStart"] >= row.organelleStart) & (minimap2_df["organelleEnd"] <= row.organelleEnd)]

