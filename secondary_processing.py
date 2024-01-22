import itertools
import json
import pandas as pd
import os
import json
import shutil

from json_to_csv import json_to_csv_processing


LOOKING_FOR = "ANCHOR_SCORE"
# LOOKING_FOR = "IUPRED_SCORE"


# Turns json_to_csv.json to a csv.

def merge_tuples(list_, MERGE_IF_CLOSER_THAN):
    merged_list = []

    if not list_:
        return merged_list

    current_tuple = list_[0]

    for next_tuple in list_[1:]:
        if next_tuple[0] - current_tuple[1] < MERGE_IF_CLOSER_THAN:
            current_tuple = (current_tuple[0], next_tuple[1])
        else:
            merged_list.append(current_tuple)
            current_tuple = next_tuple

    merged_list.append(current_tuple)

    return merged_list


def find_occurrences(main_string, substrings):
    """Finds how far through a string a substring is. E.g., "ll" is 3 through "Hello"."""
    occurrences = []

    for substring in substrings:
        start = 0
        positions = []
        while start < len(main_string):
            start = main_string.find(substring, start)
            if start == -1:
                break
            positions.append(start)
            start += 1
        if positions:
            occurrences.append((substring, positions))

    return occurrences


def count_letters(string, letters_to_count):
    """Counts the occurrence of letters in a string. Need to input which letters are being counted."""
    letter_counts = {}

    for char in letters_to_count:
        letter_counts[char] = 0

    for char in string:
        if char.isalpha():
            char = char.upper()

            if char in letter_counts:
                letter_counts[char] += 1

    sorted_letter_counts = dict(sorted(letter_counts.items()))

    return sorted_letter_counts

# IN_A_ROW_MIN = 30
# MERGE_IF_CLOSER_THAN = 10
# LOOK_FOR_ABOVE = 0.4

ALPHABET = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
ILV = ['I', 'L', 'V']
ILVM = ['I', 'L', 'V', 'M']
ILVMF = ['I', 'L', 'V', 'M', 'F']
DE = ['D', 'E']
VI = ['V', 'I']
PLVIMFY = ['P', 'L', 'V', 'I', 'M', 'F', 'Y']
PLVIM = ['P', 'L', 'V', 'I', 'M']
DES = ['D', 'E', 'S']
EDSTP = ['E', 'D', 'S', 'T', 'P']

type_1 = list(itertools.product(
    ILV, ILV, ALPHABET, ILV
))
type_2 = list(itertools.product(
    ILV, ALPHABET, ILV, ILV
))
type_3 = list(itertools.product(
    ILV, ILV, ILV, ILV
))
type_4 = list(itertools.product(
    ALPHABET, ILV, ILV, ILV
))
type_5 = list(itertools.product(
    ILV, DE, ILV, DE, ILV
))
type_alpha = list(itertools.product(
    VI, ALPHABET, VI, VI
))
type_beta = list(itertools.product(
    VI, ALPHABET, VI, ILV
))
type_a = list(itertools.product(
    PLVIM, ILVM, ALPHABET, DES, DES, DES
))
type_b = list(itertools.product(
    PLVIMFY, ILVM, "D", "L", "T"
))
type_r = list(itertools.product(
    DES, DES, DES, ILVM, ALPHABET, ILVMF, ILVMF
))
combined_list = type_1 + type_2 + type_3 + type_4 + type_5 + type_alpha + type_beta + type_a + type_b + type_r
complete_sim_list = []
for sim in combined_list:
    str_ = ''
    for item in sim:
        str_ = str_ + item
    complete_sim_list.append(str_)

def secondary_processing(iupred_number, LOOK_FOR_ABOVE, IN_A_ROW_MIN, MERGE_IF_CLOSER_THAN, IUPRED_OR_ANCHOR):

    if IUPRED_OR_ANCHOR == "IUPred":
        IUPRED_OR_ANCHOR = "IUPRED_SCORE"
    else:
        IUPRED_OR_ANCHOR = "ANCHOR_SCORE"

    try:
       shutil.rmtree("data/final_results")
    except FileNotFoundError:
        pass
    os.mkdir("data/final_results")

    directory = "data/raw_results_folder"

    final_json_to_csv = []

    count = 0

    for file in os.listdir(directory):
        #    print(file.title())

        with open(f"data/raw_results_folder/{file}") as file1:
            title = file1.readline()
        
        try:
            title = title.split("|")[1]
        except IndexError:
            print(title)

        df = pd.read_csv(
            f"data/raw_results_folder/{file}",
            sep="\t",
            header=None,
            names=["#", "AMINO_ACID", "IUPRED_SCORE", "ANCHOR_SCORE"],
            skiprows=1
        )



        '''
        The first line below gives a tag of "True" to any IUPRED score over x.
        The bottom code lines looks through rows with the tag "True" to search for consecutive ones.
        In the final one labelled "pr", the "30" defines how many in sequence it's looking for.

        Frankly, I don't really know how it works, but it does. I got it from stackoverflow.
        https://stackoverflow.com/questions/24281936/delimiting-contiguous-regions-with-values-above-a-certain-threshold-in-pandas-da
        '''

        df["tag"] = df[IUPRED_OR_ANCHOR] > LOOK_FOR_ABOVE
        fst = df.index[df["tag"] & ~ df["tag"].shift(1).fillna(False)]
        lst = df.index[df["tag"] & ~ df["tag"].shift(-1).fillna(False)]
        pr = [(i, j) for i, j in zip(fst, lst) if j > i + IN_A_ROW_MIN]

        number_before_merging = len(pr)
        pr = merge_tuples(pr, MERGE_IF_CLOSER_THAN)
        regions = []
        region_means = []
        region_mean_jsons = []
        number_regions = len(pr)
        number_merged = number_before_merging - number_regions


        '''GETTING THE AMINO ACIDS BETWEEN THE START OF THE SEQUENCE AND THE END OF IT.'''
        identifier_sims_found_dictionary_list = []
        json_to_csv_list = []
        for amino_acid_location in pr:
            aa_chain = []
            for n in range(amino_acid_location[0], amino_acid_location[1] + 1):
                row = df.iloc[n]
                amino_acid = row.iloc[1]
                aa_chain.append(amino_acid)
                number = row.iloc[0]
            amino_acid = ''.join(aa_chain)

            '''FINDING OUT IF THERE ARE ANY POTENTIAL SIMS IN CHAIN'''
            sim_dict_list = []
            sim_occurrences = find_occurrences(amino_acid, substrings=complete_sim_list)


            for occurrence in sim_occurrences:
                sim = occurrence[0]
                sim_location_in_aa = occurrence[1][0] + 1
                amino_acid_start = amino_acid_location[0]
                sim_location_in_identifier_start = sim_location_in_aa + amino_acid_start
                sim_location_in_idenfitier = (sim_location_in_identifier_start, sim_location_in_identifier_start + len(sim))

                sim_tuple = tuple(list(sim))

                sim_dict = {
                    "Amino Acid Regions where SIMs present":
                        f"{sim_location_in_identifier_start}-{sim_location_in_identifier_start + len(sim) - 1}",
                    "Sequences of the SIM": sim,
                }
                if sim_tuple in type_r:
                    sim_dict["Type of SIM"] = "Type r"
                elif sim_tuple in type_a:
                    sim_dict["Type of SIM"] = "Type a"
                elif sim_tuple in type_b:
                    sim_dict["Type of SIM"] = "Type b"
                elif sim_tuple in type_alpha:
                    sim_dict["Type of SIM"] = "Type alpha"
                elif sim_tuple in type_5:
                    sim_dict["Type of SIM"] = "Type 5"
                elif sim_tuple in type_beta:
                    sim_dict["Type of SIM"] = "Type beta"
                elif sim_tuple in type_1:
                    sim_dict["Type of SIM"] = "Type 1"
                elif sim_tuple in type_2:
                    sim_dict["Type of SIM"] = "Type 2"
                elif sim_tuple in type_3:
                    sim_dict["Type of SIM"] = "Type 3"
                elif sim_tuple in type_4:
                    sim_dict["Type of SIM"] = "Type 4"


                else:
                    sim_dict["type"] = "error"
                identifier_sims_found_dictionary_list.append(sim_dict)

                '''FINDING COUNTS E/D/S/T/P WITHIN 10 AMINO ACIDS OF A SIM'''

                search_location_start = sim_location_in_identifier_start - 10
                if search_location_start < 0:
                    search_location_start = 0
                search_location_end = sim_location_in_idenfitier[1] + 10
                if search_location_end > len(df.index):
                    search_location_end = len(df.index)
                edstp_search = []
                for n in range(search_location_start, search_location_end):
                    row = df.iloc[n]
                    aa = row.iloc[1]
                    edstp_search.append(aa)
                edstp_search = ''.join(edstp_search)
                edstp_tally = count_letters(edstp_search, EDSTP)

                sim_dict["Amino acid region of SIM"] = edstp_search

                d_e = edstp_tally["D"] + edstp_tally["E"]
                s_t = edstp_tally["S"] + edstp_tally["T"]
                p_ = edstp_tally["P"]

                sim_dict["D-E"] = d_e
                sim_dict["S-T"] = s_t
                sim_dict["P"] = p_

                if sim_dict not in sim_dict_list:
                    sim_dict_list.append(sim_dict)
                else:
                    pass


            amino_acid_letter_count = count_letters(amino_acid, ALPHABET)

            amino_acid_json = {
                "amino_acid": amino_acid,
                "amino_acid_letter_counts": amino_acid_letter_count
            }

            results = df.iloc[amino_acid_location[0]:amino_acid_location[1], :]

            region = f"AA {amino_acid_location[0] + 1}-{amino_acid_location[1]}"
            regions.append(region)

            iupred_mean = round(results[IUPRED_OR_ANCHOR].mean(), 4)
            region_mean = f"{region} = {iupred_mean}"
            region_mean_json = {
                "region": region,
                "region_start": amino_acid_location[0] + 1,
                "region_end": amino_acid_location[1],
                "mean": iupred_mean,
                "amino_acid_chain": amino_acid_json,
            }

            if len(sim_dict_list) > 0:
                json_to_csv = {
                    "Disordered region": f"{amino_acid_location[0] + 1}-{amino_acid_location[1]}",
                    "Mean disorder score": iupred_mean,
                    "Number of Sims": len(sim_dict_list),
                    "SIMs": sim_dict_list
                }
                json_to_csv_list.append(json_to_csv)

            region_means.append(region_mean)
            region_mean_jsons.append(region_mean_json)

        count += 1
        
        try:
            print(f"secondary processing... {round((count/iupred_number)*100,1)}%", end="\r")
        except TypeError:
            print(f"secondary processing... {count}", end="\r")

        if number_regions > 0:
            tuple_list = [(tuple_[0] + 1, tuple_[1]) for tuple_ in pr]

            if not identifier_sims_found_dictionary_list:
                identifier_sims_found_dictionary_list = 0

            final_json = {
                "identifier": title,
                "number_regions": number_regions,
                "sequences_merged": number_merged,
                "regions": region_mean_jsons,
                "sims_found": identifier_sims_found_dictionary_list
            }

            json_to_csv = {
                "Identifier": title,
                "Number of disordered regions": number_regions,
                "Disordered region": {
                    "disordered_region": json_to_csv_list
                }
            }
            # result_to_csv_list.append(json_to_csv)
            # with open("data/final_results/results_to_csv.json", "a") as f:
            #     json.dump(json_to_csv, f)
            #     f.write(",")
            final_json_to_csv.append(json_to_csv)


            final_entry = (
                f"Identifier: {title}\n"
                f"Number of regions: {number_regions}\n"
                f"Sequences merged: {number_merged}\n"
                f"Regions: {regions}\n"
                f"Regions tuple: {tuple_list}\n"
                f"IUPRED mean: {region_means}\n"
                f"SIMs found: {identifier_sims_found_dictionary_list}\n\n"
            )


        else:
            final_entry = f"Identifier: {title}\nNumber of regions: 0\n\n"
            final_json = {
                "identifier": title,
                "number_regions": 0
            }


    with open("data/final_results/results_to_csv.json", "w") as f:
        json.dump(final_json_to_csv, f, indent = 4)

    json_to_csv_processing()
    print("Secondary processing complete")
 

def tertiary_processing(nuclear_score):
    # Merging the iupred results and the nuclear csv files
    with open("data/nuclear_data.csv") as f:
        filedata = f.read()
    filedata = filedata.replace("Nucleus", "Nuclear Score")
    with open("data/nuclear_data.csv", "w") as f:
        f.write(filedata)

    main_df = pd.read_csv("data/output.csv")
    main_df = main_df.drop(columns=["Drop1", "Drop2"])
    additional_df = pd.read_csv("data/nuclear_data.csv")
    df = pd.merge(main_df, additional_df[["Identifier", "Nuclear Score"]], on="Identifier", how="left")

    # drop row if nuclear score less than input from Data Processing page
    df = df.drop(df[df["Nuclear Score"] < nuclear_score].index)

    # Creating the long dataset
    df.to_csv("data/final_results/long_dataset.csv")

    # Creating the short dataset. Finding the count of each identifier and then dropping all duplicate rows. 
    df = df[["Identifier"]]
    df['No. Putative SIMs'] = df.groupby('Identifier')['Identifier'].transform('count')
    df = df.drop_duplicates(subset='Identifier')
    df.to_csv("data/final_results/short_dataset.csv", index=False)

    print("Tertiary processing complete")
