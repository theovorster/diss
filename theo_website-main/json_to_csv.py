import csv
import json


def flatten_json(json_data, identifier):
    result = []
    total_disordered_regions = json_data.get("Number of disordered regions", 0)

    for region in json_data.get("Disordered region", {}).get("disordered_region", []):
        disordered_region = region.get("Disordered region", "")
        mean_disorder_score = region.get("Mean disorder score", "")
        num_sims = region.get("Number of Sims", 0)
        sims = region.get("SIMs", [])

        if num_sims == 0:
            result.append([
                identifier, total_disordered_regions, disordered_region, mean_disorder_score, num_sims, "", "", "", "",
                "", ""
            ])
        else:
            for sim in sims:
                amino_acid_regions = sim.get("Amino Acid Regions where SIMs present", "")
                sequences_of_sim = sim.get("Sequences of the SIM", "")
                type_of_sim = sim.get("Type of SIM", "")
                amino_acid_region_of_sim = sim.get("Amino acid region of SIM", "")
                d_e = sim.get("D-E", "")
                s_t = sim.get("S-T", "")
                p = sim.get("P", "")
                result.append([
                    identifier, total_disordered_regions, disordered_region, mean_disorder_score, num_sims,
                    amino_acid_regions, sequences_of_sim, type_of_sim, amino_acid_region_of_sim, d_e, s_t, p
                ])

    return result

def json_to_csv_processing():
    with open("data/final_results/results_to_csv.json", "r") as json_file:
        json_data = json.load(json_file)


    # Open CSV file for writing with comma delimiter
    with open("data/output.csv", "w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=",")

        csv_writer.writerow([
            "Identifier",
            "Drop1",
            "Disordered region",
            "Disorder score",
            "Drop2",
            "SIM Position Site",
            "SIM Sequence",
            "SIM Type",
            "SIM region sequence",
            "D/E",
            "S/T",
            "P"
        ])

        for data in json_data:
            identifier = data.get("Identifier", "")
            flattened_data = flatten_json(data, identifier)
            csv_writer.writerows(flattened_data)
