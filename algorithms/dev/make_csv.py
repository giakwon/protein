import pickle
import numpy as np
import matplotlib.pyplot as plt
import csv


def main():
    # Load data from pickle 
    with open("ga_results.pkl", "rb") as f:
        results = pickle.load(f)

    seed = results["seed"]
    seq_length = results["seq_length"]
    num_sequences = results["num_sequences"]
    trials = results["trials"]
    generations = results["generations"]
    popsize = results["popsize"]
    target_fitness = results["target_fitness"]
    all_histories_grouped = results["all_histories"]

    with open("ga_history_from_pkl.csv", "w", newline="") as csvfile:
        fieldnames = ["sequence_number", "trial", "generation",
                    "best_fit", "median_fit", "worst_fit",
                    "count_max", "time_sec"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for seq_num, sequence_histories in enumerate(all_histories_grouped):
            for trial_num, trial_history in enumerate(sequence_histories):
                for gen_data in trial_history:
                    writer.writerow({
                        "sequence_number": seq_num,
                        "trial": trial_num,
                        "generation": gen_data["generation"],
                        "best_fit": gen_data["best_fit"],
                        "median_fit": gen_data["median_fit"],
                        "worst_fit": gen_data["worst_fit"],
                        "count_max": gen_data["count_max"],
                        "time_sec": gen_data["time_sec"]
                    })

main()