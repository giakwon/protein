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

    # Plot best fits across generations
    best_fits = np.array([[gen['best_fit'] for gen in trial]
                        for sequence in all_histories_grouped
                        for trial in sequence])

    highest = np.max(best_fits, axis=0)
    median  = np.median(best_fits, axis=0)
    lowest  = np.min(best_fits, axis=0)

    fig, ax = plt.subplots()
    x = np.arange(best_fits.shape[1])
    ax.plot(x, highest, label="Highest Fitness", color="red")
    ax.plot(x, median, label="Median Fitness", color="blue")
    ax.plot(x, lowest, label="Lowest Fitness", color="green")
    ax.set_ylim(ymin=0)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Fitness")
    ax.set_title(f"GA Performance (Sequence length={seq_length}, "
                f"{num_sequences} sequences, {trials} trials each)")
    ax.legend()
    plt.show()

    # Plot count of max fitness
    count_maxes = np.array([[gen['count_max'] for gen in trial]
                            for sequence in all_histories_grouped
                            for trial in sequence])
    avg_counts = count_maxes.mean(axis=0)

    plt.figure(figsize=(8,5))
    plt.plot(range(generations), avg_counts, marker="o")
    plt.xlabel("Generation")
    plt.ylabel(f"Count of puzzles with fitness = {target_fitness}")
    plt.title(f"Count of Max-Fitness (={target_fitness}) Puzzles per Generation")
    plt.grid(True)
    plt.show()

main()