import subprocess
import numpy.random as npr
import random
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import multiprocessing as mp
import pickle
import time
import csv


fitness_cache = {}
def fitness(seq):
    """
    Parses the output from C file ('../lattice-enum/lattice-enum') and turns it 
    into a dictionary. 

    param:
        seq (str): the HP sequence of the protein

    return (int): total number of distinct number of HH contacts for the sequence
    """
    if seq in fitness_cache:
        return fitness_cache[seq]

    process = subprocess.run(['../lattice-enum/lattice-enum'] + [seq], 
                             capture_output=True, text=True, check=True)
    c_output = process.stdout.strip().split('\n')

    conformations = {}
    for line in c_output:
        kv_pair = line.strip().split(" ")
        if len(kv_pair) == 2:
            key = kv_pair[0]
            value = kv_pair[1]
            if conformations.get(key) is None:
                conformations[key] = [value]
            conformations[key].append(value)

    fit = len(conformations)
    fitness_cache[seq] = fit
    return fit


def evaluate_population(population):
    """
    Uses parllel processing to evaluate the fitness 
    of each sequence in the population.

    params:
        population (lst of str): list of HP sequences to evaluate

    return (lst of tuples): each tuple contains (fitness, sequence)
    """
    with mp.Pool(processes=mp.cpu_count()) as pool:
        fitnesses = pool.map(fitness, population)
    return list(zip(fitnesses, population))


def random_sequence(length):
    """
    Makes a random sequence with 'H' and 'P'

    params:
        length (int): lenth of sequence 

    return (str): random seqeunce 
    """
    sequence = ""
    for _ in range(length):
        sequence += random.choice(['H', 'P']) 
    return sequence 


def mutate(seq, mutation_rate=0.1):
    """
    Randomly mutuates the sequence.

    params:
        seq (str): the HP sequence of the protein
        mutation_rate (float): the probabilty of mutations

    returns (str): mutated sequence
    """
    new_seq = ""
    for res in seq:
        if random.random() < mutation_rate:
            new_seq += random.choice(['H', 'P'])
        else:
            new_seq += res
    return new_seq


def crossover(p1, p2):
    """
    Splits and combines two sequences at a random crossover point.

    params:
        p1 (str): first parent from population
        p2 (str): second parent from population

    returns:
        c1 (str): first child
        c2 (str): second child
    """
    point = random.randint(1, len(p1)-2)
    c1 = p1[:point] + p2[point:]
    c2 = p2[:point] + p1[point:]
    return c1, c2


def select(population):
    """
    Selects an element from the population from proportion based on fitness

    params:
        population (list of tuples): current population where each element is (fitness, seq).

    return: selected sequence
    """
    m = sum([c[0] for c in population])
    if(m == 0):
      selection_probs = [1/len(population) for c in population]
    else:
      selection_probs = [c[0]/m for c in population]
    return population[npr.choice(len(population), p=selection_probs)][1]


def evolve(base_seq, generations, popsize, target_fitness, mutation_rate=0.1, elitism=True):
    """
    Runs a genetic algorithm to evolve HP sequences toward higher fitness.

    Args:
        base_seq (str): initial HP sequence to define sequence length.
        generations (int, optional): number of generations to run the algorithm. 
        popsize (int, optional): Number of sequences in each generation. 
        mutation_rate (float, optional): Probability of mutating each residue in a sequence (0.1)
        elitism (bool, optional): Keeps the best individual from each generation into the next.

    Returns:
        tuple: 
            best_seq_overall (str): The highest-fitness sequence found across all generations.
            best_fitness_overall (int): Fitness value of the best sequence.
            history (list of dicts): Additional info from each generation
    """
    seq_len = len(base_seq)
    population = [random_sequence(seq_len) for _ in range(popsize)]
    history = []
    best_seq_overall = None
    best_fitness_overall = -1

    for gen in tqdm(range(generations), desc="Evolving", leave=False):
        start_time = time.time()

        scored = evaluate_population(population)
        scored.sort(reverse=True)

        best_fitness_gen, best_seq_gen = scored[0]  # best of current generation
        median_fitness_gen = np.median([f for f, _ in scored])
        worst_fitness_gen  = scored[-1][0]

        # update overall best
        if best_fitness_gen > best_fitness_overall:
            best_fitness_overall = best_fitness_gen
            best_seq_overall = best_seq_gen

        count_max = sum(1 for f, _ in scored if f == target_fitness)

        end_time = time.time() - start_time

        # history.append(best_fitness_gen)
        history.append({
            "generation": gen,
            "population": [seq for _, seq in scored],
            "fitnesses": [f for f, _ in scored],
            "best_seq": best_seq_gen,
            "best_fit": best_fitness_gen,
            "median_fit": median_fitness_gen,
            "worst_fit": worst_fitness_gen,
            "count_max": count_max,
            "time_sec": end_time, 
        })

        new_population = []
        if elitism:
            new_population.append(best_seq_gen)

        while len(new_population) < popsize:
            p1 = select(scored)
            p2 = select(scored)
            c1, c2 = crossover(p1, p2)
            c1 = mutate(c1, mutation_rate)
            c2 = mutate(c2, mutation_rate)
            new_population.extend([c1, c2])

        population = new_population[:popsize]

    return best_seq_overall, best_fitness_overall, history



if __name__ == "__main__": 

    #making a seed
    seed = 42
    random.seed(seed)
    np.random.seed(seed)
    print("Seed used: %s" % seed)

    # 1: OVERALL GA PERFORMANCE (low, med, high)
    seq_length = 8
    num_sequences = 5 
    trials = 5
    generations = 50
    popsize = 30
    target_fitness = 9

    # all_histories = []
    # for _ in tqdm(range(num_sequences), desc="sequences", leave=False):
    #     base_seq = random_sequence(seq_length)

    #     for _ in range(trials):
    #         _, _, hist = evolve(base_seq, generations=generations, popsize=popsize, target_fitness=target_fitness)
    #         all_histories.append(hist)

    all_histories_grouped = []
    for seq_idx in tqdm(range(num_sequences), desc="Sequences"):
        base_seq = random_sequence(seq_length)
        sequence_trials = []
        for trial_idx in tqdm(range(trials), desc="Trials"):
            _, _, hist = evolve(base_seq, generations=generations, popsize=popsize, target_fitness=target_fitness)
            sequence_trials.append(hist)
        all_histories_grouped.append(sequence_trials)


    results = {
        "seed": seed,
        "seq_length": seq_length,
        "num_sequences": num_sequences,
        "trials": trials,
        "generations": generations,
        "popsize": popsize,
        "target_fitness": target_fitness,
        "all_histories": all_histories_grouped
    }
    with open("ga_results.pkl", "wb") as f:
        pickle.dump(results, f)
    print("Results saved to ga_results.pkl")

    with open("ga_history.csv", "w", newline="") as csvfile:
        fieldnames = ["sequence_number", "trial", "generation", 
                    "best_fit", "median_fit", "worst_fit", "count_max", "time_sec"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        # for seq_idx, _ in enumerate(range(num_sequences)):
        #     for trial_idx in range(trials):
        #         hist = all_histories[seq_idx * trials + trial_idx]  
        #         for generation_data in hist:
        #             writer.writerow({
        #                 "sequence_number": seq_idx,
        #                 "trial": trial_idx,
        #                 "generation": generation_data["generation"],
        #                 "best_fit": generation_data["best_fit"],
        #                 "median_fit": generation_data["median_fit"],
        #                 "worst_fit": generation_data["worst_fit"],
        #                 "time_sec": generation_data["time_sec"]
        #             })

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
    print("History exported to ga_history.csv")

    # histories = np.array(all_histories)
    # best_fits = np.array([[gen['best_fit'] for gen in hist] for hist in all_histories])
    # best_fits = np.array([[gen['best_fit'] for gen in trial] 
    #                       for sequence in all_histories_grouped 
    #                       for trial in sequence])
    
    best_fits = np.array([
        [gen['best_fit'] for gen in trial[:generations]] + 
        [np.nan] * (generations - len(trial))  # pad if too short
        for sequence in all_histories_grouped
        for trial in sequence
    ])

    # Compute per-generation statistics
    # highest = np.max(best_fits, axis=0)
    # median  = np.median(best_fits, axis=0)
    # lowest  = np.min(best_fits, axis=0)
    highest = np.nanmax(best_fits, axis=0)
    median  = np.nanmedian(best_fits, axis=0)
    lowest  = np.nanmin(best_fits, axis=0)

    # Plotting GA Performance
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


    # Count of Max Fitness in each Generation
    # max_counts = np.array([[gen['count_max'] for gen in hist] for hist in all_histories])
    # count_maxes = np.array([[gen['count_max'] for gen in trial] 
    #                         for sequence in all_histories_grouped 
    #                         for trial in sequence])

    count_maxes = np.array([
        [gen['count_max'] for gen in trial[:generations]] + 
        [np.nan] * (generations - len(trial))
        for sequence in all_histories_grouped
        for trial in sequence
    ])
    # avg_counts = count_maxes.mean(axis=0)
    avg_counts = np.nanmean(count_maxes, axis=0)


    plt.figure(figsize=(8,5))
    plt.plot(range(generations), avg_counts, marker="o")
    plt.xlabel("Generation")
    plt.ylabel(f"Count of puzzles with fitness = {target_fitness}")
    plt.title(f"Count of Max-Fitness (={target_fitness}) Puzzles per Generation")
    plt.grid(True)
    plt.show()


    # 2: BEST SEQUENCE FOR ACROSS ALL TRIALS
    """    
    # base_seq = random_sequence(16)
    # trials = 5
    # best_overall = None
    # best_fit_overall = -1

    # for _ in tqdm(range(trials), desc="Trials"):
    #     seq, fit, hist = evolve(base_seq, generations=50, popsize=30)
    #     if fit > best_fit_overall:
    #         best_fit_overall = fit
    #         best_overall = seq

    # print("Best sequence across all trials:", best_overall)
    # print("Fitness:", best_fit_overall)
    # """


    # 3: GETTING SEQUENCES WITH SIMILAR FITNESS FROM BASE SEQ (PUZZLE GENERATION) 
    """    
    base_seq = "HPPHHPHP"
    best_seq, best_fitness, history = evolve(base_seq, generations=50, popsize=30)

    target_fitness = fitness(base_seq)
    final_population = [best_seq]
    final_fitness = evaluate_population(final_population)

    unique_results = {}
    for score, seq in final_fitness:
        if seq not in unique_results:
            unique_results[seq] = score

    print("Unique sequences with similar fitness:")
    for seq, score in unique_results.items():
        if abs(score - target_fitness) <= 1:
            print(seq, score)
    """