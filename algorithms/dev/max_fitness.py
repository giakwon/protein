import subprocess
import numpy.random as npr
import random
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from multiprocessing import Pool
import multiprocessing as mp


fitness_cache = {}

def fitness(seq):
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

    f = len(conformations)
    fitness_cache[seq] = f
    return f


def evaluate_population(population):
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
    seq_len = len(base_seq)
    population = [random_sequence(seq_len) for _ in range(popsize)]
    history = []  # will store counts of max-fitness individuals

    for gen in tqdm(range(generations), desc="Evolving", leave=False):
        scored = evaluate_population(population)

        # count how many reached target_fitness
        count_max = sum(1 for f, _ in scored if f == target_fitness)
        history.append(count_max)

        scored.sort(reverse=True)
        new_population = []
        if elitism:
            best_fitness, best_seq = scored[0]
            new_population.append(best_seq)

        while len(new_population) < popsize:
            p1 = select(scored)
            p2 = select(scored)
            c1, c2 = crossover(p1, p2)
            c1 = mutate(c1, mutation_rate)
            c2 = mutate(c2, mutation_rate)
            new_population.extend([c1, c2])
        population = new_population[:popsize]

    return history


if __name__ == "__main__":
    seq_length = 16
    num_sequences = 3
    trials = 5
    generations = 50
    popsize = 30
    target_fitness = 9

    all_histories = []

    for _ in tqdm(range(num_sequences), desc="Sequences", leave=False):
        base_seq = random_sequence(seq_length)
        for _ in range(trials):
            hist = evolve(base_seq, generations, popsize, target_fitness=target_fitness)
            all_histories.append(hist)

    # avg across trials
    histories = np.array(all_histories)
    avg_counts = histories.mean(axis=0)

    # plot
    plt.figure(figsize=(8,5))
    plt.plot(range(generations), avg_counts, marker="o")
    plt.xlabel("Generation")
    plt.ylabel(f"Count of puzzles with fitness = {target_fitness}")
    plt.title(f"Count of Max-Fitness (={target_fitness}) Puzzles per Generation")
    plt.grid(True)
    plt.show()
