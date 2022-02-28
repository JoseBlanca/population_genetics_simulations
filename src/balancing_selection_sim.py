from one_locus_two_alleles_simulator import (
    GenotypicFreqs,
    Fitness,
    INF,
    Population,
    AllelicFreqLogger,
)
from one_locus_two_alleles_simulator import simulate as one_locus_simulate
import plot


def simulate(
    allelic_freq_axess,
    freq_AA_pop1,
    freq_Aa_pop1,
    freq_aa_pop1,
    freq_AA_pop2,
    freq_Aa_pop2,
    freq_aa_pop2,
    pop1_size=INF,
    pop2_size=INF,
    num_generations=100,
    w11_pop1=1,
    w12_pop1=1,
    w22_pop1=1,
    w11_pop2=1,
    w12_pop2=1,
    w22_pop2=1,
    migration_rate_1_to_2=0,
    migration_rate_2_to_1=0,
):

    if (w11_pop1, w12_pop1, w22_pop1) == (1, 1, 1):
        fitness_pop1 = None
    else:
        fitness_pop1 = Fitness(w11=w11_pop1, w12=w12_pop1, w22=w22_pop1)
    if (w11_pop2, w12_pop2, w22_pop2) == (1, 1, 1):
        fitness_pop2 = None
    else:
        fitness_pop2 = Fitness(w11=w11_pop2, w12=w12_pop2, w22=w22_pop2)

    genotypic_freqs = GenotypicFreqs(
        freq_AA=freq_AA_pop1, freq_Aa=freq_Aa_pop1, freq_aa=freq_aa_pop1
    )
    pop1 = Population(
        id="pop1",
        size=pop1_size,
        genotypic_freqs=genotypic_freqs,
        fitness=fitness_pop1,
    )
    genotypic_freqs = GenotypicFreqs(
        freq_AA=freq_AA_pop2, freq_Aa=freq_Aa_pop2, freq_aa=freq_aa_pop2
    )
    pop2 = Population(
        id="pop2",
        size=pop2_size,
        genotypic_freqs=genotypic_freqs,
        fitness=fitness_pop2,
    )
    pops = {"pop1": pop1, "pop2": pop2}

    demographic_events = []
    if migration_rate_1_to_2 > 0:
        migration_1_2 = {
            "id": "migration_pop1_to_pop2",
            "type": "migration_start",
            "from_pop": pop1,
            "to_pop": pop2,
            "inmigrant_rate": migration_rate_1_to_2,
            "num_generation": 1,
        }
        demographic_events.append(migration_1_2)
    if migration_rate_2_to_1 > 0:
        migration_2_1 = {
            "id": "migration_pop2_to_pop1",
            "type": "migration_start",
            "from_pop": pop2,
            "to_pop": pop1,
            "inmigrant_rate": migration_rate_2_to_1,
            "num_generation": 1,
        }
        demographic_events.append(migration_2_1)
    allelic_freqs_logger = AllelicFreqLogger()
    one_locus_simulate(
        pops=[pop1, pop2],
        num_generations=num_generations,
        demographic_events=demographic_events,
        random_seed=None,
        loggers=[allelic_freqs_logger],
    )

    for idx, pop_id in enumerate(
        sorted(allelic_freqs_logger.values_per_generation.columns)
    ):
        axes = allelic_freq_axess[idx]
        plot.plot_allelic_freq_one_pop(
            allelic_freqs_logger,
            axes,
            pops[pop_id],
            take_color_from_color_wheel=True,
        )
    return {"allelic_freqs_logger": allelic_freqs_logger}


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    fig, axess = plt.subplots(nrows=2, sharex=True, figsize=(8, 4))

    simulate(
        allelic_freq_axess=axess,
        freq_AA=1,
        freq_Aa=0,
        freq_aa=0,
    )
