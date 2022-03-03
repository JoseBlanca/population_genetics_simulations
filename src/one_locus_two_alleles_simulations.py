from one_locus_two_alleles_simulator import (
    GenotypicFreqs,
    Fitness,
    MutRates,
    INF,
    Population,
    AllelicFreqLogger,
    GenotypicFreqsLogger,
    ExpHetLogger,
    simulate,
)
import plot


def simulate_one_locus_two_alleles_one_pop(
    allelic_freq_axes,
    genotypic_freqs_axes,
    freq_AA,
    freq_Aa,
    freq_aa,
    pop_size=INF,
    num_generations=100,
    w11=1,
    w12=1,
    w22=1,
    A2a=0,
    a2A=0,
    selfing_rate=0,
    num_populations=1,
):

    if genotypic_freqs_axes is not None and num_populations != 1:
        raise ValueError(
            "if num_populations > 1 then genotypic_freqs_axes should be None"
        )

    if (w11, w12, w22) == (1, 1, 1):
        fitness = None
    else:
        fitness = Fitness(w11=w11, w12=w12, w22=w22)

    if A2a == 0 and a2A == 0:
        mut_rates = None
    else:
        mut_rates = MutRates(A2a, a2A)

    loggers = {
        "genotypic_freqs_logger": [],
        "allelic_freqs_logger": [],
        "exp_het_logger": [],
    }
    for idx in range(num_populations):
        genotypic_freqs = GenotypicFreqs(
            freq_AA=freq_AA, freq_Aa=freq_Aa, freq_aa=freq_aa
        )
        pop = Population(
            id=f"pop{idx}",
            size=pop_size,
            genotypic_freqs=genotypic_freqs,
            fitness=fitness,
            mut_rates=mut_rates,
            selfing_rate=selfing_rate,
        )

        allelic_freqs_logger = AllelicFreqLogger()
        genotypic_freqs_logger = GenotypicFreqsLogger()
        exp_het_logger = ExpHetLogger()
        simulate(
            pops=[pop],
            num_generations=num_generations,
            demographic_events=None,
            random_seed=None,
            loggers=[allelic_freqs_logger, genotypic_freqs_logger, exp_het_logger],
        )
        loggers["genotypic_freqs_logger"].append(genotypic_freqs_logger)
        loggers["allelic_freqs_logger"].append(allelic_freqs_logger)
        loggers["exp_het_logger"].append(exp_het_logger)

        plot.plot_allelic_freq_one_pop(
            allelic_freqs_logger,
            allelic_freq_axes,
            pop,
            take_color_from_color_wheel=True,
        )
    loggers = {
        key: logger_list[0].__class__.from_loggers(logger_list)
        for key, logger_list in loggers.items()
    }

    if genotypic_freqs_axes is not None:
        plot.plot_genotypic_freqs_one_pop(
            genotypic_freqs_logger, genotypic_freqs_axes, pop
        )
        genotypic_freqs_axes.set_xlabel("Num. generations")
    return loggers


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    fig, axess = plt.subplots(nrows=2, sharex=True, figsize=(8, 4))

    simulate_one_locus_two_alleles_one_pop(
        allelic_freq_axes=axess[0],
        # genotypic_freqs_axes=axess[1],
        genotypic_freqs_axes=None,
        freq_AA=1,
        freq_Aa=0,
        freq_aa=0,
        num_populations=3,
    )
