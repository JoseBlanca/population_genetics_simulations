def plot_freqs(plot_path, allelic_freqs_logger, genotypic_freqs_logger, pop_id):
    figure, axess = plt.subplots(nrows=2, sharex=True)
    allelic_freq = allelic_freqs_logger.values_per_generation[pop_id]

    allelic_freq_axes = axess[0]
    allelic_freq_axes.plot(allelic_freq.index, allelic_freq)
    allelic_freq_axes.set_ylabel("allelic freq A")

    genotypic_freq_axes = axess[1]
    genotypic_freqs = genotypic_freqs_logger.values_per_generation
    for genotypic_freq_name, this_genotypic_freqs in genotypic_freqs.items():
        genotypic_freq_axes.plot(
            this_genotypic_freqs[pop_id].index,
            this_genotypic_freqs[pop_id],
            label=genotypic_freq_name,
        )
    genotypic_freq_axes.legend()

    figure.savefig(plot_path)


def plot_allelic_freq_one_pop(allelic_freqs_logger, axes, pop):
    freqs = allelic_freqs_logger.values_per_generation[pop.id]
    axes.plot(freqs.index, freqs)
    axes.set_ylabel("allelic freq A")


def plot_genotypic_freqs_one_pop(genotypic_freqs_logger, axes, pop):
    freqs = genotypic_freqs_logger.values_per_generation
    pop_id = pop.id

    for genotypic_freq_name, this_genotypic_freqs in freqs.items():
        axes.plot(
            this_genotypic_freqs[pop_id].index,
            this_genotypic_freqs[pop_id],
            label=genotypic_freq_name,
        )
    axes.legend()
