def plot_allelic_freq_one_pop(allelic_freqs_logger, axes, pop):
    freqs = allelic_freqs_logger.values_per_generation[pop.id]
    axes.plot(freqs.index, freqs)
    axes.set_ylabel("allelic freq A")
    axes.set_ylim([0, 1])


def plot_genotypic_freqs_one_pop(genotypic_freqs_logger, axes, pop):
    freqs = genotypic_freqs_logger.values_per_generation
    pop_id = pop.id

    labels = {"freqs_AA": "freq(AA)", "freqs_Aa": "freq(Aa)", "freqs_aa": "freq(aa)"}

    for genotypic_freq_name, this_genotypic_freqs in freqs.items():
        axes.plot(
            this_genotypic_freqs[pop_id].index,
            this_genotypic_freqs[pop_id],
            label=labels[genotypic_freq_name],
        )
    axes.set_ylim([0, 1])
    axes.set_ylabel("genotypic freqs")
    axes.legend()
