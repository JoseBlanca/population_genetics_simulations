GENTOTYPIC_FREQS_LABELS = {
    "freqs_AA": "freq(AA)",
    "freqs_Aa": "freq(Aa)",
    "freqs_aa": "freq(aa)",
    "freq_AA": "freq(AA)",
    "freq_Aa": "freq(Aa)",
    "freq_aa": "freq(aa)",
}
COLORS = {
    "freq_A": "tab:blue",
    "freq_AA": "tab:blue",
    "freqs_AA": "tab:blue",
    "freq_Aa": "tab:green",
    "freqs_Aa": "tab:green",
    "freq_aa": "tab:orange",
    "freqs_aa": "tab:orange",
}


def plot_allelic_freq_one_pop(
    allelic_freqs_logger, axes, pop, take_color_from_color_wheel=False
):
    freqs = allelic_freqs_logger.values_per_generation[pop.id]

    kwargs = {}
    if take_color_from_color_wheel:
        kwargs["color"] = COLORS["freq_A"]

    axes.plot(freqs.index, freqs, **kwargs)
    axes.set_ylabel("allelic freq A")
    axes.set_ylim([0, 1])


def plot_genotypic_freqs_one_pop(genotypic_freqs_logger, axes, pop):
    freqs = genotypic_freqs_logger.values_per_generation
    pop_id = pop.id

    labels = GENTOTYPIC_FREQS_LABELS

    for genotypic_freq_name, this_genotypic_freqs in freqs.items():
        color = COLORS[genotypic_freq_name]
        axes.plot(
            this_genotypic_freqs[pop_id].index,
            this_genotypic_freqs[pop_id],
            label=labels[genotypic_freq_name],
            color=color,
        )
    axes.set_ylim([0, 1])
    axes.set_ylabel("genotypic freqs")
    axes.legend()
