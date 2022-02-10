import msprime

import msprime_utils


def simulate_drifting_pops(
    num_indis_to_sample_per_pop,
    ancestral_pop_size,
    drifting_pops_sizes: dict[int],
    num_generations_ago_when_split_happened,
    seq_length_in_bp,
    recomb_rate=1e-8,
    add_mutations=True,
    ploidy=2,
    random_seed=None,
    sampling_times=(0,),
):
    ancestral_pop_name = "ancestral"
    demography = msprime.Demography()
    demography.add_population(name=ancestral_pop_name, initial_size=ancestral_pop_size)

    derived_pops = []
    samples = []
    for pop_name, size in drifting_pops_sizes.items():
        demography.add_population(name=pop_name, initial_size=size)
        derived_pops.append(pop_name)

        for time in sampling_times:
            sample = msprime.SampleSet(
                num_samples=num_indis_to_sample_per_pop,
                population=pop_name,
                ploidy=ploidy,
                time=time,
            )
            samples.append(sample)

    demography.add_population_split(
        num_generations_ago_when_split_happened,
        derived=derived_pops,
        ancestral=ancestral_pop_name,
    )

    return msprime_utils.simulate(
        sample_sets=samples,
        demography=demography,
        seq_length_in_bp=seq_length_in_bp,
        recomb_rate=recomb_rate,
        random_seed=random_seed,
        add_mutations=add_mutations,
        ploidy=ploidy,
    )


if __name__ == "__main__":
    num_indis_to_sample_per_pop = 20
    num_generations_ago_when_split_happened = 2000
    ancestral_pop_size = 10000
    seq_length_in_bp = 500000
    drifting_pops_sizes = {
        "pop1": ancestral_pop_size,
        "pop2": ancestral_pop_size,
        "pop3": ancestral_pop_size / 10,
    }
    sampling_times = [
        num_generations_ago_when_split_happened - 1,
        1500,
        1000,
        500,
        0,
    ]
    res = simulate_drifting_pops(
        num_indis_to_sample_per_pop=num_indis_to_sample_per_pop,
        drifting_pops_sizes=drifting_pops_sizes,
        ancestral_pop_size=ancestral_pop_size,
        num_generations_ago_when_split_happened=num_generations_ago_when_split_happened,
        seq_length_in_bp=seq_length_in_bp,
        sampling_times=sampling_times,
    )

    import pca

    for sampling_time in sampling_times:
        print(f"{sampling_time=}")
        print(res.calculate_nucleotide_diversities_per_pop(sampling_time=sampling_time))
        print(
            res.calculate_nucleotide_diversities_per_pop(
                sampling_time=sampling_time, pop_names=["pop1"]
            )
        )
        print(res.calculate_fsts(sampling_time=sampling_time))
        print(
            res.calculate_fsts(sampling_time=sampling_time, pop_names=["pop1", "pop3"])
        )
        print(
            res.calculate_allele_frequency_spectrum(sampling_time=sampling_time)[
                "pop1"
            ].num_loci_for_each_allele_freq
        )
        genotypes = res.get_genotypes(sampling_time=sampling_time).keep_only_biallelic()
        print(genotypes.mean_num_alleles_per_variant)
        pca_res = pca.do_pca(genotypes)

        import matplotlib.pyplot as plt

        fig = plt.figure()
        axes = fig.add_subplot()
        pca.plot_pca_result(pca_res, axes, classification=genotypes.classification)
        plot_path = f"/home/jose/tmp/pca_time_{sampling_time}.svg"
        fig.savefig(plot_path)
