import msprime

import msprime_utils


def simulate(
    num_indis_to_sample,
    pop_size,
    pop_size_during_bottleneck,
    bottleneck_start_time,
    bottleneck_end_time,
    seq_length_in_bp,
    recomb_rate=1e-8,
    mutation_rate=1e-8,
    add_mutations=True,
    ploidy=2,
    random_seed=None,
):
    pop_name = "pop"
    demography = msprime.Demography()
    demography.add_population(name=pop_name, initial_size=pop_size)

    demography.add_population_parameters_change(
        time=bottleneck_start_time, initial_size=pop_size_during_bottleneck
    )
    demography.add_population_parameters_change(
        time=bottleneck_end_time, initial_size=pop_size
    )

    sampling_times = [
        bottleneck_start_time + 1,
        bottleneck_start_time - 1,
        int((bottleneck_end_time + bottleneck_start_time) / 2),
        bottleneck_end_time + 1,
        bottleneck_end_time - 1,
        0,
    ]
    samples = []
    for time in sampling_times:
        sample = msprime.SampleSet(
            num_samples=num_indis_to_sample,
            population=pop_name,
            ploidy=ploidy,
            time=time,
        )
        samples.append(sample)

    sim_result = msprime_utils.simulate(
        sample_sets=samples,
        demography=demography,
        seq_length_in_bp=seq_length_in_bp,
        recomb_rate=recomb_rate,
        random_seed=random_seed,
        add_mutations=add_mutations,
        mutation_rate=mutation_rate,
        ploidy=ploidy,
    )
    sampling_times = list(reversed(sorted(sampling_times)))
    return {"sim_result": sim_result, "sampling_times": sampling_times}
