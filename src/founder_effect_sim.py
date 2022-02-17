import numpy

import msprime

import msprime_utils


def simulate(
    num_indis_to_sample,
    pop_sizes,
    num_founder_individuals,
    bottleneck_time,
    pop_founding_times,
    seq_length_in_bp,
    recomb_rate=1e-8,
    mutation_rate=1e-8,
    add_mutations=True,
    ploidy=2,
    random_seed=None,
):
    ancestral_pop_name = "pop1"
    demography = msprime.Demography()
    demography.add_population(
        name=ancestral_pop_name, initial_size=pop_sizes[0], initially_active=True
    )
    samples = {}
    msprime_utils.add_sample(
        num_samples=num_indis_to_sample,
        pop_name=ancestral_pop_name,
        ploidy=ploidy,
        time=pop_founding_times[0] - 1,
        samples=samples,
    )

    msprime_utils.add_sample(
        num_samples=num_indis_to_sample,
        pop_name=ancestral_pop_name,
        ploidy=ploidy,
        time=0,
        samples=samples,
    )

    for idx, (founding_time, founding_individuals, final_size) in enumerate(
        zip(pop_founding_times, num_founder_individuals, pop_sizes[1:])
    ):
        pop_name = f"pop{idx+2}"
        demography.add_population(
            name=pop_name, initial_size=founding_individuals, initially_active=True
        )
        demography.add_population_split(
            founding_time, derived=[pop_name], ancestral=ancestral_pop_name
        )
        demography.add_population_parameters_change(
            time=founding_time - 1,
            initial_size=founding_individuals,
            population=pop_name,
        )
        time_full_size = founding_time - bottleneck_time
        demography.add_population_parameters_change(
            time=time_full_size,
            initial_size=founding_individuals,
            population=pop_name,
        )
        demography.add_population_parameters_change(
            time=0, initial_size=final_size, population=pop_name
        )

        msprime_utils.add_sample(
            num_samples=num_indis_to_sample,
            pop_name=pop_name,
            ploidy=ploidy,
            time=time_full_size - 1,
            samples=samples,
        )
        msprime_utils.add_sample(
            num_samples=num_indis_to_sample,
            pop_name=pop_name,
            ploidy=ploidy,
            time=0,
            samples=samples,
        )

    demography.sort_events()

    sim_result = msprime_utils.simulate(
        sample_sets=list(samples.values()),
        demography=demography,
        seq_length_in_bp=seq_length_in_bp,
        recomb_rate=recomb_rate,
        random_seed=random_seed,
        add_mutations=add_mutations,
        mutation_rate=mutation_rate,
        ploidy=ploidy,
    )
    return {"sim_result": sim_result, "samples": samples}


if __name__ == "__main__":
    res = simulate(
        num_indis_to_sample=40,
        pop_sizes=[1_000_000] * 2,
        num_founder_individuals=[100],
        bottleneck_time=50,
        pop_founding_times=[100],
        seq_length_in_bp=100_00,
        recomb_rate=1e-8,
        mutation_rate=1e-8,
        add_mutations=True,
        ploidy=2,
        random_seed=None,
    )
    sim_result = res["sim_result"]
    print(sim_result.calculate_poly095_per_sample(samples=list(res["samples"])))

    if False:
        genotypes = sim_result.get_genotypes(sampling_times=[101], pop_names=None)
        print(genotypes.keep_only_biallelic().poly095)
        # print(genotypes.keep_only_biallelic().folded_sfs)
        genotypes = sim_result.get_genotypes(sampling_times=[79], pop_names=None)
        print(genotypes.keep_only_biallelic().poly095)
        # print(genotypes.keep_only_biallelic().folded_sfs)
        genotypes = sim_result.get_genotypes(sampling_times=[0], pop_names=None)
        # print(genotypes.keep_only_biallelic().folded_sfs)
        print(genotypes.keep_only_biallelic().poly095)
