import msprime

import msprime_utils


def simulate_one_simple_pop(
    num_indis_to_sample,
    pop_size,
    seq_length_in_bp,
    recomb_rate=1e-8,
    add_mutations=True,
    ploidy=2,
    random_seed=None,
):

    demography = msprime.Demography()
    demography.add_population(name="pop1", initial_size=pop_size)

    samples_pop1 = msprime.SampleSet(
        num_samples=num_indis_to_sample, population="pop1", ploidy=ploidy
    )
    samples = [samples_pop1]
    return msprime_utils.simulate(
        samples=samples,
        demography=demography,
        seq_length_in_bp=seq_length_in_bp,
        recomb_rate=recomb_rate,
        random_seed=random_seed,
        add_mutations=add_mutations,
    )


if __name__ == "__main__":
    num_individuals = 1000
    sim_res = simulate_one_simple_pop(
        num_indis_to_sample=100, pop_size=1_000_000, seq_length_in_bp=5000
    )
    afss = sim_res.allele_frequency_spectrum
    afs = next(iter(afss.values()))
    print(afs.num_loci_for_each_allele_freq)
    print(sim_res.nucleotide_diversities_per_pop)
    sim_res = simulate_one_simple_pop(
        num_indis_to_sample=100, pop_size=10_000, seq_length_in_bp=5000
    )
    afss = sim_res.allele_frequency_spectrum
    afs = next(iter(afss.values()))
    print(afs.num_loci_for_each_allele_freq)
    print(sim_res.nucleotide_diversities_per_pop)
