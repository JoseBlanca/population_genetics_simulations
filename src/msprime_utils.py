import itertools

import numpy
import pandas

import msprime


class AFS:
    def __init__(self, num_loci_for_each_allele_freq, allele_freqs, num_samples):
        self.num_loci_for_each_allele_freq = num_loci_for_each_allele_freq
        self.allele_freqs = allele_freqs
        self.num_samples = num_samples


class SeveralPopsSimulationResult:
    def __init__(self, tree_seqs, samples, demography):
        self.tree_seqs = tree_seqs
        self.samples = samples
        self.demography = demography

    @property
    def pop_names(self):
        return [sample_set.population for sample_set in self.samples]

    @property
    def sample_seq_ids_per_pop(self):

        seq_ids = [
            self.tree_seqs.samples(population=pop_idx)
            for pop_idx, _ in enumerate(self.samples)
        ]
        return seq_ids

    @property
    def nucleotide_diversities_per_pop(self):
        """Computes mean genetic diversity (also known as “pi”),
        a common citation for the definition is Nei and Li (1979) (equation 22)"""
        diversities = self.tree_seqs.diversity(sample_sets=self.sample_seq_ids_per_pop)
        diversities = pandas.Series(
            diversities,
            index=self.pop_names,
        )
        return diversities

    @property
    def _pairwise_comparison_indexes_between_pops(self):
        num_pops = len(self.pop_names)
        if num_pops == 1:
            raise RuntimeError("Only one pop, no pairwise comparison posible")
        return list(itertools.combinations(list(range(num_pops)), 2))

    @property
    def Fsts(self):
        raise NotImplementedError(
            "Almost done, just figure out how to create the square matrix"
        )
        samples = self.sample_seq_ids_per_pop
        print(samples)
        pairwise_comparisons = self._pairwise_comparison_indexes_between_pops
        print(pairwise_comparisons)
        print(
            self.tree_seqs.Fst(
                sample_sets=self.sample_seq_ids_per_pop, indexes=pairwise_comparisons
            )
        )

    @property
    def allele_frequency_spectrum(self):
        afss = {}
        for pop, samples_in_pop in zip(self.pop_names, self.sample_seq_ids_per_pop):
            afs_counts = self.tree_seqs.allele_frequency_spectrum(
                sample_sets=[samples_in_pop],
                span_normalise=False,
                polarised=False,
            )[1:-1]
            num_chroms = afs_counts.size + 1
            allele_freqs = numpy.arange(1, num_chroms) / (num_chroms)
            afss[pop] = AFS(afs_counts, allele_freqs, num_chroms)
        return afss


def simulate(
    samples,
    demography,
    seq_length_in_bp,
    recomb_rate=1e-8,
    add_mutations=True,
    random_seed=None,
):
    tree_seqs = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        recombination_rate=recomb_rate,
        sequence_length=seq_length_in_bp,
        random_seed=random_seed,
    )
    if add_mutations:
        tree_seqs = msprime.sim_mutations(tree_seqs, rate=1e-8, random_seed=54321)

    return SeveralPopsSimulationResult(
        tree_seqs=tree_seqs, samples=samples, demography=demography
    )
