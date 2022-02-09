import itertools

import numpy
import pandas

import msprime

from genotypes import Genotypes


class AFS:
    def __init__(self, num_loci_for_each_allele_freq, allele_freqs, num_samples):
        self.num_loci_for_each_allele_freq = num_loci_for_each_allele_freq
        self.allele_freqs = allele_freqs
        self.num_samples = num_samples


class SeveralPopsSimulationResult:
    def __init__(self, tree_seqs, sampling_times, demography):
        self.tree_seqs = tree_seqs
        self.sampling_times = sampling_times
        self.demography = demography

    def _get_samples(self, sampling_time=None, pop_names=None):
        tree_seqs = self.tree_seqs

        pop_ids_by_pop_name_in_tseq = {}
        pop_names_by_pop_id_in_tseq = {}
        for pop in tree_seqs.populations():
            pop_name = pop.metadata["name"]
            pop_id = pop.id
            pop_ids_by_pop_name_in_tseq[pop_name] = pop_id
            pop_names_by_pop_id_in_tseq[pop_id] = pop_name

        if pop_names:
            pop_ids = [pop_ids_by_pop_name_in_tseq[pop_name] for pop_name in pop_names]
        else:
            pop_ids = sorted(pop_ids_by_pop_name_in_tseq.values())

        if sampling_time is None:
            if len(self.sampling_times) == 1:
                sampling_time = self.sampling_times[0]
            else:
                msg = "You have sample at different times, so you have to set one sample time"
                raise ValueError(msg)
        if sampling_time not in self.sampling_times:
            msg = f"Sampling time ({sampling_time}) not in available sampling times ({self.sampling_times})"
            raise ValueError(msg)

        samples = []
        for pop_id in pop_ids:
            sample = {
                "sample_node_ids": tree_seqs.samples(
                    population=pop_id, time=sampling_time
                ),
                "sampling_time": sampling_time,
                "pop_name": pop_names_by_pop_id_in_tseq[pop_id],
            }
            if sample["sample_node_ids"].size == 0:
                continue
            samples.append(sample)
        return samples

    def _get_sample_sets_and_pop_names_from_samples(self, samples):
        sample_sets = []
        pop_names = []
        for sample in samples:
            sample_sets.append(sample["sample_node_ids"])
            pop_names.append(sample["pop_name"])
        return sample_sets, pop_names

    def calculate_nucleotide_diversities_per_pop(
        self, sampling_time=None, pop_names=None
    ):
        """Computes mean genetic diversity (also known as “pi”),

        a common citation for the definition is Nei and Li (1979) (equation 22)"""
        samples = self._get_samples(sampling_time=sampling_time, pop_names=pop_names)
        sample_sets, pop_names = self._get_sample_sets_and_pop_names_from_samples(
            samples
        )

        diversities = self.tree_seqs.diversity(sample_sets=sample_sets)
        diversities = pandas.Series(
            diversities,
            index=pop_names,
        )
        return diversities

    def _generate_combination_idxs(self, items):
        num_pops = len(items)
        if num_pops == 1:
            raise RuntimeError("Only one pop, no pairwise comparison posible")
        return list(itertools.combinations(list(range(num_pops)), 2))

    def calculate_fsts(self, sampling_time=None, pop_names=None):
        samples = self._get_samples(sampling_time=sampling_time, pop_names=pop_names)
        sample_sets, pop_names = self._get_sample_sets_and_pop_names_from_samples(
            samples
        )

        num_pops = len(pop_names)
        pairwise_comparisons = self._generate_combination_idxs(pop_names)
        fsts = self.tree_seqs.Fst(sample_sets=sample_sets, indexes=pairwise_comparisons)

        fst_matrix = numpy.zeros((num_pops, num_pops), dtype=fsts.dtype)
        for (idx1, idx2), fst in zip(pairwise_comparisons, fsts):
            fst_matrix[idx1, idx2] = fst
            fst_matrix[idx2, idx1] = fst
        fst_matrix = pandas.DataFrame(fst_matrix, index=pop_names, columns=pop_names)
        return fst_matrix

    def calculate_allele_frequency_spectrum(self, sampling_time=None, pop_names=None):
        samples = self._get_samples(sampling_time=sampling_time, pop_names=pop_names)

        afss = {}
        for sample in samples:
            afs_counts = self.tree_seqs.allele_frequency_spectrum(
                sample_sets=[sample["sample_node_ids"]],
                span_normalise=False,
                polarised=False,
            )[1:-1]
            num_chroms = afs_counts.size + 1
            allele_freqs = numpy.arange(1, num_chroms) / (num_chroms)
            afss[sample["pop_name"]] = AFS(afs_counts, allele_freqs, num_chroms)
        return afss

    def get_genotypes(self, sampling_time=None, pop_names=None):
        """Loci will be in rows and samples in columns

        Returns an m x n numpy array of the genotypes, where m is the number of sites
        and n the number of samples."""
        samples = self._get_samples(sampling_time=sampling_time, pop_names=pop_names)
        node_idxs = []
        pops = []
        for sample in samples:
            node_idxs_for_these_samples = sample["sample_node_ids"]
            pops_for_this_sample = [
                sample["pop_name"]
            ] * node_idxs_for_these_samples.size
            node_idxs.extend(node_idxs_for_these_samples)
            pops.extend(pops_for_this_sample)

        genotypes = self.tree_seqs.genotype_matrix()
        genotypes = genotypes[:, node_idxs]
        genotypes = Genotypes(genotypes, classification=pops)
        return genotypes


def simulate(
    sample_sets,
    demography,
    seq_length_in_bp,
    recomb_rate=1e-8,
    add_mutations=True,
    random_seed=None,
):
    tree_seqs = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        recombination_rate=recomb_rate,
        sequence_length=seq_length_in_bp,
        random_seed=random_seed,
    )
    if add_mutations:
        tree_seqs = msprime.sim_mutations(tree_seqs, rate=1e-8, random_seed=54321)

    sampling_times = sorted({sample_set.time for sample_set in sample_sets})

    return SeveralPopsSimulationResult(
        tree_seqs=tree_seqs, sampling_times=sampling_times, demography=demography
    )
