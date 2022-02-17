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


def add_sample(num_samples: int, ploidy: int, pop_name: str, time: int, samples: dict):
    sample = msprime.SampleSet(
        num_samples=num_samples, population=pop_name, ploidy=ploidy, time=time
    )
    samples[(pop_name, time)] = sample


class SeveralPopsSimulationResult:
    def __init__(self, tree_seqs, sampling_times, demography, ploidy):
        self.tree_seqs = tree_seqs
        self.sampling_times = sampling_times
        self.demography = demography
        self.ploidy = ploidy

    def _get_samples_from_sample_info(self, samples_info):
        samples = []
        tree_seqs = self.tree_seqs
        pop_ids_by_pop_name_in_tseq, _ = self._get_pop_ids_and_names()
        for pop_name, sampling_time in samples_info:
            pop_id = pop_ids_by_pop_name_in_tseq[pop_name]
            sample = {
                "sample_node_ids": tree_seqs.samples(
                    population=pop_id, time=sampling_time
                ),
                "sampling_time": sampling_time,
                "pop_name": pop_name,
                "sample_name": f"{pop_name}_{sampling_time}",
            }
            if sample["sample_node_ids"].size == 0:
                continue
            samples.append(sample)
        return samples

    def _get_samples(self, sampling_times=None, pop_names=None, samples=None):

        if samples is not None:
            if sampling_times is not None or pop_names is not None:
                msg = (
                    "if samples are given, sampling_times and pop_names should be None"
                )
                raise ValueError(msg)
            return self._get_samples_from_sample_info(samples)

        if sampling_times is None:
            sampling_times = [None]
        elif isinstance(sampling_times, (int, float)):
            sampling_times = [sampling_times]
        if pop_names is not None:
            pop_names = list(pop_names)

        if len(sampling_times) == 1:
            return self._get_samples_for_pops(
                sampling_time=sampling_times[0], pop_names=pop_names
            )
        elif pop_names is None or len(pop_names) == 1:
            if pop_names is None:
                pop_name = pop_names
            else:
                pop_name = pop_names[0]
            return self._get_samples_for_times(
                sampling_times=sampling_times, pop_name=pop_name
            )
        else:
            raise NotImplementedError("Fixme")

    def _get_pop_ids_and_names(self):
        tree_seqs = self.tree_seqs

        pop_ids_by_pop_name_in_tseq = {}
        pop_names_by_pop_id_in_tseq = {}
        for pop in tree_seqs.populations():
            pop_name = pop.metadata["name"]
            pop_id = pop.id
            pop_ids_by_pop_name_in_tseq[pop_name] = pop_id
            pop_names_by_pop_id_in_tseq[pop_id] = pop_name
        return pop_ids_by_pop_name_in_tseq, pop_names_by_pop_id_in_tseq

    def _get_samples_for_times(self, sampling_times, pop_name=None):
        res = self._get_pop_ids_and_names()
        pop_ids_by_pop_name_in_tseq, _ = res

        if pop_name is None:
            if len(pop_ids_by_pop_name_in_tseq) == 1:
                pop_id = list(pop_ids_by_pop_name_in_tseq.values())[0]
            else:
                msg = "You have more than one pop, so you have to set one pop"
                raise ValueError(msg)
        else:
            pop_id = pop_ids_by_pop_name_in_tseq[pop_name]

        for sampling_time in sampling_times:
            if sampling_time not in self.sampling_times:
                msg = f"Sampling time ({sampling_time}) not in available sampling times ({self.sampling_times})"
                raise ValueError(msg)

        tree_seqs = self.tree_seqs

        samples = []
        for sampling_time in sampling_times:
            sample = {
                "sample_node_ids": tree_seqs.samples(
                    population=pop_id, time=sampling_time
                ),
                "sampling_time": sampling_time,
                "pop_name": pop_name,
                "sample_name": f"{pop_name}_{sampling_time}",
            }
            if sample["sample_node_ids"].size == 0:
                continue
            samples.append(sample)
        return samples

    def _get_samples_for_pops(self, sampling_time=None, pop_names=None):

        res = self._get_pop_ids_and_names()
        pop_ids_by_pop_name_in_tseq, pop_names_by_pop_id_in_tseq = res

        if pop_names is None:
            pop_names = list(sorted(pop_ids_by_pop_name_in_tseq.keys()))

        if sampling_time is None:
            if len(self.sampling_times) == 1:
                sampling_time = self.sampling_times[0]
            else:
                msg = "You have sample at different times, so you have to set one sample time"
                raise ValueError(msg)
        if sampling_time not in self.sampling_times:
            msg = f"Sampling time ({sampling_time}) not in available sampling times ({self.sampling_times})"
            raise ValueError(msg)

        if pop_names:
            pop_ids = [pop_ids_by_pop_name_in_tseq[pop_name] for pop_name in pop_names]
        else:
            pop_ids = sorted(pop_ids_by_pop_name_in_tseq.values())

        tree_seqs = self.tree_seqs

        samples = []
        for pop_id in pop_ids:
            pop_name = (pop_names_by_pop_id_in_tseq[pop_id],)
            sample = {
                "sample_node_ids": tree_seqs.samples(
                    population=pop_id, time=sampling_time
                ),
                "sampling_time": sampling_time,
                "pop_name": pop_name,
                "sample_name": pop_name,
            }
            if sample["sample_node_ids"].size == 0:
                continue
            samples.append(sample)
        return samples

    def _get_sample_sets_and_sample_names_from_samples(self, samples):
        sample_sets = []
        pop_names = []
        for sample in samples:
            sample_sets.append(sample["sample_node_ids"])
            pop_names.append(sample["sample_name"])
        return sample_sets, pop_names

    def calculate_nucleotide_diversities_per_sample(
        self,
        sampling_times=None,
        pop_names=None,
        samples=None,
    ):
        """Computes mean genetic diversity (also known as “pi”),

        a common citation for the definition is Nei and Li (1979) (equation 22)"""
        samples = self._get_samples(sampling_times=sampling_times, pop_names=pop_names)
        sample_names = [sample["sample_name"] for sample in samples]
        sample_sets, pop_names = self._get_sample_sets_and_sample_names_from_samples(
            samples
        )

        diversities = self.tree_seqs.diversity(sample_sets=sample_sets)
        diversities = pandas.Series(
            diversities,
            index=sample_names,
        )
        return diversities

    def calculate_poly095_per_sample(
        self, sampling_times=None, pop_names=None, samples=None
    ):
        afss = self.calculate_allele_frequency_spectrum(
            sampling_times=sampling_times, pop_names=pop_names, samples=samples
        )
        poly_095 = []
        sample_names = list(afss.keys())
        for sample_name in sample_names:
            afs = afss[sample_name]
            value = numpy.sum(
                afs.num_loci_for_each_allele_freq[afs.allele_freqs > (1 - 0.95)]
            ) / numpy.sum(afs.num_loci_for_each_allele_freq)
            poly_095.append(value)
        poly_095 = pandas.Series(poly_095, index=sample_names)
        return poly_095

    def _generate_combination_idxs(self, items):
        num_pops = len(items)
        if num_pops == 1:
            raise RuntimeError("Only one pop, no pairwise comparison posible")
        return list(itertools.combinations(list(range(num_pops)), 2))

    def calculate_fsts(self, sampling_times=None, pop_names=None, samples=None):
        samples = self._get_samples(
            sampling_times=sampling_times, pop_names=pop_names, samples=samples
        )
        sample_sets, pop_names = self._get_sample_sets_and_sample_names_from_samples(
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

    def calculate_allele_frequency_spectrum(
        self, sampling_times=None, pop_names=None, samples=None
    ):
        samples = self._get_samples(
            sampling_times=sampling_times, pop_names=pop_names, samples=samples
        )
        afss = {}
        for sample in samples:
            afs_counts = self.tree_seqs.allele_frequency_spectrum(
                sample_sets=[sample["sample_node_ids"]],
                span_normalise=False,
                polarised=False,
            )[1:-1]
            num_chroms = afs_counts.size + 1
            allele_freqs = numpy.arange(1, num_chroms) / (num_chroms)
            afss[sample["sample_name"]] = AFS(afs_counts, allele_freqs, num_chroms)
        return afss

    def get_genotypes(self, sampling_times=None, pop_names=None, samples=None):
        samples = self._get_samples(
            sampling_times=sampling_times, pop_names=pop_names, samples=samples
        )
        node_idxs = []
        pops = []
        for sample in samples:
            node_idxs_for_these_samples = sample["sample_node_ids"]
            pops_for_this_sample = [
                sample["sample_name"]
            ] * node_idxs_for_these_samples.size
            node_idxs.extend(node_idxs_for_these_samples)
            pops.extend(pops_for_this_sample)

        genotypes = self.tree_seqs.genotype_matrix()
        genotypes = genotypes[:, node_idxs]
        genotypes = Genotypes.from_tree_seq_genotypes(
            genotypes, classification=pops, ploidy=self.ploidy
        )
        return genotypes


def simulate(
    sample_sets,
    demography,
    seq_length_in_bp,
    recomb_rate=1e-8,
    add_mutations=True,
    random_seed=None,
    mutation_rate=1e-8,
    ploidy=2,
):
    tree_seqs = msprime.sim_ancestry(
        samples=sample_sets,
        demography=demography,
        recombination_rate=recomb_rate,
        sequence_length=seq_length_in_bp,
        random_seed=random_seed,
    )
    if add_mutations:
        tree_seqs = msprime.sim_mutations(
            tree_seqs, rate=mutation_rate, random_seed=random_seed
        )

    sampling_times = sorted({sample_set.time for sample_set in sample_sets})

    return SeveralPopsSimulationResult(
        tree_seqs=tree_seqs,
        sampling_times=sampling_times,
        demography=demography,
        ploidy=ploidy,
    )
