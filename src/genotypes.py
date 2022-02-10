import numpy
import pandas

import allel


class Genotypes:
    def __init__(self, genotypes, classification=None, indiviual_names=None):
        self.genotypes = genotypes
        self._classifcation = None
        self.classification = classification
        self.indiviual_names = indiviual_names

    @classmethod
    def from_tree_seq_genotypes(
        cls, genotypes, ploidy, classification=None, indiviual_names=None
    ):
        genotypes = allel.HaplotypeArray(genotypes).to_genotypes(ploidy=ploidy)

        if classification is not None:
            classification = classification[::ploidy]
        return cls(
            genotypes=genotypes,
            classification=classification,
            indiviual_names=indiviual_names,
        )

    def _get_classification(self):
        return self._classfication

    def _set_classification(self, classification):
        try:
            len_class = classification.size
        except AttributeError:
            len_class = len(classification)
        if len_class != self.num_individuals:
            raise ValueError("classification length and num_individuals should match")
        self._classfication = classification

    classification = property(_get_classification, _set_classification)

    @property
    def num_individuals(self):
        return self.genotypes.n_samples

    @property
    def num_variants(self):
        return self.genotypes.n_variants

    def count_alleles_per_variant(self):
        return self.genotypes.count_alleles()

    def remove_singletons(self):
        allele_count = self.count_alleles_per_variant()
        mask = allele_count.max_allele() == 1
        genotypes = self.genotypes.compress(mask, axis=0)
        return self.__class__(genotypes=genotypes, classification=self.classification)

    def keep_only_biallelic(self):
        allele_count = self.count_alleles_per_variant()
        mask = (allele_count.max_allele() == 1) & (allele_count[:, :2].min(axis=1) > 1)
        genotypes = self.genotypes.compress(mask, axis=0)
        return self.__class__(genotypes=genotypes, classification=self.classification)

    def to_genotype_012_matrix(self):
        genotypes = pandas.DataFrame(self.genotypes.to_n_alt())
        if self.indiviual_names is not None:
            genotypes.index = self.indiviual_names
        return genotypes

    def ld_prune(self, window_size=100, step=20, threshold=0.1, n_iter=1):
        genotypes = self.to_genotype_012_matrix()
        for idx in range(n_iter):
            loc_unlinked = allel.locate_unlinked(
                genotypes, size=window_size, step=step, threshold=threshold
            )
            n = numpy.count_nonzero(loc_unlinked)
            n_remove = genotypes.shape[0] - n
            # print("iteration", i + 1, "retaining", n, "removing", n_remove, "variants")
        genotypes = self.genotypes.compress(loc_unlinked, axis=0)
        return self.__class__(genotypes=genotypes, classification=self.classification)


def calc_ld_rogers_huff_r(genotype_012_matrix):
    return allel.rogers_huff_r(genotype_012_matrix)
