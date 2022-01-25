from typing import Callable, Union, Iterable
import random
import math
from array import array
from collections import namedtuple

import numpy


class GenotypicFreqs:
    def __init__(self, freq_AA, freq_Aa, freq_aa=None):
        self.AA = freq_AA
        self.Aa = freq_Aa
        expected_freq_aa = 1 - freq_AA - freq_Aa
        if freq_aa is None:
            freq_aa = expected_freq_aa
        else:
            if not math.isclose(expected_freq_aa, freq_aa):
                raise ValueError("freq_aa should be 1 - freq_AA - freq_Aa")

        if not math.isclose(freq_aa + freq_AA + freq_Aa, 1):
            raise ValueError("Genotypic freqs should sum 1")
        self.aa = freq_aa


class AllelicFreqs:
    def __init__(self, freq_A):
        self.A = freq_A
        freq_a = 1 - freq_A
        self.a = freq_a
        if not math.isclose(freq_A + freq_a, 1):
            raise ValueError("Allelic freqs should sum 1")


Fitness = namedtuple("Fitness", ("w11", "w12", "w22"))

MutRates = namedtuple("MutRates", ("A2a", "a2A"))


def calc_allelic_freq(genotypic_freqs):
    return genotypic_freqs.AA + genotypic_freqs.Aa * 0.5


def _calc_w_avg(genotypic_freqs, fitness):
    if fitness is None:
        raise ValueError("fitness is None")
    w_avg = (
        genotypic_freqs.AA * fitness.w11
        + genotypic_freqs.Aa * fitness.w12
        + genotypic_freqs.aa * fitness.w22
    )
    return w_avg


class Population:
    def __init__(
        self,
        id: str,
        size: int,
        genotypic_freqs: GenotypicFreqs,
        fitness: Union[Fitness, None] = None,
        mut_rates: Union[MutRates, None] = None,
        selfing_coeff: float = 0,
    ):
        self.id = id
        self.size = size
        self.genotypic_freqs = genotypic_freqs
        self.selfing_coeff = selfing_coeff

        if fitness is not None:
            w_avg = _calc_w_avg(genotypic_freqs, fitness)
            if math.isclose(w_avg, 0):
                msg = "fitness is 0 because the remaining genotypes have a zero fitness, population would die"
                raise ValueError(msg)
        self.fitness = fitness

        self.mut_rates = mut_rates

    @property
    def allelic_freqs(self):
        genotypic_freqs = self.genotypic_freqs
        return AllelicFreqs(calc_allelic_freq(genotypic_freqs))

    def _choose_parent(self, freq_AA, freq_Aa):
        value = random.uniform(0, 1)
        if value < freq_AA:
            return "AA"
        elif value >= (freq_AA + freq_Aa):
            return "aa"
        else:
            return "Aa"

    def evolve(self):
        genotypic_freqs = self.genotypic_freqs

        # selection
        fitness = self.fitness
        if fitness is not None:
            w_avg = _calc_w_avg(genotypic_freqs, fitness)
            if math.isclose(w_avg, 0):
                msg = "fitness is 0 because the remaining genotypes have a zero fitness, population would die"
                raise RuntimeError(msg)

            freq_AA = genotypic_freqs.AA * fitness.w11 / w_avg
            freq_Aa = genotypic_freqs.Aa * fitness.w12 / w_avg
        else:
            freq_AA = genotypic_freqs.AA
            freq_Aa = genotypic_freqs.Aa

        # mutation
        if self.mut_rates:
            mu = self.mut_rates.A2a
            mu2 = mu ** 2
            nu = self.mut_rates.a2A
            nu2 = nu ** 2
            AA0 = freq_AA
            Aa0 = freq_Aa
            aa0 = 1 - freq_AA - freq_Aa

            new_aa = AA0 * mu2 + Aa0 * mu
            new_AA = aa0 * nu2 + Aa0 * nu
            new_Aa = AA0 * mu + aa0 * nu
            AA_removed = AA0 * mu2 + AA0 * mu
            aa_removed = aa0 * nu2 + aa0 * nu
            Aa_removed = Aa0 * mu + Aa0 * nu

            AA1 = AA0 + new_AA - AA_removed
            Aa1 = Aa0 + new_Aa - Aa_removed
            aa1 = aa0 + new_aa - aa_removed
            assert math.isclose(AA1 + Aa1 + aa1, 1)
            freq_AA = AA1
            freq_Aa = Aa1

        # drift
        selfing_coeff = self.selfing_coeff
        mendel_het_cross_segregation = ((1, 0, 0), (0, 1, 0), (0, 1, 0), (0, 0, 1))
        num_AA = 0
        num_Aa = 0
        num_aa = 0
        for _ in range(self.size):
            parent1 = self._choose_parent(freq_AA, freq_Aa)

            value = random.uniform(0, 1)
            if value < selfing_coeff:
                parent2 = parent1
            else:
                parent2 = self._choose_parent(freq_AA, freq_Aa)

            if (parent1, parent2) == ("AA", "AA"):
                num_AA += 1
            elif (parent1, parent2) == ("aa", "aa"):
                num_aa += 1
            else:
                result = random.choice(mendel_het_cross_segregation)
                num_AA += result[0]
                num_Aa += result[1]
                num_aa += result[2]
        total_indis = num_AA + num_Aa + num_aa
        assert total_indis == self.size
        freq_AA = num_AA / total_indis
        freq_Aa = num_Aa / total_indis

        self.genotypic_freqs = GenotypicFreqs(freq_AA, freq_Aa)


def simulate(
    pops: list[Population],
    num_generations: int,
    loggers: list[Callable],
    demographic_events: Union[list[dict], None] = None,
    random_seed: Union[int, None] = None,
):
    if random_seed is not None:
        numpy.random.seed(random_seed)
        random.seed(random_seed)

    if demographic_events is None:
        demographic_events = []

    for num_generation in range(1, num_generations + 1):

        for event in demographic_events:
            event_num_generation = event.get("num_generation")
            if (
                event_num_generation is not None
                and event_num_generation != num_generation
            ):
                continue

            if event["type"] == "size_change":
                event["pop"].size = event["new_size"]

        for pop in pops:
            pop.evolve()

        for logger in loggers:
            logger(pops)


class _PerPopLogger:
    def __init__(self):
        self.values_per_generation = None

    def __call__(self, pops: Iterable[Population]):
        if self.values_per_generation is None:
            values = {pop.id: array("f") for pop in pops}
            self.values_per_generation = values
        else:
            values = self.values_per_generation

        for pop in pops:
            values[pop.id].append(self._calc_value_for_pop(pop))


class AllelicFreqLogger(_PerPopLogger):
    def _calc_value_for_pop(self, pop):
        return calc_allelic_freq(pop.genotypic_freqs)


class PopSizeLogger(_PerPopLogger):
    def _calc_value_for_pop(self, pop):
        return pop.size


# TODO
# admixture
# migration

if __name__ == "__main__":

    pop1 = Population(
        id="pop1",
        size=100,
        genotypic_freqs=GenotypicFreqs(1, 0, 0),
        fitness=Fitness(w11=0.1, w12=0.5, w22=1),
        mut_rates=MutRates(a2A=0, A2a=0.1),
    )
    pop_size_change = {
        "type": "size_change",
        "pop": pop1,
        "new_size": 10000,
        "num_generation": 3,
    }
    demographic_events = [pop_size_change]
    allelic_freqs_logger = AllelicFreqLogger()
    pop_size_logger = PopSizeLogger()
    simulate(
        pops=[pop1],
        num_generations=10,
        demographic_events=demographic_events,
        random_seed=42,
        loggers=[allelic_freqs_logger, pop_size_logger],
    )
    print(allelic_freqs_logger.values_per_generation)
    print(pop_size_logger.values_per_generation)
