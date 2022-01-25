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


class Fitness:
    def __init__(self, w11, w12, w22):
        self.w11 = w11
        self.w12 = w12
        self.w22 = w22


def calc_allelic_freq(genotypic_freqs):
    return genotypic_freqs.AA + genotypic_freqs.Aa * 0.5


class Population:
    def __init__(
        self,
        id: str,
        size: int,
        genotypic_freqs: GenotypicFreqs,
        fitness: Union[Fitness, None] = None,
        selfing_coeff: float = 0,
    ):
        self.id = id
        self.size = size
        self.genotypic_freqs = genotypic_freqs
        self.selfing_coeff = selfing_coeff

        self.fitness = fitness

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
        genotipic_freqs = self.genotypic_freqs

        fitness = self.fitness
        w11 = fitness.w11
        w12 = fitness.w12
        w22 = fitness.w22
        w_avg = (
            genotipic_freqs.AA * w11
            + genotipic_freqs.Aa * w12
            + genotipic_freqs.aa * w22
        )

        # selection
        freq_AA = genotipic_freqs.AA * w11 / w_avg
        freq_Aa = genotipic_freqs.Aa * w12 / w_avg

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
# mutation
# admixture
# migration

if __name__ == "__main__":

    pop1 = Population(
        id="pop1",
        size=10,
        genotypic_freqs=GenotypicFreqs(0.5, 0, 0.5),
        fitness=Fitness(w11=1, w12=1, w22=1),
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
        num_generations=5,
        demographic_events=demographic_events,
        random_seed=42,
        loggers=[allelic_freqs_logger, pop_size_logger],
    )
    print(allelic_freqs_logger.values_per_generation)
    print(pop_size_logger.values_per_generation)
