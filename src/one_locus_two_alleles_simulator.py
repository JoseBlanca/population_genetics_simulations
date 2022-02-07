from operator import ge
from typing import Callable, Union, Iterable
import random
import math
from array import array
from collections import defaultdict, namedtuple

import numpy
import pandas

import matplotlib.pyplot as plt

INF = math.inf

MENDELIAN_SEGREGATIONS = {
    ("AA", "AA"): [(1, 0, 0)],
    ("aa", "aa"): [(0, 0, 1)],
    ("AA", "aa"): [(0, 1, 0)],
    ("aa", "AA"): [(0, 1, 0)],
    ("AA", "Aa"): ((1, 0, 0), (0, 1, 0)),
    ("Aa", "AA"): ((1, 0, 0), (0, 1, 0)),
    ("aa", "Aa"): ((0, 0, 1), (0, 1, 0)),
    ("Aa", "aa"): ((0, 0, 1), (0, 1, 0)),
    ("Aa", "Aa"): (
        (1, 0, 0),
        (0, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
    ),
}


class GenotypicFreqs:
    def __init__(self, freq_AA, freq_Aa, freq_aa=None):
        self.AA = freq_AA
        self.Aa = freq_Aa
        expected_freq_aa = 1 - freq_AA - freq_Aa
        if freq_aa is not None:
            if not math.isclose(expected_freq_aa, freq_aa, rel_tol=0.001):
                raise ValueError(
                    f"freq_aa ({freq_aa}) should be 1 - freq_AA - freq_Aa ({expected_freq_aa})"
                )
        freq_aa = expected_freq_aa

        if not math.isclose(freq_aa + freq_AA + freq_Aa, 1):
            raise ValueError("Genotypic freqs should sum 1")
        self.aa = freq_aa

    @property
    def freqs(self):
        return (self.AA, self.Aa, self.aa)


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


def calc_allelic_freqs(genotypic_freqs):
    return AllelicFreqs(calc_allelic_freq(genotypic_freqs))


def calc_hwe_genotypic_freqs(allelic_freqs):
    freq_AA = allelic_freqs.A**2
    freq_Aa = 2 * allelic_freqs.A * allelic_freqs.a
    return GenotypicFreqs(freq_AA, freq_Aa)


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
        genotypic_freqs: GenotypicFreqs,
        size: Union[int, float] = INF,
        fitness: Union[Fitness, None] = None,
        mut_rates: Union[MutRates, None] = None,
        selfing_rate: float = 0,
    ):
        self.id = id
        self.size = size
        self.genotypic_freqs = genotypic_freqs
        self.selfing_rate = selfing_rate

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

    def evolve_to_next_generation(self, migration_origins=None):
        genotypic_freqs = self.genotypic_freqs

        if migration_origins is None:
            migration_origins = []

        freq_AA = genotypic_freqs.AA
        freq_Aa = genotypic_freqs.Aa

        # migration
        this_pop_contribution = 1 - sum(
            [origin["inmigrant_rate"] for origin in migration_origins]
        )
        if this_pop_contribution < 0:
            raise ValueError(f"Too many inmigrants for pop {self.id}, more than 100%")

        freq_AA = this_pop_contribution * freq_AA + sum(
            [
                origin["from_pop"].genotypic_freqs.AA * origin["inmigrant_rate"]
                for origin in migration_origins
            ]
        )
        freq_Aa = this_pop_contribution * freq_Aa + sum(
            [
                origin["from_pop"].genotypic_freqs.Aa * origin["inmigrant_rate"]
                for origin in migration_origins
            ]
        )
        freq_aa = this_pop_contribution * genotypic_freqs.aa + sum(
            [
                origin["from_pop"].genotypic_freqs.aa * origin["inmigrant_rate"]
                for origin in migration_origins
            ]
        )
        assert math.isclose(freq_AA + freq_Aa + freq_aa, 1)

        # selection
        fitness = self.fitness
        if fitness is not None:
            w_avg = _calc_w_avg(genotypic_freqs, fitness)
            if math.isclose(w_avg, 0):
                msg = "fitness is 0 because the remaining genotypes have a zero fitness, population would die"
                raise RuntimeError(msg)

            freq_AA = freq_AA * fitness.w11 / w_avg
            freq_Aa = freq_Aa * fitness.w12 / w_avg
        else:
            freq_AA = freq_AA
            freq_Aa = freq_Aa

        # mutation
        if self.mut_rates:
            mu = self.mut_rates.A2a
            mu2 = mu**2
            nu = self.mut_rates.a2A
            nu2 = nu**2
            AA0 = freq_AA
            Aa0 = freq_Aa
            aa0 = 1 - freq_AA - freq_Aa

            new_aa = AA0 * mu2 + Aa0 * mu
            new_AA = aa0 * nu2 + Aa0 * nu
            new_Aa = 2 * AA0 * mu + 2 * aa0 * nu
            AA_removed = AA0 * mu2 + 2 * AA0 * mu
            aa_removed = aa0 * nu2 + 2 * aa0 * nu
            Aa_removed = Aa0 * mu + Aa0 * nu

            AA1 = AA0 + new_AA - AA_removed
            Aa1 = Aa0 + new_Aa - Aa_removed
            aa1 = aa0 + new_aa - aa_removed
            assert math.isclose(AA1 + Aa1 + aa1, 1)
            freq_AA = AA1
            freq_Aa = Aa1

        # drift
        selfing_rate = self.selfing_rate
        if self.size is None or math.isinf(self.size):
            # selfed indivuals
            selfed_freq_AA1 = freq_AA + freq_Aa * 0.25
            selfed_freq_Aa1 = freq_Aa * 0.5
            # non selfed individuals
            allelic_freqs = calc_allelic_freqs(GenotypicFreqs(freq_AA, freq_Aa))
            panmix_genotypic_freqs = calc_hwe_genotypic_freqs(allelic_freqs)
            # final freqs
            freq_AA = selfed_freq_AA1 * selfing_rate + panmix_genotypic_freqs.AA * (
                1 - selfing_rate
            )
            freq_Aa = selfed_freq_Aa1 * selfing_rate + panmix_genotypic_freqs.Aa * (
                1 - selfing_rate
            )
        else:
            num_AA = 0
            num_Aa = 0
            num_aa = 0
            for _ in range(self.size):
                parent1 = self._choose_parent(freq_AA, freq_Aa)

                parent2 = None
                if selfing_rate is not None:
                    value = random.uniform(0, 1)
                    if value < selfing_rate:
                        parent2 = parent1
                if parent2 is None:
                    parent2 = self._choose_parent(freq_AA, freq_Aa)

                mendelian_choices = MENDELIAN_SEGREGATIONS[(parent1, parent2)]
                if len(mendelian_choices) == 1:
                    descendants = mendelian_choices[0]
                else:
                    descendants = random.choice(mendelian_choices)
                num_AA += descendants[0]
                num_Aa += descendants[1]
                num_aa += descendants[2]

            total_indis = num_AA + num_Aa + num_aa
            assert total_indis == self.size

            freq_AA = num_AA / total_indis
            freq_Aa = num_Aa / total_indis

        self.genotypic_freqs = GenotypicFreqs(freq_AA, freq_Aa)


def _update_events(demographic_events, num_generation, active_migrations):
    for event in demographic_events:
        event_num_generation = event.get("num_generation")
        if event_num_generation is not None and event_num_generation != num_generation:
            continue

        if event["type"] == "size_change":
            event["pop"].size = event["new_size"]
        elif event["type"] == "migration_start":
            active_migrations[event["id"]] = event
        elif event["type"] == "migration_stop":
            del active_migrations[event["migration_id"]]


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

    for logger in loggers:
        logger(pops, num_generation=1)

    if demographic_events is None:
        demographic_events = []
    active_migrations = {}
    _update_events(demographic_events, 1, active_migrations)

    for num_generation in range(2, num_generations + 1):
        _update_events(demographic_events, num_generation, active_migrations)

        for pop in pops:
            migration_origin_pops = defaultdict(list)
            for migration in active_migrations.values():
                if migration["to_pop"].id == pop.id:
                    migration_origin_pops[pop.id].append(
                        {
                            "from_pop": migration["from_pop"],
                            "inmigrant_rate": migration["inmigrant_rate"],
                        }
                    )
            pop.evolve_to_next_generation(migration_origin_pops[pop.id])

        for logger in loggers:
            logger(pops, num_generation)


class _PerPopLogger:
    def __init__(self):
        self._values_per_generation = None
        self._generations = array("L")

    def __call__(self, pops: Iterable[Population], num_generation: int):
        self._generations.append(num_generation)

        if self._values_per_generation is None:
            values = {pop.id: array("f") for pop in pops}
            self._values_per_generation = values
        else:
            values = self._values_per_generation

        for pop in pops:
            values[pop.id].append(self._calc_value_for_pop(pop))

    @property
    def values_per_generation(self):
        values = self._values_per_generation
        values = pandas.DataFrame(values, index=self._generations)
        return values


class AllelicFreqLogger(_PerPopLogger):
    def _calc_value_for_pop(self, pop):
        return calc_allelic_freq(pop.genotypic_freqs)


class PopSizeLogger(_PerPopLogger):
    def _calc_value_for_pop(self, pop):
        size = pop.size
        if size is None:
            size = math.inf
        return size


GENOTYPIC_FREQS_NAMES = ["freqs_AA", "freqs_Aa", "freqs_aa"]


class GenotypicFreqsLogger:
    def __init__(self):
        self._values_per_generation = None
        self._generations = array("L")

    def __call__(self, pops: Iterable[Population], num_generation: int):
        self._generations.append(num_generation)

        if self._values_per_generation is None:
            values = {}
            for genotypic_freq_name in GENOTYPIC_FREQS_NAMES:
                values[genotypic_freq_name] = {}
                for pop in pops:
                    values[genotypic_freq_name][pop.id] = array("f")
            self._values_per_generation = values
        else:
            values = self._values_per_generation

        for pop in pops:
            for genotypic_freq_name, value in zip(
                GENOTYPIC_FREQS_NAMES, pop.genotypic_freqs.freqs
            ):
                values[genotypic_freq_name][pop.id].append(value)

    @property
    def values_per_generation(self):
        values = self._values_per_generation
        dframes = {}
        for genotypic_freq_name, freqs in values.items():
            dframes[genotypic_freq_name] = pandas.DataFrame(
                freqs, index=self._generations
            )
        return dframes


def plot_freqs(plot_path, allelic_freqs_logger, genotypic_freqs_logger, pop_id):
    figure, axess = plt.subplots(nrows=2, sharex=True)
    allelic_freq = allelic_freqs_logger.values_per_generation[pop_id]

    allelic_freq_axes = axess[0]
    allelic_freq_axes.plot(allelic_freq.index, allelic_freq)
    allelic_freq_axes.set_ylabel("allelic freq A")

    genotypic_freq_axes = axess[1]
    genotypic_freqs = genotypic_freqs_logger.values_per_generation
    for genotypic_freq_name, this_genotypic_freqs in genotypic_freqs.items():
        genotypic_freq_axes.plot(
            this_genotypic_freqs[pop_id].index,
            this_genotypic_freqs[pop_id],
            label=genotypic_freq_name,
        )
    genotypic_freq_axes.legend()

    figure.savefig(plot_path)


if __name__ == "__main__":

    pop1 = Population(
        id="pop1",
        # size=100,
        size=INF,
        genotypic_freqs=GenotypicFreqs(1, 0, 0),
        fitness=Fitness(w11=0.1, w12=0.1, w22=1),
        # mut_rates=MutRates(a2A=0.01, A2a=0.01),
        selfing_coeff=0,
    )
    pop2 = Population(
        id="pop2",
        # size=100,
        size=INF,
        genotypic_freqs=GenotypicFreqs(0, 0, 1),
        # fitness=Fitness(w11=1, w12=1, w22=0.1),
        # mut_rates=MutRates(a2A=0.01, A2a=0.01),
        selfing_coeff=0,
    )
    pop_size_change = {
        "type": "size_change",
        "pop": pop1,
        "new_size": 1000,
        "num_generation": 100,
    }
    migration_start1 = {
        "id": "migration_pop1_to_pop2",
        "type": "migration_start",
        "from_pop": pop2,
        "to_pop": pop1,
        "inmigrant_rate": 0.1,
        "num_generation": 1,
    }
    migration_start2 = {
        "id": "migration_pop2_to_pop1",
        "type": "migration_start",
        "from_pop": pop1,
        "to_pop": pop2,
        "inmigrant_rate": 0.1,
        "num_generation": 1,
    }
    migration_stop = {
        "migration_id": "migration_pop1_to_pop2",
        "num_generation": 20,
        "type": "migration_stop",
    }
    demographic_events = [
        pop_size_change,
        migration_start1,
        migration_start2,
        migration_stop,
    ]
    allelic_freqs_logger = AllelicFreqLogger()
    pop_size_logger = PopSizeLogger()
    genotypic_freqs_logger = GenotypicFreqsLogger()
    simulate(
        pops=[pop1, pop2],
        num_generations=20,
        demographic_events=demographic_events,
        random_seed=42,
        loggers=[allelic_freqs_logger, pop_size_logger, genotypic_freqs_logger],
    )

    plot_freqs("../temp_pop1.png", allelic_freqs_logger, genotypic_freqs_logger, "pop1")
    plot_freqs("../temp_pop2.png", allelic_freqs_logger, genotypic_freqs_logger, "pop2")
