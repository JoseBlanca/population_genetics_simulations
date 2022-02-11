import numpy
import pandas

from matplotlib import pyplot as plt
import seaborn

import ipywidgets as widget
from IPython.display import clear_output

import demesdraw

import drifting_pops_sim
import pca


class ThreePopDriftSimulationApp(widget.VBox):
    def __init__(
        self,
        min_pop_size=100,
        max_pop_size=100_000,
        default_pop_size=10_000,
        min_num_generations_since_split=100,
        max_num_generations_since_split=10000,
        default_num_generations_since_split=2000,
        num_indis_to_sample_per_pop=20,
        seq_length_in_bp=500_000,
        num_time_intervals_to_sample=4,
    ):

        self.seq_length_in_bp = seq_length_in_bp
        self.num_time_intervals_to_sample = num_time_intervals_to_sample
        self.num_indis_to_sample_per_pop = num_indis_to_sample_per_pop

        self.min_pop_size = min_pop_size
        self.max_pop_size = max_pop_size
        self.default_pop_size = default_pop_size

        self.min_num_generations_since_split = min_num_generations_since_split
        self.max_num_generations_since_split = max_num_generations_since_split
        self.default_num_generations_since_split = default_num_generations_since_split

        widget.VBox.__init__(self, _dom_classes=["widget-interact"])

        children_widgets = []

        self.setup_ui(children_widgets)

        self.output = widget.Output()
        children_widgets.append(self.output)

        self.children = children_widgets

    def setup_ui(self, children_widgets):
        box_border_style = "solid 1px #cccccc"

        self.ancestral_pop_size_slider = widget.IntSlider(
            min=self.min_pop_size,
            max=self.max_pop_size,
            value=self.default_pop_size,
            description="Size ancestral pop",
        )
        self.pop1_size_slider = widget.IntSlider(
            min=self.min_pop_size,
            max=self.max_pop_size,
            value=self.default_pop_size,
            description="Size pop 1",
        )
        self.pop2_size_slider = widget.IntSlider(
            min=self.min_pop_size,
            max=self.max_pop_size,
            value=self.default_pop_size,
            description="Size pop 2",
        )
        self.pop3_size_slider = widget.IntSlider(
            min=self.min_pop_size,
            max=int(self.max_pop_size / 5),
            value=self.default_pop_size,
            description="Size pop 3",
        )

        pop_sizes_box = widget.VBox(
            [
                self.ancestral_pop_size_slider,
                self.pop1_size_slider,
                self.pop2_size_slider,
                self.pop3_size_slider,
            ],
            layout=widget.Layout(border=box_border_style),
        )

        self.num_generations_since_split_slider = widget.IntSlider(
            min=self.min_num_generations_since_split,
            max=self.max_num_generations_since_split,
            value=self.default_num_generations_since_split,
        )
        label = widget.Label(value="Num. generations since split")
        generations_box = widget.HBox(
            [label, self.num_generations_since_split_slider],
            layout=widget.Layout(border=box_border_style),
        )

        children_widgets.append(pop_sizes_box)
        children_widgets.append(generations_box)

        self.run_button = widget.Button(description="Run")
        children_widgets.append(self.run_button)
        self.run_button.on_click(self.update)

    def _get_ui_simulation_parameters(self):
        kwargs = {}

        kwargs["num_indis_to_sample_per_pop"] = self.num_indis_to_sample_per_pop
        kwargs["ancestral_pop_size"] = self.ancestral_pop_size_slider.value
        kwargs["drifting_pops_sizes"] = {
            "pop_1": self.pop1_size_slider.value,
            "pop_2": self.pop2_size_slider.value,
            "pop_3": self.pop3_size_slider.value,
        }
        num_generations = self.num_generations_since_split_slider.value
        kwargs["num_generations_ago_when_split_happened"] = num_generations
        kwargs["seq_length_in_bp"] = self.seq_length_in_bp
        sampling_times = list(
            numpy.linspace(
                0, num_generations, num=self.num_time_intervals_to_sample, dtype=int
            )
        )
        sampling_times[-1] = sampling_times[-1] - 1
        kwargs["sampling_times"] = sampling_times
        return kwargs

    def generate_simulation_plots(self, sim_res, sampling_times):

        nucleotide_diversities = {}
        for sampling_time in reversed(sampling_times):
            nucleotide_diversities[
                sampling_time
            ] = sim_res.calculate_nucleotide_diversities_per_pop(
                sampling_time=sampling_time
            )
        nucleotide_diversities = pandas.DataFrame(nucleotide_diversities).T
        nucleotide_diversities.index = -numpy.array(nucleotide_diversities.index)
        fig, axess = plt.subplots(
            nrows=3 + self.num_time_intervals_to_sample, figsize=(8, 35)
        )
        axes = axess[0]
        demesdraw.tubes(sim_res.demography.to_demes(), ax=axes)

        axes = axess[1]
        seaborn.lineplot(data=nucleotide_diversities, ax=axes)
        axes.set_ylabel("Nucleotide diversity")
        axes.set_xlabel("Num. generations ago")
        axes.set_ylim((0, axes.get_ylim()[1]))

        fsts = sim_res.calculate_fsts(sampling_time=0)
        axes = axess[2]
        seaborn.heatmap(fsts, annot=True, ax=axes)
        axes.set_title("Pairwise Fsts")

        axess = axess[3:]
        for idx, sampling_time in enumerate(reversed(sampling_times)):
            genotypes = sim_res.get_genotypes(
                sampling_time=sampling_time
            ).keep_only_biallelic()
            pca_res = pca.do_pca(genotypes)

            sampling_time_str = (
                f"{sampling_time} generations ago" if sampling_time > 0 else "Now"
            )
            if False:
                print(f"PCA process ({sampling_time_str})")
                print(
                    f'Num. input variants: {pca_res["num_variants_before_ld_pruning"]}'
                )
                print(
                    f'Num. variants after LD pruning: {pca_res["num_variants_after_ld_pruning"]}'
                )

            axes = axess[idx]
            axes.set_title(sampling_time_str)
            pca.plot_pca_result(pca_res, axes, classification=genotypes.classification)

    def update(self, *_, **__):

        widget.interaction.show_inline_matplotlib_plots()
        with self.output:
            clear_output(wait=True)
            kwargs = self._get_ui_simulation_parameters()
            sim_res = drifting_pops_sim.simulate_drifting_pops(**kwargs)

            self.result = self.generate_simulation_plots(
                sim_res, kwargs["sampling_times"]
            )
            widget.interaction.show_inline_matplotlib_plots()
