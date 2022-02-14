import numpy
import pandas

from matplotlib import pyplot as plt
import seaborn

import ipywidgets as widget
from IPython.display import clear_output

import demesdraw

import bottleneck_sim
import pca


class BottleneckApp(widget.VBox):
    def __init__(
        self,
        min_pop_size=10_000,
        max_pop_size=500_000,
        default_pop_size=100_000,
        min_bottleneck_pop_size=100,
        max_bottleneck_pop_size=2_000,
        default_bottleneck_pop_size=100,
        min_num_generations_bottleneck=100,
        max_num_generations_bottleneck=2_000,
        default_num_generations_bottleneck=(800, 1000),
        num_indis_to_sample_per_pop=40,
        seq_length_in_bp=500_000,
    ):

        self.seq_length_in_bp = seq_length_in_bp
        self.num_indis_to_sample_per_pop = num_indis_to_sample_per_pop

        self.min_pop_size = min_pop_size
        self.max_pop_size = max_pop_size
        self.default_pop_size = default_pop_size

        self.min_bottleneck_pop_size = min_bottleneck_pop_size
        self.max_bottleneck_pop_size = max_bottleneck_pop_size
        self.default_bottleneck_pop_size = default_bottleneck_pop_size
        self.min_num_generations_bottleneck = min_num_generations_bottleneck
        self.max_num_generations_bottleneck = max_num_generations_bottleneck
        self.default_num_generations_bottleneck = default_num_generations_bottleneck

        widget.VBox.__init__(self, _dom_classes=["widget-interact"])

        children_widgets = []

        self.setup_ui(children_widgets)

        self.output = widget.Output()
        children_widgets.append(self.output)

        self.children = children_widgets

    def setup_ui(self, children_widgets):
        box_border_style = "solid 1px #cccccc"

        self.pop_size_slider = widget.IntSlider(
            min=self.min_pop_size,
            max=self.max_pop_size,
            value=self.default_pop_size,
            description="Pop. size",
        )
        pop_sizes_box = widget.VBox(
            [
                self.pop_size_slider,
            ],
            layout=widget.Layout(border=box_border_style),
        )

        self.bottleneck_pop_size_slider = widget.IntSlider(
            min=self.min_bottleneck_pop_size,
            max=self.max_bottleneck_pop_size,
            value=self.default_bottleneck_pop_size,
        )
        label = widget.Label(value="Pop. size during bottleneck")
        bottleneck_pop_size_box = widget.HBox(
            [label, self.bottleneck_pop_size_slider],
            layout=widget.Layout(border=box_border_style),
        )
        self.botleneck_time_range_slider = widget.IntRangeSlider(
            min=self.min_num_generations_bottleneck,
            max=self.max_num_generations_bottleneck,
            value=self.default_num_generations_bottleneck,
        )
        label = widget.Label(value="Bottleneck start and end")
        generations_box = widget.HBox(
            [label, self.botleneck_time_range_slider],
            layout=widget.Layout(border=box_border_style),
        )
        bottleneck_box = widget.VBox([bottleneck_pop_size_box, generations_box])

        children_widgets.append(pop_sizes_box)
        children_widgets.append(bottleneck_box)

        self.run_button = widget.Button(description="Run")
        children_widgets.append(self.run_button)
        self.run_button.on_click(self.update)

    def _get_ui_simulation_parameters(self):
        kwargs = {}

        kwargs["num_indis_to_sample"] = self.num_indis_to_sample_per_pop
        kwargs["pop_size"] = self.pop_size_slider.value
        kwargs["pop_size_during_bottleneck"] = self.bottleneck_pop_size_slider.value
        kwargs["bottleneck_start_time"] = self.botleneck_time_range_slider.value[0]
        kwargs["bottleneck_end_time"] = self.botleneck_time_range_slider.value[1]

        kwargs["seq_length_in_bp"] = self.seq_length_in_bp

        return kwargs

    def generate_simulation_plots(self, sim_res, sampling_times):

        nucleotide_diversities = {}
        for sampling_time in reversed(sorted(sampling_times)):
            nucleotide_diversities[
                sampling_time
            ] = sim_res.calculate_nucleotide_diversities_per_sample(
                sampling_times=[sampling_time]
            )
        nucleotide_diversities = pandas.DataFrame(nucleotide_diversities).T
        nucleotide_diversities.index = -numpy.array(nucleotide_diversities.index)
        fig, axess = plt.subplots(nrows=4, figsize=(8, 24))
        axes = axess[0]
        demesdraw.tubes(sim_res.demography.to_demes(), ax=axes)

        axes = axess[1]
        seaborn.lineplot(data=nucleotide_diversities, ax=axes)
        axes.set_ylabel("Nucleotide diversity")
        axes.set_xlabel("Num. generations ago")
        axes.set_ylim((0, axes.get_ylim()[1]))

        fsts = sim_res.calculate_fsts(sampling_times=sampling_times, pop_names=["pop"])
        axes = axess[2]
        seaborn.heatmap(fsts, annot=True, ax=axes)
        axes.set_title("Pairwise Fsts")

        sampling_time = 0
        genotypes = sim_res.get_genotypes(
            sampling_times=sampling_times, pop_names=["pop"]
        ).keep_only_biallelic()
        pca_res = pca.do_pca(genotypes)

        sampling_time_str = (
            f"{sampling_time} generations ago" if sampling_time > 0 else "Now"
        )
        axes = axess[3]
        axes.set_title(sampling_time_str)
        pca.plot_pca_result(pca_res, axes, classification=genotypes.classification)

    def update(self, *_, **__):

        widget.interaction.show_inline_matplotlib_plots()
        with self.output:
            clear_output(wait=True)
            kwargs = self._get_ui_simulation_parameters()
            res = bottleneck_sim.simulate(**kwargs)
            sim_res = res["sim_result"]

            self.result = self.generate_simulation_plots(sim_res, res["sampling_times"])
            widget.interaction.show_inline_matplotlib_plots()
