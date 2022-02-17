import numpy
import pandas

from matplotlib import pyplot as plt
import seaborn

import ipywidgets as widget
from IPython.display import clear_output

import demesdraw

import founder_effect_sim
import pca


class FounderEffectApp(widget.VBox):
    def __init__(
        self,
        min_pop_size=10_000,
        max_pop_size=500_000,
        default_pop_size=100_000,
        min_num_founders_individuals=50,
        max_num_founders_individuals=4_000,
        default_num_founders_individuals=100,
        min_founding_time=150,
        max_founding_time=500,
        default_founding_time=200,
        min_bottleneck_time=10,
        max_bottleneck_time=100,
        default_bottleneck_time=25,
        num_indis_to_sample_per_pop=50,
        seq_length_in_bp=500_000,
    ):

        if min_founding_time < max_bottleneck_time:
            raise ValueError(
                "min_founding_time should be lower than max_bottleneck_time"
            )

        self.seq_length_in_bp = seq_length_in_bp
        self.num_indis_to_sample_per_pop = num_indis_to_sample_per_pop

        self.min_pop_size = min_pop_size
        self.max_pop_size = max_pop_size
        self.default_pop_size = default_pop_size

        self.min_num_founders_individuals = min_num_founders_individuals
        self.max_num_founders_individuals = max_num_founders_individuals
        self.default_num_founders_individuals = default_num_founders_individuals
        self.min_founding_time = min_founding_time
        self.max_founding_time = max_founding_time
        self.default_founding_time = default_founding_time
        self.min_bottleneck_time = min_bottleneck_time
        self.max_bottleneck_time = max_bottleneck_time
        self.default_bottleneck_time = default_bottleneck_time

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
            description="Pop. sizes",
        )
        pop_sizes_box = widget.VBox(
            [
                self.pop_size_slider,
            ],
            layout=widget.Layout(border=box_border_style),
        )

        self.num_founders_individuals_slider = widget.IntSlider(
            min=self.min_num_founders_individuals,
            max=self.max_num_founders_individuals,
            value=self.default_num_founders_individuals,
        )
        label = widget.Label(value="Num. founding individuals")
        bottleneck_pop_size_box = widget.HBox(
            [label, self.num_founders_individuals_slider],
            layout=widget.Layout(border=box_border_style),
        )
        self.founding_time_slider = widget.IntSlider(
            min=self.min_founding_time,
            max=self.max_founding_time,
            value=self.default_founding_time,
        )
        label = widget.Label(value="Founding time")
        founding_time_box = widget.HBox(
            [label, self.founding_time_slider],
            layout=widget.Layout(border=box_border_style),
        )
        self.botleneck_time_slider = widget.IntSlider(
            min=self.min_bottleneck_time,
            max=self.max_bottleneck_time,
            value=self.default_bottleneck_time,
        )
        label = widget.Label(value="Bottleneck time duration")
        generations_box = widget.HBox(
            [label, self.botleneck_time_slider],
            layout=widget.Layout(border=box_border_style),
        )
        bottleneck_box = widget.VBox(
            [bottleneck_pop_size_box, founding_time_box, generations_box]
        )

        children_widgets.append(pop_sizes_box)
        children_widgets.append(bottleneck_box)

        self.run_button = widget.Button(description="Run")
        children_widgets.append(self.run_button)
        self.run_button.on_click(self.update)

    def _get_ui_simulation_parameters(self):
        kwargs = {}

        kwargs["num_indis_to_sample"] = self.num_indis_to_sample_per_pop
        kwargs["pop_sizes"] = [self.pop_size_slider.value] * 2
        kwargs["num_founder_individuals"] = [self.num_founders_individuals_slider.value]
        kwargs["bottleneck_time"] = self.botleneck_time_slider.value
        kwargs["pop_founding_times"] = [self.founding_time_slider.value]

        kwargs["seq_length_in_bp"] = self.seq_length_in_bp

        return kwargs

    def generate_simulation_plots(self, sim_res, samples):

        nucleotide_diversities = sim_res.calculate_nucleotide_diversities_per_sample(
            sampling_times=0
        )

        poly_095 = sim_res.calculate_poly095_per_sample(sampling_times=0)

        fig, axess = plt.subplots(nrows=5, figsize=(8, 32))
        axes = axess[0]
        demesdraw.tubes(sim_res.demography.to_demes(), ax=axes)

        axes = axess[1]
        seaborn.barplot(
            x=nucleotide_diversities.index, y=nucleotide_diversities.values, ax=axes
        )
        axes.set_ylabel("Nucleotide diversity")

        axes = axess[2]
        seaborn.barplot(x=poly_095.index, y=poly_095.values * 100, ax=axes)
        axes.set_ylabel("% polymorphic markers (95%)")

        samples = list(samples.keys())
        fsts = sim_res.calculate_fsts(samples=samples)
        axes = axess[3]
        seaborn.heatmap(fsts, annot=True, ax=axes)
        axes.set_title("Pairwise Fsts")

        genotypes = sim_res.get_genotypes(samples=samples).keep_only_biallelic()
        pca_res = pca.do_pca(genotypes)

        axes = axess[4]
        pca.plot_pca_result(pca_res, axes, classification=genotypes.classification)

    def update(self, *_, **__):

        widget.interaction.show_inline_matplotlib_plots()
        with self.output:
            clear_output(wait=True)
            kwargs = self._get_ui_simulation_parameters()
            res = founder_effect_sim.simulate(**kwargs)
            sim_res = res["sim_result"]

            self.result = self.generate_simulation_plots(sim_res, res["samples"])
            widget.interaction.show_inline_matplotlib_plots()


# LD
# number of variants
# SFS
