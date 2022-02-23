from cProfile import label
import numpy
import pandas

from matplotlib import pyplot as plt
import seaborn

import ipywidgets as widget
from IPython.display import clear_output

import demesdraw

from mutation_drift_sfs_sim import simulate
import pca


class MutationDriftApp(widget.VBox):
    def __init__(
        self,
        min_pop_size=500,
        max_pop_size=500_000,
        default_pop_size=100_000,
        ref_pop_size=100_000,
        num_indis_to_sample_per_pop=50,
        seq_length_in_bp=500_000,
    ):

        self.seq_length_in_bp = seq_length_in_bp
        self.num_indis_to_sample_per_pop = num_indis_to_sample_per_pop
        self.ref_pop_size = ref_pop_size

        self.min_pop_size = min_pop_size
        self.max_pop_size = max_pop_size
        self.default_pop_size = default_pop_size

        widget.VBox.__init__(self, _dom_classes=["widget-interact"])

        children_widgets = []

        self.setup_ui(children_widgets)

        self.output = widget.Output()
        children_widgets.append(self.output)

        self.children = children_widgets

    def setup_ui(self, children_widgets):
        box_border_style = "solid 1px #cccccc"

        self.ref_pop_size_text = widget.IntText(value=self.ref_pop_size, disabled=True)
        ref_pop_size_box = widget.HBox(
            [widget.Label("Ref. pop. size:"), self.ref_pop_size_text]
        )

        self.pop_size_slider = widget.IntSlider(
            min=self.min_pop_size,
            max=self.max_pop_size,
            value=self.default_pop_size,
            description="Pop. size",
        )
        pop_sizes_box = widget.VBox(
            [
                ref_pop_size_box,
                self.pop_size_slider,
            ],
            layout=widget.Layout(border=box_border_style),
        )

        label_col = widget.VBox(
            [
                widget.Label(""),
                widget.Label("Nucleotide diversity (π):"),
                widget.Label("% polymorphic markers (0.95):"),
                widget.Label("Num. polymorphic markers (0.95):"),
            ]
        )
        self.diversities_res = widget.FloatText(disabled=True)
        self.poly095_res = widget.FloatText(disabled=True)
        self.num_poly095_res = widget.IntText(disabled=True)
        pop_col = widget.VBox(
            [
                widget.Label("Pop."),
                self.diversities_res,
                self.poly095_res,
                self.num_poly095_res,
            ]
        )

        self.ref_diversities_res = widget.FloatText(disabled=True)
        self.ref_poly095_res = widget.FloatText(disabled=True)
        self.ref_num_poly095_res = widget.IntText(disabled=True)
        ref_pop_col = widget.VBox(
            [
                widget.Label("Ref. pop."),
                self.ref_diversities_res,
                self.ref_poly095_res,
                self.ref_num_poly095_res,
            ]
        )
        grid_box = widget.HBox([label_col, pop_col, ref_pop_col])

        results_box = widget.VBox(
            [widget.Label("Result:"), grid_box],
            layout=widget.Layout(border=box_border_style),
        )

        children_widgets.extend([pop_sizes_box, results_box])

        self.run_button = widget.Button(description="Run")
        children_widgets.append(self.run_button)
        self.run_button.on_click(self.update)

    def _get_ui_simulation_parameters(self):
        kwargs = {}

        kwargs["num_indis_to_sample"] = self.num_indis_to_sample_per_pop
        kwargs["pop_size"] = self.pop_size_slider.value
        kwargs["seq_length_in_bp"] = self.seq_length_in_bp

        return kwargs

    def _get_sfs(sefl, sim_res, normalized):
        genotypes = sim_res.get_genotypes(sampling_times=[0]).keep_only_biallelic()
        sfs = genotypes.folded_sfs
        if normalized:
            sfs = sfs / sfs.sum()
        return sfs

    def generate_simulation_plots(self, sim_res, sim_ref_res):

        self.diversities_res.value = (
            sim_res.calculate_nucleotide_diversities_per_sample(sampling_times=[0])[0]
        )
        self.poly095_res.value = sim_res.calculate_poly095_per_sample(
            sampling_times=[0]
        )[0]
        self.num_poly095_res.value = sim_res.calculate_num_poly095_per_sample(
            sampling_times=[0]
        )[0]

        self.ref_diversities_res.value = (
            sim_ref_res.calculate_nucleotide_diversities_per_sample(sampling_times=[0])[
                0
            ]
        )
        self.ref_poly095_res.value = sim_ref_res.calculate_poly095_per_sample(
            sampling_times=[0]
        )[0]
        self.ref_num_poly095_res.value = sim_ref_res.calculate_num_poly095_per_sample(
            sampling_times=[0]
        )[0]

        fig, axess = plt.subplots(nrows=2, figsize=(8, 8))

        for idx, normalized in enumerate([False, True]):
            sfs = self._get_sfs(sim_res, normalized=normalized)
            ref_sfs = self._get_sfs(sim_ref_res, normalized=normalized)
            axes = axess[idx]
            if normalized:
                title = "Allele Frequency Spectrum (normalized)"
            else:
                title = "Allele Frequency Spectrum"
            axes.set_title(title)
            axes.plot(sfs.index, sfs.values, label="pop")
            axes.plot(ref_sfs.index, ref_sfs.values, label="reference")
            axes.legend()

    def update(self, *_, **__):

        widget.interaction.show_inline_matplotlib_plots()
        with self.output:
            clear_output(wait=True)
            kwargs = self._get_ui_simulation_parameters()
            res = simulate(**kwargs)
            sim_res = res["sim_result"]
            ref_res = simulate(
                pop_size=self.ref_pop_size,
                seq_length_in_bp=self.seq_length_in_bp,
                num_indis_to_sample=self.num_indis_to_sample_per_pop,
            )
            sim_ref_res = ref_res["sim_result"]

            self.result = self.generate_simulation_plots(sim_res, sim_ref_res)
            widget.interaction.show_inline_matplotlib_plots()
