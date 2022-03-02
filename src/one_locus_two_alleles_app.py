from faulthandler import disable
import math

import ipywidgets as widget
from IPython.display import clear_output
from matplotlib import pyplot as plt

from one_locus_two_alleles_simulations import (
    simulate_one_locus_two_alleles_one_pop,
    INF,
)
from plot import GENTOTYPIC_FREQS_LABELS, COLORS


class OneLociTwoAllelesSimulationApp(widget.VBox):
    def __init__(
        self,
        allow_selection=True,
        allow_mutation=True,
        allow_selfing=True,
        asssume_hw=False,
        allow_several_populations=False,
        show_genotypic_freqs_plot=True,
    ):
        widget.VBox.__init__(self, _dom_classes=["widget-interact"])

        children_widgets = []

        if asssume_hw and allow_selfing:
            raise ValueError("assume_hw and allow_selfing options are not compatible")
        if allow_several_populations and show_genotypic_freqs_plot:
            raise ValueError(
                "show_genotypic_freqs_plot and allow_several_populations options are not compatible"
            )

        self.allow_selection = allow_selection
        self.allow_mutation = allow_mutation
        self.allow_selfing = allow_selfing
        self.assume_hw = asssume_hw
        self.allow_several_populations = allow_several_populations
        self.show_genotypic_freqs_plot = show_genotypic_freqs_plot

        self.setup_ui(children_widgets)

        self.clear_output = True
        self.output = widget.Output()
        children_widgets.append(self.output)

        self.run_button = widget.Button(description="Run")
        children_widgets.append(self.run_button)

        self.children = children_widgets

        self.run_button.on_click(self.update)
        self.show_initial_values()

    def _plot_initial_values(
        self, allelic_freq_axes, genotypic_freqs_axes, genotypic_freqs, add_label=False
    ):

        geno_labels = GENTOTYPIC_FREQS_LABELS

        freq_A = genotypic_freqs["freq_AA"] + 0.5 * genotypic_freqs["freq_Aa"]
        allelic_freq_axes.plot([1], [freq_A], marker="o", color=COLORS["freq_A"])

        if self.show_genotypic_freqs_plot:
            for freq_type in ["freq_AA", "freq_Aa", "freq_aa"]:
                kwargs = {"marker": "o", "color": COLORS[freq_type]}
                if add_label:
                    kwargs["label"] = geno_labels[freq_type]
                genotypic_freqs_axes.plot([1], [genotypic_freqs[freq_type]], **kwargs)

    def show_initial_values(self):

        widget.interaction.show_inline_matplotlib_plots()
        with self.output:
            clear_output(wait=True)
            _, axess = self._get_matplotlib_axess()

            allelic_freq_axes = axess[0]
            genotypic_freqs_axes = axess[1]
            sim_params = self._get_ui_simulation_parameters()
            sim_params = self._get_simulation_parameters(**sim_params)

            self._plot_initial_values(
                axess[0],
                axess[1],
                genotypic_freqs={
                    "freq_AA": sim_params["freq_AA"],
                    "freq_Aa": sim_params["freq_Aa"],
                    "freq_aa": sim_params["freq_aa"],
                },
                add_label=True,
            )
            if genotypic_freqs_axes is not None:
                genotypic_freqs_axes.legend()

            self._set_xy_plot_lims(
                allelic_freq_axes,
                genotypic_freqs_axes,
                sim_params["num_generations"],
            )

            widget.interaction.show_inline_matplotlib_plots()

    def _set_xy_plot_lims(
        self, allelic_freq_axes, genotypic_freqs_axes, num_generations
    ):
        allelic_freq_axes.set_ylim((0, 1))
        if genotypic_freqs_axes is None:
            allelic_freq_axes.set_xlim((1, num_generations))
        else:
            genotypic_freqs_axes.set_ylim((0, 1))
            genotypic_freqs_axes.set_xlim((1, num_generations))

    def update(self, *_):
        widget.interaction.show_inline_matplotlib_plots()
        with self.output:
            if self.clear_output:
                clear_output(wait=True)
            kwargs = self._get_ui_simulation_parameters()
            self.result = self.simulate(**kwargs)
            widget.interaction.show_inline_matplotlib_plots()

    def setup_ui(self, children_widgets):
        box_border_style = "solid 1px #cccccc"

        if self.assume_hw:
            label = widget.Label("Freq. A")
            self.freq_A_slider = widget.FloatSlider(
                min=0,
                max=1,
                value=0.5,
            )
            self.freq_A_slider.observe(self.update_freq_sliders)
            freqs_box = widget.HBox(
                [label, self.freq_A_slider],
                layout=widget.Layout(border=box_border_style),
            )
        else:
            self.AA_fraction_slider = widget.FloatSlider(
                min=0, max=1, value=0, description="AA"
            )
            self.Aa_fraction_slider = widget.FloatSlider(
                min=0, max=1, value=1, description="Aa"
            )
            self.aa_fraction_slider = widget.FloatSlider(
                min=0, max=1, value=0, description="aa"
            )
            self.AA_fraction_slider.observe(self.update_freq_sliders)
            self.Aa_fraction_slider.observe(self.update_freq_sliders)
            self.aa_fraction_slider.observe(self.update_freq_sliders)
            fraction_vbox = widget.VBox(
                [
                    widget.Label("Genotypic proportions"),
                    self.AA_fraction_slider,
                    self.Aa_fraction_slider,
                    self.aa_fraction_slider,
                ]
            )

            self.freq_AA_slider = widget.FloatSlider(
                min=0, max=1, value=0, description="AA", disabled=True
            )
            self.freq_Aa_slider = widget.FloatSlider(
                min=0, max=1, value=1, description="Aa", disabled=True
            )
            self.freq_aa_slider = widget.FloatSlider(
                min=0, max=1, value=0, description="aa", disabled=True
            )
            genotypic_freqs_vbox = widget.VBox(
                [
                    widget.Label(
                        "Initial genotypic freqs (Calculated from the proportions)"
                    ),
                    self.freq_AA_slider,
                    self.freq_Aa_slider,
                    self.freq_aa_slider,
                ]
            )
            freqs_box = widget.HBox(
                [fraction_vbox, genotypic_freqs_vbox],
                layout=widget.Layout(border=box_border_style),
            )

        label = widget.Label("Pop. size")
        self.pop_size_slider = widget.IntSlider(min=1, max=1000, value=100)
        self.pop_inf_checkbox = widget.Checkbox(description="Pop is infinite")
        pop_size_box = widget.HBox(
            [label, self.pop_size_slider, self.pop_inf_checkbox],
            layout=widget.Layout(border=box_border_style),
        )

        label = widget.Label("Num. generations")
        self.num_generations_slider = widget.IntSlider(min=10, max=1000, value=100)
        self.num_generations_slider.observe(self.update_freq_sliders)
        num_generations_box = widget.HBox(
            [label, self.num_generations_slider],
            layout=widget.Layout(border=box_border_style),
        )

        if self.allow_several_populations:
            label = widget.Label("Num. populations")
            self.num_populations_slider = widget.IntSlider(min=1, max=20, value=1)
            num_populations_box = widget.HBox(
                [label, self.num_populations_slider],
                layout=widget.Layout(border=box_border_style),
            )

        if self.allow_selection:
            label = widget.Label("Fitness")
            self.fitness_AA_slider = widget.FloatSlider(
                min=0, max=1, value=1, description="WAA"
            )
            self.fitness_Aa_slider = widget.FloatSlider(
                min=0, max=1, value=1, description="WAa"
            )
            self.fitness_aa_slider = widget.FloatSlider(
                min=0, max=1, value=1, description="Waa"
            )
            fitness_box = widget.VBox(
                [
                    label,
                    self.fitness_AA_slider,
                    self.fitness_Aa_slider,
                    self.fitness_aa_slider,
                ]
            )
        if self.allow_mutation:
            label = widget.Label("Mutation rates")
            self.mut_A2a_slider = widget.FloatSlider(
                min=0, max=0.1, value=0, description="A -> a", step=0.01
            )
            self.mut_a2A_slider = widget.FloatSlider(
                min=0, max=0.1, value=0, description="a -> A", step=0.01
            )
            mut_box = widget.VBox([label, self.mut_A2a_slider, self.mut_a2A_slider])
        if self.allow_selfing:
            label = widget.Label("Selfing rate")
            self.self_rate_slider = widget.FloatSlider(
                min=0, max=1, value=0, description="selfing"
            )

        accordion_tabs = []
        if self.allow_selection:
            accordion_tabs.append(("Selection", fitness_box))
        if self.allow_mutation:
            accordion_tabs.append(("Mutation", mut_box))
        if self.allow_selfing:
            accordion_tabs.append(("Selfing", self.self_rate_slider))

        if accordion_tabs:
            accordions = []
            for idx, tab_info in enumerate(accordion_tabs):
                accordion = widget.Accordion(
                    children=[tab_info[1]], selected_index=None
                )
                accordion.set_title(0, tab_info[0])
                accordions.append(accordion)
            accordion_box = widget.VBox(accordions)

        self.exp_het_text = widget.FloatText(disable=True)
        if self.allow_several_populations:
            text = "Mean final Exp. Het."
        else:
            text = "Final Exp. Het."
        exp_het_box = widget.HBox([widget.Label(text), self.exp_het_text])
        self.freq_A_text = widget.FloatText(disable=True)
        if self.allow_several_populations:
            text = "Mean final freq. A"
        else:
            text = "Final freq. A"
        freq_A_box = widget.HBox([widget.Label(text), self.freq_A_text])
        result_box = widget.HBox([exp_het_box, freq_A_box])

        children_widgets.append(freqs_box)
        children_widgets.append(pop_size_box)
        children_widgets.append(num_generations_box)
        if self.allow_several_populations:
            children_widgets.append(num_populations_box)
        if accordion_tabs:
            children_widgets.append(accordion_box)
        children_widgets.append(result_box)

    def update_freq_sliders(self, change):
        if not self.assume_hw:
            fraction_AA = self.AA_fraction_slider.value
            fraction_Aa = self.Aa_fraction_slider.value
            fraction_aa = self.aa_fraction_slider.value

            frac_sum = fraction_AA + fraction_Aa + fraction_aa
            if math.isclose(frac_sum, 0):
                self.AA_fraction_slider.value = 1
                self.Aa_fraction_slider.value = 1
                self.aa_fraction_slider.value = 1
                fraction_AA = self.AA_fraction_slider.value
                fraction_Aa = self.Aa_fraction_slider.value
                fraction_aa = self.aa_fraction_slider.value

            frac_sum = fraction_AA + fraction_Aa + fraction_aa
            freq_AA = fraction_AA / frac_sum
            freq_Aa = fraction_Aa / frac_sum
            freq_aa = fraction_aa / frac_sum

            self.freq_AA_slider.value = freq_AA
            self.freq_Aa_slider.value = freq_Aa
            self.freq_aa_slider.value = freq_aa

        self.show_initial_values()

    def _get_ui_simulation_parameters(self):
        kwargs = {}
        if self.assume_hw:
            freq_A = self.freq_A_slider.value
            freq_AA = freq_A**2
            freq_Aa = 2 * freq_A * (1 - freq_A)
            freq_aa = 1 - freq_AA - freq_Aa
        else:
            freq_AA = self.freq_AA_slider.value
            freq_Aa = self.freq_Aa_slider.value
            freq_aa = self.freq_aa_slider.value
        kwargs["freq_AA"] = freq_AA
        kwargs["freq_Aa"] = freq_Aa
        kwargs["freq_aa"] = freq_aa

        pop_size = self.pop_size_slider.value
        if self.pop_inf_checkbox.value:
            pop_size = INF
        kwargs["pop_size"] = pop_size

        if self.allow_selection:
            wAA = self.fitness_AA_slider.value
            wAa = self.fitness_Aa_slider.value
            waa = self.fitness_aa_slider.value
            if math.isclose(wAA + wAa + waa, 3):
                wAA = 1
                wAa = 1
                waa = 1
            kwargs["wAA"] = wAA
            kwargs["wAa"] = wAa
            kwargs["waa"] = waa
        if self.allow_mutation:
            mut_a2A = self.mut_a2A_slider.value
            mut_A2a = self.mut_A2a_slider.value
            if math.isclose(mut_a2A + mut_A2a, 0):
                mut_a2A = 0
                mut_A2a = 0
            kwargs["mut_a2A"] = mut_a2A
            kwargs["mut_A2a"] = mut_A2a
        if self.allow_selfing:
            kwargs["selfing_rate"] = self.self_rate_slider.value

        kwargs["num_generations"] = self.num_generations_slider.value

        if self.allow_several_populations:
            kwargs["num_populations"] = self.num_populations_slider.value

        return kwargs

    def _get_matplotlib_axess(self):
        if self.show_genotypic_freqs_plot:
            fig, axess = plt.subplots(nrows=2, sharex=True, figsize=(8, 8))
        else:
            fig, axes = plt.subplots(figsize=(8, 4))
            axess = (axes, None)

        return fig, axess

    def _get_simulation_parameters(self, **kwargs):
        sim_kwargs = {}
        sim_kwargs["freq_AA"] = kwargs["freq_AA"]
        sim_kwargs["freq_Aa"] = kwargs["freq_Aa"]
        sim_kwargs["freq_aa"] = kwargs["freq_aa"]
        sim_kwargs["pop_size"] = kwargs["pop_size"]
        sim_kwargs["num_generations"] = kwargs["num_generations"]
        if self.allow_selection:
            sim_kwargs["w11"] = kwargs["wAA"]
            sim_kwargs["w12"] = kwargs["wAa"]
            sim_kwargs["w22"] = kwargs["waa"]
        if self.allow_mutation:
            sim_kwargs["A2a"] = kwargs["mut_a2A"]
            sim_kwargs["a2A"] = kwargs["mut_a2A"]
        if self.allow_selfing:
            sim_kwargs["selfing_rate"] = kwargs["selfing_rate"]
        if self.allow_several_populations:
            sim_kwargs["num_populations"] = kwargs["num_populations"]

        return sim_kwargs

    def simulate(self, **kwargs):
        _, axess = self._get_matplotlib_axess()

        sim_kwargs = {}
        sim_kwargs["allelic_freq_axes"] = axess[0]
        sim_kwargs["genotypic_freqs_axes"] = axess[1]
        for key, value in self._get_simulation_parameters(**kwargs).items():
            sim_kwargs[key] = value

        self._plot_initial_values(
            axess[0],
            axess[1],
            genotypic_freqs={
                "freq_AA": sim_kwargs["freq_AA"],
                "freq_Aa": sim_kwargs["freq_Aa"],
                "freq_aa": sim_kwargs["freq_aa"],
            },
        )

        res = simulate_one_locus_two_alleles_one_pop(**sim_kwargs)
        mean_exp_het = res["exp_het_logger"].values_per_generation.iloc[-1, :].mean()
        self.exp_het_text.value = mean_exp_het
        mean_freq_A = (
            res["allelic_freqs_logger"].values_per_generation.iloc[-1, :].mean()
        )
        self.freq_A_text.value = mean_freq_A
        self._set_xy_plot_lims(axess[0], axess[1], sim_kwargs["num_generations"])


def get_app(app_name):
    if app_name == "allelic_freqs_hw_assumed_app":
        return OneLociTwoAllelesSimulationApp(
            allow_selection=False,
            allow_mutation=False,
            allow_selfing=False,
            asssume_hw=True,
            allow_several_populations=True,
            show_genotypic_freqs_plot=False,
        )
    elif app_name == "simple_one_pop_app":
        return OneLociTwoAllelesSimulationApp(
            allow_selection=False,
            allow_mutation=False,
            allow_selfing=False,
            asssume_hw=False,
            allow_several_populations=False,
            show_genotypic_freqs_plot=True,
        )
    elif app_name == "full_one_pop_app":
        return OneLociTwoAllelesSimulationApp(
            allow_selection=True,
            allow_mutation=True,
            allow_selfing=True,
            asssume_hw=False,
            allow_several_populations=False,
            show_genotypic_freqs_plot=True,
        )
    else:
        raise ValueError(f"Unknown app: {app_name}")
