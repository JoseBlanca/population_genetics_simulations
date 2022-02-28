from faulthandler import disable
from matplotlib import pyplot as plt

import ipywidgets as widget
from IPython.display import clear_output

from one_locus_two_alleles_simulator import AllelicFreqs
from balancing_selection_sim import simulate, INF


class App(widget.VBox):
    def __init__(
        self,
        freq_A_pop1=0.5,
        freq_A_pop2=0.5,
        pop_size_range=(100, 10_000),
        pop_size=1000,
        num_generations=100,
        w11_pop1=1,
        w12_pop1=1,
        w22_pop1=1,
        w11_pop2=1,
        w12_pop2=1,
        w22_pop2=1,
        migration_rate_range=(0, 0.2),
        migration_rate_1_to_2=0,
        migration_rate_2_to_1=0,
    ):
        self.allelic_freqs_pop1 = AllelicFreqs(freq_A_pop1)
        self.allelic_freqs_pop2 = AllelicFreqs(freq_A_pop2)
        self.pop_size_range = pop_size_range
        self.pop_size = pop_size
        self.num_generations = num_generations
        self.w11_pop1 = w11_pop1
        self.w12_pop1 = w12_pop1
        self.w22_pop1 = w22_pop1
        self.w11_pop2 = w11_pop2
        self.w12_pop2 = w12_pop2
        self.w22_pop2 = w22_pop2
        self.migration_rate_range = migration_rate_range
        self.migration_rate_1_to_2 = migration_rate_1_to_2
        self.migration_rate_2_to_1 = migration_rate_2_to_1

        widget.VBox.__init__(self, _dom_classes=["widget-interact"])

        children_widgets = []

        self.setup_ui(children_widgets)

        self.output = widget.Output()
        children_widgets.append(self.output)

        self.children = children_widgets

    def setup_ui(self, children_widgets):
        box_border_style = "solid 1px #cccccc"

        labels = [
            widget.Label(""),
            widget.Label("Freq. A"),
            widget.Label("wAA fitness"),
            widget.Label("wAa fitness"),
            widget.Label("waa fitness"),
            widget.Label("Inmigration rate"),
        ]
        labels_box = widget.VBox(labels)
        widget_cols = [labels_box]

        mig_values = {1: self.migration_rate_2_to_1, 2: self.migration_rate_1_to_2}

        for idx in range(1, 3):
            widgets = [widget.Label(f"Pop. {idx}")]
            slider = widget.FloatSlider(
                min=0, max=1, value=getattr(self, f"allelic_freqs_pop{idx}").A
            )
            setattr(self, f"pop{idx}_freq_A_slider", slider)
            widgets.append(slider)

            AA_slider = widget.FloatSlider(min=0, max=1, value=self.w11_pop1)
            setattr(self, f"pop{idx}_fitness_AA_slider", AA_slider)
            widgets.append(AA_slider)
            Aa_slider = widget.FloatSlider(min=0, max=1, value=self.w12_pop1)
            setattr(self, f"pop{idx}_fitness_Aa_slider", Aa_slider)
            widgets.append(Aa_slider)
            aa_slider = widget.FloatSlider(min=0, max=1, value=self.w22_pop1)
            setattr(self, f"pop{idx}_fitness_aa_slider", aa_slider)
            widgets.append(aa_slider)

            mig_rate_range = self.migration_rate_range
            mig_rate_value = mig_values[idx]
            mig_rate_slider = widget.FloatSlider(
                min=mig_rate_range[0],
                max=mig_rate_range[1],
                value=mig_rate_value,
                step=(mig_rate_range[1] - mig_rate_range[0]) / 10,
            )
            setattr(self, f"pop{idx}_migration_rate", mig_rate_slider)
            widgets.append(mig_rate_slider)

            pop_col = widget.VBox(widgets)
            widget_cols.append(pop_col)

        parameter_set_widget_matrix = widget.HBox(
            widget_cols, layout=widget.Layout(border=box_border_style)
        )
        children_widgets.append(parameter_set_widget_matrix)

        pop_size_range = self.pop_size_range
        self.pop_size_slider = widget.IntSlider(
            min=pop_size_range[0], max=pop_size_range[1], value=self.pop_size
        )
        self.pop_is_inf_checkbox = widget.Checkbox(description="is inf.")
        pop_size_box = widget.HBox(
            [
                widget.Label("Pop. sizes"),
                self.pop_size_slider,
                self.pop_is_inf_checkbox,
            ],
            layout=widget.Layout(border=box_border_style),
        )
        children_widgets.append(pop_size_box)

        label = widget.Label("Num. generations")
        self.num_generations_slider = widget.IntSlider(
            min=10, max=1000, value=self.num_generations
        )
        gen_box = widget.HBox(
            [label, self.num_generations_slider],
            layout=widget.Layout(border=box_border_style),
        )
        children_widgets.append(gen_box)

        self.pop1_exp_het = widget.FloatText(disable=True)
        self.pop2_exp_het = widget.FloatText(disable=True)
        self.total_exp_het = widget.FloatText(disable=True)
        exp_het_box = widget.VBox(
            [
                widget.Label("Exp. het."),
                widget.HBox([widget.Label("Pop. 1"), self.pop1_exp_het]),
                widget.HBox([widget.Label("Pop. 2"), self.pop2_exp_het]),
                widget.HBox([widget.Label("Total"), self.total_exp_het]),
            ]
        )
        children_widgets.append(exp_het_box)

        self.run_button = widget.Button(description="Run")
        children_widgets.append(self.run_button)
        self.run_button.on_click(self.update)

    def _calc_exp_het(self, freq_A):
        freq_AA = freq_A**2
        freq_aa = (1 - freq_A) ** 2
        return 1 - freq_AA - freq_aa

    def _write_final_diversity(self, allelic_freqs_logger):
        values = allelic_freqs_logger.values_per_generation
        pop1_freq_A = values.iloc[-1, 0]
        pop2_freq_A = values.iloc[-1, 1]
        freq_A = (pop1_freq_A + pop2_freq_A) / 2
        self.pop1_exp_het.value = self._calc_exp_het(pop1_freq_A)
        self.pop2_exp_het.value = self._calc_exp_het(pop2_freq_A)
        self.total_exp_het.value = self._calc_exp_het(freq_A)

    def _get_ui_simulation_parameters(self):
        kwargs = {}

        for idx in range(1, 3):
            freq_A = getattr(self, f"pop{idx}_freq_A_slider").value
            freq_a = 1 - freq_A
            freq_AA = freq_A**2
            freq_aa = freq_a**2
            kwargs[f"freq_AA_pop{idx}"] = freq_AA
            kwargs[f"freq_aa_pop{idx}"] = freq_aa
            kwargs[f"freq_Aa_pop{idx}"] = 1 - freq_AA - freq_aa

            pop_size = self.pop_size_slider.value
            if self.pop_is_inf_checkbox.value:
                pop_size = INF
            kwargs[f"pop{idx}_size"] = pop_size

            kwargs[f"w11_pop{idx}"] = getattr(self, f"pop{idx}_fitness_AA_slider").value
            kwargs[f"w12_pop{idx}"] = getattr(self, f"pop{idx}_fitness_Aa_slider").value
            kwargs[f"w22_pop{idx}"] = getattr(self, f"pop{idx}_fitness_aa_slider").value

        kwargs["migration_rate_1_to_2"] = getattr(self, "pop2_migration_rate").value
        kwargs["migration_rate_2_to_1"] = getattr(self, "pop1_migration_rate").value

        kwargs["num_generations"] = self.num_generations_slider.value

        return kwargs

    def update(self, *_, **__):

        widget.interaction.show_inline_matplotlib_plots()
        _, axess = plt.subplots(nrows=2, sharex=True, figsize=(8, 8))
        with self.output:
            clear_output(wait=True)
            kwargs = self._get_ui_simulation_parameters()
            kwargs["allelic_freq_axess"] = axess
            res = simulate(**kwargs)
            widget.interaction.show_inline_matplotlib_plots()

            self._write_final_diversity(res["allelic_freqs_logger"])


# hacer una gráfica con la freq A conjunta
