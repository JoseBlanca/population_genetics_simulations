import ipywidgets as widget
from IPython.display import display, clear_output
from matplotlib import pyplot as plt


class OneLociTwoAllelesSimulationApp(widget.VBox):
    def __init__(self):
        widget.VBox.__init__(self, _dom_classes=["widget-interact"])

        children_widgets = []

        self.setup_ui(children_widgets)

        self.clear_output = True
        self.output = widget.Output()
        children_widgets.append(self.output)

        self.run_button = widget.Button(description="Run")
        children_widgets.append(self.run_button)

        self.children = children_widgets

        self.run_button.on_click(self.update)
        self.update()

    def update(self, *_):
        widget.interaction.show_inline_matplotlib_plots()
        with self.output:
            if self.clear_output:
                clear_output(wait=True)
            # for widget in self.kwargs_widgets:
            #    value = widget.get_interact_value()
            #    self.kwargs[widget._kwarg] = value
            kwargs = self.prepare_update_kwargs()
            self.result = self.generate_simulation_plot(**kwargs)
            widget.interaction.show_inline_matplotlib_plots()
            # if self.auto_display and self.result is not None:
            #    display(self.result)

    def setup_ui(self, children_widgets):
        box_border_style = "solid 1px #cccccc"

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

        label = widget.Label("Mutation rates")
        self.mut_A2a_slider = widget.FloatSlider(
            min=0, max=0.1, value=0, description="A -> a", step=0.01
        )
        self.mut_a2A_slider = widget.FloatSlider(
            min=0, max=0.1, value=0, description="a -> A", step=0.01
        )
        mut_box = widget.VBox([label, self.mut_A2a_slider, self.mut_a2A_slider])

        label = widget.Label("Selfing rate")
        self.self_rate_slider = widget.FloatSlider(
            min=0, max=1, value=0, description="selfing"
        )

        accordion_tabs = [
            ("Selection", fitness_box),
            ("Mutation", mut_box),
            ("Selfing", self.self_rate_slider),
        ]

        accordions = []
        for idx, tab_info in enumerate(accordion_tabs):
            accordion = widget.Accordion(children=[tab_info[1]], selected_index=None)
            accordion.set_title(0, tab_info[0])
            accordions.append(accordion)
        accordion_box = widget.VBox(accordions)

        children_widgets.append(freqs_box)
        children_widgets.append(pop_size_box)
        children_widgets.append(accordion_box)

    def update_freq_sliders(self, change):
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

    def prepare_update_kwargs(self):
        kwargs = {}
        kwargs["freq_AA"] = self.freq_AA_slider.value
        kwargs["freq_Aa"] = self.freq_Aa_slider.value
        kwargs["freq_aa"] = self.freq_aa_slider.value

        pop_size = self.pop_size_slider.value
        if self.pop_inf_checkbox.value:
            pop_size = None
        kwargs["pop_size"] = pop_size

        kwargs["W_AA"] = self.fitness_AA_slider.value
        kwargs["W_Aa"] = self.fitness_Aa_slider.value
        kwargs["W_aa"] = self.fitness_aa_slider.value

        kwargs["mut_a2A"] = self.mut_a2A_slider.value
        kwargs["mut_A2a"] = self.mut_A2a_slider.value

        kwargs["selfing_rate"] = self.self_rate_slider.value
        return kwargs

    def generate_simulation_plot(self, **kwargs):
        print(kwargs)
        fig, axes = plt.subplots(figsize=(8, 4))
        axes.plot([kwargs["freq_AA"]] * 10)
