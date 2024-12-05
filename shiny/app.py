from dropsolver import calculate
from dropsolver.progress import Progress
from dropsolver.util import inclusive_range
from matplotlib import pyplot
from pandas import DataFrame
from pathlib import Path
from scipy.interpolate import interp1d
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.ui import HTML
import dropsolver.parameters
import numpy

def create_numeric_input(parameter):
    if 'dimension' in parameter:
        parameter['description'] += ' [{}]'.format(parameter['dimension'])
    if 'html' in parameter:
        parameter['description'] += ', {}'.format(parameter['html'])
    if 'min' not in parameter:
        parameter['min'] = None
    if 'max' not in parameter:
        parameter['max'] = None
    return ui.input_numeric(parameter['parameter'],
                            HTML(parameter['description']),
                            parameter['default_value'],
                            min=parameter['min'],
                            max=parameter['max'])

measurements = []
disperse_phase = []
continuous_phase = []
other = []

for parameter in dropsolver.parameters.parameters():
    if 'default_value' not in parameter or 'display' in parameter:
        pass
    elif 'dimension' in parameter and parameter['dimension'] == 'μm':
        measurements.append(parameter)
    elif parameter['description'].endswith('(disperse phase)'):
        disperse_phase.append(parameter)
    elif parameter['description'].endswith('(continuous phase)'):
        continuous_phase.append(parameter)
    else:
        other.append(parameter)

app_ui = ui.page_auto(
    ui.include_css(Path(__file__).parent / "dropsolver.css"),
    ui.card(
        ui.card_header("Surfactant type"),
        ui.input_radio_buttons(
            "is_ionic",
            "",
            { 1: "Ionic", 0: "Non-ionic" },
        ),
    ),
    ui.card(
        ui.card_header("Disperse phase viscosity"),
        ui.input_radio_buttons(
            "is_disperse_newtonian",
            "",
            { 1: "Newtonian", 0: "Non-Newtonian" },
        ),
    ),
    ui.card(
        ui.card_header("Continuous phase viscosity"),
        ui.input_radio_buttons(
            "is_continuous_newtonian",
            "",
            { 1: "Newtonian", 0: "Non-Newtonian" },
        ),
    ),
    ui.card(
        ui.card_header("Geometry"),
        ui.layout_columns(
            ui.card([create_numeric_input(p) for p in measurements]),
            ui.card({"class": "center"}, ui.output_image("junction")),
        ),
    ),
    ui.card(
        ui.card_header("Disperse phase"),
        ui.layout_columns(
            ui.card(
                create_numeric_input(disperse_phase[0]),
                ui.panel_conditional(
                    "input.is_disperse_newtonian == 0",
                    [create_numeric_input(p) for p in disperse_phase[1:]],
                ),
            ),
            ui.card({"class": "center"}, ui.HTML("""
<p>
    <math>
        <mrow>
            <msub>
                <mi>&eta;</mi>
                <mi>d</mi>
            </msub>
            <mo>=</mo>
            <msubsup>
                <mi>&eta;</mi>
                <mi>&infin;</mi>
                <mi>(d)</mi>
            </msubsup>
            <mo>+</mo>
            <mfrac>
                <mrow>
                    <msub>
                        <mi>K</mi>
                        <mi>d</mi>
                    </msub>
                    <mo>-</mo>
                    <msubsup>
                        <mi>&eta;</mi>
                        <mi>&infin;</mi>
                        <mi>(d)</mi>
                    </msubsup>
                </mrow>
                <mrow>
                    <mn>1</mn>
                    <mo>+</mo>
                    <msup>
                        <mrow>
                            <mo>(</mo>
                            <msub>
                                <mi>B</mi>
                                <mn>1</mn>
                            </msub>
                            <mo>&InvisibleTimes;</mo>
                            <mover>
                                <mi>&gamma;</mi>
                                <mo>.</mo>
                            </mover>
                            <mo>)</mo>
                        </mrow>
                        <mi>p</mi>
                    </msup>
                </mrow>
            </mfrac>
        </mrow>
    </math>
</p>
            """)),
        ),
    ),
    ui.card(
        ui.card_header("Continuous phase"),
        ui.layout_columns(
            ui.card(
                create_numeric_input(continuous_phase[0]),
                ui.panel_conditional(
                    "input.is_continuous_newtonian == 0",
                    [create_numeric_input(p) for p in continuous_phase[1:]],
                ),
            ),
            ui.card({"class": "center"}, ui.HTML("""
<p>
    <math>
        <mrow>
            <msub>
                <mi>&eta;</mi>
                <mi>c</mi>
            </msub>
            <mo>=</mo>
            <msubsup>
                <mi>&eta;</mi>
                <mi>&infin;</mi>
                <mi>(c)</mi>
            </msubsup>
            <mo>+</mo>
            <mfrac>
                <mrow>
                    <msub>
                        <mi>&eta;</mi>
                        <mi>0</mi>
                    </msub>
                    <mo>-</mo>
                    <msubsup>
                        <mi>&eta;</mi>
                        <mi>&infin;</mi>
                        <mi>(c)</mi>
                    </msubsup>
                </mrow>
                <mrow>
                    <mn>1</mn>
                    <mo>+</mo>
                    <msup>
                        <mrow>
                            <mo>(</mo>
                            <msub>
                                <mi>B</mi>
                                <mn>2</mn>
                            </msub>
                            <mo>&InvisibleTimes;</mo>
                            <mover>
                                <mi>&gamma;</mi>
                                <mo>.</mo>
                            </mover>
                            <mo>)</mo>
                        </mrow>
                        <mi>n</mi>
                    </msup>
                </mrow>
            </mfrac>
        </mrow>
    </math>
</p>
            """)),
        ),
    ),
    ui.card([create_numeric_input(p) for p in other]),
    ui.input_action_button("calculate", "Calculate"),
    ui.output_data_frame("result_dataframe"),
    ui.output_plot("result_plot"),
)

def server(input: Inputs, output: Outputs, session: Session):
    df: reactive.Value[pd.DataFrame] = reactive.Value()

    @render.image
    def junction():
        return {"src": "shiny/junction.svg", "width": "400em"}

    @render.data_frame
    def result_dataframe():
        return render.DataGrid(df())

    @render.plot
    def result_plot():
        figure, axes = pyplot.subplots()
        axes.scatter(list(df()["Qoil [μL/hr]"]), list(df()["Vdrop [pL]"]))
        cubic = interp1d(list(df()["Qoil [μL/hr]"]), list(df()["Vdrop [pL]"]), kind="cubic")
        x = numpy.arange(min(list(df()["Qoil [μL/hr]"])), max(list(df()["Qoil [μL/hr]"])), 1)
        axes.plot(x, cubic(x))
        axes.set_xlabel("Qoil [μL/hr]")
        axes.set_ylabel("Vdrop [pL]")
        return figure

    @reactive.Effect
    @reactive.event(input.is_disperse_newtonian)
    def _():
        if bool(int(input.is_disperse_newtonian())):
            # This value has to be set to equal
            ui.update_numeric("Kd", value=input.etaINF1())

    @reactive.Effect
    @reactive.event(input.is_continuous_newtonian)
    def _():
        if bool(int(input.is_continuous_newtonian())):
            # This value has to be set to equal
            ui.update_numeric("EtaZero", value=input.EtaInf())

    @reactive.Effect
    @reactive.event(input.calculate)
    def _():
        npoints = len(inclusive_range(input.QoilStart(), input.QoilEnd(), input.QoilStep()))
        if npoints > 10:
            msg = ui.modal("Web application can process only up to 10 continuous-phase flow rate values.",
                           title="Too many data points",
                           easy_close=True,
                           footer=None)
            ui.modal_show(msg)
        else:
            parameters = dropsolver.parameters.parameters()
            args = { 'is_ionic': bool(int(input.is_ionic())) }
            for parameter in parameters:
                if parameter['parameter'] not in input:
                    continue
                args[parameter['parameter']] = input[parameter['parameter']]()
                if parameter['parameter'].startswith('is_'):
                    args[parameter['parameter']] = bool(args[parameter['parameter']])
                else:
                    args[parameter['parameter']] = float(args[parameter['parameter']])
                if 'SI_multiplier' in parameter:
                    args[parameter['parameter']] *= parameter['SI_multiplier']
            with ui.Progress(min=1, max=npoints) as progress_bar:
                progress_bar.set(message="Preparing")
                progress = Progress(progress_bar, "Calculating")
                table = calculate(**args, reporter=progress)
            df.set(DataFrame(table, columns=["Qoil [μL/hr]", "Vdrop [pL]"]))

app = App(app_ui, server)
