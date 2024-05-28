from pandas import DataFrame
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from shiny.ui import HTML
from dropsolver import calculate
from dropsolver.progress import Progress
from dropsolver.util import inclusive_range
import pandas
import dropsolver.parameters

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
    ui.card(
        ui.input_switch("is_ionic", "Ionic", True),
        ui.input_switch("is_newtonian", "Newtonian", False),
    ),
    ui.card(
        ui.layout_columns(
            ui.page_auto([create_numeric_input(p) for p in measurements]),
            ui.page_auto(ui.output_image("junction")),
        ),
    ),
    ui.card(
        ui.layout_columns(
            ui.page_auto([create_numeric_input(p) for p in disperse_phase]),
            ui.page_auto(ui.output_image("dispersed_phase")),
        ),
    ),
    ui.card(
        ui.layout_columns(
            ui.page_auto([create_numeric_input(p) for p in continuous_phase]),
            ui.page_auto(ui.output_image("continuous_phase")),
        ),
    ),
    ui.card([create_numeric_input(p) for p in other]),
    ui.input_action_button("calculate", "Calculate"),
    ui.output_data_frame("result_dataframe"),
)

def server(input: Inputs, output: Outputs, session: Session):
    df: reactive.Value[pd.DataFrame] = reactive.Value()

    @render.image
    def junction():
        return {"src": "shiny/junction.svg", "width": "300px"}

    @render.image
    def dispersed_phase():
        return {"src": "shiny/dispersed-phase.svg", "width": "300px"}

    @render.image
    def continuous_phase():
        return {"src": "shiny/continuous-phase.svg", "width": "300px"}

    @render.data_frame
    def result_dataframe():
        return render.DataGrid(df())

    @reactive.Effect
    @reactive.event(input.is_newtonian)
    def _():
        if input.is_newtonian():
            # These two values have to be equal
            ui.update_numeric("Kd", value=input.etaINF1())
            ui.update_numeric("EtaZero", value=input.EtaInf())
        else:
            pass

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
            args = { 'is_ionic': bool(input.is_ionic()) }
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
