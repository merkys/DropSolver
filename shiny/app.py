from pandas import DataFrame
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from surfactant.progress import Progress
from surfactant.surfactant import calculate
import numpy
import pandas
import surfactant.parameters

def create_numeric_input(parameter):
    if 'dimension' in parameter:
        parameter['description'] += ' [{}]'.format(parameter['dimension'])
    if 'min' not in parameter:
        parameter['min'] = None
    if 'max' not in parameter:
        parameter['max'] = None
    return ui.input_numeric(parameter['parameter'],
                            parameter['description'],
                            parameter['default_value'],
                            min=parameter['min'],
                            max=parameter['max'])

numeric_inputs = [create_numeric_input(p) for p in filter(lambda p: 'default_value' in p, surfactant.parameters.parameters())]

app_ui = ui.page_fluid(
    ui.output_image("junction"),
    ui.input_switch("is_ionic", "Ionic", True),
    ui.input_switch("is_newtonian", "Newtonian", False),
    *numeric_inputs,
    ui.input_action_button("calculate", "Calculate"),
    ui.output_data_frame("result_dataframe"),
)

def server(input: Inputs, output: Outputs, session: Session):
    df: reactive.Value[pd.DataFrame] = reactive.Value()

    @render.image
    def junction():
        return {"src": "shiny/junction.svg", "width": "300px"}

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
        npoints = len(numpy.arange(input.QoilStart(), input.QoilEnd(), input.QoilStep()))
        if npoints > 10:
            msg = ui.modal("Web application can process only up to 10 continuous-phase flow rate values.",
                           title="Too many data points",
                           easy_close=True,
                           footer=None)
            ui.modal_show(msg)
        else:
            parameters = surfactant.parameters.parameters()
            args = {}
            for parameter in parameters:
                if parameter['parameter'] not in input:
                    continue
                args[parameter['parameter']] = input[parameter['parameter']]()
            with ui.Progress(min=1, max=npoints) as progress_bar:
                progress_bar.set(message="Preparing")
                progress = Progress(progress_bar, "Calculating")
                table = calculate(**args, reporter=progress)
            df.set(DataFrame(table, columns=["Qoil [μl/hr]", "Tdrop [pl]"]))

app = App(app_ui, server)
