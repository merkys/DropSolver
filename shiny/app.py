from pandas import DataFrame
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
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
    *numeric_inputs,
    ui.input_switch("is_ionic", "Ionic", True),
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
    @reactive.event(input.calculate)
    def _():
        parameters = surfactant.parameters.parameters()
        table = calculate(is_ionic=input.is_ionic(), omega=input.omega())
        df.set(DataFrame(table, columns=["Qoil [μl/hr]", "Tdrop [pl]"]))

app = App(app_ui, server)
