from pandas import DataFrame
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from surfactant.surfactant import calculate
import numpy
import pandas

app_ui = ui.page_fluid(
    ui.input_numeric("omega", "Surfactant concentration wt% (dissolved in oil)", 0.006, min=0, max=1, step=0.001),
    ui.input_switch("is_ionic", "Ionic", True),
    ui.input_action_button("calculate", "Calculate"),
    ui.output_data_frame("result_dataframe"),
)

def server(input: Inputs, output: Outputs, session: Session):
    df: reactive.Value[pd.DataFrame] = reactive.Value()

    @render.data_frame
    def result_dataframe():
        return render.DataGrid(df())

    @reactive.Effect
    @reactive.event(input.calculate)
    def _():
        table = calculate(is_ionic=input.is_ionic(), omega=input.omega())
        df.set(DataFrame(table, columns=["Qoil", "Tdrop"]))

app = App(app_ui, server)
