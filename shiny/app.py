from pandas import DataFrame
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from surfactant.surfactant import calculate
import numpy
import pandas

app_ui = ui.page_fluid(
    ui.output_image("junction"),
    ui.input_numeric("wn", "width of the focussing nozzle (wn), [m]", 7*10 ** -5),
    ui.input_numeric("Ln", "length of the focussing nozzle (Ln), [m]", 7*10 ** -5),
    ui.input_numeric("wcont", "Width of the Continuous-phase channel (wcont)", 6*10 ** -5),
    ui.input_numeric("wdisp", "Width of the Dispersed-phase channel (wdisp)", 7*10 ** -5),
    ui.input_numeric("wout", "Width of the outlet channel (wout)", 11*10 ** -5),
    ui.input_numeric("omega", "Surfactant concentration wt% (dissolved in oil)", 0.006, min=0, max=1, step=0.001),
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
        table = calculate(is_ionic=input.is_ionic(), omega=input.omega())
        df.set(DataFrame(table, columns=["Qoil", "Tdrop"]))

app = App(app_ui, server)
