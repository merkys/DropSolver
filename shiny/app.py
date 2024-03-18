from pandas import DataFrame
from shiny import App, Inputs, Outputs, Session, reactive, render, ui
from surfactant.surfactant import calculate
import numpy
import pandas

app_ui = ui.page_fluid(
    ui.output_image("junction"),

    # NON-NEWTONIAN VISCOSITY PARAMETERS
    ui.input_numeric("Kd", "Zero-shear viscosity (disperse phase) [Pa*s]", 0.001, min=0.0005, max=5),
    ui.input_numeric("etaINF1", "Infinite-shear viscosity (disperse phase) [Pa*s]", 0.001, min=0.0002, max=0.05),
    ui.input_numeric("B1", "Consistency (Cross time constant), (disperse phase) [s]", 4.691, min=0.1, max=2),
    ui.input_numeric("p", "Cross rate constant (disperse phase)", 1.0, min=0.1, max=2),
    ui.input_numeric("Kvisc", "Consistency index (continuous phase) [Pa*s]", 0.0014),

    # CHIP GEOMETRIC PARAMETERS
    ui.input_numeric("wn", "Nozzle width (wn) [μm]", 70, min=5, max=200), # * 10 ** -6
    ui.input_numeric("Ln", "Nozzle length (Ln), [μm]", 70, min=5, max=200),
    ui.input_numeric("H", "Channel height [μm]", 80, min=5, max=200),
    ui.input_numeric("wcont", "Width of the Continuous-phase channel (wcont) [μm]", 60, min=5, max=300),
    ui.input_numeric("wdisp", "Width of the Dispersed-phase channel (wdisp) [μm]", 70, min=5, max=300),
    ui.input_numeric("wout", "Width of the outlet channel (wout) [μm]", 110, min=5, max=300),

    ui.input_numeric("omega", "Fraction of dissolved surfactant", 0.006, min=0, max=1, step=0.001),
    ui.input_numeric("sigmaEQ", "Interface tension without any surfactant [mN/m]", 52, min=0.1, max=72), # * 0.001
    ui.input_numeric("Qw", "Dispersed-phase flow rate [μl/hr]", 220 * 2.78 * 10 ** -13), # FIXME: Is this OK?
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
        df.set(DataFrame(table, columns=["Qoil [μl/hr]", "Tdrop [pl]"]))

app = App(app_ui, server)
