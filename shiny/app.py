from shiny import App, Inputs, Outputs, Session, reactive, ui
from surfactant.surfactant import calculate

app_ui = ui.page_fluid(
    ui.input_numeric("omega", "Surfactant concentration wt% (dissolved in oil)", 0.006, min=0, max=1, step=0.001),
    ui.input_switch("is_ionic", "Ionic", True),
    ui.input_action_button("calculate", "Calculate"),
)

def server(input: Inputs, output: Outputs, session: Session):

    @reactive.Effect
    @reactive.event(input.calculate)
    def _():
        calculate(is_ionic=input.is_ionic(), omega=input.omega())

app = App(app_ui, server)
