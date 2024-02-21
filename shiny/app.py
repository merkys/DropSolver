from shiny import App, Inputs, Outputs, Session, reactive, ui
from surfactant.surfactant import calculate

app_ui = ui.page_fluid(
    ui.input_switch("is_ionic", "Ionic", True),
    ui.input_action_button("calculate", "Calculate"),
)

def server(input: Inputs, output: Outputs, session: Session):

    @reactive.Effect
    @reactive.event(input.calculate)
    def _():
        calculate(is_ionic=input.is_ionic())

app = App(app_ui, server)
