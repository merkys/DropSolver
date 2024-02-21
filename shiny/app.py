from shiny import App, Inputs, Outputs, Session, reactive, ui

app_ui = ui.page_fluid(
    ui.input_switch("ionic", "Ionic", True),
    ui.input_action_button("calculate", "Calculate"),
)

def server(input: Inputs, output: Outputs, session: Session):

    @reactive.Effect
    @reactive.event(input.calculate)
    def _():
        pass

app = App(app_ui, server)
