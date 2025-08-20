import param
import panel as pn

pn.extension()


class Carrousel(param.Parameterized):
    # List of SVGs to display
    svgs = param.List([])
    current_index = param.Integer(0)
    svg_viewer = pn.pane.SVG(object=None, sizing_mode="stretch_both", margin=0)
    prev_button = param.Action(
        default=lambda x: x.param.trigger("prev_button"), label="◀️"
    )
    next_button = param.Action(
        default=lambda x: x.param.trigger("next_button"), label="▶️"
    )
    beginning_button = param.Action(
        default=lambda x: x.param.trigger("beginning_button"), label="⏮️"
    )
    end_button = param.Action(
        default=lambda x: x.param.trigger("end_button"), label="⏭️"
    )
    reset_button = param.Action(
        default=lambda x: x.param.trigger("reset_button"), label="Reset"
    )

    def __init__(self, **params):
        super().__init__(**params)
        self.buttons = pn.Row(
            self.param["beginning_button"],
            self.param["prev_button"],
            self.param["next_button"],
            self.param["end_button"],
            align="center",
            sizing_mode="fixed",
            width=150,
            height=35,
            margin=(1, 1),
        )
        self.update_carrousel()

    @param.depends("svgs", watch=True)
    def update_carrousel(self):
        # If the list is empty, hide buttons and show no SVG
        if not self.svgs:
            self.svg_viewer.object = None
            self.buttons.visible = False
        else:
            self.beginning_carrousel()
            # NOTE: below needs to be explicitly called since calling begining_carrousel
            # may not necessarily trigger the update_image if the index was already 0.
            self.update_image()
            if len(self.svgs) == 1:
                self.buttons.visible = False
            else:
                self.buttons.visible = True

    @param.depends("current_index", watch=True)
    def update_image(self):
        svg = self.svgs[self.current_index]
        self.svg_viewer.object = svg

    @param.depends("prev_button", watch=True)
    def prev_image(self):
        if self.current_index > 0:
            self.current_index -= 1

    @param.depends("next_button", watch=True)
    def next_image(self):
        if self.current_index < len(self.svgs) - 1:
            self.current_index += 1

    @param.depends("beginning_button", watch=True)
    def beginning_carrousel(self):
        self.current_index = 0

    @param.depends("end_button", watch=True)
    def end_carrousel(self):
        self.current_index = len(self.svgs) - 1

    def view(self):
        return pn.Column(
            pn.Row(self.svg_viewer, height=650, margin=(0, 0)),
            self.buttons,
            sizing_mode="stretch_both",
            margin=0,
        )


# Example usage:
# carrousel = Carrousel(svgs=['<svg>...</svg>', '<svg>...</svg>'])
# carrousel_view = carrousel.panel()
# carrousel_view.servable()
