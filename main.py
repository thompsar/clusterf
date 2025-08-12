from clusterf.app import ClusterFApp
import panel as pn

pn.extension()

app = ClusterFApp()
dashboard = app.serve()  # returns a layout

pn.panel(dashboard).servable()
