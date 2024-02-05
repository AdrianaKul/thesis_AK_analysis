---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.0
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

```python
import sys

# # uncomment the following lines to install the required packages
# !{sys.executable} -m pip install gmsh
# !{sys.executable} -m pip install scipy
# !{sys.executable} -m pip install matplotlib
# !{sys.executable} -m pip install pandas
# !{sys.executable} -m pip install numpy
# !{sys.executable} -m pip install pyvista
# !{sys.executable} -m pip install ipywidgets
# !{sys.executable} -m pip install pyvirtualdisplay

# print("All packages installed!")
```

```python
# maths
import numpy as np
import pandas as pd

# plotting graphs settings
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [12, 9]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = "Serif"
plt.rcParams['font.size'] = 15

# create mesh with gmsh
import gmsh

# plotting vtks
import pyvista as pv
pv.set_plot_theme("document")
import ipywidgets as widgets
from pyvirtualdisplay import Display

print("Imports finished! :D")
```

```python
# locations of the executables

# directories
um_view = "$HOME/um_view_release" # change this to the location of your um_view_release or um_view with mofem_data_driven_finite_elements installed
tools_dir = um_view + "/bin"
dd_dir = um_view + "/mofem_data_driven_finite_elements"

# tools
read_med = tools_dir + "/read_med"
mofem_part = tools_dir + "/mofem_part"

# create datasets
create_csv_dataset = dd_dir + "/create_csv_dataset"

# classic diffusion
classic_diffusion = dd_dir + "/classic_diffusion"
hdiv_diffusion = dd_dir + "/hdiv_diffusion"

# data driven diffusion
data_driven_diffusion = dd_dir + "/data_driven_diffusion"
data_driven_diffusion_snes = dd_dir + "/data_driven_diffusion_snes"
hdiv_data_driven_diffusion = dd_dir + "/hdiv_data_driven_diffusion"
hdiv_data_driven_diffusion_snes = dd_dir + "/hdiv_data_driven_diffusion_snes"

# nonlinear diffusion
diffusion_nonlinear_graphite = dd_dir + "/diffusion_nonlinear_graphite"
diffusion_nonlinear_hdiv_graphite = dd_dir + "/diffusion_nonlinear_hdiv_graphite"
diffusion_nonlinear_snes = dd_dir + "/diffusion_nonlinear_snes"
```

```python
class AttrDict(dict):
    def __getattr__(self, attr):
        if attr in self:
            return self[attr]
        raise AttributeError(f"'AttrDict' object has no attribute '{attr}'")
    
params = AttrDict() # Attribute dictionary for storing the parameters
```

```python
def show_results(params):
    out_to_vtk = !ls -c1 {params.show_file}*vtk
    last_file=out_to_vtk[0]

    mesh = pv.read(last_file[:-3] + "vtk")
    if params.warp_field:
        mesh = mesh.warp_by_vector(vectors=params.warp_field, factor=params.warp_factor)

    if params.show_edges:
        mesh=mesh.shrink(0.95)
    
    jupyter_backend='ipygany'
    cmap = "turbo"

    p = pv.Plotter(notebook=True)
    p.add_mesh(mesh, scalars=params.show_field, component=None, smooth_shading=True, cmap=cmap)
    p.camera_position = "xy"
    
    p.enable_parallel_projection()
    p.enable_image_style()
    
    p.show(jupyter_backend=jupyter_backend)

params.warp_field = ""    # no warp field
params.warp_factor = 1.0  # warp factor
params.show_edges = True # show edges
```
