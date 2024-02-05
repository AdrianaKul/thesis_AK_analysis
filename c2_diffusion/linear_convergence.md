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
%run ./../imports.ipynb

```

```python
# Convergence analysis parameters
order_list = [1, 2, 3, 4] # approximation order p
elem_size_list = [0.5, 0.2, 0.1, 0.05] # element size h

ana_name = "ana_square_top"
sumanalys = "FEM_errors.csv"

run_test = True
run_analysis = True

naming = ["order", "gaussnum", "iterations","volume", "datanum","rmsPoiErr", "errorIndicator",
          "L2norm", "H1seminorm","fluxErr", "orderRefinementCounter"]

error_name_list = ["L2norm", "H1seminorm", "fluxErr"]
error_label_list = [(r'Global error $L^2$-norm'),
               (r'Global error $H^1$-seminorm'), (r'Global Flux error')]
```

```python
params.conductivity = 1.0 # linear conductivity
params.element_size = 0.3 # element size in the regular mesh
params.order = 1 # approximation order for displacements

# Pre-processing parameters
params.mesh_file = "square_top"
params.length_x = 1
params.length_y = 1
params.length_z = 0
params.show_mesh = True

# boundary condition configuration
params.config_file = "bc.cfg"

# solution parameters
params.log_file = "log" # log file name 
params.nproc = 1 # number of processors

```

```python

display = Display(backend="xvfb", visible=False, size=(1024, 768))
display.start()

```

```python

# generation of a config file - what attributes should the blocksets have
def generate_config(params):
    # Open the file for writing
    with open(params.config_file, 'w') as f:
        # FLUX_SQUARE_TOP boundary condition
        data = ['[SET_ATTR_FLUX_SQUARE_TOP]', 'number_of_attributes=1', 'user1='+str(params.conductivity), ' ']
        # Use a for loop to write each line of data to the file
        for line in data:
            f.write(line + '\n')
            # print the data as it is written to the file
            print(line)
        # PRESSURE_UNIFORM boundary condition set to 0
        data = ['[SET_ATTR_PRESSURE_UNIFORM_0]', 'number_of_attributes=1', 'user1=0.0', ' ']
        # Use a for loop to write each line of data to the file
        for line in data:
            f.write(line + '\n')
            # print the data as it is written to the file
            print(line)
```

```python
# gmsh creation of a 3D beam + visualisation of it
def generate_mesh_box(params):
    # Initialize gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 3)

    square1 = gmsh.model.occ.add_rectangle(0, 0, 0, params.length_x, params.length_y)

    # Create the relevant Gmsh data structures from Gmsh model.
    gmsh.model.occ.synchronize()

    # # ensuring a structured mesh with required element size 
    N = int(params.length_x / params.element_size) + 1

    for n in range(len(gmsh.model.getEntities(1))):
        gmsh.model.mesh.setTransfiniteCurve(n+1, N,'Progression', 1.0)

    gmsh.model.mesh.setTransfiniteSurface(square1)

    # gmsh.model.addPhysicalGroup(dimention, [number of element], name="name")
    gmsh.model.addPhysicalGroup(1, [3], name="FLUX_SQUARE_TOP")
    gmsh.model.addPhysicalGroup(1, [1,2,4], name="PRESSURE_UNIFORM_0")
    gmsh.model.addPhysicalGroup(2, [square1], name="DOMAIN")

    # generate a 3D mesh
    gmsh.model.mesh.generate(2)
    
    # save as a .med file
    med_file = params.mesh_file + ".med"
    gmsh.write(med_file)
    
    # close gmsh
    gmsh.finalize()
    
    # translate .med file to a format readable by MoFEM and assign values to physical groups
    h5m_file=params.mesh_file + ".h5m"    
    !{read_med} -med_file {med_file} -output_file {h5m_file} -meshsets_config {params.config_file} -dim 2 -adj_dim 1 -log_sl error
    
    # visualise the mesh
    if params.show_mesh:
        vtk_file=params.mesh_file + ".vtk"
        !mbconvert {h5m_file} {vtk_file}

        mesh = pv.read(vtk_file)
        mesh = mesh.shrink(0.98)

        p = pv.Plotter(notebook=True)
        p.add_mesh(mesh, smooth_shading=False)

        p.camera_position = "xy"
        p.show(jupyter_backend='ipygany')
    
    return

```

```python
# Testing mesh generation
if run_test:
    params.show_mesh = True
    generate_config(params)
    generate_mesh_box(params)
```

```python
# Testing running analysis
if run_test:
    params.part_file = params.mesh_file + "_" + str(params.nproc) + "p.h5m"
    !{mofem_part} -my_file {params.mesh_file + ".h5m"} -my_nparts {params.nproc} -output_file {params.part_file} -dim 2 -adj_dim 1
    !{classic_diffusion} -file_name {params.part_file} -my_order {params.order} -ana_square_top

    !convert.py out*

```

```python
if run_test:
    params.show_file = "out_ori_result"
    params.show_field = "P_reference"
    params.show_edges = True
    show_results(params)
```

# Convergence analysis

```python
if run_analysis:
    !rm {sumanalys}
    !rm ./output_files/out*
    for elem_size in elem_size_list:
        params.element_size = elem_size
        params.show_mesh = False
        generate_mesh_box(params)
        params.part_file = params.mesh_file + "_" + str(params.nproc) + "p.h5m"
        !{mofem_part} -my_file {params.mesh_file + ".h5m"} -my_nparts {params.nproc} -output_file {params.part_file} -dim 2 -adj_dim 1
        for order in order_list:
            params.order = order
            !{classic_diffusion} -file_name {params.part_file} -my_order {params.order} -ana_square_top
    !mv {sumanalys} {ana_name}.csv

```

```python
# read analysis results from ana_name
analysis = []
data = pd.read_csv(ana_name+'.csv', header=None,
                   names=naming,  index_col=False)
analysis.append(data)
analysis[0] = data
```

```python
# separate by order
classic = []

for i in order_list:
    classic.append(analysis[0].query('order == ' + str(i)))
```

```python
for (error_name, error_label) in zip(error_name_list,error_label_list):
    fig = plt.figure()
    ax = plt.axes()
    for i in range(len(order_list)):
        ax.plot(classic[i]["gaussnum"], classic[i][str(error_name)], '-x',
                label=('order = ' + str(order_list[i])))
        ax.set_xscale('log')
        ax.legend(loc='best')
        ax.set_yscale('log')
        ax.set_ylabel(error_label)
        ax.set_xlabel("Number of integration points")
        ax.grid(True, ls=':')

    plt.tight_layout()
    plt.savefig('c2_linear_gauss_' + error_name + '.svg')
    plt.savefig('c2_linear_gauss_' + error_name + '.png')
```

```python
for (error_name, error_label) in zip(error_name_list,error_label_list):
    fig = plt.figure()
    ax = plt.axes()
    for i in range(len(order_list)):
        ax.plot(elem_size_list, classic[i][str(error_name)], '-x',
                label=('order = ' + str(order_list[i])))
        ax.set_xscale('log')
        ax.legend(loc='best')
        ax.set_yscale('log')
        ax.set_ylabel(error_label)
        ax.set_xlabel("Element size")
        ax.grid(True, ls=':')

    plt.tight_layout()
    plt.savefig('c2_linear_lenghtEle_' + error_name + '.svg')
    plt.savefig('c2_linear_lenghtEle_' + error_name + '.png')
```

```python

```
