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
# import settings and functions
%run ./../imports.ipynb

```

## What mesh?

Copy your choice to the next cell

for SquareTop:
```
analytical_solution_tag = "-ana_square_top"
generate_config = generateConfig_squareTop
generate_mesh = generateMesh_squareTop
```

for SquareSinCos:
```
analytical_solution_tag = "-ana_square_sincos"
generate_config = generateConfig_squareSinCos
generate_mesh = generateMesh_squareSinCos
```

```python
# Change according to instruction above
analytical_solution_tag = "-ana_square_sincos"
generate_config = generateConfig_squareSinCos
generate_mesh = generateMesh_squareSinCos
```

## Analysis setup

```python
# Convergence analysis parameters
order_list = [1, 2, 3, 4] # approximation order p
elem_size_list = [0.2, 0.1, 0.05, 0.02, 0.01] # element size h
# elem_size_list = [0.5, 0.2, 0.1] # element size h

ana_name = "ana_square_top"
sumanalys = "FEM_errors.csv"

run_test = False
run_analysis = False

naming = ["order", "gaussnum", "iterations","volume", "datanum","rmsPoiErr", "errorIndicator",
          "L2norm", "H1seminorm","fluxErr", "orderRefinementCounter"]

error_name_list = ["L2norm", "H1seminorm", "fluxErr"]
error_label_list = [(r'Global error $L^2$-norm'),
               (r'Global error $H^1$-seminorm'), (r'Global Flux error')]
```

```python
params.conductivity = 1.0 # linear conductivity
params.element_size = 0.1 # element size in the regular mesh
params.order = 1 # approximation order for displacements

# params.triangle_mesh = False # use triangular mesh

# Pre-processing parameters
params.mesh_file = "square_top"
params.length_x = 1
params.length_y = 1
params.length_z = 0
params.show_mesh = True


# solution parameters
params.log_file = "log" # log file name 
params.nproc = 1 # number of processors

```

## Run test

```python
# start display for showing results
display = Display(backend="xvfb", visible=False, size=(1024, 768))
display.start()
```

```python
# Testing mesh generation
if run_test:
    params.show_mesh = True
    generate_config(params)
    generate_mesh(params)
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

## Convergence analysis

```python
if run_analysis:
    !rm {sumanalys}
    !rm ./output_files/out*
    for elem_size in elem_size_list:
        params.element_size = elem_size
        params.show_mesh = False
        generate_mesh(params)
        params.part_file = params.mesh_file + "_" + str(params.nproc) + "p.h5m"
        !{mofem_part} -my_file {params.mesh_file + ".h5m"} -my_nparts {params.nproc} -output_file {params.part_file} -dim 2 -adj_dim 1
        for order in order_list:
            params.order = order
            !{classic_diffusion} -file_name {params.part_file} -my_order {params.order} {analytical_solution_tag}
    !mv {sumanalys} {ana_name}.csv

```

## Read and organise results

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

## Plot results

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
    plt.savefig('c2_linear_gauss_' + error_name + '.pdf')
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
    plt.savefig('c2_linear_lenghtEle_' + error_name + '.pdf')
```
