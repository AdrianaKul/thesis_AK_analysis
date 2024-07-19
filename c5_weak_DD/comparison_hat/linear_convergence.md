```python
# import settings and functions
%run ./../../imports.ipynb

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
analytical_solution_tag = "-ana_mexi_hat"
generate_config = generateConfig_squareMexiHat
generate_mesh = generateMesh_squareMexiHat

# analytical_solution_tag = "-ana_square_top"
# generate_config = generateConfig_squareTop
# generate_mesh = generateMesh_squareTop

# analytical_solution_tag = "-ana_square_sincos"
# generate_config = generateConfig_squareSinCos
# generate_mesh = generateMesh_squareSinCos
```

## Analysis setup


```python
# which executable?

exe = hdiv_data_driven_diffusion_snes
sumanalys = "sumanalys.csv"
ana_name = "ana_square_mexi_mixed_order"

ana_compare_exe = [hdiv_data_driven_diffusion_snes, hdiv_diffusion, classic_diffusion, data_driven_diffusion_snes]
ana_compare_name = ["ana_square_mexi_dd_weak", "ana_square_mexi_mixed", "ana_square_mexi_classic", "ana_square_mexi_dd"]
# ana_compare_name = ["ana_square_mexi_mixed"]
ana_compare_sum = ["sumanalys.csv", "sumanalys.csv", "FEM_errors.csv", "sumanalys.csv"]

# Convergence analysis parameters
order_list = [1, 2, 3] # approximation order p
elem_size_list = [0.1, 0.05, 0.02] # element size h
# order_list = [2, 3] # approximation order p
# elem_size_list = [0.5, 0.2, 0.1] # element size h
params.triangle_mesh = True
params.nproc = 1 # number of processors
jumps = ""
if params.nproc == 1:
    jumps = "-get_jumps"
# jumps = "-get_jumps"

use_line = "-use_line"
# use_line = ""

run_test = True
run_analysis = True
run_refinement_analysis = True
run_refinement_mesh_analysis = True
run_refinement_hp_analysis = True

run_test = False
run_analysis = False
run_refinement_analysis = False
run_refinement_mesh_analysis = False
# run_refinement_hp_analysis = False

naming = ["order", "gaussnum", "iterations","volume", "datanum","rmsPoiErr", "errorEstimator",
          "L2norm", "H1seminorm","fluxErr", "orderRefinementCounter", "errorIndicatorGrad", "errorIndicatorDiv", "jumpL2", "jumpHdiv", "eleNum"]
# naming = ["order", "gaussnum", "iterations","volume", "datanum","rmsPoiErr", "errorEstimator",
#           "L2norm", "H1seminorm","fluxErr", "orderRefinementCounter"]

error_name_list = ["L2norm", "H1seminorm", "fluxErr"]
error_label_list = [(r'Global error $L^2$-norm'),
               (r'Global error $H^1$-seminorm'), (r'Global Flux error')]
```


```python
params.conductivity = 1.0 # linear conductivity
params.element_size = elem_size_list[0] # element size in the regular mesh
params.order = 1 # approximation order for displacements

# params.triangle_mesh = False # use triangular mesh

# Pre-processing parameters
params.mesh_file = "square_mexi"
params.length_x = 1
params.length_y = 1
params.length_z = 0
params.show_mesh = True


# solution parameters
params.log_file = "log" # log file name 

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
    !rm out*
    params.part_file = params.mesh_file + "_" + str(params.nproc) + "p.h5m"
    !{mofem_part} -my_file {params.mesh_file + ".h5m"} -nparts {params.nproc} -output_file {params.part_file} -dim 2 -adj_dim 1
    !mpirun -np {params.nproc} {exe} -file_name {params.part_file} -my_order {params.order} {analytical_solution_tag} {jumps} {use_line}

    !convert.py out*

```


```python
if run_test:
    params.field_part = -1
    params.show_file = "out_ite"
    params.show_field = "T"
    params.show_edges = True
    params.p_cmap = color_temperature
    params.p_save = "c5_T.pdf"
    show_results(params)
```


```python
if run_test:
    params.field_part = 10
    params.show_file = "out_ite"
    params.show_field = "T"
    params.show_edges = True
    params.clim = [0, 12]
    params.p_cmap = color_gradient
    params.p_save = "c5_grad_T.pdf"
    show_results(params)

```


```python
if run_test:
    params.field_part = -1
    params.show_file = "out_ite"
    params.show_field = "G"
    params.show_edges = True
    params.p_cmap = color_gradient
    params.p_save = "c5_G.pdf"
    show_results(params)
```


```python
if run_test:
    params.show_file = "out_error"
    params.show_field = "ERROR_FLUX"
    params.clim = None
    params.show_edges = True
    params.p_cmap = "jet"
    params.p_save = "c5_err_flux.pdf"
    show_results(params)
```


```python
if run_test:
    if jumps:
        params.show_file = "out_error"
        params.show_field = "JUMP_L2"
        params.show_edges = True
        params.p_cmap = "jet"
        params.p_save = "c5_err_ind_jump.pdf"
        show_results(params)
```


```python
if run_test:
    params.show_file = "out_error"
    params.show_field = "ERROR_ESTIMATOR"
    params.show_edges = True
    params.p_cmap = "jet"
    params.p_save = "c5_err_est.pdf"
    show_results(params)
```


```python

if run_test:
    params.show_file = "out_error"
    params.show_field = "ERROR_H1_SEMINORM"
    params.show_edges = True
    params.p_save = "c5_err_H1_.pdf"
    show_results(params)
```


```python
# if run_test:
#     # params.show_file = "out_ori_result"
#     params.show_file = "out_result"
#     params.show_field = "P_reference"
#     params.warp_field_scalar = "P_reference"
#     params.warp_factor = 0.4  # warp factor
#     params.show_edges = True
#     params.p_save = "run_test_p.pdf"
#     show_results(params)
```


```python
if run_test:
    params.show_file = "out_error"
    params.show_field = "ERROR_INDICATOR_DIV"
    params.show_edges = True
    params.warp_field_scalar = ""
    # params.warp_factor = 0.4  # warp factor
    params.p_save = "run_test_err_ind_div.pdf"
    show_results(params)
```


```python
if run_test:
    params.show_file = "out_error"
    params.show_field = "ERROR_INDICATOR_GRAD"
    params.show_edges = True
    params.warp_field_scalar = ""
    # params.warp_factor = 0.4  # warp factor
    params.p_save = "run_test_err_ind_grad.pdf"
    show_results(params)
```

## Comparison between standard and mixed


```python
if run_analysis:    
    for i in range(len(ana_compare_name)):
        !rm {ana_compare_sum[i]}
        !rm ./out_*
        for elem_size in elem_size_list:
            params.element_size = elem_size
            params.show_mesh = False
            generate_mesh(params)
            params.part_file = params.mesh_file + "_" + str(params.nproc) + "p.h5m"
            !{mofem_part} -my_file {params.mesh_file + ".h5m"} -my_nparts {params.nproc} -output_file {params.part_file} -dim 2 -adj_dim 1
            for order in order_list:
                params.order = order
                !mpirun -np {params.nproc} {ana_compare_exe[i]} -file_name {params.part_file} -my_order {params.order} {analytical_solution_tag} {jumps} {use_line}
        !mv {ana_compare_sum[i]} {ana_compare_name[i]}.csv
    

```


```python
# exe = hdiv_diffusion
sumanalys = "sumanalys.csv"
ana_ref_ord_name = "ana_square_mexi_mixed_order"

refinement_style = 1
ref_iter_num = 7
ref_control = 1.0
params.nproc = 1

if run_refinement_analysis:    
    !rm {sumanalys}
    !rm ./out_*
    elem_size =  elem_size_list[0]
    params.element_size = elem_size
    params.show_mesh = True
    generate_mesh(params)
    params.part_file = params.mesh_file + "_" + str(params.nproc) + "p.h5m"
    !{mofem_part} -my_file {params.mesh_file + ".h5m"} -my_nparts {params.nproc} -output_file {params.part_file} -dim 2 -adj_dim 1
    order = order_list[0]
    params.order = order
    !mpirun -np {params.nproc} {exe} -file_name {params.part_file} -my_order {params.order} {analytical_solution_tag} -refinement_style {refinement_style} -ref_iter_num {ref_iter_num} -ref_control {ref_control} {jumps} {use_line}
    !mv {sumanalys} {ana_ref_ord_name}.csv
```


```python
!convert.py out*

if run_test:
    params.show_file = "out_error"
    params.show_field = "ORDER"
    params.show_edges = True
    params.p_cmap = "rainbow"
    # params.p_save = "run_test_err_ind_grad.pdf"
    show_results(params)

    params.show_file = "out_error"
    params.show_field = "ERROR_INDICATOR_GRAD"
    params.show_edges = True
    params.p_cmap = "jet"
    # params.p_save = "run_test_err_ind_grad.pdf"
    show_results(params)
```


```python
# exe = hdiv_diffusion
sumanalys = "sumanalys.csv"
ana_ref_mesh_name = "ana_square_mexi_mixed_mesh"

if run_refinement_mesh_analysis:   
    refinement_style = 2
    ref_iter_num = 5
    ref_control = 1.0
    params.nproc = 1

    !rm {sumanalys}
    !rm ./out_*
    elem_size =  elem_size_list[0]
    params.element_size = elem_size
    params.show_mesh = True
    generate_mesh(params)
    params.part_file = params.mesh_file + "_" + str(params.nproc) + "p.h5m"
    !{mofem_part} -my_file {params.mesh_file + ".h5m"} -my_nparts {params.nproc} -output_file {params.part_file} -dim 2 -adj_dim 1
    order = order_list[0]
    params.order = order
    !mpirun -np {params.nproc} {exe} -file_name {params.part_file} -my_order {params.order} {analytical_solution_tag} -refinement_style {refinement_style} -ref_iter_num {ref_iter_num} -ref_control {ref_control} {jumps} {use_line}
    !mv {sumanalys} {ana_ref_mesh_name}.csv
```


```python
!convert.py out*

if run_refinement_mesh_analysis:
    params.show_file = "out_error"
    params.show_field = "ERROR_INDICATOR_GRAD"
    params.show_edges = True
    params.p_cmap = "jet"
    # params.p_save = "run_test_err_ind_grad.pdf"
    show_results(params)
```


```python

ana_ref_hp_name = "ana_hat_hp"

if run_refinement_hp_analysis:   
    refinement_style = 3
    ref_iter_num = 7
    ref_control = 1.0
    params.nproc = 1

    !rm {sumanalys}
    !rm ./out_*
    elem_size =  elem_size_list[1]
    order = order_list[1]

    params.element_size = elem_size
    params.show_mesh = True
    generate_mesh(params)
    params.part_file = params.mesh_file + "_" + str(params.nproc) + "p.h5m"
    !{mofem_part} -my_file {params.mesh_file + ".h5m"} -my_nparts {params.nproc} -output_file {params.part_file} -dim 2 -adj_dim 1
    
    params.order = order
    !mpirun -np {params.nproc} {exe} -file_name {params.part_file} -my_order {params.order} {analytical_solution_tag} -refinement_style {refinement_style} -ref_iter_num {ref_iter_num} -ref_control {ref_control} {jumps} {use_line}
    !mv {sumanalys} {ana_ref_hp_name}.csv
```


```python
!convert.py out*

if run_refinement_hp_analysis:
    params.show_file = "out_error"
    params.show_field = "ORDER"
    params.show_edges = True
    # params.p_cmap = "rainbow"
    # params.p_save = "run_test_err_ind_grad.pdf"
    show_results(params)
```


```python
if run_refinement_hp_analysis:
    params.show_file = "out_error"
    params.show_field = "ERROR_INDICATOR_GRAD"
    params.show_edges = True
    # params.p_cmap = "jet"
    # params.p_save = "run_test_err_ind_grad.pdf"
    show_results(params)
```


```python
if run_refinement_hp_analysis:
    params.show_file = "out_result"
    params.show_field = "T"
    params.show_edges = True
    # params.p_save = "run_test_err_ind_grad.pdf"
    show_results(params)
```

### Load analysis


```python
print(naming)
```


```python
error_name_list = []
error_label_list = []

error_name_list.append("L2norm")
error_label_list.append(r'Global error $L^2$-norm')
error_name_list.append("H1seminorm")
error_label_list.append(r'Global error $H^1$-seminorm')
error_name_list.append("fluxErr")
error_label_list.append(r'Global Flux error')

error_name_list.append("errorEstimator")
error_label_list.append(r'Global error estimator')
error_name_list.append("errorIndicatorGrad")
error_label_list.append(r'Global error indicator grad')
error_name_list.append("errorIndicatorDiv")
error_label_list.append(r'Global error indicator div')

if jumps:
    error_name_list.append("jumpL2")
    error_label_list.append(r'Global jump L2')
    # error_name_list.append("jumpHdiv")
    # error_label_list.append(r'Global jump Hdiv')


```


```python
filename_prefix = "c5_mixed_hat_"

mixed_DD_ana = Analysis(ana_compare_name[0], naming, order_list, error_name_list, error_label_list, filename_prefix, elem_size_list,  marker='*', linestyle='-', plot_gradients=True, label="Weak DD")
mixed_ana = Analysis(ana_compare_name[1], naming, order_list, error_name_list, error_label_list, filename_prefix, elem_size_list,  marker='x', linestyle='--', plot_gradients=True, label="Mixed")
classic_ana = Analysis(ana_compare_name[2], naming, order_list, error_name_list, error_label_list, filename_prefix, elem_size_list,  marker='o', linestyle=';', plot_gradients=False, label="Standard")
classic_DD_ana = Analysis(ana_compare_name[3], naming, order_list, error_name_list, error_label_list, filename_prefix, elem_size_list,  marker='v', linestyle=':', plot_gradients=True, label="Standard DD")

order_ref_ana = Analysis(ana_ref_ord_name, naming, order_list, error_name_list, error_label_list, filename_prefix, elem_size_list,  marker='*', linestyle=':', plot_gradients=False, label="Order Refinement", color = 'black')

mesh_ref_ana = Analysis(ana_ref_mesh_name, naming, order_list, error_name_list, error_label_list, filename_prefix, elem_size_list,  marker='v', linestyle=':', plot_gradients=False, label="Mesh Refinement", color = 'black')


hp_ref_ana = Analysis(ana_ref_hp_name, naming, order_list, error_name_list, error_label_list, filename_prefix, elem_size_list,  marker='*', linestyle=':', plot_gradients=False, label="HP Refinement", color = 'red')

# ana_ref_ord_name
```

### Plot results


```python
# mixed_ana.plot_both_analyses_by_elem_size([classic_ana])
mixed_DD_ana.plot_both_analyses_by_elem_size([])
# mixed_ana.plot_both_analyses_by_gaussnum([classic_ana], [order_ref_ana, mesh_ref_ana])
mixed_DD_ana.plot_both_analyses_by_gaussnum([])
```


```python
mixed_DD_ana.plot_gradients = False
mixed_DD_ana.legend_fond_size = 12

mixed_DD_ana.filename_prefix = "c5_mixed_hat_mixed_"
mixed_DD_ana.plot_both_analyses_by_gaussnum([mixed_ana])
mixed_DD_ana.plot_both_analyses_by_elem_size([mixed_ana])

mixed_DD_ana.filename_prefix = "c5_mixed_hat_DD_"
mixed_DD_ana.plot_both_analyses_by_gaussnum([classic_DD_ana])
mixed_DD_ana.plot_both_analyses_by_elem_size([classic_DD_ana])

mixed_DD_ana.filename_prefix = "c5_mixed_hat_ref_order_"
mixed_DD_ana.plot_both_analyses_by_gaussnum([], [order_ref_ana])
mixed_DD_ana.filename_prefix = "c5_mixed_hat_ref_mesh_"
mixed_DD_ana.plot_both_analyses_by_gaussnum([], [mesh_ref_ana])
mixed_DD_ana.filename_prefix = "c5_mixed_hat_ref_order_mesh_"
mixed_DD_ana.plot_both_analyses_by_gaussnum([], [order_ref_ana, mesh_ref_ana])

mixed_DD_ana.filename_prefix = "c5_mixed_hat_ref_hp_"
mixed_DD_ana.plot_both_analyses_by_gaussnum([], [order_ref_ana, mesh_ref_ana, hp_ref_ana])
```


```python

mixed_DD_ana.legend_fond_size = 12
# mixed_ana.plot_gradients = True
mixed_DD_ana.filename_prefix = "c5_mixed_hat_compare_"
mixed_DD_ana.plot_both_analyses_by_gaussnum([classic_DD_ana], [order_ref_ana])
```


```python

```


```python

```
