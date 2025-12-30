# ADMESH+: An automatic unstructured simplex mesh generator

ADMESH+ is a two-dimensional, automatic unstructured mesh generator for shallow water models. Starting with only target minimum and maximum element sizes and points defining the boundary and bathymetry/topography of the domain, the goal of the mesh generator is to automatically produce a high-quality mesh from this minimal set of input. From the geometry provided, properties such as local feature sizes, curvature of the boundary, bathymetric/topographic gradients, and approximate flow characteristics can be extracted, which are then used to determine local element sizes. The result is a high-quality mesh, with the correct amount of refinement where it is needed to resolve all the geometry and flow characteristics of the domain.

<img width="1392" height="915" alt="ADMESH_GUI" src="https://github.com/user-attachments/assets/de6f7981-4c3d-4a1f-9339-3ea887fbdfd4" />

## Quick Start: Basic Usage Example
1. From `File` on menubar, click `Open File`.
2. Open `examples/1. Basic examples/Baranja_Hill.mat`.
3. (Optional) Change mesh generation parameters such as min/max element size in `Settings` pannel. For details of the parameters, see [[2]](#conroy2012).
4. (Optional) Import elevation data using `Import > Elevation File`. The elevation data can be imported after mesh is generated (after Step 5). In this example, skip this step.
5. Click `Run ADMESH` below the `Settings` pannel to start mesh generation.
6. After mesh generation, the results are shown in the `ADMESH Results` pannel.
7. To identify low quality elements, `Set element quality threshold` and  `Zoom to low quality element` on toolbar can be used.

### Notes
The `.mat` is an ADMESH+ input file. ADMESH+ uses its own input file format to save and load ADMESH+ inputs (e.g., boundaries of mesh, elevation data, and user-specified parameters). Users can use a shapefile (`.shp`) to define boundaries of mesh and start mesh generation from the scratch. The shapefile must be closed lines, where the first attribute line is the external boundary and the other lines are internal boundaries.

## Additional Features
ADMESH+ has additional features over the basic mesh generation. In current version, two features can be found in `Tools`:
* `Extract channels (land)`: Extract channels from given DEM using TopoToolbox. It is designed to generate meshes for watersheds, and it is assumed that the domain do not contain water portion (e.g., lake or ocean). After identifying channels, 1D meshes are generated along the channels first and then 2D meshes on the entire domain are generated.
* `Extract channels (water)`: It is designed to generated meshes for coastal regions, especially with vast small-scale channels. The small-scale channels are identified (based on user-specified minimum channel width) and handled as 1D domain (see [[1]](#kang2024) for details).

## Notes
ADMESH+ GUI was originally developed with GUI Layout Toolbox and is now being migrated to MATLAB App Designer. As the migration has not been completed yet, a few features are not operational in this release (for example, editing nodes with GUI to fix low quality elements). These features will be restored in future updates.

## References
<a id="kang2024">[1]</a> Y. Kang and E. J. Kubatko, “An automatic mesh generator for coupled 1D–2D hydrodynamic models,” Geoscientific Model Development, vol. 17, no. 4, pp. 1603–1625, Feb. 2024, doi: 10.5194/gmd-17-1603-2024.

<a id="conroy2012">[2]</a> C. J. Conroy, E. J. Kubatko, and D. W. West, “ADMESH: An advanced, automatic unstructured mesh generator for shallow water models,” Ocean Dynamics, vol. 62, no. 10–12, pp. 1503–1517, Dec. 2012, doi: 10.1007/s10236-012-0574-0.

<a id="persson2004">[3]</a> P.-O. Persson and G. Strang, “A Simple Mesh Generator in MATLAB,” SIAM Rev., vol. 46, no. 2, pp. 329–345, Jan. 2004, doi: 10.1137/S0036144503429121.
<!--
## Input File
1. ADMESH+ input file (`.mat`): ADMESH+ uses its own input file format to save and load ADMESH+ inputs (e.g., boundaries of mesh, elevation data, and user-specified parameters)
2. Shapefile (`.shp`): A shapefile can be opened to define boundaries of mesh and start mesh generation from the scratch. The shapefile must be closed lines, where the first attribute line is the external boundary and the other lines are internal boundaries.
3. Meshes (`.mat, .14, .grd, .2dm`): A given mesh can be loaded to investigate or modify the mesh.

## Menu
### Menubar
1. File: Open the following files to start mesh generation or investigate/modify given mesh.
   * Shapefile (`.shp`): A shapefile can be opened to define boundaries of mesh and start mesh generation from the scratch. The shapefile must be closed lines, where the first attribute line is the external boundary and the other lines are internal boundaries.
   * ADMESH+ input file (`.mat`): ADMESH+ uses its own input file format to save and load ADMESH+ inputs (e.g., boundaries of mesh, elevation data, and user-specified parameters).
   * Meshes (`.mat, .14, .grd, .2dm`): A given mesh can be loaded to investigate or modify the mesh.
2. Edit
   * Clear: Clear all loaded data
3. Import: Import mesh, NOAA coastline, elevation data, or KML file.
4. Tools: Additional features over the basic mesh generation. In current version, it supports two features:
   * Extract channels (land): Extract channels from given DEM using TopoToolbox. It is designed to generate meshes for watersheds, and it is assumed that the domain do not contain water portion (e.g., lake or ocean). After identifying channels, 1D meshes are generated along the channels first and then 2D meshes on the entire domain are generated.
   * Extract channels (water). It is designed to generated meshes for coastal regions, especially with vast small-scale channels. The small-scale channels are identified (based on user-specified minimum channel width) and handled as 1D domain (see [[1]](#kang2024) for details).

### Toolbar
1. Zoom
2. Reset to original view
3. Pan
4. ~~Edit nodes~~
6. Set element quality threshold: Users can set element quality threshold (between 0 to 1). The elements below this threshold are marked as red. 
7. Zoom to low quality element: Zoom to the element below the threshold. Note that it shows only the first element (not lowest quality) below the threshold.
8. ~~Flag location~~

### Main Figure

### Setting Pannel

### Result Pannel

### Run Button
-->
