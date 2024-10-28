# ADMESH+
An automatic unstructured simplex mesh generator.

ADMESH+ is a two-dimensional, automatic unstructured mesh generator for shallow water models. Starting with only target minimum and maximum element sizes and points defining the boundary and bathymetry/topography of the domain, the goal of the mesh generator is to automatically produce a high-quality mesh from this minimal set of input. From the geometry provided, properties such as local feature sizes, curvature of the boundary, bathymetric/topographic gradients, and approximate flow characteristics can be extracted, which are then used to determine local element sizes. The result is a high-quality mesh, with the correct amount of refinement where it is needed to resolve all the geometry and flow characteristics of the domain.

## NOTES:
1. ADMESH+ GUI is originally developed with GUI Layout Toolbox, and now is migrating to MATLAB App Designer. As the migration has not completed yet, a few features are not operational in this release. These features include:
   1. **Fix low quality elements with GUI.** The user has the ability to ``fix" low quality elements by manually moving a node after the mesh is generated. Click on the red triangle to define the lower limit of an acceptable mesh quality, and then click the red triangle with a "+" symbol to locate each element below the limit. Once a "bad" element is located, click the blue triangle button to manually adjust the node placement.
   2. **Subdomain mesh feature.**
2. 

## References:
1. Y. Kang and E. J. Kubatko, “An automatic mesh generator for coupled 1D–2D hydrodynamic models,” Geoscientific Model Development, vol. 17, no. 4, pp. 1603–1625, Feb. 2024, doi: 10.5194/gmd-17-1603-2024.
2. C. J. Conroy, E. J. Kubatko, and D. W. West, “ADMESH: An advanced, automatic unstructured mesh generator for shallow water models,” Ocean Dynamics, vol. 62, no. 10–12, pp. 1503–1517, Dec. 2012, doi: 10.1007/s10236-012-0574-0.
3. P.-O. Persson and G. Strang, “A Simple Mesh Generator in MATLAB,” SIAM Rev., vol. 46, no. 2, pp. 329–345, Jan. 2004, doi: 10.1137/S0036144503429121.
