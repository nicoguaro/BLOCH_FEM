# -------------------------------------------------------------------- #

    Bloch BC Imposition Routine - DAMIAN-BLOCH
    ------------------------------------------

Compute dispersion relations for a composite material using a single cell
and Bloch Theorem (see [1,3]). The idea is to find the bulk behaviour
for a material made with this periodic microstructure (see [2,4,5]).

# -------------------------------------------------------------------- #

Author: Nicolas Guarin Z.
------

Date: September 19, 2012
-----

# -------------------------------------------------------------------- #

Input:
------
        - "Mass" matrix M (NDOF, NDOF).
        - "Stiffness" matrices K (NDOF, NDOF).
        - Coordinates (DIM, NNODE).
        - Image nodes list, reference nodes list (See [2,6]).
        - kxmin, kymin, kxmax, kymax: minimum/maximum wavenumber in x/y
                                      direction-
        - nkx, nky: number of wavenumber steps for x/y direction.
        - Number of Degrees of Freedom per node.


Output:
-------
        - Dispersion relation list [kx, ky, omega^2].


Pseudocode:
-----------

Start
    Associate node numbering with DOF numbering
    
    For i=1 to nkx
        For j=1 to nky
            Apply Bloch_BC to K --> KR
            Apply Bloch_BC to M --> MR
            Solve generalized eigenproblem [KR]{uR} = omega^2[MR]{uR}
            Store kx, ky, omega^2
        EndFor
    EndFor
    
    Return list [kx, ky, omega^2]
End

# -------------------------------------------------------------------- #

References
----------

[1] Brillouin, León (1946). Wave Propagation in Periodic
Structures. McGraw-Hill, 1946.

[2] Guarín Z., Nicolás (2012). Simulación Numérica de Problemas de 
Propagación de Ondas: Dominios Infinitos y Semi-infinitos. 
MSc in Engineering Thesis, Universidad EAFIT.

[3] Charles Kittel (1996). Introduction to Solid State Physics.
Wiley; 7 edition, 1996.

[4] Ruzzene, Massimo, Scarpa, Fabrizio and Soranna, Francesco
(2003). Wave beaming effects in two-dimensional cellular
structures. Smart Mat. Struct. 12, 363–372.

[5] Srikantha, A., J. Woodhouse y N.A. Flecka: Wave propagation in 
two-dimensional periodic lattices. Journal of the
Acoustical Society of America, 119:1995–2005, 2006.

[6] N. Sukumar and J.E. Pask (2009). Classical and enriched
finite element formulations for Bloch-periodic boundary
conditions. Int. J. Numer. Meth. Engng. 7, 8, 1121-1138,
2009.



