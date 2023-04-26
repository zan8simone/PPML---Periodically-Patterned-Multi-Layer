# PPML---Periodically-Patterned-Multi-Layer
## Electromagnetic properties of patterned multilayers based on RCWA (Rigorous Coupled Wave Analysis)

Rigorous coupled wave analysis (RCWA) based on the scattering matrix (SM) algorithm is one of the most powerful tools for the electromagnetic simulation of patterned multilayer structures. 

PPML - RCWA is a project which implements the SM-RCWA, based on the formalisms of [a-d]. 
Three groups of functions are available: one is for 1-d patterns under TM polarization, another is for 1-d anisotropic (biaxial) patterns, the third for certain 2-d patterns.
For 1-d TM patterns, currently available are functions for the calculation of
- intensity reflectance, transmittance, and layer-by-layer absorbance 
- the full 2x2 scattering matrix 
- the E and S fields inside the structure
 
The proper factorization rules [c,d] make the code extremely performing, and fully suitable for the simulation of metal components (plasmonic gratings).
Out-of-plane uniaxial materials can be treated.
For 1-d biaxial media, the following function is available:
-  amplitude, phase and polarization of transmitted diffracted waves

For 2-d patterns, currently available are functions for the calculation of
-  zero-diffraction order scattering matrix (i.e, the 4x4 matrix whose submatrices are the transmission and reflection Jones matrices.)

Allowed unit cell geometries are rectangular and L-shaped inclusions (see the Manual for details). The proper factorization rules are implemented [c,d], thus allowing for a fast convergence even in presence of metallic inclusions.
Conducting interfaces can be seamlessly included in all the geometries by specifying their in-plane conductivity.
The present code is distributed for free, but we kindly ask you to cite its source and, if applies, the publications below.
Several publications are based on PPML (see below). For some of them, you can find the corresponding tutorial in the software package.


## Publications based on PPML (not exhaustive list)
1. S. Zanotto, ... A. Pitanti, "Photonic bands, superchirality, and inverse design of a chiral minimal metasurface", Nanophotonics (2019), DOI: 10.1515/nanoph-2019-0321
2. S. Zanotto, ... A. Pitanti, "Optomechanics of Chiral Dielectric Metasurfaces", Advanced Optical Materials (2020), DOI: 10.1002/adom.201901507
3. S. Zanotto, ... D. S. Wiersma, "Multichannel remote polarization control enabled by nanostructured Liquid Crystalline Networks", Applied Physics Letters (2019), DOI 10.1063/1.5096648
4. S. Zanotto, G. C. La Rocca, and A. Tredicucci, “Understanding and overcoming fundamental limits of asymmetric light-light switches”, Optics Express 26, 3, 3618 (2018).
5. S. Zanotto, ..., A. Melloni, "Metasurface reconfiguration through lithium ion intercalation in a transition metal oxide", Advanced Optical Materials 2017, 5, 1600732 (2017).
6. S. Zanotto and A. Tredicucci, "Universal lineshapes at the crossover between weak and strong critical coupling in Fano-resonant coupled oscillators", Scientific Reports 6, 24592 (2016).
7. L. Baldacci, ..., A. Tredicucci, “Interferometric control of absorption in thin plasmonic metamaterials: general two port theory and broadband operation”, Optics Express 23, 9202 (2015).
8. S. Zanotto, ... A. Tredicucci, “Perfect energy feeding into strongly coupled systems and interferometric control of polariton absorption”, Nature Physics 10, 830 (2014).
9. J.-M. Manceau, ..., R. Colombelli, “Mid-infrared intersubband polaritons in dispersive metal-insulator-metal resonators”, Appl. Phys. Lett. 105, 081105 (2014).
10. J.-M. Manceau, ..., R. Colombelli, “Optical critical coupling into highly confining metal-insulator-metal resonators”, Appl. Phys. Lett. 103, 091110 (2013).
11. S. Zanotto, ..., A. Tredicucci, “Ultrafast optical bleaching of intersubband cavity polaritons”, Phys. Rev. B 86, 201302(R) (2012).
12. S. Zanotto, ... A. Tredicucci, “Analysis of line shapes and strong coupling with intersubband transitions in one-dimensional metallodielectric photonic crystal slabs”, Phys. Rev. B 85, 035307 (2012).
13. R. Degl'Innocenti, ... A. Tredicucci, “One-dimensional surface-plasmon gratings for the excitation of intersubband polaritons in suspended membranes”, Solid State Comm. 151, 1725-1727 (2011).
14. S. Zanotto, ... A. Tredicucci, “Intersubband polaritons in a one-dimensional surface plasmon photonic crystal”, Appl. Phys. Lett. 97, 231123 (2010).

## References for the method
a. 	D. M. Whittaker & I. S. Culshaw, "Scattering-matrix treatment of patterned multilayer photonic structures",
Phys. Rev. B 60, 2610 (1999).

b.	M. Liscidini, D. Gerace, L. C. Andreani & J. E. Sipe, "Scattering-matrix analysis of periodically patterned multilayers with asymmetric unit cells and birefringent media", Phys. Rev. B 77, 035324 (2008).

c. 	L. Li. "Use of Fourier series in the analysis of discontinuous periodic structures". J. Opt. Soc. Am. A 13, 1870 (1996).

d.	Lalanne, Philippe, and G. Michael Morris, "Highly improved convergence of the coupled-wave method for TM polarization", JOSA A 13, 779 (1996).
