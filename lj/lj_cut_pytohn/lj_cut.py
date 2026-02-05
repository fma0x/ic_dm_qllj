"""
LAMMPS Python Pair Potential - Lennard-Jones Cut
=================================================

This module implements the Lennard-Jones potential with cutoff
to be used with LAMMPS pair_style python command.

This implementation reproduces the behavior of the built-in
lj/cut pair style in LAMMPS.

Usage in LAMMPS:
    pair_style python 2.5
    pair_coeff * * lj_cut_python.LJCut Ar

For multiple atom types:
    pair_coeff * * lj_cut_python.LJCut Ar Ar

Author: Paulo Giovani (Created based on LAMMPS pair_lj_cut.cpp implementation)
"""

from __future__ import print_function


class LAMMPSPairPotential(object):
    """Base class for LAMMPS Python pair potentials."""
    
    def __init__(self):
        self.pmap = dict()
        """self.units = 'lj'"""
        self.units = 'real'
    
    def map_coeff(self, name, ltype):
        """Map LAMMPS atom type to internal parameter name."""
        self.pmap[ltype] = name
    
    def check_units(self, units):
        """Verify units consistency."""
        if units != self.units:
            raise Exception("Conflicting units: %s vs. %s" % (self.units, units))


class LJCut(LAMMPSPairPotential):
    """
    Lennard-Jones potential with cutoff.
    
    U(r) = 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]  for r < r_cut
    U(r) = 0                                             for r >= r_cut
    
    The force is:
    F(r) = 48 * epsilon * sigma^12 / r^13 - 24 * epsilon * sigma^6 / r^7
    
    Parameters stored as pre-computed coefficients:
    - lj1 = 48 * epsilon * sigma^12
    - lj2 = 24 * epsilon * sigma^6
    - lj3 = 4 * epsilon * sigma^12
    - lj4 = 4 * epsilon * sigma^6
    """
    
    def __init__(self):
        super(LJCut, self).__init__()
        
        # Use LJ reduced units (same as LAMMPS 'lj' units)
        self.units = 'real'
        
        # Define coefficients for Argon in LJ reduced units
        # epsilon = 1.0, sigma = 1.0 (reduced units)
        epsilon = 1.0
        sigma = 1.0
        
        # Pre-compute LJ coefficients for efficiency
        # These match the LAMMPS implementation in pair_lj_cut.cpp
        sigma6 = sigma**6
        sigma12 = sigma6 * sigma6
        
        lj1 = 48.0 * epsilon * sigma12  # For force calculation
        lj2 = 24.0 * epsilon * sigma6   # For force calculation
        lj3 = 4.0 * epsilon * sigma12   # For energy calculation
        lj4 = 4.0 * epsilon * sigma6    # For energy calculation
        
        # Store coefficients in nested dictionary
        # Format: self.coeff[itype][jtype] = (lj1, lj2, lj3, lj4)
        self.coeff = {
            'Ar': {
                'Ar': (lj1, lj2, lj3, lj4)
            }
        }
    
    def compute_force(self, rsq, itype, jtype):
        """
        Compute the pairwise force divided by distance.
        
        Args:
            rsq: Square of the distance between atoms (r^2)
            itype: LAMMPS atom type of first atom
            jtype: LAMMPS atom type of second atom
        
        Returns:
            Force divided by distance (F/r)
            
        Note:
            Returns F/r, not F, following LAMMPS convention.
            To get force vector: Fx = (F/r) * dx, Fy = (F/r) * dy, Fz = (F/r) * dz
        """
        # Get coefficients for this atom type pair
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        
        # Extract pre-computed coefficients
        lj1 = coeff[0]  # 48 * epsilon * sigma^12
        lj2 = coeff[1]  # 24 * epsilon * sigma^6
        
        # Compute inverse powers of r
        r2inv = 1.0 / rsq           # 1/r^2
        r6inv = r2inv * r2inv * r2inv  # 1/r^6
        
        # Force calculation from LAMMPS pair_lj_cut.cpp:
        # fpair = (r6inv * (lj1*r6inv - lj2)) * r2inv
        # This gives F/r (force divided by distance)
        force_over_r = r6inv * (lj1 * r6inv - lj2) * r2inv
        
        return force_over_r
    
    def compute_energy(self, rsq, itype, jtype):
        """
        Compute the pairwise potential energy.
        
        Args:
            rsq: Square of the distance between atoms (r^2)
            itype: LAMMPS atom type of first atom
            jtype: LAMMPS atom type of second atom
        
        Returns:
            Potential energy for this pair
        """
        # Get coefficients for this atom type pair
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        
        # Extract pre-computed coefficients
        lj3 = coeff[2]  # 4 * epsilon * sigma^12
        lj4 = coeff[3]  # 4 * epsilon * sigma^6
        
        # Compute inverse powers of r
        r2inv = 1.0 / rsq           # 1/r^2
        r6inv = r2inv * r2inv * r2inv  # 1/r^6
        
        # Energy calculation from LAMMPS pair_lj_cut.cpp:
        # evdwl = r6inv * (lj3*r6inv - lj4)
        energy = r6inv * (lj3 * r6inv - lj4)
        
        return energy


class LJCutReal(LAMMPSPairPotential):
    """
    Lennard-Jones potential with cutoff using LAMMPS 'real' units.
    
    Suitable for simulations with physical units:
    - Distance: Angstroms
    - Energy: kcal/mol
    - Time: femtoseconds
    
    Parameters for Argon in real units:
    - epsilon = 0.238 kcal/mol
    - sigma = 3.405 Angstroms
    """
    
    def __init__(self):
        super(LJCutReal, self).__init__()
        
        # Use LAMMPS 'real' units
        self.units = 'real'
        
        # Argon parameters in real units
        epsilon = 0.238  # kcal/mol
        sigma = 3.405    # Angstroms
        
        # Pre-compute LJ coefficients
        sigma6 = sigma**6
        sigma12 = sigma6 * sigma6
        
        lj1 = 48.0 * epsilon * sigma12
        lj2 = 24.0 * epsilon * sigma6
        lj3 = 4.0 * epsilon * sigma12
        lj4 = 4.0 * epsilon * sigma6
        
        self.coeff = {
            'Ar': {
                'Ar': (lj1, lj2, lj3, lj4)
            }
        }
    
    def compute_force(self, rsq, itype, jtype):
        """Compute force/distance for real units."""
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        lj1 = coeff[0]
        lj2 = coeff[1]
        
        r2inv = 1.0 / rsq
        r6inv = r2inv * r2inv * r2inv
        
        return r6inv * (lj1 * r6inv - lj2) * r2inv
    
    def compute_energy(self, rsq, itype, jtype):
        """Compute energy for real units."""
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        lj3 = coeff[2]
        lj4 = coeff[3]
        
        r2inv = 1.0 / rsq
        r6inv = r2inv * r2inv * r2inv
        
        return r6inv * (lj3 * r6inv - lj4)


class LJCutMultiType(LAMMPSPairPotential):
    """
    Lennard-Jones potential supporting multiple atom types with mixing rules.
    
    Implements geometric mixing for different atom type pairs:
    - epsilon_ij = sqrt(epsilon_i * epsilon_j)
    - sigma_ij = sqrt(sigma_i * sigma_j)
    
    This class can be extended to support arithmetic or other mixing rules.
    """
    
    def __init__(self):
        super(LJCutMultiType, self).__init__()
        
        self.units = 'lj'
        
        # Define parameters for different atom types
        # Example: Ar-like and Ne-like particles in reduced units
        params = {
            'Ar': {'epsilon': 1.0, 'sigma': 1.0},
            'Ne': {'epsilon': 0.7, 'sigma': 0.85}
        }
        
        # Build coefficient matrix with geometric mixing
        self.coeff = {}
        for type1_name, type1_params in params.items():
            self.coeff[type1_name] = {}
            for type2_name, type2_params in params.items():
                # Geometric mixing
                eps_ij = (type1_params['epsilon'] * type2_params['epsilon'])**0.5
                sig_ij = (type1_params['sigma'] * type2_params['sigma'])**0.5
                
                # Compute LJ coefficients
                sig6 = sig_ij**6
                sig12 = sig6 * sig6
                
                lj1 = 48.0 * eps_ij * sig12
                lj2 = 24.0 * eps_ij * sig6
                lj3 = 4.0 * eps_ij * sig12
                lj4 = 4.0 * eps_ij * sig6
                
                self.coeff[type1_name][type2_name] = (lj1, lj2, lj3, lj4)
    
    def compute_force(self, rsq, itype, jtype):
        """Compute force/distance with mixing rules."""
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        lj1 = coeff[0]
        lj2 = coeff[1]
        
        r2inv = 1.0 / rsq
        r6inv = r2inv * r2inv * r2inv
        
        return r6inv * (lj1 * r6inv - lj2) * r2inv
    
    def compute_energy(self, rsq, itype, jtype):
        """Compute energy with mixing rules."""
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        lj3 = coeff[2]
        lj4 = coeff[3]
        
        r2inv = 1.0 / rsq
        r6inv = r2inv * r2inv * r2inv
        
        return r6inv * (lj3 * r6inv - lj4)


# Convenience aliases
LJCutMelt = LJCut  # Commonly used name in examples
