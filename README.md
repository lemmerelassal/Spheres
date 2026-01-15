# Molecular Dynamics and PDB Visualization System

A comprehensive molecular visualization and simulation platform built in Unreal Engine 5, featuring real-time molecular dynamics, binding affinity calculations, and interactive 3D visualization of protein structures.

## Overview

This system provides a complete toolkit for molecular modeling and analysis, including:

- PDB file loading and visualization
- Real-time molecular dynamics simulation
- MM/GBSA binding affinity calculations
- Molecular interaction analysis (hydrogen bonds, salt bridges, pi-stacking, hydrophobic contacts)
- Interactive 3D camera controls
- Comprehensive UI for structure navigation and control

## Core Components

### PDBViewer
The central visualization engine that handles molecular structure display and interaction detection.

**Key Features:**
- Load and parse PDB files
- Render atoms and bonds with CPK coloring
- Calculate molecular interactions
- Manage multiple ligand structures
- Dynamic visibility control

**Main Classes:**
- `APDBViewer` - Primary viewer actor
- `UPDBTreeNode` - Hierarchical structure representation
- `UPDBMoleculeNode` - Individual molecule data

### Molecular Dynamics Engine
High-performance O(N) molecular dynamics simulation with spatial hashing optimization.

**Key Features:**
- Multiple integration methods (Velocity Verlet, Leap-Frog, Runge-Kutta 4)
- Lennard-Jones and electrostatic forces
- Bond constraint enforcement (RATTLE algorithm)
- Temperature control via thermostat
- Real-time energy monitoring

**Main Classes:**
- `AMolecularDynamics` - Simulation engine
- `FAtomState` - Per-atom state tracking
- `FBondConstraint` - Bond constraint data

**Parameters:**
- Time step: 0.0001 - 0.01 fs
- Temperature: 0 - 1000 K
- Damping factor: 0.0 - 1.0
- Cutoff distance: Configurable for force calculations

### MM/GBSA Binding Affinity Calculator
Molecular Mechanics / Generalized Born Surface Area method for calculating protein-ligand binding free energies.

**Key Features:**
- Complete energy decomposition (electrostatic, VDW, solvation, surface area)
- Structure quality assessment
- Clash detection and severity analysis
- Automatic charge assignment
- Ki/Kd dissociation constant calculation

**Main Classes:**
- `AMMGBSA` - Main calculation engine
- `FBindingAffinityResult` - Comprehensive results structure
- `FEnergyComponents` - Energy breakdown
- `FStructureQuality` - Quality assessment

**Energy Components:**
- Electrostatic interactions (Coulomb's law with distance-dependent dielectric)
- Van der Waals forces (Lennard-Jones 6-12 potential with soft-core)
- GB solvation (Generalized Born model)
- Nonpolar solvation (SASA-based)
- Conformational entropy penalty

**Output:**
- Binding free energy (ΔG in kcal/mol)
- Dissociation constant (Ki in nM)
- Affinity classification (Very Strong, Strong, Moderate, Weak)
- Detailed energy breakdown
- Structure quality report

### Hydrogen Generator
Automated system for adding explicit hydrogens to molecular structures.

**Key Features:**
- Valence-aware hydrogen placement
- Support for sp3, sp2, and sp hybridization
- Geometry-specific positioning (tetrahedral, trigonal planar, linear)
- Bond order consideration

**Main Classes:**
- `FHydrogenGenerator` - Static hydrogen generation methods
- `FAtomBondInfo` - Bond connectivity tracking

### Interaction Analysis
Comprehensive molecular interaction detection and visualization.

**Interaction Types:**
- Hydrogen bonds (distance and angle criteria)
- Salt bridges (charged residue pairs)
- Pi-stacking (aromatic ring interactions)
- Hydrophobic contacts

**Main Classes:**
- `FMolecularInteraction` - Interaction data structure
- Implemented in `PDBViewer_Interactions.cpp`

### User Interface Widgets

#### MoleculeListWidget
List view of all loaded molecules with visibility controls and binding affinity display.

**Features:**
- Per-molecule visibility toggle
- One-click binding affinity calculation
- Color-coded affinity display
- Real-time updates

#### InteractionControlWidget
Control panel for molecular interaction analysis.

**Features:**
- Calculate interactions on demand
- Toggle interaction type visibility
- Protein-protein and protein-ligand filtering
- Interaction count statistics
- Detailed interaction list

#### MDControlWidget
Molecular dynamics simulation control panel.

**Features:**
- Start/stop/reset simulation
- Real-time energy monitoring
- Temperature control
- Integration method selection
- Force field parameter adjustment

#### PDBStructureWidget
Hierarchical structure browser with load/save functionality.

**Features:**
- Tree view of chains, residues, and atoms
- Load PDB files
- Save modified structures
- Clear current structure

### Camera System
Advanced orbit camera with intuitive controls.

**Controls:**
- Right mouse button: Rotate around pivot
- Middle mouse button: Pan view
- Mouse wheel: Zoom in/out
- Spacebar: Fit molecule to screen

**Main Classes:**
- `APDBCameraComponent` - Camera controller

## Technical Implementation

### Performance Optimizations

**Spatial Hashing:**
- O(N) neighbor search for force calculations
- Dynamic cell size based on cutoff distance
- Efficient hash function for 3D grid

**Verlet Neighbor Lists:**
- Reduced force calculation frequency
- Skin distance for list validity
- Automatic rebuild when needed

**Lazy Updates:**
- On-demand interaction calculation
- Cached binding affinity results
- Conditional mesh updates

### Physics Models

**Force Fields:**
- AMBER-compatible atom typing
- Partial charge assignment based on electronegativity
- Van der Waals radii from experimental data

**Integration Methods:**
- Velocity Verlet (default, 2nd order symplectic)
- Leap-Frog (energy-conserving, 2nd order)
- Runge-Kutta 4 (4th order accuracy, slower)

**Constraints:**
- RATTLE algorithm for bond length constraints
- Iterative satisfaction with configurable tolerance
- Maintains detailed balance for thermodynamics

### Thermodynamics

**Energy Conservation:**
- Total energy monitoring (kinetic + potential)
- Configurable damping for stability
- Temperature scaling for NVT ensemble

**Entropy Calculations:**
- Base conformational entropy penalty
- Heavy atom count scaling
- Ligand-specific corrections

## Usage Examples

### Loading a Structure
```cpp
APDBViewer* Viewer = GetWorld()->SpawnActor<APDBViewer>();
Viewer->LoadPDBFromPath("/Path/To/Structure.pdb");
```

### Running Molecular Dynamics
```cpp
AMolecularDynamics* MD = GetWorld()->SpawnActor<AMolecularDynamics>();
MD->InitializeFromViewer(Viewer);
MD->SetTemperature(300.0f);
MD->SetTimeStep(0.002f);
MD->StartSimulation();
```

### Calculating Binding Affinity
```cpp
AMMGBSA* Calculator = GetWorld()->SpawnActor<AMMGBSA>();
Calculator->InitializeFromViewer(Viewer);
FBindingAffinityResult Result = Calculator->CalculateBindingAffinity("LigandKey");
UE_LOG(LogTemp, Log, TEXT("ΔG = %.2f kcal/mol, Ki = %.2f nM"), 
    Result.DeltaG_Binding, Result.Ki_nM);
```

### Analyzing Interactions
```cpp
Viewer->CalculateAllInteractions(true, true); // Protein-protein and protein-ligand
TArray<FMolecularInteraction> HBonds = Viewer->GetInteractionsByType(EInteractionType::HydrogenBond);
```

## Configuration

### MM/GBSA Parameters
- **Temperature:** 298.15 K (physiological)
- **Interior Dielectric:** 10.0 (protein interior)
- **Exterior Dielectric:** 78.5 (water)
- **Salt Concentration:** 0.15 M (physiological)
- **Charge Scaling:** 0.95 (for unminimized structures)
- **GB Scale Factor:** 0.15 (empirical fit)

### MD Parameters
- **Time Step:** 0.002 fs (default, adjustable)
- **Cutoff Distance:** 1200 Å (12 Å in real units)
- **Neighbor List Skin:** 200 Å (2 Å buffer)
- **Constraint Iterations:** 3 (RATTLE)

## File Structure

```
├── PDBViewer.h/cpp                    # Main visualization engine
├── PDBViewer_Interactions.cpp         # Interaction detection
├── MolecularDynamics.h/cpp            # MD simulation engine
├── MMGBSA.h/cpp                       # Binding affinity calculator
├── HydrogenGenerator.h/cpp            # Hydrogen addition
├── PDBCameraComponent.cpp             # Camera controls
├── MoleculeListWidget.h/cpp           # Molecule list UI
├── MoleculeListEntry.h/cpp            # List entry widget
├── InteractionControlWidget.h/cpp     # Interaction UI
├── MDControlWidget.h/cpp              # MD control UI
├── PDBStructureWidget.h/cpp           # Structure browser UI
└── README.md                          # This file
```

## Dependencies

- Unreal Engine 5.6+
- UMG (Unreal Motion Graphics) for UI
- Standard C++ library

## Future Enhancements

- GPU acceleration for force calculations
- MMFF94 force field implementation
- Enhanced solvation models (Poisson-Boltzmann)
- Trajectory recording and playback
- Energy minimization algorithms (steepest descent, conjugate gradient)
- Multi-molecule docking
- Pharmacophore mapping
- QSAR model integration

## License

This project is provided as-is for educational and research purposes.

## References

- Molecular Mechanics: Leach, A.R. "Molecular Modelling: Principles and Applications"
- MM/GBSA: Genheden, S. & Ryde, U. "The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities"
- Generalized Born: Onufriev, A. et al. "Exploring protein native states and large-scale conformational changes"
- RATTLE Algorithm: Andersen, H.C. "Rattle: A 'velocity' version of the shake algorithm"

## Contact

For questions, issues, or contributions, please contact:

Email: elassallemmer@gmail.com
