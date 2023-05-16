# Python wrapper implementation

To compile you may need to change the LIBS line in the Makefile.

An example script is [tddft_test.py](tddft_test.py).

For different examples from the paper, there are separate scripts
Example_1.py 1 particle, 1 dimension separating wavepack on harmonic oscillator
Example_2.py 1 particle, 3 dimensions bilinearly coupled harmonic oscillator potential
Example_3.py 2 particles, 1 dimension, driven wavepacket from ground state
Example_4.py 5 partiles, 1 dimension, spin system with 15 lattice points

## Potential Energy Surfaces from adiabatic density

Relevant data is in python/data folder.

### Singlet
Run `python lih_along.py`
Analyse data with `plot_densities.ipynb`

### Triplet
Run `python lih_along_t.py`
Analyse data with `plot_densities_t.ipynb`

### Noisy Singlet
Run `lih_noisy_start.py`
Analyse data with `plot_densities_noisy_start.ipynb`

### From quantum algorithm
Run `pes_without_energy_qc.ipynb` to calculate measured densities
Run `lih_from_qc.ipynb` to run TDDFTinversion
Analyse data with `plot_densities_qc.ipynb`

### Three and Four particle systems
`3_part.ipynb`
`4_part.ipynb`
