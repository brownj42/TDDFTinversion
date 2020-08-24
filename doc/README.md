# Tutorial

This is a basic introduction to using the program.

## Initialization

The initialization mainly uses derived types and subroutines located in the [derivedtypes module](derivedtypes.md).  

There are three different mandatory derived types that are used in the program.

* [systemparameters](systemparameters.md): This contains all the information that defines the system
being studied
* [ksvalues](ksvalues.md): This contains all the information about the Kohn-Sham system
* [sharedvalues](sharedvalues.md): This contains the one-body potential, time-independant 1-body
potential, the 2-body interaction potential.

Additionally, there is the optional

* [fullvalues](fullvalues.md): This contains the full n-body wavefunction and the n-body
potential

The program uses a direct-product of grids. To initialize the systemparameters
derived types, you need to define the number of grid points in one-dimension (L)
and then  

* Fortran:

```f90
   use derivedtypes
   integer, parameter :: number_of_1d_points=L
   type(systemparameters) :: sysparams
   call init_systemparameters(sysparams,number_of_1d_points)
```

* Python:

```python
   import TDDFTinversion as td
   sysparams=td.derivedtypes.init_systemparameters(number_of_1d_points)
```

You now have access to the systemparameters derived type that you *must* fill
out a few values with.

* Fortran

```f90
  sysparams%np1=L ![*integer*]
  sysparams%nd= !number of dimensions [*integer*]
  sysparams%npart = !number of particles [*integer*] 
  sysparams%xmin = !minimum grid point position [*integer*]
  sysparams%xmax = !maximum grid point value [*integer*]
  sysparams%ct = !current time [*real(8)*]
  sysparams%dth = !goal time step size [*real(8)*]
  sysparams%dvksmax = !maximum allowable potential derivative with respect to
                      !time (1000 is generally a reasonable choice) [*real(8)*]
  sysparams%pinv0minresqlp1 = !use pseudoinverse, set to 0 (only for 1D
                              !problems); use MINRES-QLP, set to 1 (The preferred choice) [*integer*]
  sysparams%quantization = !First Quantization, set to 1; Second quantization,
                           !set to 2 [*integer*]
  sysparams%energy = !KS system energy, can be anthing [*real(8)*]
  sysparams%T = !1D Laplacian d^2/dx^2 (can use ``call buildkeo(sysparams)" 
                !for using sinc basis) [*real(8)*, *dimension(L,L)*]
  sysparams%xlattice = !positions of lattice points [*real(8), dimension(L)*]
```

* Python

```python
  sysparams.np1=L #(*np.int32*)
  sysparams.nd= #number of dimensions (*np.int32*)
  sysparams.npart = #number of particles (*np.int32*)
  sysparams.xmin = #minimum grid point position (*np.int32*)
  sysparams.xmax = #maximum grid point value (*np.int32*)
  sysparams.ct = #current time (can be anything) (*np.float64*)
  sysparams.dth = #goal time step size (*np.float64*)
  sysparams.dvksmax = #maximum allowable potential derivative with respect to
                      #time (1000 is generally a reasonable choice) (*np.float64*)
  sysparams.pinv0minresqlp1 = #use pseudoinverse, set to 0 (only for 1D
                              #problems); use MINRES-QLP, set to 1 (The preferred choice) (*np.int32*)
  sysparams.quantization = #First Quantization, set to 1; Second quantization,
                           #set to 2 (*np.int32*)
  sysparams.energy = #KS system energy, can be anthing (*np.float64*)
  sysparams.T = #1D Laplacian d^2/dx^2 (can use ``call buildkeo(sysparams)"
                #for using sinc basis) (*np.float64, shape((L,L))*)
  sysparams.xlattice = #position of lattice points [*np.float64,shape(L)*]
```

## Fill out system parameters

There are some parts of the sysparams derived type that need to be calculated
using the above information

* Fortran:

```f90
  call fill_systemparameters(sysparams)
```

* Python:

```python
  td.derivedtypes.fill_systemparameters(sysparams)
```

## Generate sharedvalues derived type

* Fortran:

```f90
  type(sharedvalues):: sharedvals
  call init_sharedvalues(sysparams,sharedvals)
```

* Python:

```python
  sharedvals=td.derivedtypes.init_sharedvalues(sysparams)
```

## Add information to sharedvalues

There are two portions of sharedvalues that need to be defined

* The 1-body potential of length L^(number of dimensions) sharedvals%v1
* The 1-body time-independant potential of length L^(number of dimensions)
sharedvals%vin
* The 2-body potential of length L^(number of dimensions)^2

For certain parameters, can obtain the examples in the paper using

* Fortran:

```f90
  use potential
  call generate_1bodypot(sysparams,sharedvals)
  call generate_2bodypot(sysparams,sharedvals)
```

* Python:

```python
  td.potential.generate_1bodypot(sysparams,sharedvals)
  td.potential.generate_2bodypot(sysparams,sharedvals)
```

## Initialize the fullvalues derived type

* Fortran

```f90
  type(fullvalues) :: fullvals
  call init_fullvals(sysparams,fullvals)
```

* Python

```python
  fullvals=td.derivedtypes.init_fullvals(sysparams)
```

## Generate the full n-body potential

More information about the potential functions can be found in [potentials.md](potentials.md) 

* Fortran

```f90
  call generate_nbodypot(sysparams,sharedvals,fullvals)
```

* Python

```python
  td.potential.generate_nbodypot(sysparams,sharedvals,fullvals)
```

## Initialize the full system

For certain example parameters, you can call

* Fortran:

```f90
  call initializefullsystem(sysparams,fullvals)
```

* Python

```python
  td.initial_states.initializedfullsystem(sysparams,fullvals)
```

You can also define your own initial fullvals%psi of length
sysparams%ntot (sysparams.ntot). If the number of particles (sysparams%npart)
is 1 then you *must* define fullvals%psi of length sysparams%ntot1=L^(number
of dimensions)

## Create vector of density

Create a double precision or numpy.float64 vector of density of length
sysparams%ntot1. Put the initial density of the system in this vector. If
examining one of the examples from the paper

* Fortran:

```f90
  use density
  real(8), allocatable :: dpe(:)
  allocate(dpe(sysparams%ntot1))
  call fullwf_density(sysparams,fullvals.psi,dpe)
```

* Python

```python
  dpe=np.zeros(sysparams.ntot1,dtype=np.float64)
  td.density.fullwf_density(sysparams,fullvals%psi,dpe)
```

## Initialize the KS system.

The program attemps to initialize the Kohn-Sham system using

* Fortran:

```f90
  type(ksvalues) :: ksvals
  call initializekssystem(sysparams,sharedvals,dpe,fullvals,ksvals)
```

* Python:

```python
  ksvals=td.initial_states.initializekssytem(sysparams,sharedvals,dpe,fullvals)
```

## Create other four double precision or numpy.float64 vectors needed to

   propagate system

* Fortran:

```f90
  real(8), allocatable :: dpenew(:),dnx(:),ddnx(:),ddnxnew(:)
  allocate(dpenew(sysparams%ntot1)) !density at time step sysparams%ct+sysparams%dt
  allocate(dnx(sysparams%ntot1)) !time derivative of density at time sysparams%ct
  allocate(ddnx(sysparams%ntot1)) !second time derivative of density at
    !time sysparams%ct
  allocate(ddnxnew(sysparams%ntot1)) !second time derivative of density at
    !time sysparams%ct+sysparams%dt
```

* Python

```python
  dpenew=np.zeros(sysparams.ntot1,dtype=np.float64) #density after time step
  dnx=np.zeros(sysparams.ntot1,dtype=np.float64) #current time derivative of density
  ddnx=np.zeros(sysparams.ntot1,dtype=np.float64) #current second time
    #derivative of density
  ddnxnew=np.zeros(sysparams.ntot1,dtype=np.float64) #second time derivative
    #of density after time step sysparams.dt
```



# Main propagation step

## Obtain necessary vectors

To obtain one time step, you need to calculate the five vectors,
dpe,dpenew,dnx,ddnx,and ddnxnew. For the examples in the paper, this can be
done by

* Fortran:
```f90
  use propagate
  complex(8), allocatable :: psinew(:)
  call fullwf_density(sysparams,fullvals%psi,dpe)
  call calcdnx(sysparams,sharedvals,sysparams%ntot1,fullvals%psi,fullvals%v,dnx)
  call calcddnx(sysparams,sharedvals,sysparams%ntot1,fullvals%psi,fullvals%v,ddnx)
  !Propagate full system and calculate new density and second derivative of density
  allocate(psinew(sysparams%ntot))
  call advancewf(sysparams,sharedvals,25,fullvals%v,fullvals%psi,psinew)
  call fullwf_density(sysparams,psinew,dpenew)
  call calcddnx(sysparams,sharedvals,sysparams%ntot1,psinew,fullvals%v,ddnxnew)
```
* Python
```python
  td.density.fullwf_density(sysparams,fullvals.psi,dpe)
  td.density.calcdnx(sysparams,sharedvals,sysparams.ntot1,fullvals.psi,fullvals.v,dnx)
  td.density.calcddnx(sysparams,sharedvals,sysparams.ntot1,fullvals.psi,fullvals.v,ddnx)
  #propagate full system
  td.propagate.advancewf(sysparams,sharedvals,25,fullvals.v,fullvals.psi,psinew)
  #Calculate density and second derivative of density at time ct+dt using psinew
  td.density.fullwf_density(sysparams,psinew,dpenew)
  td.density.calcddnx(sysparams,sharedvals,sysparams.ntot1,psinew,fullvals.v,ddnxnew)
```

If one is obtaining these values by some other method, (like using a quantum
computer. dpe,dnx,ddnx are the values at time sysparams%ct and dpenew,ddnxnew
are the values at time sysparams%ct+sysparams%dt.

## Attempt to propagate system

* Fortran:
```f90
  call advancekssystem(dpe,dpenew,dnx,ddnx,ddnxnew,sysparams,KSvals,sharedvals,info)
```
* Python
```python
  info = td.ksadvance_mod.advancekssystem(dpe,dpenew,dnx,ddnx,ddnxnew,sysparams,KSvals,sharedvals)
```

### If info=1
* The orbitals advanced successfully.
* Data has been output to
  * times.dat has the current time sysparams%dt
  * wave.dat has the current KS orbitals at times in times.dat
  * pots.dat has the potential values at times in times.dat
  * dense.dat has the full density and KS density at times in times.dat

wave.dat, pots.dat and dense.dat all separate different times by a double
space in accordance with the gnuplot format index style. A sample output for
each with titles for each column is given in
* [times_sample.dat](./times_sample.dat)
* [wave_sample.dat](./wave_sample.dat)
* [pots_sample.dat](./pots_sample.dat)
* [dense_sample.dat](dense_sample.dat)


### If info=0


* The time step attempt has failed and (sysparams%dt) has been halved. A new ddnxnew and dpenew will need to
be calculated at time sysparams%ct+sysparams%dt.


## Obtaining information about KS system
Apart from the information being outputed to a data file. All the information
about the KS system is included in the ksvals derived type. More information
about that derived type can be found at.
