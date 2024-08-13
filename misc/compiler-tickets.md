# Help tickets 

__Authors:__ Spencer Bryngleson (Georgia Tech), Henry Le Berre (Georgia Tech), Steve Abbott (HPE/Cray), Reuben Budiardja (OLCF)

MFC has been used to isolate many compiler bugs. A short summary is below.

## OLCF (Frontier)

* OLCFHELP-15413 'Performance on Frontier'
* OLCFHELP-14881 'Frontier Hackathon'
* OLCFDEV-1416 'CCE ACC routine seq'
* OLCFHELP-19298 'CCE OpenACC doesn't respect `acc declare create` statements when an abstract interface is in the same compilation unit'

## CCE (Cray, internal)

* HPECOE-28 'Compiler performance bug on autogen kernels - LLVM IR'
* CAST-31898, Case 5370334261
* PE-47678
* PE-47679, 'CCE-15.0.1,16.0.0 Unsupported OpenACC construct '!$acc routine seq''
* PE-53844, 'CCE >17 crashes on host_data use_device for module variables'
* PE-55153, 'xor crashes logical comparison for Fortran extensions'
* HPECOE-71, associated with `OLCFHELP-19298`

## NVHPC (NVIDIA, internal)

* NV-33317 'NVHPC 22.5 fort2 TERMINATED by signal 11.'

## GNU - gfortran

* GNU-106643 '[gfortran + OpenACC] Allocate in module causes refcount error.'
* GNU-108895 '[13.0.1 (exp)] Fortran + gfx90a !$acc update device produces a segfault.'
* GNU-107227 'Compiler bug in private allocatable array in OpenACC compute statement'
