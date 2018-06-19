# qgraf 

## quick guide
- run the script: ``install``
- specify the process at ``qgraf.dat``
- tune the SM file to avoid particles in the diagrams
- run the ``qgraf.py``

## What does qgraf.py do? 
- a folder called ``data``  will be created
- inside there would be two files: amps.qgraf and qgraf.obj
- amps.qgraf: this is the output of qgraf.exe. It is a raw
  representatin of the feynman amps. This file is processed 
  by qgraf.exe to compute the squared amps.
- qgraf.obj: this is the output of qgraf.py. It is not human
  readable, but loadable from python. It contains string 
  represtantions of amps2 to be processed by form. 

## What to do afer this?
- the qgraf.obj is the final product of this repository.
  The user needs to mv the file to the next part which is
  the form manipulation. See the repo form for more instructions




