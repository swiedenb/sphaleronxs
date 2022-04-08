# sphaleronxs
Script to calculate the sphaleron cross section for different CoM and PDF sets

It uses numpy, scipy, optparse and lhapdf for python 

It can be simply run by:

python xs.py 

Possible options:

-e: Energy of the Sphaleron (in GeV) (default:9000)
-p: Name of the PDFSet (must be available in LHAPDF) (default:"CT10")
-f: Bool to switch between the integration, including the exponential for energies from 0 to Esph (default: 1)
-c: CoM Energy (default:13000)

The output follows: Esp xs p/m integration error p/m PDF err (if possible)

