#! /usr/bin/env python
import lhapdf
import numpy as np
from scipy.integrate import quad,romberg
import optparse
def calc_xs(energy,p,full, ecm):
    #Center-of-mass energy
    E_CM = ecm
    #quark list
    quarks = [-1,-2,-3,-4,1,2,3,4]
    #alpha_w taken from paper
    alpha_w = 1./30.
    # assume c = 2 for S(E) approximation
    c = 2.
    # initialize the sphaleron threshold energy
    E_sph = energy
    # assume a for S(E) approximation
    a = -0.005
    # W mass (GeV)
    m_W = 80.73
    # conversion factor from 1/GeV^2 to fb
    tobarn = 389379290.5*1e3
    # initialize error and result variables
    err = 0.
    result = 0.
    # loop through all possible quark combinations
    for i in quarks:
        for j in quarks:
            # definition of the parton luminosity function
            def function(y,E):
                x1 = np.sqrt(E**2/E_CM**2) * np.exp(y)
                x2 = np.sqrt(E**2/E_CM**2) * np.exp(-y)
                return 2*E/(E_CM**2) * p.xfxQ2(i,x1, x1*x2*E_CM**2)*p.xfxQ2(j,x2,x1*x2*E_CM**2)
            # definition of the cross section integral
            def y_integral(E):
                if(E/E_sph >= 1.):
                    return quad(function,np.log(np.sqrt(E**2/E_CM**2)),-np.log(np.sqrt(E**2/E_CM**2)),args=(E))[0]
                else:
                    return np.exp(c*(4*np.pi/alpha_w)*((1-a)*(E/E_sph)+a*(E/E_sph)**2-1))*quad(function,np.log(np.sqrt(E**2/E_CM**2)),-np.log(np.sqrt(E**2/E_CM**2)),args=(E))[0]
            # perform the integral
            if full:
                int_res = quad(y_integral,0.,E_CM)
            else:
                int_res = quad(y_integral,E_sph,E_CM)
            # sum over all possible quark combinations (for integral result and the integral error)
            # factor if 2/3 if quarks are the same
            a=1.
            if i == j:
                a = 2./3.
            result += int_res[0]*a
            err += ((int_res[1]*a)**2)

    # final calculation of the squared errors and converting into fb
    errfinal = np.sqrt(err)*1./(m_W**2) * 1./4. * tobarn
    xs = result*1./(m_W**2) * 1./4. * tobarn
    xs_vec = [xs,errfinal]
    return xs_vec

def pdf4lhc_formula(weights):
    '''Calculate and return mean + uncertainty of PDF shifts

    The calculation follows the latest implementation of PDF4LHC. Available at:
    http://arxiv.org/pdf/1510.03865v1.pdf
    Fallback method described in:
    https://indico.cern.ch/event/459797/contribution/2/attachments/1181555/1800214/mcaod-Feb15-2016.pdf
    '''

    # first weight is from production PDF
    mean_production = weights.pop(0)
    weights.sort()
    # use weights marking the beginning of the 16% and 84% quantile
    lower_index = int(round(len(weights) * 0.16));
    upper_index = int(round(len(weights) * 0.84));

    uncertainty = (weights[upper_index] - weights[lower_index]) / 2.0;
    return mean_production + uncertainty;
def main():
    parser = optparse.OptionParser( description = 'Calculate sphaleron xs for different energies and pdf sets', usage = 'usage: %prog energy pdfset')
    parser.add_option("-e",action="store",type="float",dest="energy",default=9000)
    parser.add_option("-p",action="store",type="str",dest="pdfset",default="CT10")
    parser.add_option("-f",action="store",type="int",dest="full",default=1)
    parser.add_option("-c",action="store",type="float",dest="ecm",default=13000)
    (options, args ) = parser.parse_args()
    pdf = lhapdf.mkPDFs(options.pdfset)
    for e in range(8000,10500,500):
        xs = calc_xs(e, pdf[0], options.full, options.ecm)
        xs_pdf_diff = []
        if len(pdf)>1:
            for i in range(len(pdf)):
                xs_pdf_diff.append(calc_xs(e, pdf[i], options.full, options.ecm)[0])
            xs_pdf_err = abs(xs[0] - pdf4lhc_formula(xs_pdf_diff))
        else: xs_pdf_err = 0.
        print("Mass: ", e, " " + str(xs[0]) + " +/- " + str(xs[1])  + " (Int. Err) " + "+/- " + str(xs_pdf_err) + "(PDF Err) [fb]")
if __name__ == '__main__':
   main()

