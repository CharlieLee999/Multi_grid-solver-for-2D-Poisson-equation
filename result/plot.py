

import numpy as np
import matplotlib.pyplot as plt

errRatio_1 = np.genfromtxt("errRatio_n_4_nu1_1_nu2_1.csv", delimiter=',')

errRatio_2 = np.genfromtxt("errRatio_n_4_nu1_2_nu2_1.csv", delimiter=',')
errRatio_3 = np.genfromtxt("errRatio_n_7_nu1_1_nu2_1.csv", delimiter=',')
errRatio_4 = np.genfromtxt("errRatio_n_7_nu1_2_nu2_1.csv", delimiter=',')

iterVec = errRatio_1[:, 1]

errRatioVec_1 = errRatio_1[:, 0]
errRatioVec_2 = errRatio_2[:, 0]
errRatioVec_3 = errRatio_3[:, 0]
errRatioVec_4 = errRatio_4[:, 0]

plt.figure()
plt.plot(iterVec, np.log10(errRatioVec_1), label=r'$n=4, \nu_1=1, \nu_2 = 1$')
plt.plot(iterVec, np.log10(errRatioVec_2), label=r'$n=4, \nu_1=2, \nu_2 = 1$')
plt.plot(iterVec, np.log10(errRatioVec_3), label=r"$n=7, \nu_1=1, \nu_2 = 1$")
plt.plot(iterVec, np.log10(errRatioVec_4), label=r"$n=7, \nu_1=2, \nu_2 = 1$")

plt.legend()
plt.title(r"iteration - $log_{10}(\frac{||e^{(m)}||_{\infty}}{||e^{(0)}||_{0}})$")
plt.xlabel('n')
plt.ylabel(r'$log_{10}(\frac{||e^{(m)}||_{\infty}}{||e^{(0)}||_{0}})$')
plt.savefig("errRatio.pdf", format='pdf')


resRatio_1 = np.genfromtxt("resRatio_n_4_nu1_1_nu2_1.csv", delimiter=',')

resRatio_2 = np.genfromtxt("resRatio_n_4_nu1_2_nu2_1.csv", delimiter=',')
resRatio_3 = np.genfromtxt("resRatio_n_7_nu1_1_nu2_1.csv", delimiter=',')
resRatio_4 = np.genfromtxt("resRatio_n_7_nu1_2_nu2_1.csv", delimiter=',')

iterVec = errRatio_1[:, 1]

resRatioVec_1 = resRatio_1[:, 0]
resRatioVec_2 = resRatio_2[:, 0]
resRatioVec_3 = resRatio_3[:, 0]
resRatioVec_4 = resRatio_4[:, 0]


plt.figure()
plt.plot(iterVec, np.log10(resRatioVec_1), label=r'$n=4, \nu_1=1, \nu_2 = 1$')
plt.plot(iterVec, np.log10(resRatioVec_2), label=r'$n=4, \nu_1=2, \nu_2 = 1$')
plt.plot(iterVec, np.log10(resRatioVec_3), label=r"$n=7, \nu_1=1, \nu_2 = 1$")
plt.plot(iterVec, np.log10(resRatioVec_4), label=r"$n=7, \nu_1=2, \nu_2 = 1$")

plt.legend()
plt.xlabel('n')
plt.ylabel(r'$log_{10}(\frac{||r^{(m)}||_{\infty}}{||r^{(0)}||_{0}})$')
plt.title(r"iteration - $log_{10}(\frac{||r^{(m)}||_{\infty}}{||r^{(0)}||_{0}})$")
plt.savefig("resRatio.pdf", format='pdf')