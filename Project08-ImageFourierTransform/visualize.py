import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit

sns.set_theme()


def num_2_prime_factors(num):
    pf = 0
    cur = int(num) 

    while(cur % 2 == 0):
        cur /= 2
        pf += 1

    return pf


n2pf = np.vectorize(num_2_prime_factors)

lens, means, stds = np.genfromtxt("build/output/times_lin.npy", unpack=True)
# plt.errorbar(lens, means, yerr=stds)
n2pf_lens = n2pf(lens)
cmap = plt.get_cmap(plt.rcParams["image.cmap"]+"_r", np.max(n2pf_lens) - np.min(n2pf_lens) +1)
sc = plt.scatter(lens, means, s=2, c=n2pf_lens, cmap=cmap)
plt.colorbar(sc, ticks=np.arange(np.min(n2pf_lens), np.max(n2pf_lens)+1)+.5, format="%i", label="\# 2-Prime Factors")

plt.xlabel("FFT length $N$")
plt.ylabel("FFT Avg. Time $t/ms$")

plt.tight_layout()
plt.savefig("build/plots/times_lin_c.pdf")
plt.clf()

sc = plt.scatter(lens, means, s=2)
plt.xlabel("FFT length $N$ ")
plt.ylabel("FFT Avg. Time $t/ms$")

plt.tight_layout()
plt.savefig("build/plots/times_lin.pdf")
plt.clf()

lens, means, stds = np.genfromtxt("build/output/times_p2.npy", unpack=True)


def scaling(x, a):
    return x*np.log2(x)*a


x = np.linspace(np.min(lens), np.max(lens), 1000)
popt, pcov = curve_fit(scaling, lens, means, p0=(means[-1]/lens[-1]))
plt.plot(np.log2(x), scaling(x, *popt), label=r"$c\cdot N\log\! N$-fit")
plt.plot(np.log2(lens), means, "x", label="Data")

plt.xticks(np.log2(lens))

plt.xlabel(r"FFT length $\log_2 N$")
plt.yscale("log")
plt.ylabel("FFT Avg. Time $t/ms$")

plt.legend()

plt.tight_layout()
plt.savefig("build/plots/times_p2.pdf")
plt.cla()

f, fft_r, fft_imag = np.genfromtxt("build/output/plot.npy", unpack=True)

x = np.linspace(-1, 1, f.shape[0])
plt.plot(x, f)

plt.xlabel("$x$")
plt.ylabel("$f(x)$")

plt.tight_layout()
plt.savefig("build/plots/fft_example.pdf")
plt.cla()

plt.xlabel(r"$\nu$")
plt.ylabel(r"$\hat f(\nu x)$")

x = np.linspace(-2, 2, f.shape[0])
plt.plot(x, fft_r, label="Real")
plt.plot(x, fft_imag, label="Imaginary")

plt.xlim(-.3, .3)
plt.ylim(-1e-1, 1e-1)

plt.legend()
plt.tight_layout()
plt.savefig("build/plots/fft_example_fft.pdf")
plt.cla()