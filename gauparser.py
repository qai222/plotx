import re
import numpy as np
import matplotlib.pyplot as plt
from math import exp, log
'''
plot uv-vis spectrum
'''


plt.switch_backend('agg')


class ES:
    def __init__(self, os, pe, transmat, number, multi, sym):
        self.sym = sym
        self.multi = multi
        self.number = number  # start from 1
        self.os = os
        self.pe = pe
        self.transmat = transmat

    def __repr__(self):
        return '%s\t%s\t%s\t%f\t%f\n' % (self.number, self.sym, self.multi, self.nm, self.os)

    @property
    def nm(self):
        return 1239.84193/self.pe

    @staticmethod
    def eslst2eigen(eslist, unit):
        eigens = np.zeros((len(eslist), 2))
        for i in range(len(eslist)):
            if unit == 'ev':
                eigens[i] = [eslist[i].pe, eslist[i].os]
            elif unit == 'nm':
                eigens[i] = [eslist[i].nm, eslist[i].os]
        return eigens

    @staticmethod
    def splites(eslist):
        tes = [es for es in eslist if es.multi == 3]
        ses = [es for es in eslist if es.multi == 1]
        return ses, tes


def parse(log):
    alles = []
    with open(log, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if re.search(r"^\sExcited State\s*\d", lines[i]):
                items = lines[i].strip().split()
                multi = 0
                if items[3].split('-')[0] == 'Singlet':
                    multi = 1
                elif items[3].split('-')[0] == 'Triplet':
                    multi = 3
                transitions = []
                j = i + 1
                while len(lines[j].strip().split()) == 4 and '->' in lines[j]:
                    frommo, arrow, tomo, coe = lines[j].strip().split()
                    transitions.append([int(frommo), int(tomo), float(coe)])
                    j += 1
                es = ES(float(items[8][2:]), float(items[4]), transitions, int(items[2][:-1]), multi, items[3].split('-')[1])
                alles.append(es)
    return alles

def gsmear(peaks, xs, xe, npts, fwhm):
    d = np.zeros((npts, 2))
    d[:, 0] = np.linspace(xs, xe, npts)
    for i in range(npts):
        yi = 0
        for p in peaks:
            px, py = p
            yi += py * exp(-4 * log (2) * (px - d[i][0])**2/fwhm**2)
        d[i][1] = yi
    return d

class ESplt:
    def __init__(self, loglst, xsnm, xenm, fwhm, npts):
        self.namelst = [log.split('.')[-2] for log in loglst]
        self.xsnm = xsnm
        self.xenm = xenm
        self.xsev = self.nm2ev(xsnm)
        self.xeev = self.nm2ev(xenm)
        self.npts = npts
        self.fwhm = fwhm
        self.loglst = loglst  # [mol1.log, mo2.log, ...]
        self.fig, self.ax1 = plt.subplots()
        self.ax1.set_ylabel('Oscillator Strength')
        self.ax1.set_xlabel('Wavelength (nm)')
        self.ax1.set_xlim([self.xsnm, self.xenm])
        self.ax2 = self.ax1.twiny()
        self.colors = ['r', 'b', 'g', 'y', 'm', 'y', 'c']
        self.markers= ['.', 'v', '1', '3', '+', 'x', 'd']

    @staticmethod
    def nm2ev(nm):
        return 1239.84193/nm

    def plot(self, filename):
        for i in range(len(self.loglst)):
            eslst = parse(self.loglst[i])
            eigensev = ES.eslst2eigen(eslst, 'ev')
            dev = gsmear(eigensev, xs=self.xsev, xe=self.xeev, npts=self.npts, fwhm=self.fwhm)
            eigensnm = ES.eslst2eigen(eslst, 'nm')

            brighteigens = np.array([eigen for eigen in eigensnm if eigen[1] > 1e-5])
            if len(brighteigens) > 0:
                self.ax1.vlines(brighteigens[:, 0], [0], brighteigens[:, 1], linestyles='solid', colors=self.colors[i], alpha=0.1) #, label = 'Singlet Eigenvalues ' + self.namelst[i])
            darkeigens = np.array([eigen for eigen in eigensnm if eigen[1] <= 1e-5])
            if len(darkeigens) > 0:
                darkeigens[:, 1] = darkeigens[:, 1] + -0.1 * (i+1)
                self.ax1.scatter(darkeigens[:, 0], darkeigens[:,1], c=self.colors[i], marker=self.markers[i])

            dnm = dev
            dnm[:, 0] = [self.nm2ev(x) for x in dnm[:, 0]]

            self.ax1.plot(dnm[:, 0], dnm[:, 1], ls='dotted', label = self.namelst[i], color=self.colors[i])
            tick_loc = self.ax1.get_xticks()
            top_tiks = np.linspace(self.nm2ev(tick_loc[0]), self.nm2ev(tick_loc[-1]), 5)
            top_ticks = ["%.2f" % z for z in top_tiks]
            top_tick_loc = np.array( [self.nm2ev(float(i)) for i in top_ticks] )
            top_tick_loc = (top_tick_loc - min(top_tick_loc)) / (max(top_tick_loc) - min(top_tick_loc))
            self.ax2.set_xticks(top_tick_loc)
            self.ax2.set_xticklabels(top_ticks)
            self.ax2.set_xlabel('Photon Energy (eV)')
        self.ax1.legend()
        if filename.split('.')[-1] != 'eps':
            self.fig.savefig(filename, dpi=800)
        else:
            self.fig.savefig(filename)
