from __future__ import annotations

import re
from collections import defaultdict
from math import log, exp
from typing import Union

import numpy as np
import seaborn as sns
from loguru import logger
from monty.json import MSONable
from pandas._typing import FilePath

'''
extract excitation info (td) from gaussian log file
'''
sns.set_style('darkgrid')


class Transition(MSONable):
    def __init__(self, from_orb: str, to_orb: str, frac: float, index: int = None):
        self.index = index
        self.frac = frac
        self.to_orb = to_orb
        self.from_orb = from_orb

    def __gt__(self, other: Transition):
        if self.index != other.index:
            logger.warning("comparing transitions belong to different excited states")
        return self.frac > other.frac

    def __lt__(self, other):
        return not self.__gt__(other)

    def __repr__(self):
        return f"{self.__class__.__name__}: {self.from_orb} --> {self.to_orb} // {self.frac}"


class ExcitedState:
    def __init__(
            self,
            oscillator_strength: float,
            excitation_energy: float,
            transitions: list[Transition],
            index: int,
            multiplicity: float,
            symmetry: str,
            manifold_index: int = None,
    ):
        self.manifold_index = manifold_index
        self.symmetry = symmetry
        self.multiplicity = multiplicity
        self.index = index
        self.transitions = sorted(transitions)
        self.excitation_energy = excitation_energy
        self.oscillator_strength = oscillator_strength

    def __repr__(self):
        s = f"Excited State #{self.index}: {self.excitation_energy} eV / {self.nm} nm @ {self.oscillator_strength} @ {self.multiplicity}-{self.manifold_index} @ {self.symmetry}"
        for t in self.transitions:
            s += "\n" + t.__repr__()
        return s

    @property
    def nm(self):
        return ev2nm(self.excitation_energy)

    def __gt__(self, other: ExcitedState):
        return self.excitation_energy > other.excitation_energy

    def __lt__(self, other):
        return not self.__gt__(other)

    @property
    def is_bright(self):
        return self.oscillator_strength > 1e-5


class TdCalc(MSONable):

    def __init__(self, excited_states: list[ExcitedState], properties: dict = None):
        if properties is None:
            properties = dict()
        self.properties = properties
        self.excited_states = excited_states

    @property
    def possible_multiplicities(self) -> list[float]:
        return sorted(set([es.multiplicity for es in self.excited_states]))

    def get_excited_states(
            self, multiplicity: float = None, x_lim_nm: list[float] = None, f_lim: str = "all",
    ) -> list[ExcitedState]:
        states = []
        for es in self.excited_states:
            if multiplicity is not None and es.multiplicity != multiplicity:
                continue
            if x_lim_nm is not None and not x_lim_nm[0] < es.nm < x_lim_nm[1]:
                continue

            if f_lim == "bright" and not es.is_bright:
                continue
            elif f_lim == "dark" and es.is_bright:
                continue
            elif f_lim not in ("all", 'bright', 'dark'):
                raise ValueError(f"f_lim should be one of bright, dark, all, but we have: {f_lim}")

            states.append(es)
        return states

    def get_spectrum_data(self, use_ev=True, multi: float = None, x_lim_nm: list[float] = None,
                          f_lim='all') -> np.ndarray:
        states = self.get_excited_states(multi, x_lim_nm, f_lim=f_lim)
        d = np.zeros((len(states), 2))
        for i, es in enumerate(states):
            if use_ev:
                d[i] = [es.excitation_energy, es.oscillator_strength]
            else:
                d[i] = [es.nm, es.oscillator_strength]
        return d

    @classmethod
    def from_logfile(cls, logfile: FilePath, **kwargs):
        states = []
        manifolds = defaultdict(list)
        with open(logfile, 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if re.search(r"^\sExcited State\s*\d", lines[i]):
                    items = lines[i].strip().split()
                    multi_str, sym = items[3].split("-")
                    try:
                        multi = float(multi_str)
                        multi = int(multi)
                    except ValueError:
                        if multi_str.lower() == "singlet":
                            multi = 1
                        elif multi_str.lower() == 'triplet':
                            multi = 3
                        else:
                            multi = None
                    assert multi is not None, f"cannot resolve multiplicity from {items[3]}"
                    es_index = int(items[2][:-1])
                    manifolds[multi].append(es_index)

                    transitions = []
                    j = i + 1
                    while len(lines[j].strip().split()) == 4 and '->' in lines[j]:
                        frommo, arrow, tomo, coe = lines[j].strip().split()
                        transition = Transition(frommo, tomo, float(coe), index=len(states))
                        transitions.append(transition)
                        j += 1
                    es = ExcitedState(
                        oscillator_strength=float(items[8][2:]),
                        excitation_energy=float(items[4]),
                        transitions=transitions,
                        index=es_index,
                        multiplicity=multi,
                        symmetry=sym,
                        manifold_index=len(manifolds[multi])
                    )
                    states.append(es)
        return cls(states, dict(**kwargs))


def gau_smear(peaks, xs, xe, npts, fwhm):
    d = np.zeros((npts, 2))
    d[:, 0] = np.linspace(xs, xe, npts)
    for i in range(npts):
        yi = 0
        for p in peaks:
            px, py = p
            yi += py * exp(-4 * log(2) * (px - d[i][0]) ** 2 / fwhm ** 2)
        d[i][1] = yi
    return d


def ev2nm(x: Union[np.ndarray, float]):
    if isinstance(x, float):
        return 1239.84193 / x
    elif isinstance(x, np.ndarray):
        assert x.ndim == 1
        return np.array([1239.84193 / xx for xx in x])
    raise TypeError(f"expecting a float or np.ndarray, got: {x.__class__.__name__}")


def nm2ev(x):
    return ev2nm(x)
