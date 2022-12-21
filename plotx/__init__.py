import matplotlib.pyplot as plt

from .schema import *

_DarkSingletSymbol = "$S$"
_DarkTripletSymbol = "$T$"


def plotx_one_trace(
        td_calc: TdCalc,
        ax: plt.Axes,
        trace_index: int,
        label: str = "Singlet",
        color: str = "k",
        line_style_vline: str = "solid",
        line_style_smear: str = "solid",
        multiplicity: int = 1,
        show_dark=True,
        smear_fwhm=None,
        smear_label="",
        dark_offset=-0.1,
):
    data_bright_nm = td_calc.get_spectrum_data(use_ev=False, multi=multiplicity, x_lim_nm=None, f_lim='bright')
    data_dark_nm = td_calc.get_spectrum_data(use_ev=False, multi=multiplicity, x_lim_nm=None, f_lim='dark')

    p_vline = ax.vlines(
        data_bright_nm[:, 0], [0], data_bright_nm[:, 1], linestyles=line_style_vline, colors=color, alpha=0.5,
        label=label, linewidth=2
    )
    if multiplicity == 1:
        marker = _DarkSingletSymbol
    else:
        marker = _DarkTripletSymbol
    dark_nm = data_dark_nm.copy()
    dark_nm[:, 1] = dark_nm[:, 1] + dark_offset * (trace_index + 1)

    if show_dark:
        p_dark = ax.scatter(dark_nm[:, 0], dark_nm[:, 1], c=color, marker=marker, alpha=0.5)
        trace_xs = data_dark_nm[:, 0].tolist() + data_bright_nm[:, 0].tolist()
    else:
        trace_xs = data_bright_nm[:, 0].tolist()

    trace_xmin = min(trace_xs)
    trace_xmax = max(trace_xs)
    logger.info(f"trace_index {trace_index}: min={trace_xmin} nm, max={trace_xmax} nm")

    data_bright_ev = td_calc.get_spectrum_data(use_ev=True, multi=multiplicity, x_lim_nm=None, f_lim='bright')
    if smear_fwhm is not None and len(data_bright_ev) > 0:
        smear_ends_multiplier = 3
        smear_start = min(data_bright_ev[:, 0]) - smear_ends_multiplier * smear_fwhm
        smear_end = max(data_bright_ev[:, 0]) + smear_ends_multiplier * smear_fwhm

        smear_start = min([smear_start, nm2ev(trace_xmax)])
        smear_end = max([smear_end, nm2ev(trace_xmin)])

        smear_data = gau_smear(data_bright_ev, smear_start, smear_end, 500, smear_fwhm)
        smear_data[:, 0] = ev2nm(smear_data[:, 0])
        p_smear = ax.plot(smear_data[:, 0], smear_data[:, 1], ls=line_style_smear, label=smear_label, color=color,
                          alpha=0.2)
    # TODO shared legend as in https://stackoverflow.com/questions/31478077/
    return data_bright_nm, data_dark_nm


def plotx_set_ax(ax: plt.Axes):
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Oscillator strength")
    yticks = []
    yticklabels = []
    for t, l in zip(ax.get_yticks(), ax.get_yticklabels()):
        if t >= 0:
            yticks.append(t)
            yticklabels.append(l)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.legend()
