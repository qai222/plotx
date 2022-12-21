import os

from plotx import *


def compare_pa_and_apa(
        pa_logfile: FilePath,
        apa_logfile: FilePath,
        ax: plt.Axes,
        color: str,
        xlim: list[float] = None,
):
    for logfile, line_style in zip([pa_logfile, apa_logfile], ['solid', 'dotted']):
        plotx_one_trace(
            TdCalc.from_logfile(logfile),
            ax, trace_index=1, label=os.path.basename(logfile).replace(".log", ""), color=color,
            line_style_vline=line_style, multiplicity=1, show_dark=False, smear_fwhm=None
        )
    ax.set_xlim(xlim)


fig, ax = plt.subplots()
colors = sns.color_palette(n_colors=10)
i = 0
for pa_log, apa_log in [
    ["apa/p3.log", "apa/3p3.log", ],
    ["apa/p4.log", "apa/4p4.log", ],
    ["apa/p5.log", "apa/5p5.log", ],
]:
    compare_pa_and_apa(
        pa_log, apa_log, ax, colors[i],
        xlim=[365, 423]
    )
    i += 1
ax.set_ylim([-0.02, 0.15])
ax.legend()
ax.set_xlabel("Wavelength (nm)")
ax.set_ylabel("Oscillator strength")
ax.set_title("Singlet (td uLC-w*HPBE/def2tzvp)")
fig.savefig("apa.png", dpi=600)
