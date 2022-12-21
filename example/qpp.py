import matplotlib.pyplot as plt

from plotx import *

sns.set_style('darkgrid')

if __name__ == '__main__':
    fig, ax = plt.subplots()
    qpp1_half = TdCalc(
        TdCalc.from_logfile("qpp/qpp1_h_half_tuned_td_singlet.log").excited_states +
        TdCalc.from_logfile("qpp/qpp1_h_half_tuned_td_triplet.log").excited_states
    )

    qpp1 = TdCalc(
        TdCalc.from_logfile("qpp/qpp1_h_tuned_td_singlet.log").excited_states +
        TdCalc.from_logfile("qpp/qpp1_h_tuned_td_triplet.log").excited_states
    )

    plotx_one_trace(
        qpp1_half, ax, trace_index=1, label="QPP1-H Singlet", color="blue", line_style_vline="solid", multiplicity=1,
        show_dark=True,
        smear_fwhm=0.1, dark_offset=-0.2
    )
    plotx_one_trace(
        qpp1_half, ax, trace_index=1, label="QPP1-H Triplet", color="navy", line_style_vline="solid", multiplicity=3,
        show_dark=True,
        smear_fwhm=0.1, dark_offset=-0.2
    )
    plotx_one_trace(
        qpp1, ax, trace_index=2, label="QPP1-H Singlet", color="red", line_style_vline="solid", multiplicity=1,
        show_dark=True,
        smear_fwhm=0.1, dark_offset=-0.2
    )
    plotx_one_trace(
        qpp1, ax, trace_index=3, label="QPP1-H Triplet", color="darkred", line_style_vline="solid", multiplicity=3,
        show_dark=True,
        smear_fwhm=0.1, dark_offset=-0.2
    )

    plotx_set_ax(ax)
    ax2 = ax.secondary_xaxis('top', functions=(lambda x: 1239.84193/x, lambda x: 1239.84193/x))
    ax2.set_xlabel('Excitation energy (eV)')

    fig.savefig("qpp.png", dpi=600)
