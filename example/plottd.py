import gauparser as gp


logs = [
    'qpp1_h_tuned_td_singlet.log',
    'qpp1_h_half_tuned_td_singlet.log',
    'qpp1_h_tuned_td_triplet.log',
    'qpp1_h_half_tuned_td_triplet.log',
]

p = gp.ESplt(logs, xsnm=200, xenm=600, fwhm=0.1, npts=500)
p.ax1.set_ylim([-0.5, 3.5])
p.plot('qpp1_h.png')

# logs = [
#     'qpp1_tbu_tuned_td_singlet.log',
#     'qpp1_tbu_half_tuned_td_singlet.log',
#     'qpp1_tbu_tuned_td_triplet.log',
#     'qpp1_tbu_half_tuned_td_triplet.log',
# ]
#
# p = gp.ESplt(logs, xsnm=200, xenm=600, fwhm=0.1, npts=500)
# p.ax1.set_ylim([-0.5, 3.5])
# p.plot('qpp1_tbu.png')
#
# logs = [
#     'qpp2_h_tuned_td_singlet.log',
#     'qpp2_h_half_tuned_td_singlet.log',
#     'qpp2_h_tuned_td_triplet.log',
#     'qpp2_h_half_tuned_td_triplet.log',
# ]
#
# p = gp.ESplt(logs, xsnm=200, xenm=600, fwhm=0.1, npts=500)
# p.ax1.set_ylim([-0.5, 3.5])
# p.plot('qpp2_h.png')
#
# logs = [
#     'qpp2_tbu_tuned_td_singlet.log',
#     'qpp2_tbu_half_tuned_td_singlet.log',
#     'qpp2_tbu_tuned_td_triplet.log',
#     'qpp2_tbu_half_tuned_td_triplet.log',
# ]
#
# p = gp.ESplt(logs, xsnm=200, xenm=600, fwhm=0.1, npts=500)
# p.ax1.set_ylim([-0.5, 3.5])
# p.plot('qpp2_tbu.png')
