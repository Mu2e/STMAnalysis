import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.stats import norm, lognorm, chisquare, ks_2samp

debug=False
save_plots=True
show_conservative=True # set to False if you don't want to see the conservative calculation

#
# Parameters
# - proton beam
running_time = 2*60*60 # second # np.linspace(10,2400,240)
beam_off_time = 1*60*60 # beam off at this time
time_unit = 1 # second
#N_POT_per_second = 5.71e12 # 2BB
N_POT_per_second = 3e12 # 1BB
on_spill_duty_factor = 0.323

# - muons
mu_stop_rate = 0.0007 # mu-stop / POT
f_cap = 0.609

# - efficiencies / acceptances
geom_acceptance = 1e-9
ge_eff = 0.30

# - photons
I_1809 = 0.51
lifetime_1809 = 864e-9 # s
I_844 = 0.13*0.718
half_life_844 = 9.5*60 # s
lifetime_844 = half_life_844 / np.log(2) # log = natural log

def plot_time_lines(ax):
    ypos = (ax.get_ylim()[1]-ax.get_ylim()[0])/2
    ax.axvline(60, color='black', linestyle='dashed')
    ax.text(x=60+2, y=ypos, s="1 min", rotation=-90, fontsize=12)
    ax.axvline(10*60, color='black', linestyle='dashed')
    ax.text(x=10*60+2, y=ypos, s="10 min", rotation=-90, fontsize=12)
    ax.axvline(30*60, color='black', linestyle='dashed')
    ax.text(x=30*60+2, y=ypos, s="30 min", rotation=-90, fontsize=12)
    ax.axvline(60*60, color='black', linestyle='dashed')
    ax.text(x=60*60+2, y=ypos, s="60 min", rotation=-90, fontsize=12)


#
# Generate data
#
n_times = np.round(running_time / time_unit)+1
times = np.linspace(0, running_time, num=int(n_times))

n_times_beam_on = np.round((beam_off_time)/time_unit)+1
beam_on_off = np.ones(int(n_times_beam_on))
beam_on_off = np.append(beam_on_off, np.zeros(int(n_times - n_times_beam_on)))

mu_captures_per_second = N_POT_per_second*mu_stop_rate*f_cap*beam_on_off*time_unit
N_1809_states_per_second =  mu_captures_per_second*I_1809
N_844_states_per_second = mu_captures_per_second*I_844


# We will fill these
total_mu_captures = np.zeros(len(times))

N_1809_emitted_per_second = np.zeros(len(times))
N_1809_detected_per_second = np.zeros(len(times))
N_1809_emitted = np.zeros(int(n_times))
N_1809_detected = np.zeros(int(n_times))

N_844_emitted_per_second = np.zeros(len(times))
N_844_detected_per_second = np.zeros(len(times))
N_844_emitted = np.zeros(int(n_times))
N_844_detected = np.zeros(int(n_times))

#
# First go through and calculate how many photons will be emitted each second
# Note: we need to know the number of muons that were captured in previous second
#
for i,time in enumerate(times):
    if debug:
        print(f"t = {time} s: ")
    for j,prev_time in enumerate(times[times<=time]):
        total_mu_captures[i] = total_mu_captures[i] + mu_captures_per_second[j]
        if debug:
            print(f"\tPrev 844 keV emitted / second = {N_844_emitted_per_second[i]}")
            print(f"\tStates created at t = {prev_time}: {N_844_states_per_second[j]}")
            print(f"\t\tFraction still surviving: {np.exp(-(time-prev_time)/lifetime_844)}")
            print(f"\t\tFraction still sirviving one second later: {np.exp(-(time-prev_time+1)/lifetime_844)}")
            print(f"\t\t844 keV emitted this second = {N_844_states_per_second[j]*(np.exp(-(time-prev_time)/lifetime_844) - np.exp(-(time-prev_time+time_unit)/lifetime_844))}")
        N_1809_emitted_per_second[i] = N_1809_emitted_per_second[i] + N_1809_states_per_second[j]*(np.exp(-(time-prev_time)/lifetime_1809) - np.exp(-(time-prev_time+time_unit)/lifetime_1809))
        N_844_emitted_per_second[i] = N_844_emitted_per_second[i] + N_844_states_per_second[j]*(np.exp(-(time-prev_time)/lifetime_844) - np.exp(-(time-prev_time+time_unit)/lifetime_844))


    N_1809_detected_per_second[i] = N_1809_emitted_per_second[i]*geom_acceptance*ge_eff
    N_844_detected_per_second[i] = N_844_emitted_per_second[i]*geom_acceptance*ge_eff

#
# Now go through and sum up the total we will have emitted and detected
# Note: consider three scenarios for the 844 keV line
#
N_1809_emitted_total = np.zeros(len(times))
N_844_emitted_total = np.zeros(len(times))
N_1809_detected_total = np.zeros(len(times))
N_844_detected_total = np.zeros(len(times))
N_844_conservative_total = np.zeros(len(times))
for i,time in enumerate(times):
    for j,prev_time in enumerate(times[times<=time]):
        N_1809_emitted_total[i] = N_1809_emitted_total[i] + N_1809_emitted_per_second[j]
        N_844_emitted_total[i] = N_844_emitted_total[i] + N_844_emitted_per_second[j]
        N_1809_detected_total[i] = N_1809_detected_total[i] + N_1809_detected_per_second[j]
        N_844_detected_total[i] = N_844_detected_total[i] + N_844_detected_per_second[j]
        if (prev_time > beam_off_time):
            N_844_conservative_total[i] = N_844_conservative_total[i] + N_844_detected_per_second[j]
        else:
            N_844_conservative_total[i] = N_844_conservative_total[i] + (N_844_detected_per_second[j])*(1-on_spill_duty_factor)


#
# Finally, calculate the statistical uncertianties
#
N_1809_detected_error = np.sqrt(N_1809_detected_total)
N_1809_stat_error_pct = (N_1809_detected_error / N_1809_detected_total) * 100

N_844_detected_error = np.sqrt(N_844_detected_total)
N_844_stat_error_pct = (N_844_detected_error / N_844_detected_total) * 100

N_844_conservative_error = np.sqrt(N_844_conservative_total)
N_844_conservative_stat_error_pct = (N_844_conservative_error / N_844_conservative_total) * 100


# Fig.1 - muon captures
fig, axs = plt.subplots(2, 1, figsize=(16,9), tight_layout=True, sharex=True)
fig.suptitle("Muon Captures", fontsize=16)
axs[0].plot(times, mu_captures_per_second, label=f"$N_{{\\mu-cap}}$ = {N_POT_per_second:.2e} POT / s * {mu_stop_rate:.3e} $\\mu_{{stop}}$ / POT * {f_cap:.3f} $\\mu_{{cap}}$ / $\\mu_{{stop}}$")
axs[1].plot(times, total_mu_captures, label='Muon Captures')
axs[0].grid(True)
axs[1].grid(True)
axs[1].set_xlabel("Time [s]", fontsize=14)
axs[0].set_ylabel("Muon Captures per Second", fontsize=14)
axs[1].set_ylabel("Total Muon Captures", fontsize=14)
axs[0].set_ylim(bottom=0)
axs[1].set_ylim(bottom=0)
ypos = (axs[0].get_ylim()[1]-axs[0].get_ylim()[0])*0.8
axs[0].hlines(y=ypos, xmin=0, xmax=beam_off_time, color='red')
axs[0].text((beam_off_time/2.), ypos*0.95, 'beam on', color='red', fontsize=14)
for ax in axs:
    plot_time_lines(ax)
    ax.legend(fontsize=14)

if save_plots:
    fig.savefig('muon_captures.pdf')


# Fig.2 - update with the emitted photons
fig, axs = plt.subplots(2, 1, figsize=(16,9), tight_layout=True, sharex=True)
fig.suptitle("Muon Captures and Emitted Photons", fontsize=16)
axs[0].plot(times, mu_captures_per_second, label=f"$N_{{\\mu-cap}}$ = {N_POT_per_second:.2e} POT / s * {mu_stop_rate:.3e} $\\mu_{{stop}}$ / POT * {f_cap:.3f} $\\mu_{{cap}}$ / $\\mu_{{stop}}$")
axs[0].plot(times, N_1809_emitted_per_second, label=f"$N_{{1809}}$ = $N_{{\\mu-cap}}$ * {I_1809*100:.0f}% (intensity) ($\\tau$ = {lifetime_1809*1e9:.0f} ns)")
axs[0].plot(times, N_844_emitted_per_second, label=f"$N_{{844}}$ = $N_{{\\mu-cap}}$ * {I_844*100:.2f}% (intensity) ($\\tau$ = {lifetime_844/60:.1f} min)")

axs[1].plot(times, total_mu_captures, label='Muon Captures')
axs[1].plot(times, N_1809_emitted_total, label=f"1809 keV")
axs[1].plot(times, N_844_emitted_total, label=f"844 keV")

axs[0].grid(True)
axs[1].grid(True)
axs[1].set_xlabel("Time [s]", fontsize=14)
axs[0].set_ylabel("N / Second", fontsize=14)
axs[1].set_ylabel("Total N", fontsize=14)
axs[0].set_ylim(bottom=0)
axs[1].set_ylim(bottom=0)
for ax in axs:
    plot_time_lines(ax)
    ax.legend(fontsize=14)

xpos = running_time*0.75 # find a good position for the text
ypos = N_844_emitted_per_second[(times>xpos-time_unit) & (times<xpos+time_unit)][0] # [0] because could get a few times
axs[0].text(xpos, ypos*2, '844 keV photons still being emitted', color='red', fontsize=14, horizontalalignment='center')
if save_plots:
    fig.savefig('muon_captures_and_photons.pdf')

# Fig.3 - photons detected
fig, axs = plt.subplots(2, 1, figsize=(16,9), tight_layout=True, sharex=True)
fig.suptitle("Detected Photons", fontsize=16)
axs[0].plot(times, N_1809_detected_per_second, label=f"$N_{{1809-det}}$ = $N_{{1809}}$ * {geom_acceptance:.1e} (geom. acc) * {ge_eff:.2f} (det. eff))  ($\\tau$ = {lifetime_1809*1e9:.0f} ns)")
axs[0].plot(times, N_844_detected_per_second, label=f"$N_{{844-det}}$ = $N_{{844}}$ * {geom_acceptance:.1e} (geom. acc) * {ge_eff:.2f} (det. eff))  ($\\tau$ = {lifetime_844/60:.1f} min)")

axs[1].plot(times, N_1809_detected_total, label=f"1809 keV")
axs[1].plot(times, N_844_detected_total, label=f"844 keV")

axs[0].grid(True)
axs[1].grid(True)
axs[1].set_xlabel("Time [s]", fontsize=14)
axs[0].set_ylabel("N / Second", fontsize=14)
axs[1].set_ylabel("Total N", fontsize=14)
axs[0].set_ylim(bottom=0)
axs[1].set_ylim(bottom=0)
for ax in axs:
    plot_time_lines(ax)
    ax.legend(fontsize=14)

if save_plots:
    fig.savefig('photons_detected.pdf')


# Fig.4 - photons detected (scenarios)
fig, axs = plt.subplots(2, 1, figsize=(16,9), tight_layout=True, sharex=True)
fig.suptitle("Statistical Uncertainties", fontsize=16)
axs[0].plot(times, N_1809_detected_total, label=f"1809 keV")
if show_conservative:
    axs[0].plot(times, N_844_detected_total, label=f"844 keV (optimistic = OnSpill+OffSpill+NoBeam)")
    axs[0].plot(times, N_844_conservative_total, label=f"844 keV (conservative = OffSpill+NoBeam)")
else:
    axs[0].plot(times, N_844_detected_total, label=f"844 keV")

axs[1].plot(times, N_1809_stat_error_pct, label='1809 keV')
axs[1].plot(times, N_844_stat_error_pct, label='844 keV')
axs[1].set_xlabel("Time [s]", fontsize=14)
axs[0].set_ylabel("Total N", fontsize=14)
axs[1].set_ylabel("Stat. Uncertainty [%]", fontsize=14)
axs[0].grid(True)
axs[1].grid(True)
axs[0].set_ylim(bottom=0)
axs[1].set_ylim(0, 25)
axs[1].axhline(10, color='red', linestyle='dashed', label='10% stat. uncertainty')

reaches_10pct_1809 = N_1809_stat_error_pct[-1]<10.0
if reaches_10pct_1809:
    time_to_10pct_1809=times[(N_1809_stat_error_pct>9.99) & (N_1809_stat_error_pct<10.01)][0]
    axs[1].annotate(f"1809 keV reaches\n10% stat. uncertainty\nin $\\sim{time_to_10pct_1809/60:.0f}$ min", xy=(time_to_10pct_1809, 10.0), xytext=(time_to_10pct_1809, 8.0), arrowprops=dict(color='red', arrowstyle='->', lw=2), color='red', horizontalalignment='center', verticalalignment='top', fontsize=14)
else:
    plateau_1809=N_1809_stat_error_pct[-1]*np.ones(int(n_times))
    difference_1809=N_1809_stat_error_pct - plateau_1809
    time_to_plateau_1809=times[(difference_1809>0.1)][-1]
    axs[1].axhline(plateau_1809[0], color='gray', linestyle='dashed', label=f"1809 keV, best uncertainty = {plateau_1809[0]:.1f}%")
    axs[1].annotate(f"844 keV reaches\nplateau\nin $\\sim{time_to_plateau_1809/60:.0f}$ min", xy=(time_to_plateau_1809, plateau_1809[0]), xytext=(time_to_plateau_1809, 7.5), arrowprops=dict(color='red', arrowstyle='->', lw=2), color='red', horizontalalignment='center', verticalalignment='top', fontsize=14)

reaches_10pct_844 = N_844_stat_error_pct[-1]<10.0
if reaches_10pct_844:
    time_to_10pct_844=times[(N_844_stat_error_pct>9.99) & (N_844_stat_error_pct<10.01)][0]
    axs[1].annotate(f"844 keV reaches\n10% stat. uncertainty\nin $\\sim{time_to_10pct_844/60:.0f}$ min", xy=(time_to_10pct_844, 10.0), xytext=(time_to_10pct_844, 8.0), arrowprops=dict(color='red', arrowstyle='->', lw=2), color='red', horizontalalignment='center', verticalalignment='top', fontsize=14)
else:
    plateau_844=N_844_stat_error_pct[-1]*np.ones(int(n_times))
    difference_844=N_844_stat_error_pct - plateau_844
    time_to_plateau_844=times[(difference_844>0.1)][-1]
    axs[1].axhline(plateau_844[0], color='gray', linestyle='dashed', label=f"844 keV, best uncertainty = {plateau_844[0]:.1f}%")
    axs[1].annotate(f"844 keV reaches\nplateau\nin $\\sim{time_to_plateau_844/60:.0f}$ min", xy=(time_to_plateau_844, plateau_844[0]), xytext=(time_to_plateau_844, 7.5), arrowprops=dict(color='red', arrowstyle='->', lw=2), color='red', horizontalalignment='center', verticalalignment='top', fontsize=14)


if show_conservative:
    axs[1].plot(times, N_844_conservative_stat_error_pct, label='844 keV (conservative)')
    reaches_10pct_844_conservative = N_844_conservative_stat_error_pct[-1]<10.0
    if reaches_10pct_844_conservative:
        time_to_10pct_844_conservative=times[(N_844_conservative_stat_error_pct>9.99) & (N_844_conservative_stat_error_pct<10.01)][0]
        axs[1].annotate(f"844 keV reaches\n10% stat. uncertainty\nin $\\sim{time_to_10pct_844_conservative/60:.0f}$ min", xy=(time_to_10pct_844_conservative, 10.0), xytext=(time_to_10pct_844_conservative, 8.0), arrowprops=dict(color='red', arrowstyle='->', lw=2), color='red', horizontalalignment='center', verticalalignment='top', fontsize=14)
    else:
        plateau_844_conservative=N_844_conservative_stat_error_pct[-1]*np.ones(int(n_times))
        difference_844_conservative=N_844_conservative_stat_error_pct - plateau_844_conservative
        time_to_plateau_844_conservative=times[(difference_844_conservative>0.1)][-1]
        axs[1].axhline(plateau_844_conservative[0], color='gray', linestyle='dashed', label=f"844_conservative keV, best uncertainty = {plateau_844_conservative[0]:.1f}%")
        axs[1].annotate(f"844 keV reaches\nplateau\nin $\\sim{time_to_plateau_844_conservative/60:.0f}$ min\n(conservative)", xy=(time_to_plateau_844_conservative, plateau_844_conservative[0]), xytext=(time_to_plateau_844_conservative, 7.5), arrowprops=dict(color='red', arrowstyle='->', lw=2), color='red', horizontalalignment='center', verticalalignment='top', fontsize=14)


axs[1].legend(fontsize=12)
for ax in axs:
    plot_time_lines(ax)
    ax.legend(fontsize=14)

if save_plots:
    fig.savefig('stat_uncertainty.pdf')


plt.show()
exit(0)
