# Libraries
import numpy as np, pandas as pd
import matplotlib.pyplot as plot

F0_K = 666.7
F0_24 = 7.17

R = F0_24/F0_K

cols = [
    "objid", "ra", "dec", "crowded", "badbkgmatch",
    "m1_fluxflag", "m1_dflag", "m1_extfrac", "m1_brtfrac",
    "m1_snr", "m1_f_psf", "m1_df_psf", "k", "dk",
    "Plx", "RPlx", "PM", "Gmag", "BP-RP", "RUWE", "angDist"
]

DataFrame = pd.read_csv("C:/Users/TaMaM/Downloads/SEIP_data.csv", usecols=cols)

print(len(DataFrame))

DataFrame["Kmag"] = -2.5 * np.log10(
    (DataFrame["k"] * 1e-6) / F0_K
)
DataFrame["MG"] = np.nan
valid_plx = DataFrame["Plx"] > 0

DataFrame.loc[valid_plx, "MG"] = (
    DataFrame.loc[valid_plx, "Gmag"]
    + 5 * np.log10(DataFrame.loc[valid_plx, "Plx"])
    - 10
)

OkData = (DataFrame.k.notna() & (DataFrame.dk < DataFrame.k) & (DataFrame.Kmag > 5)
          & (DataFrame.m1_fluxflag == 0)
          & (DataFrame.m1_extfrac == 0) & (DataFrame.m1_brtfrac == 0)
          & (DataFrame.m1_dflag == "N")
          & (DataFrame.badbkgmatch == 0) & (DataFrame.crowded == 0)
          & (DataFrame.m1_df_psf < 1e6)
          & DataFrame["BP-RP"].notna()
          )

CleanData = DataFrame[OkData].copy()
print(len(CleanData))

CleanData["phot_snr"] = R*CleanData.k / CleanData.m1_df_psf
testable = CleanData[CleanData.phot_snr >= 5]
print(len(testable))

dwarf = testable[(testable.MG > 4.5 ) & (testable.RPlx > 5)].copy()
print(len(dwarf))

dwarf["c"] = 2.5 * np.log10(
    dwarf.m1_f_psf / (R * dwarf.k)
)

dwarf["sc"] = 1.0857 * np.sqrt(
    (dwarf.m1_df_psf / dwarf.m1_f_psf)**2
    + (dwarf.dk / dwarf.k)**2
)

x = dwarf["BP-RP"].values
c = dwarf.c.values
sc = dwarf.sc.values

keep = np.ones(len(dwarf), bool)

for _ in range(10):
    p = np.polyfit(x[keep],c[keep], 3)
    res = c - np.polyval(p,x)
    sig = 1.4826 * np.median(np.abs(res[keep] - np.median(res[keep])))
    keep = (res < 3*sig) & (res > -4*sig)

sig_int = np.sqrt(max(sig**2 - np.median(sc[keep]**2), 0.01*2))
dwarf["locus"] =  np.polyval(p,x)

dwarf["chi"] = (dwarf.c - dwarf.locus) / np.sqrt(dwarf.sc**2 + sig_int**2)
cand = dwarf[dwarf.chi > 3].sort_values("chi", ascending=False)
print(len(cand), (dwarf.chi > 5).sum())     # I got 22 and 13

print((dwarf.chi < -3).sum())   # stars with a significant *negative* excess

phot_flux = R*dwarf.k*10**(0.4*dwarf.locus)           # predicted photosphere flux at 24 µm
base = dwarf.chi <= 3                                  # only test stars not already flagged
for f in [0.1, 0.2, 0.3, 0.5, 1.0, 2.0]:               # excess as a fraction of the photosphere
    F = dwarf.m1_f_psf + f*phot_flux
    chi = (2.5*np.log10(F/(R*dwarf.k)) - dwarf.locus) / np.sqrt(dwarf.sc**2 + sig_int**2)
    print(f"{f:4.0%} excess -> recovered {np.mean(chi[base] > 3):.0%}")


# Vibrant astronomy-style graphs
plot.style.use("dark_background")

# Graph 1: Chi distribution
fig, ax = plot.subplots(figsize=(10, 6))

ax.hist(
    dwarf.chi.clip(-8, 15),
    bins=60,
    color="#38BDF8",
    edgecolor="#0F172A",
    alpha=0.9,
    log=True
)

ax.axvline(
    3, color="#FF4D6D", linestyle="--",
    linewidth=2.5, label="Candidate threshold (χ = 3)"
)

ax.set_title("24 µm Excess Significance Distribution",
             fontsize=16, fontweight="bold", pad=14)
ax.set_xlabel("Excess significance (χ)", fontsize=12)
ax.set_ylabel("Number of stars (log scale)", fontsize=12)
ax.grid(True, alpha=0.18, linestyle="--")
ax.legend(frameon=False)
fig.tight_layout()


# Graph 2: Color vs. BP-RP, highlighting candidates
fig, ax = plot.subplots(figsize=(11, 7))

ax.scatter(
    dwarf["BP-RP"], dwarf.c,
    s=13, color="#38BDF8", alpha=0.55,
    edgecolors="none", label="Tested stars", rasterized=True
)

ax.scatter(
    cand["BP-RP"], cand.c,
    s=55, color="#FF4D6D", alpha=0.95,
    edgecolors="white", linewidths=0.5,
    label=f"Excess candidates (n={len(cand)})",
    zorder=3
)

# Photospheric color-locus fit
x_line = np.linspace(dwarf["BP-RP"].min(),
                     dwarf["BP-RP"].max(), 300)
ax.plot(
    x_line, np.polyval(p, x_line),
    color="#FACC15", linewidth=2.5,
    label="Fitted photospheric locus", zorder=4
)

ax.set_ylim(-1, 3)
ax.set_title("Stellar Color vs. Mid-Infrared Excess",
             fontsize=16, fontweight="bold", pad=14)
ax.set_xlabel("Gaia BP − RP color (mag)", fontsize=12)
ax.set_ylabel("24 µm color (mag)", fontsize=12)
ax.grid(True, alpha=0.18, linestyle="--")
ax.legend(frameon=False, loc="best")
fig.tight_layout()

plot.show()
