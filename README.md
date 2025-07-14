# PFDHA Code

This repository includes the PFDHA program used in Liou and Abrahamson (2025), Non-Ergodic Probabilistic Fault-Displacement Hazard
Analysis.

line 1: Flag for P(surface rupture |M) model
0: P(SR) = 1
1: simple model
2: based on thickness (not working yet)
3: Wells and Coppersmith (1993)

line 2: Flag for rupture length model
1: Lavrentidias et al (2023)

line 3: file name for the SSC inputs
This is the out2 file from Haz45 for the rates of magnitudes

line 4: output file 1 name: mean hazard

line 5: output file 2 name: epistemic fractiles

line 6: output file 3 name: disaggregation

line 7: N_z for hazard calc
number of displacements for computing the hazard

line 8: displacements (m) for hazard calc

line 9: X/L and Rx (km)
X/L defines the location of the site relative to the end of the rupture
Rx is the distance off the fault (measured perp to strike).
Set Rx=0 for primary rupture

line 10: Number of SOF for the fault (aleatory)
This is a number between 1-3.

line 11: SOF values
0 = SS
1 = RV
-1 = NML

line 12: SOF wts (aleatory)
wts for the SOF

line 13: flt length (km), seismothickness (km), dip (degrees)
the seismogenic thickness is the vertical thickness (not down-dip width)

line 14: MagBinStep, number of mag bins, Mag_0
This sets the mag bin size for the disaggregation
the Mag_0 is the starting magnitude for the first bin

line 15: EpsBinStep, number of esp bins, eps_0
This sets the epsilon bin size for the disaggregation
the eps_0 is the starting epsilon for the first bin

line 16: X/L BinStep, number of X/L bins, X/L_0
This sets the X/L bin size for the disaggregation
the X/L_0 is the starting epsilon for the first bin

line 17: number of branches for the scaled backbone FDMs

line 18: list of the adjustment terms for the scaled backbone FDM
set of delta D values (in m^0.3 units) that are added to the FDM
- currently, only the LA23 FDM is included
  
line 19: wts for the branches of the backbone FDM

line 20: sigma_flag
sigma_flag = 0: use the sigma from the FDM
sigma_flag = 1: use user defined sigma in place of the sigma from the FDM

line 21: user defined sigma values
These are in m^(0.3) units
There are 4 sigma values to enter: fix_sigma_Dagg, fix_sigma_DP,
fix_sigma_Dagg_prime, fix_sigma_DP_prime
Currently, the program only uses the &quot;fix_sigma_DP_prime&quot;
For now, Enter the same values for all four values
