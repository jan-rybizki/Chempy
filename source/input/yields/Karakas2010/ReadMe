J/MNRAS/403/1413     Updated stellar yields from AGB models      (Karakas, 2010)
================================================================================
Updated stellar yields from asymptotic giant branch models.
    Karakas A.I.
   <Mon. Not. R. Astron. Soc., 403, 1413-1425 (2010)>
   =2010MNRAS.403.1413K
================================================================================
ADC_Keywords: Models ; Stars, giant ; Abundances
Keywords: nuclear reactions, nucleosynthesis, abundances -
          stars: AGB and post-AGB - stars: Population II - ISM: abundances

Abstract:
    An updated grid of stellar yields for low- to intermediate-mass
    thermally pulsing asymptotic giant branch (AGB) stars is presented.
    The models cover a range in metallicity Z=0.02, 0.008, 0.004 and
    0.0001, and masses between 1 and 6M_{sun}_. New intermediate-mass
    (M>=3M_{sun}_) Z=0.0001 AGB models are also presented, along with a
    finer mass grid than used in previous studies. The yields are computed
    using an updated reaction rate network that includes the latest NeNa
    and MgAl proton capture rates.

File Summary:
--------------------------------------------------------------------------------
 FileName   Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe         80        .   This file
tablea1.dat   241      727   Structural information from the new AGB models
tablea2.dat   123     1232   The stellar yields for the Z = 0.02 model
tablea3.dat   123     1155   The stellar yields for the Z = 0.008 model
tablea4.dat   123     1155   The stellar yields for the Z = 0.004 model
tablea5.dat   123     1232   The stellar yields for the Z = 0.0001 model
tablea6.dat   123      308   The stellar yields for the models with partial
                             mixing zones
--------------------------------------------------------------------------------

See also:
    J/A+AS/107/445 : Envelopes of oxygen-rich AGB stars (Hashimoto, 1994)
    J/ApJS/92/125  : Post-AGB evolution (Vassiliadis+, 1994)
    J/A+AS/126/39  : Models of circumstellar dust shells (Steffen+ 1997)
    J/A+A/512/A10  : Evolution of massive AGB stars. III. (Siess, 2010)

Byte-by-byte Description of file: tablea1.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label    Explanations
--------------------------------------------------------------------------------
   1-  4  F4.2  solMass M0       [1/6] Initial mass
   6- 11  F6.4  ---     Z0       Initial metallicity
                                  (0.0200, 0.0080, 0.0040, 0.0001)
  13- 16  I4    ---     Pulse    Pulse number
  20- 31  E12.7 solMass Mcore    Core mass
  34- 45  E12.7 solMass Mcsh     Maximum mass of the intershell convection zone
  48- 59  E12.7 yr      tcsh     Duration of intershell convection
  62- 73  E12.7 solMass Ddredge  Mass dredged into the envelope
  76- 87  E12.7 ---     lam      Extent of the partial mixing zone, lambda
  90-101  E12.7 ---     lam.dup  {lambda}_dup_ defined by Goriely & Mowlavi
                                  (2000A&A...362..599G) (1)
 104-115  E12.7 K       T(He-sh) Maximum temperature in the He-shell
 118-129  E12.7 K       Tbce     Maximum temperature at the base of the
                                  convective envelope during the previous
                                  interpulse period
 132-143  E12.7 K       T(H-sh)  Maximum temperature in the H-shell during
                                  the previous interpulse period
 146-157  E12.7 yr      iPulse   Interpulse period (note that the first entry
                                  is set to zero)
 160-171  E12.7 solMass Mtot     Total mass
 174-185  E12.7 solLum  Lmax     Maximum radiated luminosity during the
                                  previous interpulse period
 188-199  E12.7 solLum  LHe.max  Maximum He-luminosity during a thermal pulse
 202-213  E12.7 solRad  Rmax     Maximum stellar radius during the previous
                                  interpulse period
 215-227  E13.8 mag     Mbol     Bolometric magnitude
 230-241  E12.7 K       Teff     Effective temperature
--------------------------------------------------------------------------------
Note (1): {lambda}_dup_={lambda}*(DMh/Mcsh) where Mcsh is the maximum mass
     of the intershell convection zone.
--------------------------------------------------------------------------------

Byte-by-byte Description of file: tablea[23456].dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1-  4  F4.2  solMass M0        [1/6.5] Initial mass
   6- 11  F6.4  ---     Z0        Initial metallicity
  13- 17  F5.3  solMass M1        Final mass (1)
  22- 25  A4    ---     El        Species i (2)
  31- 32  I2    ---     A         Atomic mass of the species i
  35- 48  E14.9 solMass Yield     Net yield
  50- 63  E14.9 solMass M(i)lost  Mass of species i lost in the wind
  65- 78  E14.9 solMass M(i)0     Mass of species i initially present in the
                                   wind
  80- 93  E14.9 ---     <X(i)>    Average abundance of species i in the wind
                                   (mass fraction)
  95-108  E14.9 ---     X0(i)     initial mass fraction of species i
                                   (mass fraction)
 110-123  E14.9 ---     f         Production factor, f=log10(<X(i)>/X0(i))
--------------------------------------------------------------------------------
Note (1): For (M0, Z0, M1):
   * 6.00 0.0001 1.008: Reimer's mass loss
   * 6.00 0.0001 1.006: VW93 mass loss
   * 3.00 0.0200 0.682: partial mixing zone = 2x10^-3^M_{sun}_
   * 3.00 0.0080 0.694: partial mixing zone = 2x10^-3^M_{sun}_
   * 3.00 0.0040 0.731: partial mixing zone = 2x10^-3^M_{sun}_
   * 2.00 0.0001 0.685: partial mixing zone = 2x10^-3^M_{sun}_
Note (2): g represents the sum of abundances from ^64^Ni to Bi; an increase
     in g indicates that neutron captures have occurred beyond the end of
     the network.
--------------------------------------------------------------------------------

History:
    From electronic version of the journal

================================================================================
(End)                                      Patricia Vannier [CDS]    04-Oct-2010
