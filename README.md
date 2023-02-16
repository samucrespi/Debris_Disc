# Debris Disc
Study of how the post-collisional debris distribute in a forming planetary system.

## Phase 1: Select Collisions

To investigate the distribution of debris in a forming system, I will simulate post-collisional debris and analyze their evolution. Several factors can influence the debris population, including the number and mass distribution of major bodies in the vicinity of the collision, as well as resonances with Jupiter and Saturn.

The number and mass distribution of protoplanets surrounding the collision depend strongly on the system's age. After the gas disappears, the disk consists of a bimodal population of planetary embryos and planetesimals. Approximately one billion years later, the system comprises several fully-formed terrestrial planets.

I will analyze collisions occurring in four temporal bins: before $10^6$ years, between $10^6$ and $10^7$ years, between $10^7$ and $10^8$ years, and after $10^9$ years. For each temporal bin, I will consider one collision in the inner part of the system (inside 1 AU), where no strong resonances with Jupiter and Saturn are expected, and one collision in the outer part of the system (outside 1 AU), where strong resonances are present. Additionally, I will only consider collisions that are as close as possible to those present in my SPH catalog, to minimize interpolation effects. To keep the integration time feasible while maintaining statistical consistency, I will select collisions that produce approximately 1000 debris.

The selected collsions, which are obtained from the FwoD simulations (see Crespi et al. 2023, in preparation), are the following:
- FwoD2/9-35: Ndeb=3626 rcoll=0.72  tcol=9.78e+05 v=3.44  a=19.3 mtot=2.40e-02 g=1.00
- FwoD2/7-28: Ndeb=268 rcoll=2.90  tcol=6.01e+05 v=1.32  a=39.4 mtot=2.33e-02 g=1.00
- FwoD2/9-43: Ndeb=3626 rcoll=0.92  tcol=1.82e+06 v=3.29  a=20.9 mtot=2.41e-02 g=0.99
- FwoD2/0-36: Ndeb=4618 rcoll=2.46  tcol=1.33e+06 v=4.97  a=39.0 mtot=2.48e-02 g=1.00
- FwoD1/7-120: Ndeb=4748 rcoll=0.67  tcol=1.07e+07 v=4.58  a=39.1 mtot=1.76e-02 g=0.50
- FwoD2/5-133: Ndeb=115 rcoll=1.96  tcol=5.21e+07 v=1.56  a=58.5 mtot=2.28e-01 g=0.06
- FwoD2/1-191: Ndeb=2130 rcoll=0.72  tcol=4.21e+08 v=3.14  a=42.4 mtot=2.59e-01 g=0.08
