# Searching for 24 μm Excess Around Low-Luminosity Dwarfs

## Research question
Can we detect debris disk around nearby low-luminosity dwarfs (likely older) dwarf stars, and can we apply a neural network to improve on a simple color-excess search?

## Backgroung
Debris disks are rings/belts of dust, ice, rock, and small organic objects like asteroids and comets orbiting mature stars, for that are formed as leftover material from the planet formation process and other things like collisions. These objects constantly collide with each other, which grinds them down into smaller and smaller pieces of debris, eventually becoming a fine dust.

In our context, Low-luminosity stars are hard to test: their photosphere are dim at 24 μm, so only bright excess stands out, and because of this their ages are hard to determine, so young stars with brighter primordial disks can be mistaken for old debris disks.

## Method (baseline)
1. Use Spitzer SEIP sources matched to 2MASS and Gaia DR3 (17,464 rows)
2. Apply quality cuts for needed columns from the SEIP DATA (flags, crowding, saturation), which left us with 5,544 stars.
3. Keep stars whose photosphere would be detectable at 24 μm (S/N ≥ 5).
4. Select low-luminosity dwarfs (M_G > 4.5, parallax S/N > 5): (Found: 288 stars).
5. Compute Ks - [24], fit the normal-star color locus against Gaia BP-RP, and flag stars with excess significance χ > 3.
6. Measure detection limits by injecting fake excesses into real fluxes.

## Results
- 22 of 288 dwarfs (7.6%) show χ > 3; 13 show χ > 5.
- 9 of the 22 lie in your star-forming regions and are likely primordial disks.
- Detection Limit: ~50% at + 50% excess, excess, ~98% at +100%.

  ## Limitations
  - The input table only contains 24 μm detection, so only stars with detectable photospheres were tested.
  - No stellar ages; low luminosity is only one factor to determining an "old" dwarf.
  - Canditates have not been checked in image cutouts for potential background galaxies.


## TODO
-- Create a neural network to see if a NN can recover fainter excess than the χ cut.

## DATA
SEIP Source List from IRSA

## Run it 
[todo]

