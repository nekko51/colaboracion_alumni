# colaboracion_alumni
The aim of this GitHub is to delve deeper into biophysics, using simulated annealing and science data to humanize murine antibodies.
---
167 days without incident (merge conflicts)

## Important things to take into account:
- We **CANNOT** have more than one "main", or else our beautiful .json won't work, and we'd have to add each .c file manually and all that (maybe we should consider making a "makefile")
- When using AAcid indexes & properties, **MAKE SURE you are using the correct order** (maybe we should consider streamlining the idx finding and properties handling when building new functions that need these informations)
- When using Random(), ini_ran should be called first or else you'll always get a **0**
- Compiler will warn of truncation in string format since it doesn't know string length at compile time in ``batch_dir``, ``seq_dir``, ``filepath``

## File formatting
- **Chain**: each row is a position in the chain with the frequencies for every possible aminoacid and property. Aminoacids and properties are, individually, fundamental famillies of the sample space (aminoacids' add up to 1, and properties' also add up to 1)
    - ``idx`` - position index in the chain
    - ``-AC...Y`` - frequencies for each possible amoniacid in said position
    - ``-(PROP)...CHARGEDMINUS`` - frequencies for each possible property in said position

- **Entropies**: each row is a position in the chain with every entropy (for aminoacids and for properties)

- **Metropolis**:
Organised in batch folders named by date and time, for every seed a folder is generated, inside of which a .txt is generated for each metropolis run, which contains a report summarizing the simulation results, consisting of:
    1.  **Header** - s line indicating the report version.
        ```
        ***************Metropolis data report version X.Y.Z***************
        ```

    2.  **Global Parameters**: Parameters of the entire simulation run.
        *   ``weight_log``, ``weight_properties`` - weights used for energy calculation.
        *   ``n_betas``- total number of beta values simulated.
        *   ``n_sweeps``- number of Metropolis sweeps performed for each beta.

    3.  **Overall Statistics**:
        *   ``overall acceptance med`` - mean of the mean acceptance rates from each beta simulation.
        *   ``overall acceptance var`` - variance of the mean acceptance rates.

    4.  **Initial State**: Information about the starting sequence before any simulation steps.
        *   ``log humanness``, ``prop humanness``, ``wanderer penalty``, ``total energy`` - energy components of the initial sequence.
        *   ``seed`` - initial murine sequence seed.

    5.  **Per-Beta Results**: This section is repeated for each beta value.
        *   A header line showing the current ``beta`` value, its index, the mean acceptance rate for that beta (``acceptance med``), and the variance of the acceptance rate (``acceptance var``).
        *   ``delta_e``-  change in total energy between the result of the previous beta and the current one.
        *   ``log humanness``, ``prop humanness``, ``wanderer penalty``, ``total energy``- energy components of the final sequence for this beta.
        *   ``resulting chain`` - the final, (hopefully) more humane sequence after all sweeps at this beta are complete.

    6.  **Footer** - a line marking the end of the report.
        ```
        ***************Metropolis data report end***************
        ```

## Formulas
### Entropies:
- Shannon entropy: $$S=-\sum_ip_i\log p_i$$
- Linear entropy: $$S=\sum_ip_i(1-p_i)$$
- $q$-order Renyi entropy, $q\in[0,1)$: $$S=\frac{1}{1-q}\log\sum_ip_i^q$$
- $q$-order Tsallis entropy, $q\in[0,1)$: $$\frac{1}{q-1}\left(1-\sum_ip_i^q\right)$$


## Found errors
-
