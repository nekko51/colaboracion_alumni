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
Organised in batch folders named by date and time. For every seed, a folder is generated containing the results for each metropolis run. Each run produces two files:
    1. **Main Report (`run_N.txt`)**: A summary of the simulation results.
    2. **Beta Matrix File (`run_N_betas.txt`)**: A detailed dump of the position-specific beta values used in the simulation.

    **Main Report File Formatting:**
    The main report file consists of:
    1.  **Header**: A line indicating the report version.
        ```
        ***************Metropolis data report version X.Y.Z***************
        ```

    2.  **Global Parameters**: Parameters for the entire simulation run.
        *   ``weight_log``, ``weight_properties``, ``weight_penalty`` - weights used for energy calculation.
        *   ``n_betas``- total number of beta steps simulated.
        *   ``n_sweeps``- number of Metropolis sweeps performed for each beta step.

    3.  **Overall Statistics**:
        *   ``overall acceptance med`` - mean of the median acceptance rates from each beta step.
        *   ``overall acceptance var`` - variance of the median acceptance rates.

    4.  **Initial State**: Information about the starting sequence before any simulation steps.
        *   ``log humanness``, ``prop humanness``, ``wanderer penalty``, ``total energy`` - energy components of the initial sequence.
        *   ``seed`` - initial murine seed.

    5.  **Per-Beta Results**: This section is repeated for each beta step.
        *   A header line showing a summary of the beta values for the current step (``beta (min/avg/max)``), its index, the median acceptance rate for that beta (``acceptance med``), and the variance of the acceptance rate (``acceptance var``).
        *   ``delta_e`` -  change in total energy between the result of the previous beta step and the current one.
        *   ``log humanness``, ``prop humanness``, ``wanderer penalty``, ``total energy``- energy components of the final sequence for this beta step.
        *   ``resulting chain`` - the final, (hopefully) more humane sequence after all sweeps at this beta are complete.

    6.  **Footer**: A line marking the end of the report.
        ```
        ***************Metropolis data report end***************
        ```
    
    **Beta Matrix File Formatting:**
    This file (`run_N_betas.txt`) contains the raw `n_betas` x `CHAINLEN` matrix of beta values used, allowing for detailed analysis and debugging.

## Formulas
### Entropies:
- Shannon entropy: $$S=-\sum_ip_i\log p_i$$
- Linear entropy: $$S=\sum_ip_i(1-p_i)$$
- $q$-order Renyi entropy, $q\in[0,1)$: $$S=\frac{1}{1-q}\log\sum_ip_i^q$$
- $q$-order Tsallis entropy, $q\in[0,1)$: $$\frac{1}{q-1}\left(1-\sum_ip_i^q\right)$$

### Correlation:
Mutual information is given by:
$$
I(X;Y) = \sum_{x\in\mathcal X}\sum_{y\in\mathcal Y}p(x, y)\log_2\frac{p(x,y)}{p(x)p(y)},
$$
where $X$ and $Y$ are random variables that take values, $x$ and $y$, in $\mathcal X$ and $\mathcal Y$; and $p(x)=P(\{X=x\}),\ p(y)=P(\{Y=y\}),\ p(x,y)=P(\{X=x\}\cap \{Y=y\})$.

In ``correlations.c`` it is normalized as follows:
$$
i(X;Y) = \frac{I(X;Y)}{\min\{H(X),H(Y)\}},
$$
where $H(X)$ and $H(Y)$ are the Shannon entropies of the variables $X$ and $Y$.

When $\min\{H(X), H(Y)\} = 0$, the position is fully conserved and $I(X;Y)=0$ as well. In this case, the code outputs $i(X;Y) \equiv 0$, by convention.

### Metropolis:
- $\beta$-values are defined to be variable, through two different mechanisms (what I call the engine and the damper)
    - **Engine**: $$\displaystyle\beta_i^j = k_i \cdot \frac{\alpha}{S_j + \varepsilon}, \quad k_i = \left(c_r\right)^i\ \forall i \in [1,N_\beta], j \in [1, L_\text{chain}]$$ where $c_r$ is the cooling rate, $\varepsilon$ avoids division by 0, $\alpha$ is a scale factor, and $S_j$ is entropy for position $j$ in Schrödinger chain.
    - **Damper**: $$\varphi(\sigma, \lambda) = \frac{L}{\displaystyle 1+(1-L)e^{-\lambda \sigma}}$$ where $\lambda$ is a global variable (the dampening coefficient), $L$ is a global variable too (it's the maximum value our $\beta$ will be multiplied by), and $\sigma$ is the deviation of the last sweep acceptance from the target acceptance. Since we want acceptance to decay as our system explores the energy landscape, we begin with a ``STARTING_TARGET_ACCEPTANCE`` value, $\alpha_0$, (defined in head.h), and slowly decrease it using the following equation: $$\alpha_i = \alpha_0 \left(1-\frac{i}{N_\beta}\right)$$ where $N_\beta$ is the number of betas we have, and $i$ is the current iteration with our pre-defined $\beta$ given by the "engine"
    $\newline$ Lastly, $\sigma = \alpha_{ij}-\alpha_i$, where $\alpha_{ij}$ is acceptance from last sweep iteration $(j\in[1,\text{n\_sweeps}])$

- The total energy, $E_\text{total}$, is a weighted sum of three components, which we aim to minimize:
    $$E_\text{total} = w_\text{log} E_\text{log} + w_\text{prop} E_\text{prop} + w_\text{penalty} E_\text{penalty}$$
    The weights for each component (`WEIGHT_LOG`, `WEIGHT_PROP`, `WEIGHT_PENALTY`) are defined in `head.h`.

    1.  **Log Humanness Energy ($E_\text{log}$)**: This term quantifies how "human-like" a sequence is based on AA frequencies from a reference set of human sequences (the Schrödinger chain):
        $$E_\text{log} = -\sum_{i=1}^{L_\text{chain}} \log(p_i(a_i))$$
        where $p_i(a_i)$ is the frequency of the amino acid $a_i$ at position $i$ in the human reference. A large penalty (`ZERO_FREQ_PENALTY_LOG`) is added for amino acids with a frequency below a small threshold (`EPSILON`).

    2.  **Property Distance Energy ($E_\text{prop}$)**: This term measures the similarity between the properties of the sequence and the human reference; it's calculated as the squared Euclidean distance in the property space for each position:
        $$E_\text{prop} = \sum_{i=1}^{L_\text{chain}} \sum_{j=1}^{N_\text{props}} \left( \text{prop}_j(a_i) - \overline{\text{prop}_j(i)} \right)^2$$
        where $\text{prop}_j(a_i)$ is the value of property $j$ for amino acid $a_i$, and $\overline{\text{prop}_j(i)}$ is the average value of property $j$ at position $i$ in the Schrödinger chain.

    3.  **Wander(er/ing) Penalty Energy ($E_\text{penalty}$)**: This term discourages the sequence from deviating too far from the original murine sequence; we thus only incentivize mutation if there's a relatively big gain in humanness. It consists of a penalty based on property distance and a (for now) constant penalty for any mutation:
        $$E_\text{penalty} = \sum_{i=1}^{L_\text{chain}} \left[ \left( \sum_{j=1}^{N_\text{props}} (\text{prop}_j(a_i) - \text{prop}_j(a_i^\text{original}))^2 \right) + P \cdot \delta_{\displaystyle  a_i}^{\displaystyle a_i^{\text{original}}}\right]$$
        where $a_i^\text{original}$ is the amino acid at position $i$ in the original murine sequence, $P$ is the `AA_MUTATION_PENALTY`, and $\delta$ is the "generalized Kronecker delta" (logical if gate). This approach ensures each mutation has a minimum energy cost, while heavily penalizing big property changes.


## Found errors
-
