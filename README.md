# colaboracion_alumni
The aim of this GitHub is to toy around and familiarize ourselves with the notions of antibody humanization, 
while learning programming and the foundations of biophysics.

Important things to take into account:
- We CANNOT have for than one "main", or else our beautiful .json won't work, and we'd have to do like in Ising, adding each .c file manually and all that (maybe we should consider making a "makefile")
- When using AAcid indexes & properties, MAKE SURE you are using the correct order (maybe we should consider streamlining the idx finding and properties handling when building new functions that need these informations)
- generate_murine_seed is not finished (maybe we could take the schrödinger mouse chain to begin with? yeah that seems like a good idea)
- When using Random, ini_ran should be called first
- Compiler will warn of truncation in string format since it doesn't know string length at compile time in batch_dir, seq_dir, filepath

## Found errors
-

## Notes
### Entropy formulas
- Shannon entropy: $$S=-\sum_ip_i\log p_i$$
- Linear entropy: $$S=\sum_ip_i(1-p_i)$$
- $q$-order Renyi entropy, $q\in[0,1)$: $$S=\frac{1}{1-q}\log\sum_ip_i^q$$
- $q$-order Tsallis entropy, $q\in[0,1)$: $$\frac{1}{q-1}\left(1-\sum_ip_i^q\right)$$

### File formatting
- Chain: each row is a position in the chain with the frequencies for every possible aminoacid and property. Aminoacids and properties are, individually, fundamental famillies of the sample space (aminoacids' add up to 1, and properties' also add up to 1)
    - ``idx`` - position index in the chain
    - ``-AC...Y`` - frequencies for each possible amoniacid in said position
    - ``-(PROP)...CHARGEDMINUS`` - frequencies for each possible property in said position

- Entropies: each row is a position in the chain with every entropy (for aminoacids and for properties)

# 166 without incident (merge conflicts)
0 days without branches