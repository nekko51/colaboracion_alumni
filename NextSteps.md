# Next steps for C code:

## General:
### Optimization;
- In principle, the code could be **heavily** optimized if instead of returning Chain, Energy as outputs of functions (or anything bigger than 8 bytes, which is what a pointer occupies in memory), we pass them a pointer and make them fill in relevant data (basically turning all functions into void or int)
- For functions that take a const datatype as an input, we can just pass them by reference (remember pointers can also be made of const type) to prevent memory duplicates
- Furthermore, we should delve in the differences between heap and stack memory (this last one seems to be pretty limited, and we might run out of it in execution) (https://www.reddit.com/r/cs2a/comments/1egymrr/stack_vs_heap_memory/)
- Implement parallel threading

## Functions:
### mega_metropolis;
- The start could be optimized with realloc (instead of opening file twice)
- Prevent gap changing??

### giga_metropolis conceptualization;
- Variable betas for each position implementation (betas would now be a matrix of dimensions `n_betas*CHAINLEN`) where $\beta^i_j = k_i \cdot \frac{c}{S_j+\varepsilon}$, where $S_j$ is entropy for each chain position (ranging from `0` to `CHAINLEN-1`) and $k_i$ is a monotonous sequence (sucesión monótona, no sé cómo lo dicen los ingleses) which is function of (¿?), that dictates how `betas[i][]` should change as thermalization is reached for `betas[i-1][]`
- Prevent gap changing??


---
## Completed next steps:
- Add energy penalization for each AA change, that way we make the global minimum different for each murine seed (heavily penalizing complete AA changes and marginally penalizing AA change that's close in the properties sample space) -- PROBLEM: given this approach, energy will be heavily linked to n_steps -- is that a problem though?
- Finish work on mega_metropolis

















<!-- # Next steps for Python code (LP - Low priority; MP - Medium Priority; HP - High Priority; UHP - Critical Priority)
- [UHP] Extract the position of frequencies that are over a certain value (e.g. 0,4) to compare the binary results; maybe go back to the idea of making subplot every 25-50 positions that was discarded?
Maybe, to compare, parallel to making the binary plot, it'd be nice to graph JUST the (aromatics, for example), to see which change.

- [VISUAL - LP] Make it so the name of the aminoacid shows up in its correspondig point on the graph if its frequency is higher than some value (e.g. 0.5)

- DISCARDED - [DATA - MP] Maybe make sub-plots (every 25-50 positions, make a new graph) to compare more easily

- [HP] Come up with more ways to compare the two plots.

- [TOYING AROUND - LP] Try poking around and tweaking the sequences so that there exists an offset between the data (don't know what that'd be good for, but I'll just throw it in the list i guess)

- [MP] Maybe rethink the functions so that we can work with the "frequency" and "appearance" arrays (e.g. to generate various plots with different colormaps efficiently, without calculating them for every plot)




# Ways to implement them:
- Implement a generic filename_output directory and then append the images output names

 -->

