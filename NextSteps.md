<!-- # Next steps (LP - Low priority; MP - Medium Priority; HP - High Priority; UHP - Critical Priority)
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


# Next steps for C code:

- Add energy penalization for each AA change, that way we make the global minimum different for each murine seed (heavily penalizing complete AA changes and marginally penalizing AA change that's close in the properties sample space) -- PROBLEM: given this approach, energy will be heavily linked to n_steps -- is that a problem though?
- Finish work on mega_metropolis
- Start work on giga_metropolis
- Add parallel threading
- The code can be optimized HEAVILY (in principle) if instead of returning a Chain, etc. in structs we pass a pointer and fill it out with our functions; when passing Chains/big structs to functions (bigger than 8 bytes, which is what a pointer occupies in memory), we can also just pass them by reference to improve execution speed (remember we can also make const pointers).
Furthermore, we should delve in the differences between heap and stack memory (this last one seems to be pretty limited, and we might run out of it in execution) (https://www.reddit.com/r/cs2a/comments/1egymrr/stack_vs_heap_memory/)
- The start of mega_metropolis could be optimized with realloc, but it's 00 already
           