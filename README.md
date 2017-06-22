# MarkovChains_monotonicities
Module containing functions for checking various monotonicities in Markov chains.
Details are available in the article

Pawel Lorek and Piotr Markowski (2017).  "Monotonicity requirements for efficient exact sampling with Markov chains". To appear in Markov Processes and Related Fields. Available at http://www.math.uni.wroc.pl/~lorek/papers/

## Overview
`MarkovChains_monotonicity.jl` contains functions for checking several monotonicities of Markov chains. All the monotonicities are defined with respect to the partial ordering of the state space. All the functions require transition matrix `P` and 0-1 valued ordering matrix `C`. For any two states **e** and **e'**  we have `C`(**e**,**e'**)=1 if **e** is less or equal to **e'** and 0 otherwise.


| Monotonicity| Function |
| ------ | ------ |
|  Weak up monotonicity  |  ````is_weak_up(C,P)  ```` |
|  Weak up monotonicity  |  ````is_weak_down(C,P) ```` |
|  Mobius up monotonicity  |  ````is_mobius_up(C,P)```` |
|  Mobius down monotonicity  |  ````is_mobius_down(C,P)```` |
|  Stochastic monotonicity  |  ````is_stochastic(C,P)```` |


## Usage

Loading the module (from current directory)

````julia
# Julia 0.5.0 does not search for modules 
# in current path by default. You may be forced to use first the following:
julia> push!(LOAD_PATH,pwd());  

julia> using MarkovChains_monotonicity
````

**Sample usgage**

Symmetric random walk on 2-dimensional cube (states (0,0), (1,0), (0,1), (1,1),
coordinate-wise ordering)

````julia
julia> julia> P=[1//3 1//3 1//3 0; 1//3 1//3 0 1//3; 1//3 0 1//3 1//3; 0 1//3 1//3 1//3];
julia> C=[1 1 1 1; 0 1 0 1; 0 0 1 1; 0 0 0 1];
julia> is_stochastic(C,P2)
true                                                                                                                        
julia> is_mobius_up(C,P2)
false                                                                                                                        
````
### Examples of chains and ordering from the article
A file `MarkovChains_monotonicity.jl` checks all the implemented monotonicities of 19 chains from 
[the article](http://www.math.uni.wroc.pl/~lorek/papers/)


