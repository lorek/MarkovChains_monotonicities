# MarkovChains_monotonicities
# Module contains functions for checking various monotonicities in Markov chains.
# Details available in article

Pawel Lorek and Piotr Markowski (2017).  "Monotonicity requirements for efficient exact sampling with Markov chains". To appear in Markov Processes and Related Fields. Available at http://www.math.uni.wroc.pl/~lorek/papers


MarkovChains_monotonicity.jl contains function for checking the following monotonicity of a Markov chain with transition matrix P with respect to the partial ordering defined in matrix C:

-- Weak up monotonicity [is_weak_up(C,P)]

-- Weak down monotonicity [is_weak_down(C,P)]
-- Mobius up monotonicity [is_mobius_up(C,P)]
-- Mobius down monotonicity [is_mobius_down(C,P)]
-- Stochastic monotonicity [is_stochastic(C,P)]


