
# Sample file for module  MarkovChains_monotonicity
# We check monotonicites of matrices P1 ... P20 from article
# Lorek, P., Markowski, P. "Monotonicity requirements for exact sampling with Markov chains" 


# Remark
# Julia 0.5.0 does not search for modules 
# in current path by default. You may be forced to use first the following:
push!(LOAD_PATH, pwd())


using MarkovChains_monotonicity;



# Check weak_up, weak_down, mobius_up, mobius_down and stochastic monotonicity for given C and P
function check_all_monotonicities(C::Array{Int64,2},P::Union{Array{Rational{Int64},2},Array{Float64,2}})
   println("\t Weak up: \t ", is_weak_up(C,P));
   println("\t Weak down: \t ", is_weak_down(C,P));
   println("\t Mobius up: \t ", is_mobius_up(C,P));
   println("\t Mobius down: \t ", is_mobius_down(C,P));
   println("\t Stochastic: \t ", is_stochastic(C,P)); 
end




# general partial orderings
C1=[1 1 1 1 1 1; 0 1 0 1 1 1; 0 0 1 1 1 1; 0 0 0 1 0 1; 0 0 0 0 1 1; 0 0 0 0 0 1]
C2=[1 1 1 1 1 1 1 1; 0 1 0 0 1 1 0 1; 0 0 1 0 1 0 1 1; 0 0 0 1 0 1 1 1; 0 0 0 0 1 0 0 1; 0 0 0 0 0 1 0 1; 0 0 0 0 0 0 1 1; 0 0 0 0 0 0 0 1]
C3=[1 1 1 1 1; 0 1 0 0 1; 0 0 1 0 1; 0 0 0 1 1; 0 0 0 0 1]

P1=[3//8 1//8 1//8 1//8 0 1//8 1//8 0; 1//4 1//4 0 1//8 1//8 0 1//8 1//8; 1//4 0 3//8 0 0 1//4 0 1//8; 1//4 0 1//8 0 1//8 1//4 1//8 1//8; 0 1//8 1//8 0 0 1//4 1//4 1//4; 0 0 1//4 1//4 0 1//8 0 3//8; 0 0 1//8 1//8 1//4 1//4 1//8 1//8; 0 0 1//8 0 1//8 1//8 0 5//8]
P2=[5//8 0 1//8 1//8 0 1//8 0 0; 1//8 1//8 1//4 1//4 1//8 1//8 0 0; 3//8 0 1//8 0 1//4 1//4 0 0; 1//4 1//4 1//4 0 0 1//8 1//8 0; 1//8 1//8 1//4 1//8 0 1//8 0 1//4; 1//8 0 1//4 0 0 3//8 0 1//4; 1//8 1//8 0 1//8 1//8 0 1//4 1//4; 0 1//8 1//8 0 1//8 1//8 1//8 3//8]
P3=[1//2 1//6 0 1//3 0 0; 1//3 1//6 1//6 1//3 0 0; 1//3 1//6 0 1//3 1//6 0; 1//6 1//6 1//6 0 1//6 1//3; 1//6 0 1//6 1//6 1//3 1//6; 0 1//6 1//3 0 1//6 1//3]
P4=[2//5 1//5 1//5 1//5 0; 2//5 1//5 1//5 1//5 0; 0 2//5 2//5 1//5 0; 0 2//5 1//5 2//5 0; 0 1//5 2//5 0 2//5]
P5=[2//5 0 2//5 1//5 0; 0 2//5 1//5 2//5 0; 0 1//5 2//5 2//5 0; 0 1//5 1//5 1//5 2//5; 0 1//5 1//5 1//5 2//5]
P6=[17//24 0 0 1//8 1//8 1//24; 1//8 5//16 5//16 1//12 1//12 1//12; 1//8 5//16 5//16 1//12 1//12 1//12; 1//12 1//12 1//12 5//16 5//16 1//8; 1//12 1//12 1//12 5//16 5//16 1//8; 1//24 1//8 1//8 0 0 17//24]
P7=[1//3 1//3 1//3 0 0 0; 1//3 1//3 0 1//3 0 0; 1//3 0 1//3 1//3 0 0; 0 1//3 1//3 1//3 0 0; 0 1//6 1//6 1//6 1//6 1//3; 0 1//6 1//6 1//6 1//6 1//3]
P8=[2//5 1//5 1//5 1//5 0; 2//5 1//5 1//5 1//5 0; 2//5 0 1//5 1//5 1//5; 1//5 1//5 2//5 1//5 0; 0 2//5 1//5 0 2//5]
P9=[2//5 0 1//5 2//5 0; 0 1//5 2//5 1//5 1//5; 1//5 1//5 1//5 0 2//5; 0 1//5 1//5 1//5 2//5; 0 1//5 1//5 1//5 2//5]
P11=[1//3 1//6 1//6 1//6 1//6 0; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 0 1//6 1//6 1//6 1//6 1//3]
P12=[1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 0 1//3 1//6 1//6 1//6 1//6]
P13=[1//6 1//6 1//6 1//6 1//3 0; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6]
P14=[1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6; 1//6 1//6 1//6 1//6 1//6 1//6]
P15=[17//24 0 0 1//8 1//8 1//24; 1//8 5//16 5//16 1//12 1//12 1//12; 1//8 5//16 5//16 1//12 1//12 1//12; 1//12 1//12 1//12 5//16 5//16 1//8; 1//12 1//12 1//12 5//16 5//16 1//8; 1//24 1//16 1//16 1//16 1//16 17//24]
P16=[17//24 1//16 1//16 1//16 1//16 1//24; 1//8 5//16 5//16 1//12 1//12 1//12; 1//8 5//16 5//16 1//12 1//12 1//12; 1//12 1//12 1//12 5//16 5//16 1//8; 1//12 1//12 1//12 5//16 5//16 1//8; 1//24 1//8 1//8 0 0 17//24]


# tree-like ordering 

C4=[1 0 0 0 1 0 1; 0 1 0 0 1 0 1; 0 0 1 0 0 1 1; 0 0 0 1 0 1 1; 0 0 0 0 1 0 1; 0 0 0 0 0 1 1 ; 0 0 0 0 0 0 1]

P17=[1//2 1//2 0 0 0 0 0; 0 0 1//2 1//2 0 0 0; 1//2 1//2 0 0 0 0 0; 0 0 1//2 1//2 0 0 0;   0 0 0 0 1//2 1//2 0;  0 0 0 0 1//2 1//2 0;  0 0 0 0 1//2 1//2 0]
P18=[1//7 1//7 1//7 1//7 1//7 1//7 1//7 ;1//7 1//7 1//7 1//7 1//7 1//7 1//7 ;1//7 1//7 1//7 1//7 1//7 1//7 1//7 ;1//7 1//7 1//7 1//7 1//7 1//7 1//7 ;1//7 1//7 1//7 1//7 1//7 1//7 1//7 ;1//7 1//7 1//7 1//7 1//7 1//7 1//7 ;1//7 1//7 1//7 1//7 1//7 1//7 1//7 ]
P19=eye(7)
P20=[1//2 1//2 0 0 0 0 0; 0 0 1//2 1//2 0 0 0; 1//2 1//2 0 0 0 0 0; 0 0 1//2 1//2 0 0 0;   0 0 0 0 1//2 1//2 0;  0 0 0 0 1//2 1//2 0;  0 0 0 0 0 0 1]


 
 
println("\nChecking Example 1 (C2,P1) ..");
check_all_monotonicities(C2,P1)

println("\nChecking Example 2 (C2,P2) ..");
check_all_monotonicities(C2,P2)


println("\nChecking Example 3 (C1,P3) ..");
check_all_monotonicities(C1,P3)

println("\nChecking Example 4 (C3,P4) ..");
check_all_monotonicities(C3,P4)

println("\nChecking Example 5 (C3,P5) ..");
check_all_monotonicities(C3,P5)

println("\nChecking Example 6 (C1,P6) ..");
check_all_monotonicities(C1,P6)


println("\nChecking Example 7 (C1,P7) ..");
check_all_monotonicities(C1,P7)


println("\nChecking Example 8 (C3,P8) ..");
check_all_monotonicities(C3,P8)

println("\nChecking Example 9 (C3,P9) ..");
check_all_monotonicities(C3,P9)



println("\nChecking Example 11 (C1,P11) ..");
check_all_monotonicities(C1,P11)

println("\nChecking Example 12 (C1,P12) ..");
check_all_monotonicities(C1,P12)

println("\nChecking Example 13 (C1,P13) ..");
check_all_monotonicities(C1,P13)

println("\nChecking Example 14 (C1,P14) ..");
check_all_monotonicities(C1,P14)

println("\nChecking Example 15 (C1,P15) ..");
check_all_monotonicities(C1,P15)

println("\nChecking Example 16 (C1,P16) ..");
check_all_monotonicities(C1,P16)

println("\nChecking Example 17 (C4,P17) ..");
check_all_monotonicities(C4,P17)

println("\nChecking Example 18 (C4,P18) ..");
check_all_monotonicities(C4,P18)

println("\nChecking Example 19 (C4,P19) ..");
check_all_monotonicities(C4,P19)

println("\nChecking Example 20 (C4,P20) ..");
check_all_monotonicities(C4,P20)


