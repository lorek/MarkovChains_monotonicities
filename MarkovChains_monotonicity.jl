module MarkovChains_monotonicity
# Module for checking monotonicities in Markov chains.
# 2017: Pawel Lorek, Piotr Markowski


#
# All functions require order matrix C and transition matrix P
# Consult Lorek, P., Markowski, P. "Monotonicity requirements for efficient exact sampling with Markov chains" 
# for details concerning monotonicities.
#
# Tested in Julia 0.5.0






export 

  is_mobius_down,
  is_mobius_up,
  is_weak_down,
  is_weak_up,
  is_stochastic,
  which_is_not_real




function rationalize_matrix(P::Array{Float64,2})
    s=size(P)[1]
    PP=zeros(Rational{Int64},s,s)    
    for i=1:s
      for j=1:s      
        PP[i,j]=rationalize(P[i,j])
      end
    end
    return PP
end

# multiplying elements of matrix by lcm of their denominators
function convert_to_min_Int(P::Array{Rational{Int64},2})
   return convert(Array{Int64,2},lcm(map(den,P))*P)
end



# convert C into upper_triangular. Adjust P accordingly.
function upper_triangular(C::Array{Int64,2},P::Array{Rational{Int64},2})
  s=size(C)[1]
  c=zeros(Int64,s,2)
  for i=1:s
    c[i,:]=[i,sum(C[i,:])]
  end
  c=sortrows(c,by=x -> x[2],rev=true)
  CC=zeros(Int64,s,s)
  PP=zeros(P)
  for i=1:s
    for j=1:s
      CC[i,j]=C[c[i,1],c[j,1]]
      PP[i,j]=P[c[i,1],c[j,1]]
    end
  end
  return CC,PP
end


# create offspring matrix (consult the article)
function offspring_matrix(C::Array{Int64,2})
  D=deepcopy(C)
  s=size(D)[1]
  PP=zeros(Int64,0,s)
  for i=1:s
    for k=i+1:s
      if D[i,k]==1
        D[i,:]-=D[k,:]
        PP=[PP;sparse([1,1,1],[i,k,s],[1,-1,0])]
      end
    end
  end
  return PP
end

# create all upsets
function up_sets(C::Array{Int64,2})
  s=size(C)[1]
  D=zeros(Int,0,s)
  for i=1:s
    d=deepcopy(C[i,:])'
    for k=(i+1):s if d[1,k]==0
        for j=1:size(d)[1] if d[j,k]==0
            d=[d;d[j,:]'+C[k,:]']
        end end
    end end
    D=[D;d]
  end
  return sign(D)'
end

# check mobius_down monotonicity
function is_mobius_down(C::Array{Int64,2},P::Array{Rational{Int64},2})
  (C,P)=upper_triangular(C,P)
  R=convert_to_min_Int(P)
  return minimum(inv(C)*R*C .>= 0)
end

function is_mobius_down(C::Array{Int64,2},P::Array{Float64,2})
   return is_mobius_down(C,rationalize_matrix(P))
end

# check mobius_up monotonicity
function is_mobius_up(C::Array{Int64,2},P::Array{Rational{Int64},2})
  (C,P)=upper_triangular(C,P)
  R=convert_to_min_Int(P)
  return minimum(inv(C')*R*C' .>= 0)
end

function is_mobius_up(C::Array{Int64,2},P::Array{Float64,2})
   return is_mobius_up(C,rationalize_matrix(P))
end


# check weak_down monotonicity
function is_weak_down(C::Array{Int64,2},P::Array{Rational{Int64},2})
  (C,P)=upper_triangular(C,P)
  R=convert_to_min_Int(P)
  return minimum(offspring_matrix(C)*R*C .>= 0)
end

function is_weak_down(C::Array{Int64,2},P::Array{Float64,2})
   return is_weak_down(C,rationalize_matrix(P))
end



# check weak_up monotonicity
function is_weak_up(C::Array{Int64,2},P::Array{Rational{Int64},2})
  #(C,P)=upper_triangular(C,P)
  R=convert_to_min_Int(P)
  return minimum(offspring_matrix(C)*R*C' .<= 0)
end

function is_weak_up(C::Array{Int64,2},P::Array{Float64,2})
   return is_weak_up(C,rationalize_matrix(P))
end


# check stochastic monotonicity
function is_stochastic(C::Array{Int64,2},P::Array{Rational{Int64},2})
  (C,P)=upper_triangular(C,P)
  R=convert_to_min_Int(P)
  return minimum(offspring_matrix(C)*R*up_sets(C) .<= 0)
end

function is_stochastic(C::Array{Int64,2},P::Array{Float64,2})
   return is_stochastic(C,rationalize_matrix(P))
end



 


 

end
