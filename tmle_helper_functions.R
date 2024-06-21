# function to compute gradient for finding mles in 1df MOR test
grad_func_MOR = function(u, var0, var1, var2, m0, m1, m2)
{
  u0 = u[1]; u1 = u[2]
  partial_u0 = (1/var0)*(m0 - u0) + (1/var1)*(m1 - u0*u1)*u1 + (1/var2)*(m2 - u0*u1^2)*u1^2
  partial_u1 = (1/var1)*(m1 - u0*u1)*u0 + (1/var2)*(m2 - u0*u1^2)*(2*u0*u1)
  
  res = c(partial_u0, partial_u1)
  return(res)
}

# function for hessian for 1df MOR MLE solutions
MOR_mle_fisher_info = function(u, var0, var1, var2, m0, m1, m2)
{
  u0 = u[1]; u1 = u[2]
  
  # first generate hessian
  partial_u0_2 = -(1/var0) -(1/var1)*u1^2 - (1/var2)*u1^4
  partial_u1_2 = -(1/var1)*u0^2 + (2/var2)*(m2*u0 - 3*u0^2*u1^2)
  partial_u0_u1 = (1/var1)*(m1 - 2*u0*u1) + (1/var2)* (2*m2*u1 - 4*u0*u1^3)
  
  hessian = cbind(c(partial_u0_2,partial_u0_u1), c(partial_u0_u1, partial_u1_2)) 
  fisher_info = -1*solve(hessian)
  return(fisher_info)
}

MOR_mle_jacobian = function(u, var0, var1, var2, m0, m1, m2)
{
  u0 = u[1]; u1 = u[2]
  
  # first generate hessian
  partial_u0_2 = -(1/var0) -(1/var1)*u1^2 - (1/var2)*u1^4
  partial_u1_2 = -(1/var1)*u0^2 + (2/var2)*(m2*u0 - 3*u0^2*u1^2)
  partial_u0_u1 = (1/var1)*(m1 - 2*u0*u1) + (1/var2)* (2*m2*u1 - 4*u0*u1^3)
  
  jacobian = cbind(c(partial_u0_2,partial_u0_u1), c(partial_u0_u1, partial_u1_2)) 
  return(jacobian)
}

# rebound a set of values such that a <= values <= b
rebound = function(values,a,b)
{
  c1 = min(values)
  c2 = max(values)
  
  rebounded_values = ((values - c1) / (c2 - c1))*(b-a) + a 
  return(rebounded_values)
}

# rebound a set of values such that anything less than a is mapped to a 
# and anything greater than b is mapped to b 
bound_limits = function(values,a,b)
{
  values[values<a] = a
  values[values>b] = b
  
  return(values)
}

