
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

