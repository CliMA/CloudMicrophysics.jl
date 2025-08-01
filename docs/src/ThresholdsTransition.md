# Smooth transition at thresholds

To avoid abrupt change of values at thresholds we can use the logistic function to make the transition smooth. The logistic function is defined as
```math
f(t, k) = \frac{1}{1+e^{-kt}}
```
where
- ``t`` is the independent parameter,
- and ``k`` is the growth rate of the curve characterizing the steepness of the transition.

The value of the logistic function changes from zero at ``t \rightarrow -\infty`` to one at ``t \rightarrow \infty`` with the majority of the change happening around ``t = 0``. In our microphysics applications, the independent parameter (typically specific humidities) takes only non-negative values, and when this parameter is zero the function should return zero. Thus, we use the change of variable ``t = x/x_0 - x_0/x``, where ``x_0`` is the threshold value of ``x``. Therefore, the logistic function for smooth transitioning from ``f(x) = 0``, for ``x < x_0``, to ``f(x)=1``, for ``x > x_0``, is defined as
```math
f(x, x_0, k) = \frac{1}{1+e^{-k(x/x_0 - x_0/x)}}
```
Note that when both ``x`` and ``x_0`` are zero the value given by the above equation is undefined. In this case we return a zero value instead.

## Smooth transition of derivative at thresholds

When the function itself is continuous but its derivative changes abruptly at a threshold we can use the indefinite integral of the logistic function for smooth transitioning. In this case, the following function can be used

```math
f(x, x_0, k) = \frac{x_0}{k} \ln\left(1+e^{k(x/x_0 - 1 + a_{trnslt})} \right) - a_{trnslt} x_0
```
where ``a`` is a fixed value that is defined to enforce zero at ``x=0``:
```math
a_{trnslt} = -\frac{1}{k}\ln\left(1 - e^{-k}\right)
```
The curve defined by the above equation smoothly transition from ``f(x) = 0``, for ``x < x_0``, to ``f(x)=x-x_0``, for ``x > x_0``. Note that when both ``x`` and ``x_0`` are zero the value of the function is undefined. In this case we return a zero value instead.

## Example figures

```@example
include("plots/Thersholds_transitions.jl")
```
![](q_lcl_K1969.svg)
![](q_lcl_TC1980_LD2004.svg)
![](N_d_B1994.svg)
