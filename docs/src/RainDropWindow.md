# The raindrop mean-mass window

The rain size distribution is the exponential

```math
n_r(D) = N_0 \exp(-λ_r D), \qquad λ_r ≡ 1/\overline{D}_r,
```

whose two prognostic moments fix its two parameters exactly:

```math
N_r = N_0 \overline{D}_r, \qquad
L_r = π ρ_w N_0 \overline{D}_r^4, \qquad
\overline{x}_r = \frac{L_r}{N_r} = π ρ_w \overline{D}_r^3 .
```

[WackerSeifert2001](@cite) showed that a one- or two-moment scheme develops numerical
artifacts as ``N_r → 0`` and ``L_r → 0``, where ``\overline{x}_r = L_r/N_r`` is ill-defined.
The quantity whose degeneracy is the problem is the mean mass, so that is the quantity to
bound - once, and in one place.
`RainParticlePDF_SB2006_windowed` does exactly that: it evaluates the inversion at the
mean-mass-bounded number, and lets ``λ_r`` and ``N_0`` inherit their ranges through the identities
above. The alternative `RainParticlePDF_SB2006_limited`
clamps ``\overline{x}_r``, ``N_0`` and ``λ_r`` to three separately prescribed ranges, which
cannot in general be satisfied at once: where any of them binds, the returned triple describes
no exponential distribution consistent with ``(L_r, N_r)``.

The window is what `is_limited = true` builds. The cascade is retired: it remains constructible
as `RainParticlePDF_SB2006_limited` so that the comparison arms and the tests characterising it
keep working, and it should not be selected for a run.

This page derives the two ends of the window. Neither is a tuned number: each follows from a
stated physical criterion, and changing the criterion moves the number through the tables below
rather than by editing the parameter.

## What the window is load-bearing for, beyond mean size

The rain-number Jacobian entry depends on it. The self-collection/breakup pair is written as a
relaxation toward the collisional equilibrium number ``n_\text{eq} = L/x_\text{eq}``, and its
frozen-shape diagonal ``-1/τ_\text{eff}`` is non-positive only when the breakup efficiency's sign
follows ``n_r - n_\text{eq}``. That holds when the mean size is an honest function of ``(L, N_r)``,
and it survives the mean-mass clamp for one specific reason: clamping ``\overline{x}_r`` into
``[x_\text{min}, x_\text{max}]`` cannot move it across ``x_\text{eq}``, because

```math
x_\text{min} < x_\text{eq} = \tfrac{π}{6} ρ_w D_\text{eq}^3 < x_\text{max},
```

an inequality between parameters set independently: ``D_\text{eq}`` from the SB2006 breakup fit,
the window ends from the separation mass and the tail-mass criterion. It is asserted in the test
suite rather than left implicit. Measured, 850 states at both precisions: with the shipping
parameters the correspondence never breaks; move the window off ``x_\text{eq}`` and it breaks on
half of them, while a window only 10 % wide that still contains ``x_\text{eq}`` stays clean - so
it is containment that matters, not width.

Under the cascade the correspondence is broken outright, and a sign floor on ``1/τ_\text{eff}``
keeps the Jacobian entry from becoming a growth direction there. The floor costs nothing here: it
never binds on the window or on the unbounded inversion. Where it did bind - on the cascade - the
relaxation contributed no damping at all and the rain-number row's damping was the number
adjustment's alone, a regime that disappears with the cascade.

## The lower end: the cloud/rain separation mass

``\overline{x}_{r,\text{min}} = x_* = 2.6 \times 10^{-10}`` kg is the mass at which
[SeifertBeheng2006](@cite) stops calling a drop a cloud droplet and starts calling it a
raindrop. A raindrop is by definition at least that heavy, so the bound is inherited from the
autoconversion closure rather than chosen here, and fresh rain born at ``x_*`` sits exactly on
it.

## The upper end: a tail-mass criterion at the drop stability limit

Drops larger than ``D_\text{stable} ≈ 8`` mm are hydrodynamically unstable and shatter
spontaneously, so an exponential distribution claiming a substantial share of its mass beyond
``D_\text{stable}`` describes a population that cannot exist. For the exponential PSD the mass
beyond a diameter ``D`` is the Erlang-4 (Gamma(4)) upper tail,

```math
F(z) = e^{-z}\left(1 + z + \tfrac{z^2}{2} + \tfrac{z^3}{6}\right),
\qquad z = \frac{D}{\overline{D}_r},
```

because the mass density is ``∝ D^3 n_r(D)``. The criterion is a choice of how much mass may
sit beyond the stability limit; the mean mass follows.

```@example erlang4
using Printf

const D_stable = 8e-3   # m, hydrodynamic drop stability limit
const ρ_w = 1000.0      # kg/m³

erlang4_tail(z) = exp(-z) * (1 + z + z^2 / 2 + z^3 / 6)
x̄_from_Dm(Dm) = π * ρ_w * Dm^3

println(" D̄ᵣ [mm]   x̄ [kg]      z = D_stable/D̄ᵣ   mass fraction beyond 8 mm")
for Dm_mm in [0.6, 0.8, 1.0, 1.1685, 1.3, 1.45, 1.6, 2.0]
    Dm = Dm_mm * 1e-3
    z = D_stable / Dm
    @printf(" %7.4f   %.4e   %8.3f          %8.4f\n",
        Dm_mm, x̄_from_Dm(Dm), z, erlang4_tail(z))
end
```

Inverting the criterion: for a target tail fraction ``F``, solve ``F(z) = F`` for ``z``, then
``\overline{D}_r = D_\text{stable}/z`` and ``\overline{x}_r = π ρ_w \overline{D}_r^3``.

```@example erlang4
function x̄_max_from_criterion(F; z_lo = 1.0, z_hi = 40.0)
    for _ in 1:200                     # bisection; erlang4_tail is monotone decreasing
        z = (z_lo + z_hi) / 2
        erlang4_tail(z) > F ? (z_lo = z) : (z_hi = z)
    end
    z = (z_lo + z_hi) / 2
    Dm = D_stable / z
    return (; F, z, Dm, x̄ = x̄_from_Dm(Dm))
end

println(" criterion   z       D̄ᵣ [mm]   x̄_max [kg]")
for F in [0.20, 0.10, 0.05, 0.01, 0.001]
    r = x̄_max_from_criterion(F)
    @printf(" %7.3f   %6.3f   %7.4f   %.4e\n", r.F, r.z, r.Dm * 1e3, r.x̄)
end
```

The value `SB2006_raindrops_max_mass` already carries, ``5 \times 10^{-6}`` kg, sits at a tail
fraction of 9.0 %; a round 10 % criterion would put it at ``5.4 \times 10^{-6}`` kg. So the
prescribed number is the "at most about a tenth of the rain mass beyond the stability limit"
criterion to within the precision such a criterion has, and the percentage rather than the mass
becomes the honest knob. A stricter 1 % criterion would tighten the window to about
``1.6 \times 10^{-6}`` kg (mean-volume drop 1.45 mm), which visibly caps heavy-rain fall speeds
and reflectivity - a physical choice, not a numerical one.

## What the retired windows implied

The clamp cascade carries ``λ_r ∈ [10^3, 10^4]`` m``^{-1}`` and
``N_0 ∈ [2.5 \times 10^5, 2 \times 10^7]`` m``^{-4}`` in addition to the mass window. Through
``\overline{x}_r = π ρ_w / λ_r^3`` the ``λ_r`` range is itself a mean-mass window, and it is not
the same one:

```@example erlang4
for (name, λ) in [("λ_max = 1e4", 1e4), ("λ_min = 1e3", 1e3)]
    @printf("%-12s  D̄ᵣ = %6.3f mm   implies x̄ = %.4e kg\n", name, 1e3 / λ, π * ρ_w / λ^3)
end
@printf("%-12s                   window     x̄ ∈ [%.2e, %.2e] kg\n", "mass window", 2.6e-10, 5e-6)
```

Two things travelling with the cascade are not limiting at all, and are named here so they do
not persist silently once the cascade retires.

The first is the intercept range, which has no mean-mass reading, so it is not a third window
on the same quantity; it is an observational plausibility statement about ``N_0``. Under the
mean-mass window it is retired as a clamp and kept as `CloudDiagnostics.rain_intercept_plausibility`,
which reports whether a state's implied ``N_0`` leaves the range and modifies nothing. An
out-of-range intercept is not an error: ``N_0 = λ_r N_r`` grows with drop number at fixed mean
size, so ordinary heavy rain leaves the upper end and sparse large-drop populations leave the
lower end. What the flag identifies is where the cascade *would* have intervened.

The second is a velocity integral. For `SB2006VelType` the individual-drop fit
``v = a_R - b_R e^{-c_R D}`` is negative below the diameter where it crosses zero, so the bulk
moments should be integrated only over the range where it is positive. The unbounded variant does
this; the cascade variant returns untruncated factors instead. The truncated integral is the
implementation, and the windowed variant uses it. The consequence for reading any comparison
between the two variants: a golden on an `SB2006VelType` path has two reasons to move, the
limiting and the truncation, so `Chen2022VelTypeRain` - which production sedimentation uses - is
the single-knob column.

So the scheme carried two floors on the mean mass differing by a factor of 12, and two ceilings
differing by a factor of 1.6, with different consumers reading different ones: breakup,
evaporation and the freezing timescales read ``\overline{x}_r`` while terminal velocity, the
freezing moments and the P3 collision channels read ``\overline{D}_r = 1/λ_r``. The ``N_0``
window is an observational plausibility range (the Marshall-Palmer intercept sits inside it) and
has no mean-mass reading at all; it is what corrupts states that are entirely inside the
sanctioned mass range, because ``N_0`` scales with ``N_r`` at fixed mean mass and so pins on
ordinary heavy rain.
