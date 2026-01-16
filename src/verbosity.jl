"""
    DDEVerbosity <: AbstractVerbositySpecifier

Verbosity configuration for DelayDiffEq.jl solvers, providing fine-grained control over
diagnostic messages, warnings, and errors during DDE solution.

# Fields

## ODE Verbosity
- `ode_verbosity`: Verbosity configuration for the underlying ODE solver

## Delay-Specific Group
- `discontinuity_tracking`: Messages about discontinuity propagation tracking
- `delay_evaluation`: Messages about delay term evaluation
- `constrained_step`: Messages when step size is constrained by discontinuities
- `residual_control`: Messages about residual control in implicit methods
- `neutral_delay`: Messages specific to neutral delay equations
- `state_dependent_delay`: Messages when state-dependent delays are detected/evaluated

# Constructors

    DDEVerbosity(preset::AbstractVerbosityPreset)

Create a `DDEVerbosity` using a preset configuration:
- `SciMLLogging.None()`: All messages disabled
- `SciMLLogging.Minimal()`: Only critical errors and fatal issues
- `SciMLLogging.Standard()`: Balanced verbosity (default)
- `SciMLLogging.Detailed()`: Comprehensive debugging information
- `SciMLLogging.All()`: Maximum verbosity

    DDEVerbosity(; ode_verbosity=nothing, delay_specific=nothing, kwargs...)

Create a `DDEVerbosity` with group-level or individual field control.

# Examples

```julia
# Use a preset
verbose = DDEVerbosity(SciMLLogging.Standard())

# Set ODE verbosity and delay-specific group
verbose = DDEVerbosity(
    ode_verbosity = ODEVerbosity(SciMLLogging.Detailed()),
    delay_specific = SciMLLogging.InfoLevel()
)

# Set individual fields
verbose = DDEVerbosity(
    discontinuity_tracking = SciMLLogging.InfoLevel(),
    delay_evaluation = SciMLLogging.WarnLevel()
)

# Mix group and individual settings
verbose = DDEVerbosity(
    delay_specific = SciMLLogging.InfoLevel(),  # Set all delay-specific to InfoLevel
    state_dependent_delay = SciMLLogging.WarnLevel()  # Override specific field
)
```
"""
function DDEVerbosity end 

@verbosity_specifier DDEVerbosity begin
    toggles = (
        :ode_verbosity,
        :discontinuity_tracking,
        :delay_evaluation,
        :constrained_step,
        :residual_control,
        :neutral_delay,
        :state_dependent_delay
    )

    presets = (
        None = (
            ode_verbosity = ODEVerbosity(None()),
            discontinuity_tracking = Silent(),
            delay_evaluation = Silent(),
            constrained_step = Silent(),
            residual_control = Silent(),
            neutral_delay = Silent(),
            state_dependent_delay = Silent()
        ),
        Minimal = (
            ode_verbosity = ODEVerbosity(Minimal()),
            discontinuity_tracking = Silent(),
            delay_evaluation = Silent(),
            constrained_step = Silent(),
            residual_control = Silent(),
            neutral_delay = Silent(),
            state_dependent_delay = WarnLevel()
        ),
        Standard = (
            ode_verbosity = ODEVerbosity(),
            discontinuity_tracking = Silent(),
            delay_evaluation = Silent(),
            constrained_step = Silent(),
            residual_control = Silent(),
            neutral_delay = Silent(),
            state_dependent_delay = Silent()
        ),
        Detailed = (
            ode_verbosity = ODEVerbosity(Detailed()),
            discontinuity_tracking = InfoLevel(),
            delay_evaluation = InfoLevel(),
            constrained_step = InfoLevel(),
            residual_control = InfoLevel(),
            neutral_delay = InfoLevel(),
            state_dependent_delay = WarnLevel()
        ),
        All = (
            ode_verbosity = ODEVerbosity(All()),
            discontinuity_tracking = InfoLevel(),
            delay_evaluation = InfoLevel(),
            constrained_step = InfoLevel(),
            residual_control = InfoLevel(),
            neutral_delay = InfoLevel(),
            state_dependent_delay = InfoLevel()
        )
    )

    groups = (
        delay_specific = (
            :discontinuity_tracking,
            :delay_evaluation,
            :constrained_step,
            :residual_control,
            :neutral_delay,
            :state_dependent_delay
        ),
    )
end

function Base.getproperty(v::DDEVerbosity, s::Symbol)
    # For parametric types, we need to use the base type
    DDE_fields = fieldnames(typeof(v).name.wrapper)
    if s in DDE_fields
        return getfield(v, s)
    else
        # Try to delegate to ODE verbosity
        ode_v = getfield(v, :ode_verbosity)
        ODE_fields = fieldnames(typeof(ode_v).name.wrapper)
        if s in ODE_fields
            return getfield(ode_v, s)
        else
            return error("type DDEVerbosity has no field ", s)
        end
    end
end