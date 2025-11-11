"""
    DDEVerbosity <: AbstractVerbositySpecifier

Verbosity configuration for DelayDiffEq.jl solvers, providing fine-grained control over
diagnostic messages, warnings, and errors during DDE solution.

# Fields

## ODE Verbosity
- `ode_verbosity`: Verbosity configuration for the underlying ODE solver

## Delay-Specific Group
- `discontinuity_tracking`: Messages about discontinuity propagation tracking
- `history_interpolation`: Messages about history function interpolation
- `history_extrapolation`: Messages when history function extrapolates beyond bounds
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
    history_interpolation = SciMLLogging.WarnLevel()
)

# Mix group and individual settings
verbose = DDEVerbosity(
    delay_specific = SciMLLogging.InfoLevel(),  # Set all delay-specific to InfoLevel
    state_dependent_delay = SciMLLogging.WarnLevel()  # Override specific field
)
```
"""
@concrete struct DDEVerbosity <: AbstractVerbositySpecifier
    # ODE solver verbosity
    ode_verbosity
    # Delay-specific options
    discontinuity_tracking
    history_interpolation
    history_extrapolation
    delay_evaluation
    constrained_step
    residual_control
    neutral_delay
    state_dependent_delay
end

# Group classifications
const delay_specific_options = (:discontinuity_tracking, :history_interpolation,
    :history_extrapolation, :delay_evaluation, :constrained_step, :residual_control,
    :neutral_delay, :state_dependent_delay)

function option_group(option::Symbol)
    if option in delay_specific_options
        return :delay_specific
    else
        error("Unknown verbosity option: $option")
    end
end

# Get all options in a group
function group_options(verbosity::DDEVerbosity, group::Symbol)
    if group === :delay_specific
        return NamedTuple{delay_specific_options}(getproperty(verbosity, opt)
                                                  for opt in delay_specific_options)
    else
        error("Unknown group: $group")
    end
end

function DDEVerbosity(;
        ode_verbosity = nothing, delay_specific = nothing, kwargs...)
    # Validate group arguments
    if delay_specific !== nothing && !(delay_specific isa AbstractMessageLevel)
        throw(ArgumentError("delay_specific must be a SciMLLogging.AbstractMessageLevel, got $(typeof(delay_specific))"))
    end

    if ode_verbosity !== nothing && !(ode_verbosity isa ODEVerbosity)
        throw(ArgumentError("ode_verbosity must be an ODEVerbosity, got $(typeof(ode_verbosity))"))
    end

    # Validate individual kwargs
    for (key, value) in kwargs
        if !(key in delay_specific_options)
            throw(ArgumentError("Unknown verbosity option: $key. Valid options are: $(delay_specific_options)"))
        end
        if !(value isa AbstractMessageLevel)
            throw(ArgumentError("$key must be a SciMLLogging.AbstractMessageLevel, got $(typeof(value))"))
        end
    end

    # Build arguments using NamedTuple for type stability
    default_args = (
        ode_verbosity = ode_verbosity === nothing ? ODEVerbosity() : ode_verbosity,
        discontinuity_tracking = Silent(),
        history_interpolation = Silent(),
        history_extrapolation = WarnLevel(),
        delay_evaluation = Silent(),
        constrained_step = Silent(),
        residual_control = Silent(),
        neutral_delay = Silent(),
        state_dependent_delay = Silent()
    )

    # Apply group-level settings
    final_args = if delay_specific !== nothing
        NamedTuple{keys(default_args)}(
            _resolve_arg_value(key, default_args[key], delay_specific)
            for key in keys(default_args)
        )
    else
        default_args
    end

    # Apply individual overrides
    if !isempty(kwargs)
        final_args = merge(final_args, NamedTuple(kwargs))
    end

    DDEVerbosity(values(final_args)...)
end

# Constructor for verbosity presets following the hierarchical levels:
# None < Minimal < Standard < Detailed < All
# Each level includes all messages from levels below it plus additional ones
function DDEVerbosity(verbose::AbstractVerbosityPreset)
    if verbose isa Minimal
        # Minimal: Only fatal errors and critical warnings
        DDEVerbosity(
            ode_verbosity = ODEVerbosity(Minimal()),
            discontinuity_tracking = Silent(),
            history_interpolation = Silent(),
            history_extrapolation = WarnLevel(),
            delay_evaluation = Silent(),
            constrained_step = Silent(),
            residual_control = Silent(),
            neutral_delay = Silent(),
            state_dependent_delay = WarnLevel()
        )
    elseif verbose isa Standard
        # Standard: Everything from Minimal + non-fatal warnings
        DDEVerbosity()
    elseif verbose isa Detailed
        # Detailed: Everything from Standard + debugging/solver behavior
        DDEVerbosity(
            ode_verbosity = ODEVerbosity(Detailed()),
            discontinuity_tracking = InfoLevel(),
            history_interpolation = InfoLevel(),
            history_extrapolation = WarnLevel(),
            delay_evaluation = InfoLevel(),
            constrained_step = InfoLevel(),
            residual_control = InfoLevel(),
            neutral_delay = InfoLevel(),
            state_dependent_delay = WarnLevel()
        )
    elseif verbose isa All
        # All: Maximum verbosity - every possible logging message at InfoLevel
        DDEVerbosity(
            ode_verbosity = ODEVerbosity(All()),
            discontinuity_tracking = InfoLevel(),
            history_interpolation = InfoLevel(),
            history_extrapolation = InfoLevel(),
            delay_evaluation = InfoLevel(),
            constrained_step = InfoLevel(),
            residual_control = InfoLevel(),
            neutral_delay = InfoLevel(),
            state_dependent_delay = InfoLevel()
        )
    end
end

@inline function DDEVerbosity(verbose::None)
    DDEVerbosity(
        ODEVerbosity(None()),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent()
    )
end

# Helper function to resolve argument values based on group membership
@inline function _resolve_arg_value(key::Symbol, default_val, delay_specific)
    if key === :ode_verbosity
        return default_val
    elseif key in delay_specific_options && delay_specific !== nothing
        return delay_specific
    else
        return default_val
    end
end

function Base.getproperty(v::DDEVerbosity, s::Symbol)
    if s in fieldnames(DDEVerbosity)
        return getfield(v, s)
    elseif s in fieldnames(ODEVerbosity)
        return getfield(getfield(v, :ode_verbosity), s)
    else
        return error("type DDEVerbosity has no field ", s)
    end
end