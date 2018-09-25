using Base.tail

# Outer constructor for DEBStructure
DEBStructure(name::Symbol, param_specs::ParamSpecs, params, flags::DEBFlags,
             functions, settings) = begin
    (tmin, tmax) = ceil.(Int, settings[:tspan])
    timerange = tmin:tmax
    arraylength = tmax - tmin + 1

    rates = fill(0.0, arraylength)
    param_ids = [p.id for (k, p) in param_specs]

    A = 0.0
    init_params = deepcopy(params)
    u = init_state(typeof(settings[:u0][1]), settings[:state_type])
    num_structures = floor(Int64, length(settings[:u0])/length(u))
    assim_state = uzeros(PhotoIntermediate)

    state_names = fieldnames(u)
    state1_names = get_state1_names(u)
    Jbase = build_axis(state_names, TRANS, arraylength, timerange)
    J1base = build_axis(state1_names, TRANS1, arraylength, timerange)
    J = split_flux(Jbase, 1)
    J1 = split_flux(J1base, 1)
    structure = DEBStructure(name, param_specs, param_ids, init_params, params, flags, 
                             functions, u, A, Jbase, J1base, J, J1, rates, assim_state)
    return structure
end


function split_flux(base::AbstractArray, t::Int)
    view(base, :, :, t)
end

function unitize(a, t)
    if t.name.name == :Quantity
        t(a)
    else
        a
    end
end

function uconstruct(datatype::DataType, args)
    fts = fieldtype.(datatype, fieldnames(datatype))
    fields = map(unitize, args, fts)
    datatype(fields...)
end

function uzeros(datatype::DataType)
    args = zeros(Float64, length(fieldnames(datatype)))
    uconstruct(datatype, args)
end

function uones(datatype::DataType)
    args = ones(Float64, length(fieldnames(datatype)))
    uconstruct(datatype, args)
end


# -----------------------------------------------------------------------------
# Structure update functions. 
# Recursion is used instead of loops for type stability.

split_state!(ss::Tuple{<:DEBStructure,Vararg}, offset::Int, u::A where A)::Nothing = begin
    states = length(fieldnames(ss[1].u))
    for i = 1:states
        ss[1].u[i] = u[i + offset]
    end
    split_state!(tail(ss), offset + states, u)
end
split_state!(::Tuple{}, ::Int, _)::Nothing = nothing


initialise_params!(s) = begin
    for name in fieldnames(s.params)
        setfield!(s.params, name, getfield(s.init_params, name))
    end
end


set_current_flux!(ss::Tuple{<:DEBStructure,Vararg}, t::Int)::Nothing = begin
  s = ss[1]
  s.J = split_flux(s.Jbase, t + 1)
  s.J1 = split_flux(s.J1base, t + 1)

  set_current_flux!(tail(ss), t)
end
set_current_flux!(::Tuple{}, ::Int)::Nothing = nothing


sum_flux!(du, p::P where P<:AbstractStructuredSettings)::Nothing = begin
    ss = p.structures
    num_structures = length(p.structures)
    trans = length(TRANS)
    states = length(fieldnames(ss[1].u))
    sum_flux!(du, ss, num_structures, 0, states, trans)
    return nothing
end

sum_flux!(du, ss::Tuple{<:DEBStructure,Vararg},
          num_structures::Int, offset::Int, states::Int, trans::Int)::Nothing = begin
    J = ss[1].J
    for x = 1:states
        sum = 0.0unit(J[1,1])
        for y = 1:trans
            sum += J[x, y]
        end
        du[x + offset] = sum * 1.0u"hr"
    end
    sum_flux!(du, tail(ss), num_structures, offset + states, states, trans)
end
sum_flux!(du, ::Tuple{}, ::Int, ::Int, ::Int, ::Int)::Nothing = nothing
