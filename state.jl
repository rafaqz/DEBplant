abstract type AbstractState{N,T} <: FieldVector{N,T} end
abstract type AbstractStateE{N,T} <: AbstractState{N,T} end
abstract type AbstractStateCN{N,T} <: AbstractState{N,T} end
abstract type AbstractStateCNE{N,T} <: AbstractState{N,T} end

@mix @label @units @default_kw struct P{T} P::T | 0.0  | u"mol" | "Production"       end
@mix @label @units @default_kw struct V{T} V::T | 1e-4 | u"mol" | "Structure"        end
@mix @label @units @default_kw struct M{T} M::T | 0.0  | u"mol" | "Maturity"         end
@mix @label @units @default_kw struct C{T} C::T | 1e-4 | u"mol" | "Carbon Reserve"   end
@mix @label @units @default_kw struct N{T} N::T | 1e-4 | u"mol" | "Nitrogen Reserve" end
@mix @label @units @default_kw struct E{T} E::T | 10.0 | u"mol" | "General Reserve"  end

@E mutable struct StateE{} <: AbstractStateE{1,T} end
@C @N mutable struct StateCN{} <: AbstractStateCN{2,T} end
@C @N @E mutable struct StateCNE{} <: AbstractStateCNE{3,T} end
@V mutable struct StateV{} <: AbstractStateE{1,T}  end
@V @E mutable struct StateVE{} <: AbstractStateE{2,T} end
@V @C @N mutable struct StateVCN{} <: AbstractStateCN{3,T} end
@V @C @N @E mutable struct StateVCNE{} <: AbstractStateCN{3,T} end
@V @M @E mutable struct StateVME{} <: AbstractStateE{3,T} end
@V @M @C @N mutable struct StateVMCN{} <: AbstractStateE{3,T} end
@V @M @C @N @E mutable struct StateVMCNE{} <: AbstractStateE{3,T} end
@P @V mutable struct StatePV{} <: AbstractState{2,T} end
@P @V @E mutable struct StatePVE{} <: AbstractStateE{3,T} end
@P @V @C @N mutable struct StatePVCN{} <: AbstractStateCN{4,T} end
@P @V @C @N @E mutable struct StatePVCNE{} <: AbstractStateCNE{5,T} end
@P @V @M @E mutable struct StatePVME{} <: AbstractStateE{4,T} end
@P @V @M @C @N mutable struct StatePVMCN{} <: AbstractStateCN{5,T} end
@P @V @M @C @N @E mutable struct StatePVMCNE{} <: AbstractStateCNE{6,T} end

type HasM end
type NoM end
type HasP end
type NoP end

has_M(::Type{T}) where T <: AbstractState = :M in fieldnames(T) ? HasM() : NoM()
has_P(::Type{T}) where T <: AbstractState = :P in fieldnames(T) ? HasP() : NoP()

get_state1_names(state::AbstractStateE)   = [:E]
get_state1_names(state::AbstractStateCN)  = [:C, :N, :E]
get_state1_names(state::AbstractStateCNE) = [:EE, :CN, :C, :N, :E]

init_state(field, statetype) = statetype([0.0unit(field) for n in fieldnames(statetype)]...)
