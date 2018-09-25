struct IsIntermediate end
struct NotIntermediate end

intermediate(o) = typeof(o.params.transprop) == Nothing ? NotIntermediate() : IsIntermediate()
