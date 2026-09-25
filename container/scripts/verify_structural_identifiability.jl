#!/usr/bin/env julia

using Pkg
using StructuralIdentifiability

VERSION == v"1.10.12" || error("Unexpected Julia version: $(VERSION)")

package = only(
    entry for entry in values(Pkg.dependencies())
    if entry.name == "StructuralIdentifiability"
)
package.version == v"0.5.33" ||
    error("Unexpected StructuralIdentifiability version: $(package.version)")

# Published package example with known individual and combination results.
model = @ODEmodel(
    x1'(t) = -(a01 + a21) * x1(t) + a12 * x2(t) + u(t),
    x2'(t) = a21 * x1(t) - a12 * x2(t) - x3(t) / b,
    x3'(t) = x3(t),
    y(t) = x2(t),
)

status = assess_identifiability(model)
status[a21] == :globally || error("a21 should be globally identifiable")
status[a01] == :locally || error("a01 should be locally identifiable")
status[a12] == :locally || error("a12 should be locally identifiable")
status[b] == :nonidentifiable || error("b should be nonidentifiable")

combination_status = assess_identifiability(
    model;
    funcs_to_check = [a01 + a12, a01 * a12],
)
all(value == :globally for value in values(combination_status)) ||
    error("Expected parameter sum and product to be globally identifiable")

functions_found = replace.(string.(find_identifiable_functions(model)), " " => "")
for expected in ("a21", "a01+a12", "a01*a12")
    expected in functions_found || error("Missing identifiable function: $expected")
end

println("Julia structural-identifiability verification: PASS")
println("julia_version=$(VERSION)")
println("StructuralIdentifiability_version=$(package.version)")
