using Test, DelimitedFiles, Logging
include(joinpath(@__DIR__, "..", "RAD.jl"))

series = Dict{String,Vector{Float64}}()
for line in eachline(joinpath(@__DIR__, "series.csv"))
    parts = split(line, ",")
    series[parts[1]] = parse.(Float64, parts[2:end])
end
expected, _ = readdlm(joinpath(@__DIR__, "expected.csv"), ',', Any; header = true)

@testset "RAD" begin
    @testset "matches reference" begin
        for r in eachrow(expected)
            f = with_logger(NullLogger()) do
                RAD(series[r[1]], Int(r[2]), r[3] == 1)
            end
            e = Float64(r[4] isa AbstractString ? parse(Float64, r[4]) : r[4])
            if isnan(e)
                @test isnan(f)
            else
                @test f ≈ e rtol = 1e-12
            end
        end
    end
    @test RAD(series["ar1"]) == RAD(series["ar1"], 1, true)
    @test_logs (:warn,) isnan(RAD(series["constant"]))
    @test_throws ArgumentError RAD(series["ar1"], 0)
end
