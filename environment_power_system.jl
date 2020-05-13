
Pkg.add("DataFrames");
Pkg.add("XLSX");
Pkg.add("GLPK");
Pkg.add("Gurobi");
Pkg.add("MathOptInterface");
Pkg.add("MathOptFormat");
Pkg.add("Statistics")

using DataFrames
using ExcelReaders
#using GLPK
using MathOptInterface, MathOptFormat
const MOI = MathOptInterface
using XLSX
using Gurobi
using JuMP
using Statistics
