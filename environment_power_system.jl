# the new one 
Pkg.add("DataFrames");
Pkg.add("XLSX");
Pkg.add("GLPK");
Pkg.add("Gurobi");
Pkg.add("MathOptInterface");
Pkg.add("MathOptFormat");
Pkg.add("Statistics")
Pkg.add("ExcelReaders")

using DataFrames
using ExcelReaders
#using GLPK
using Gurobi
using MathOptInterface
using MathOptFormat
const MOI = MathOptInterface
using XLSX

using JuMP
using Statistics
# the old one 


# Pkg.add("DataFrames");
# Pkg.add("XLSX");
# Pkg.add("GLPK");
# Pkg.add("Gurobi");
# Pkg.add("MathOptInterface");
# Pkg.add("MathOptFormat");
# Pkg.add("Statistics")

# using DataFrames
# using ExcelReaders
# #using GLPK
# using MathOptInterface, MathOptFormat
# const MOI = MathOptInterface
# using XLSX
# using Gurobi
# using JuMP
using Statistics
