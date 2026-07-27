"""
Export reconciliation solutions to Excel files.

Uses XLSX.jl to write variable values, uncertainties, and solution metadata into a
spreadsheet, with one sheet per reconciliation problem.
"""
module ExcelModule

using XLSX

using ..ReconciliationModule

export export_excel

"Export the solution to an excel file with the given filename (full path)."
function export_excel(sol::ReconciliationSolution, filename::String)
  XLSX.openxlsx(filename, mode = "w") do xf
    sheet = XLSX.addsheet!(xf, sol.problem.name)

    sheet["A1"] = "Variable Name"
    sheet["B1"] = "Value"

    row = 2
    for (var, value) in pairs(sol.variables_value)
      sheet["A$row"] = string(var)
      sheet["B$row"] = isnothing(value) || ismissing(value) ? "?" : value
      row += 1
    end
  end
end

end
