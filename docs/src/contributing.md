# Contributing to OpenSEFA.jl

We welcome contributions, feedback, and suggestions for improvement!

## Ways to Contribute

- Report bugs or request features via [GitHub Issues](https://github.com/MatplusGmbH/opensefa-reconciliation-solver/issues)
- Submit pull requests with bug fixes or new features
- Improve documentation
- Share your use cases and examples
- Help answer questions from other users

## Development Workflow

### Branch Strategy

- **`main`** is the default branch
- Create feature branches from `main` for your work
- All pull requests should target the `main` branch

### Setting Up Your Development Environment

1. Fork or clone the repository:
   ```bash
   git clone https://github.com/MatplusGmbH/opensefa-reconciliation-solver.git
   cd opensefa-reconciliation-solver
   ```

2. Install dependencies:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

3. Run tests to verify your setup:
   ```bash
   julia --project=. test/runtests.jl
   ```

### Making Changes

1. Create a feature branch:
   ```bash
   git checkout main
   git pull origin main
   git checkout -b feature/your-feature-name
   ```

2. Make your changes and add tests if applicable

3. Run the test suite to ensure nothing breaks:
   ```bash
   julia --project=. test/runtests.jl
   ```

4. Commit your changes with clear, descriptive messages:
   ```bash
   git add .
   git commit -m "Add feature: description of your changes"
   ```

5. Push to your branch:
   ```bash
   git push origin feature/your-feature-name
   ```

6. Create a pull request targeting the `main` branch

## Contributing to Documentation

Documentation is a critical part of this project. Here's how to contribute:

### Building Documentation Locally

#### Prerequisites

- Julia 1.10 or later
- Git

#### Setup and Build

1. Navigate to the documentation directory:
   ```bash
   cd docs
   ```

2. Install documentation dependencies:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

3. Build the documentation:
   ```bash
   julia --project=. make.jl
   ```

4. Preview in your browser:
   ```bash
   # On macOS
   open build/index.html

   # On Linux
   xdg-open build/index.html

   # On Windows
   start build/index.html
   ```

### Troubleshooting Documentation Builds

#### Issue: Package installation fails

**Symptoms:**
- Error messages during `Pkg.instantiate()`
- Missing package errors

**Solutions:**
- Ensure you're in the `docs` directory when running commands
- Try removing `Manifest.toml` and running `Pkg.instantiate()` again:
  ```bash
  rm Manifest.toml
  julia --project=. -e 'using Pkg; Pkg.instantiate()'
  ```
- Verify the parent `OpenSEFA` package builds successfully first
- Check that your Julia version is up to date

#### Issue: Documentation build fails

**Symptoms:**
- Errors during `julia --project=. make.jl`
- Missing cross-references
- Failed doctests

**Solutions:**
- Check that you've successfully built the main package first
- Verify all code examples in markdown files are valid Julia code
- Look for broken cross-references in the build output (links with `@ref`)
- Check for syntax errors in markdown files
- Ensure all referenced functions/types exist and are exported
- Run with verbose output to see detailed errors:
  ```bash
  julia --project=. --color=yes make.jl
  ```

#### Issue: Changes not reflecting

**Symptoms:**
- Updated markdown files don't show changes in built HTML
- Old content still appears

**Solutions:**
- Delete the `build/` directory and rebuild:
  ```bash
  rm -rf build
  julia --project=. make.jl
  ```
- Clear your browser cache when viewing updated docs
- Do a hard refresh in your browser (Cmd+Shift+R on macOS, Ctrl+Shift+R on Linux/Windows)
- Check that you saved your markdown files before rebuilding

#### Issue: Mathematical formulas not rendering

**Symptoms:**
- LaTeX equations show as raw code
- Math symbols display incorrectly

**Solutions:**
- Verify you're using the correct syntax:
  ```markdown
  # Inline math
  This is inline math: $x^2 + y^2 = z^2$

  # Block math
  ```math
  \begin{aligned}
  x &= y + z \\
  a &= b + c
  \end{aligned}
  \```
  ```
- Check for unclosed braces or brackets in LaTeX code
- Escape special characters if needed

### Documentation Style Guide

When contributing documentation:

1. **Write clear, concise content**
   - Use active voice
   - Keep sentences short
   - Avoid jargon when possible

2. **Use proper formatting**
   - Use code blocks with language hints: ` ```julia`
   - Use inline code for function names: `` `function_name()` ``
   - Use bold for emphasis: `**important**`
   - Use italic sparingly: `*optional*`

3. **Include examples**
   - Every new feature should have at least one example
   - Examples should be self-contained when possible
   - Test all code examples before committing

4. **Cross-reference properly**
   - Link to other sections: `[Getting Started](@ref getting_started)`
   - Link to functions: `` [`solve`](@ref) ``
   - Use section IDs consistently: `# [Section Title](@id section_id)`

5. **Verify cross-references work**
   - All `@ref` links should resolve without warnings
   - Check the build output for "missing docs" warnings

6. **Check mathematical formulas**
   - Test LaTeX rendering in the built documentation
   - Use proper alignment in multi-line equations
   - Define variables clearly

### Continuous Documentation Development

For iterative documentation work:

```bash
# Make changes to markdown files in docs/src/
# Then rebuild:
julia --project=. make.jl

# Preview changes:
open build/index.html  # or xdg-open/start depending on OS
```

Consider using a file watcher for automatic rebuilds (external tools like `entr` or `watchexec`):

```bash
# Example with entr (install separately)
ls src/*.md | entr julia --project=. make.jl
```

## Code Style Guidelines

- Follow [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- Use 4 spaces for indentation (no tabs)
- Keep lines under 92 characters when reasonable
- Add docstrings to all public functions
- Use meaningful variable names

### Docstring Format

Use the standard Julia docstring format:

```julia
"""
    function_name(arg1, arg2; kwarg1=default)

Brief description of what the function does.

# Arguments
- `arg1::Type`: Description of arg1
- `arg2::Type`: Description of arg2
- `kwarg1::Type`: Description of kwarg1 (default: `default`)

# Returns
- `ReturnType`: Description of return value

# Examples
```julia
julia> function_name(1, 2)
3
\```

# See also
- `related_function`: Brief description of a related function
"""
function function_name(arg1, arg2; kwarg1=default)
    # implementation
end
```

## Testing Guidelines

- Add tests for all new features
- Ensure tests are deterministic and reproducible
- Test edge cases and error conditions
- Use descriptive test names with `@testset`

Example:
```julia
@testset "Feature Name" begin
    @testset "Normal case" begin
        result = my_function(input)
        @test result == expected_value
    end

    @testset "Edge case" begin
        @test_throws ErrorType my_function(bad_input)
    end
end
```

## Updating Release Notes

When making significant changes, update `NEWS.md`:

1. Add a new version section at the top if starting a new release
2. Use appropriate sections:
   - `### Breaking Changes` - Incompatible API changes
   - `### New Features` - New functionality
   - `### Improvements` - Enhancements to existing features
   - `### Bug Fixes` - Bug fixes
   - `### Testing` - Test improvements
3. Format entries as bullet points with backticks for code:
   ```markdown
   - Added `new_function` for improved performance [#XXX](https://github.com/MatplusGmbH/opensefa-reconciliation-solver/issues/XXX)
   ```
4. Link to relevant issues or pull requests

## Submitting Pull Requests

1. Ensure your code follows the style guidelines
2. Add/update tests as needed
3. Update documentation if you've changed APIs
4. Update `NEWS.md` if appropriate
5. Ensure all tests pass locally
6. Write a clear pull request description:
   - What problem does this solve?
   - What changes were made?
   - Any breaking changes?
   - Related issues?
7. Confirm that you have read and accept the [Contributor License Agreement](https://github.com/MatplusGmbH/opensefa-reconciliation-solver/blob/main/CLA.md)

## Getting Help

- Check existing [Issues](https://github.com/MatplusGmbH/opensefa-reconciliation-solver/issues) for similar problems
- Review the [documentation](https://matplusgmbh.github.io/opensefa-reconciliation-solver/stable/)
- Contact the maintainers (see Project.toml for contact details)

## License

OpenSEFA.jl is licensed under the [PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/). The full license text is in the [`LICENSE`](https://github.com/MatplusGmbH/opensefa-reconciliation-solver/blob/main/LICENSE) file at the package root.

By contributing to OpenSEFA.jl, you accept the [Contributor License Agreement](https://github.com/MatplusGmbH/opensefa-reconciliation-solver/blob/main/CLA.md), which provides for two grants:

- **To the community:** your contributions are licensed under the PolyForm Noncommercial License 1.0.0 when distributed as part of the project.
- **To Matplus:** you grant Matplus GmbH a perpetual, irrevocable license to use your contributions for any purpose, including commercial products and services.

Submitting a pull request constitutes acceptance of the Contributor License Agreement.

Contributors acting on behalf of a company must have authority to grant the agreement; contact [Matplus GmbH](https://www.matplus.eu/) for corporate contributor arrangements if needed.

Thank you for contributing! 🎉
