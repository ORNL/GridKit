# Developer Guidelines

## Git workflows

### Keeping your local repository updated

Main development branch is `develop`. The `main` branch is used for stable,
ready-to-release code. Use `git pull` to keep `develop` branch in your local
repo updated. If you get conflicts when pulling from upstream reset your
local `develop` to the upstream:
```
git fetch
git checkout develop
git reset --hard origin/develop
```
Here `origin` is the name of the remote upstream repository in your local
repository. If you set custom name for your remote, replace `origin` with
that name.

### Feature branches

For each new feature feature create a new branch from `develop`. Name your
feature branch `<developer name>/<short_feature_description>_dev. Bug fix
branches should follow similar pattern, only ending with `_fix`. Always
use underscores to separate different words.
```
git checkout -b shaked/sparse_matrix_transpose_dev
```
In some instances you also might want to branch off a feature branch, for
example to fix a bug or suggest significant changes.

Before engaging in extensive development, check the existing issues and pull 
requests for known issues and fixes. Also consider creating an issue to get 
feedback on your plans.

### Merging your feature branch

Feature branches should be merged back to the branch from where they were
branched off (typically `develop`). All merging should be done through
GitHub pull request.

Once your branch is ready, create pull request using the template provided
in the GridKit™ repository. Make sure that:
- All tests pass.
- Code compiles cleanly with flags `-Wall -Wpedantic -Wconversion -Wextra`.
- The new code follows GridKit™ style guidelines.
- There are unit tests for the new code.
- The new code is documented.
- The feature branch is rebased with respect to the target branch.
- The [CHANGELOG.md](/CHANGELOG.md) has been updated to reflect the changes.

Ensure that the proposed changes have a clear scope, and create separate pull 
requests for unrelated changes. We recommend for each developer to have a 
finite (typically less than 5) number of pull requests open at once. 
Check with the maintainers if you are waiting for an older pull request to 
be reviewed and merged.

It is the responsibility of the developer to ensure their pull request is reviewed in 
timely fashion. For pull requests affecting large number of files, please 
coordinate with reviewers (offline) to make sure they will be available to review. 
To rebase your feature branch, first ensure the target branch is up-to-date.
Then use
```
git checkout <feature branch>
git rebase -i <target branch>
```
to rebase your branch to the target branch. Follow the instructions on the
screen. You may need to resolve rebase conflicts.

There has to be at least one approval before the pull request can be merged.

### Using draft pull requests
Use draft pull requests to test your proposed changes in continuous integration 
(CI) and/or to get early feedback on proposed changes. 
Some reviewers will not look at the code until a pull request is marked as 
ready for review, so make sure to communicate the intended purpose of your draft
pull request.

Draft pull requests are not meant for casual collaboration. 
You can simply point other developers to your feature branch if you need to 
discuss your proposed changes without running CI or requesting a formal review.

Draft pull requests should be in situations where developer wants to confirm their
approach is acceptable before making a large scale change. For example, if a 
developer wants to make a change that affects all component models, the developer
is advised to make such change to only one model (or small subset of models) and
make it available for review in a draft pull request. Once the development direction is
approved, the developer can proceed and make the wholesale change.

### Deleting your feature branch

Once your branch has been merged or is otherwise no longer needed, delete it 
from the GridKit™ GitHub repository. We recommend cleaning up soon after 
the merge and periodically going through your existing branches to cleanup.

## Documenting Code

### Doxygen
All comments in the code should follow [Doxygen](https://www.doxygen.nl/manual/index.html)
markup. For uniformity, we recommend C-style Doxygen comments starting with
two `*`.

### Documenting functions
Functions should be documented in source files. The rationale is to have
the documentation near the implementation, so that is handy to a developer and
can be updated quickly when the function modified.

### Minimal function documentation
Function documentation should include specification of function parameters,
template parameters (if any), return value, preconditions, postconditions,
and invariants.

```c++
/**
 * @brief <BRIEF DESCRIPTION>
 *
 * @param <NAME> <DESCRIPTION>
 * @tparam <NAME> <DESCRIPTION>
 * @return <DESCRIPTION OF RETURN VALUE>
 *
 * @pre <PRECONDITION>
 * @post <POSTCONDITION>
 *
 * LONGER DESCRIPTION, RUNTIME, EXAMPLES, ETC
 */
```

* `@brief` marks the text that will be displayed in summaries and index lists.
Typically you would put here a few words description of your function.
* `@param` describes a function parameter and takes an optional direction:
`@param[in]` means the parameter's value is only read and not modified within
the function, `@param[out]` means the parameter is not read and is only
modified, and `@param[in,out]` means the parameter is both read and modified.
* `@pre` and `@post` define the pre- and postconditions, which should be
precise but brief. When in doubt, attempt rigorous conditions but keep in mind
that some concepts such as "validity" may be difficult or impossible to define
precisely. Specifications are primarily for human consumption.
* `@param`, `@pre`, and `@post` sections should be repeated as many time as
required.

### Doxygen and Markdown

Doxygen supports Markdown markup and it should be used to make documentation
more clear. For example,
```text
 * @return The size of `a`
```
is clearer than
```text
 * @return The size of a
```
when read in plain text and in formatted documentation.

## Code Style

### Existing non-compliant code
The code that does not comply with these guidelines will be fixed in separate
pull request(s). Any new contribution should follow closely these guidelines
but should not change style in the existing GridKit code.

### Error handling
Return values of member functions should be of type `int` and used for error
handling. Functions return 0 if no error is encounter, return positive value
for warnings and recoverable error, and negative value for irrecoverable
errors.

### Output
If an output is needed (for example, a warning needs to be displayed), use 
the logger. Use `"\n"` rather than `"std::endl"`, unless flushing really is needed.
There should be a space before and after each `<<`.
If the line needs to be broken, the `<<` operators should be aligned:
```c++
GridKit::Utilities::Logger::warning() << "Index out of bounds."
                                      << " Row " << i
                                      << " starts at: " << start 
                                      << " and ends at " << end << "\n";
```

### File names

Each class should be implemented in a `*.cpp` and a `*.hpp` files with
the same name. Files containing collection of standalone functions should
have a descriptive name starting with lowercase character.

### Using `clang-format`

Coding style rules are specified in `.clang-format` file in GridKit root
directory. Use
```shell
clang-format -i <filename>
```
to make sure your code follows the style guidelines.

### Line length

Keep line length below 80 characters unless a longer line is needed for better
clarity of the code. If your line is longer than 90 characters, you are
probably doing something wrong.

`clang-format` is configured so it does not enforce line length; that is
developer's responsibility. Break the long lines at places where it will
improve readability of the code. Running `clang-format` after that will help
align the code properly in accordance with GridKit coding style.

### Local variables naming

Local variable names should use C-style name format.
```c++
double local_variable; // Yes
double localVariable;  // No, using lowercase camel instead of C-style name format
```

### Member variables naming

Member variable names should use C-style name format and end with trailing underscore `_`.
```c++
double member_variable_; // Yes
double another_member;   // No, there is no trailing underscore to distinguish it from nonmember variables
double memberVariable_;  // No, using lowercase camel instead of C-style name format
```

#### Exceptions

Public member variables that are accessed directly do not need trailing
underscores. For example, consider this code:
```c++
struct ModelData
{
  int id;
  double value;
};

ModelData data;
data.id = 1;
data.value = 2.0;
```
Member variables of struct `data` are accessed directly outside the struct
and do not need to be denoted with trailing underscores `_`.

### Function names

Use lowercase camel format for function names.
```c++
int myFunction(double x); // Yes
int another_function();   // No, using C-style name format
int YetAnotherFunction(); // No, using uppercase camel name format
```

### Class names

Class names should use uppercase camel name format.
```c++
class MyClass // Yes
{
  ...
}

class Matrix // Yes
{
  ...
}

class My_Class // No, using underscore in class name
{
  ...
}
```

### Type declarations and template parameters

Always declare type aliases with the `using` keyword rather than `typedef typename`. 

```c++
  using RealT = typename GridKit::ScalarTraits<ScalarT>::RealT; // Yes
  typedef typename GridKit::ScalarTraits<ScalarT>::RealT RealT; // No
```

Types that are used as template parameters should use `snake_case` format, with
the `_type` suffix to indicate that it is a type. Use the keyword `typename`
rather than `class` for template type parameters.

```c++
  template <typename RealT>; // No, using format reserved for interface types
  template <class real_type>; // No, using `class` keyword

  template <typename real_type> ...; // Yes
```
Types used within class implementation code and provided in the public class
interface should use `UpperCamelT` format, similar to class names, but with a
`T` suffix.  For consistency, create aliases for template type parameters
(*e.g.,* `using ScalarT = scalar_type;`) and use the same name everywhere the
same type is used. This applies across multiple classes in which the type has
the same meaning.

```c++
  template <typename scalar_type, typename index_type>
  class ExTemp
  {
    // Convention used in all applicable classes for scalar, index, and
    // primitive floating-point types
    using ScalarT = scalar_type;
    using IdxT    = index_type;
    using RealT   = GridKit::ScalarTraits<ScalarT>::RealT;

    // void foo(scalar_type value); // No, using template type parameter
    void foo(ScalarT value); // Yes, using interface types in code
  };
```

### Constants

If a constant is used in more than one file, define it in `Constants.hpp`. For
constants with long names, use underscores to separate words in the constant
name. Use all caps (screaming snake case).

```c++
   constexpr double Pi = 3.1415;      // No, use all caps for the constant name
   constexpr double SQRT_TWO = 1.4142 // Yes
   constexpr double SQRTTWO = 1.4142  // No, the two words not separated by "_"
   constexpr double EXP = 2.7183      // Yes
```

### Macros

Macros should use all caps, same as constants.

### Enums (enumerated types)

Always define `enum`s inside `GridKit` namespace. The `enum` name should
be upper camel case, same as class names. The `enum` element names should
match symbol for physics quantity they represent or they should be uppercase
(same as names for constants) if they do not represent a physics quantity.
For example, `enum` element for real component of voltage $V_r$ should be `Vr`.
A name for `enum` element for a "fast mode", for example, should be something
like `FAST_MODE`, capitalized with underscores separating words (but no
underscore at the end!).

```c++
  enum ExampleEnum { Vr, // Yes, it matches symbol for real voltage component
                     VR, // No, the element name should match the physics symbol
                     FAST_MODE}; // Yes, element name is all caps.
```

### Pointers and references

The pointer `*` or reference `&` belong to the type and there should be no 
space between them and the type name.
```c++
double* x;     // Yes
int& n;        // Yes
double *x, *y; // No, the pointer symbol is a part of `double*` type
int & n;       // No, the reference symbol is a part of `int&` type
```

### Indentation
Use only spaces for indentation, not tabs. Indent size is 2 spaces.

When defining a class, the code blocks after `private`, `public` and `protected`
should be aligned with opening/closing braces. There should be an empty line
before each definition (except the first one). See example below.
```c++
class SomeClass
{
public:
  SomeClass();
  ~SomeClass();

private:
  int some_variable_;

protected:
  void someFunction();
};
```

### Braces
All braces should follow Allman style:
```c++
namespace SomeNamespace
{
  //some code
}
```
For short functions (i.e., empty constructor), do not inline braces.
```c++
ClassA::ClassA()
{
}
```
Have opening brace at the next line following  `for`, `if`, or `while`
statement. When using `else`, follow the example below.
```c++
if (cond == true)
{
  // some code
}
else
{
  // some other code
}
 ```
Have a space between keywords `for`, `while` and `if` and the parenthesis as
shown here:
```c++
for (int i = 0; i < n; ++i)
{
  // some code
}
```

Do not use one-line `if`s and `for`s. Always use braces.

### Use of spaces and newlines
There should be spaces between arithmetic operators.
```c++
x = c * (a + b);  // Yes
x = c*(a+b).      // No, the clarity is better if there are spaces between
                  // binary operators and operands.
```
When defining member functions, use one empty line between the functions.
```c++
struct MyStruct
{
  int memberFunction()
  {
    // some code
  }

  int anotherMemberFunction()
  {
    // some other code
  }
};
```

### Include files

Leave one empty line between all the includes and the first line of the actual code.
```c++
#include <iostream>

int main()
{
  std::cout
}
```

Header files should be included in 3 separate blocks: standard libraries,
GridKit external dependencies, and GridKit header files. There should be an
empty line between the blocks. External libraries should always use `<...>`.
GridKit headers should be included with `<GridKit/...>` (using full path to
file) with one exception: headers that are local to a compilation unit should
use `"..."`. That is, only in a `.cpp` file and only those headers that are
local to the component (headers from other project components should still use
`<...>`).

```c++
#include "Ida.hpp"     // GridKit local internal header

#include <iostream>    // Standard libs headers
#include <cmath>

#include <sundials.h>  // GridKit dependencies
#include <idas.h>

#include <GridKit/Model/Evaluator.hpp> // GridKit header from another component

```

```c++
#include <cstring>

#include <SparseMatrix/CooMatrix.hpp>

int main()
{
  //some code
  return 0;
}
```
The `system` includes should always be listed first.

### Long function names

Function declarations and function calls should be written on a single line.
If a function has large number of parameters, the function should be broken
over multiple lines so that each function parameter is on a separate line.
All parameters should be aligned with the first parameter in the function.

```c++
int myFunction(double x, double y); // Yes, function on a single line

int myFunction(double x1,  // Yes, function parameters on a separate lines.
               double x2,
               double x3,
               double x4,
               double x5,
               double x6,
               double x7);

int myFunction(double x1,  // No, function broken inconsistently.
               double x2, double x3,
               double x4,
               double x5, double x6, double x7);

int myFunction(
  double x1,  // No, first parameter should follow the parenthesis.
  double x2,
  double x3,
  double x4,
  double x5,
  double x6,
  double x7
);

```

### Initializer list in class constructor

Short initializer list should follow the constructor on the same line with
space around `:`. Long initializer lists should begin on the next line,
indented and starting with `:`. Initializers should be aligned.
```c++
// Short initializer list
MyClass(n, m) : n_(n), m_(m)
{
  // ...
}

// Long initializer list
MyClass(n, m)
  : n_(n),
    m_(m),
    pX_(nullptr),
    pY_(nullptr),
    pZ_(nullptr),
    size_(0)
{
  // ...
}
```

### Using namespaces
All classes should be in namespace `GridKit`. If needed, define additional
namespaces inside `GridKit`.
```c++
namespace GridKit
{
  class Solver  // Yes, class defined inside GridKit namespace
  {
    // some code;
  };

  namespace LinearAlgebra
  {
    class Vector  // Yes, class defined inside GridKit namespace
    {
      // vector code
    };
  }
}

class Matrix   // No, class is outside GridKit namespace
{
  // matrix code
};
```

## Development Container
A development container is available for all developers using VS Code to develop. 
This will automatically install all pre-requisite software you need to develop 
in GridKit. Any developer who wishes to use this setup can follow 
[this tutorial](https://code.visualstudio.com/docs/devcontainers/tutorial) and 
simply use the option "Reopen Folder in Container" rather than 
"New Dev Container...", which will automatically build the included container.

## Electric Grid Test Cases

When adding a new test case to to the repository using the GridKit input file 
format, you should follow the following guidelines to remain consistent with 
existing GridKit cases. For each test case, the associated README.md should 
contain the following:
- A high resolution oneline diagram ($\geq$ 600 dpi)
- No overlapping labels
- Use common electrical symbols for components (transformers, generators, loads, etc.)
- Use a calm color pallet

Within the README file for the test case, specify characteristics such as:
- Which component models are used
- Quantity of each component model used
- Types of events well-suited for the case
- Any other relevant case characteristics 

Additionally, specify existing multiple resolutions of the case, i.e., if 
there is a high-fidelity EMT network model along with a phasor-domain network 
model that is well known and available, we should indicate that both these 
model resolutions are available within GridKit. 
