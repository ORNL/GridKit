# Developer Guidelines

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
`std::cout` and not `printf` as shown below. There should be a space before
and after each `<<`. If the line needs to be broken, the `<<` operators should
be aligned:

```c++
std::cout << "index out of bounds. Row " << i << " starts at: " << start 
          << " and ends at " << end << std::endl;
```

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

### Function names

Use lowercase camel format for function names.
```c++
int myFunction(double x); // Yes
int another_function();   // No, using C-style name format
int YetAnotherFunction(); // No, using uppercase camel name format
```

### Class names

Class names should start with a capital letter. For instance, `Vector` and
`Matrix` are valid class names, while `point` is not.

### Enums (enumerated types)

Always define `enum`s inside `GridKit` namespace. Type names should be
capitalized and the constant names should be uppercase with underscores
(but there is no underscore at the end!).

```c++
  enum ExampleEnum { CONST_ONE = 0,
                     CONST_TWO = 8, 
                     YET_ANOTHER_CONST = 17 };
``` 

### Constants

If a constant is used in more than one file, define it in `Common.h`. For
constants with long names, use underscores to separate words in the constant
name. Use all caps (screaming snake case).

```c++
   constexpr double Pi = 3.1415; // No, the constant name should be all caps
   constexpr double SQRT_TWO = 1.4142 // Yes
   constexpr double SQRTTWO_ = 1.4142 // No, there is a trailing underscore
   constexpr double EXP = 2.7183 // Yes   
```

### Pointers and references

The pointer `*` or reference `&` belong to the type and there should be no space between them and the type name.
```c++
double* x;     // Yes
int& n;        // Yes
double *x, *y; // No, the pointer symbol is a part of `double*` type
int & n;       // No, the reference symbol is a part of `int&` type
```

### Indentation
Use only spaces for indentation, not tabs. Indent size is 4 spaces.

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
namespace someNamespace
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
x = c*(a+b).      // No, the clarity is better if there are spaces between binary operators and operands.
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
empty line between the blocks. External libraries should use `<...>`, while
GridKit headers should be included with `"..."`.

```c++
#include <iostream>    // Standard libs headers
#include <cmath>

#include <sundials.h>  // GridKit dependencies
#include <idas.h>

#include "Ida.hpp"     // GridKit internal header
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

