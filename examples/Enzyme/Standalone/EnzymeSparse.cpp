#include <assert.h>
#include <cmath>
#include <stdio.h>
#include <vector>

template <typename T>
struct Triple
{
  size_t row;
  size_t col;
  T      val;
  Triple(Triple&&) = default;

  Triple(size_t row, size_t col, T val)
    : row(row),
      col(col),
      val(val)
  {
  }
};

__attribute__((enzyme_sparse_accumulate)) static void inner_storeflt(int64_t row, int64_t col, float val, std::vector<Triple<float>>& triplets)
{
  triplets.emplace_back(row, col, val);
}

__attribute__((enzyme_sparse_accumulate)) static void inner_storedbl(int64_t row, int64_t col, double val, std::vector<Triple<double>>& triplets)
{
  triplets.emplace_back(row, col, val);
}

template <typename T>
__attribute__((always_inline)) static void sparse_store(T val, int64_t idx, size_t i, std::vector<Triple<T>>& triplets)
{
  if (val == 0.0)
    return;
  idx /= sizeof(T);
  if constexpr (sizeof(T) == 4)
    inner_storeflt(i, idx, val, triplets);
  else
    inner_storedbl(i, idx, val, triplets);
}

template <typename T>
__attribute__((always_inline)) static T sparse_load(int64_t idx, size_t i, std::vector<Triple<T>>& triplets)
{
  return 0.0;
}

template <typename T>
__attribute__((always_inline)) static void ident_store(T, int64_t idx, size_t i)
{
  assert(0 && "should never load");
}

template <typename T>
__attribute__((always_inline)) static T ident_load(int64_t idx, size_t i)
{
  idx /= sizeof(T);
  return (T) (idx == i);
}

extern int enzyme_dup;
extern int enzyme_const;
extern int enzyme_dupnoneed;

template <typename T, typename... Tys>
extern T __enzyme_fwddiff(void*, Tys...) noexcept;

template <typename T, typename... Tys>
extern T __enzyme_todense(Tys...) noexcept;

template <typename T>
__attribute__((always_inline)) static void f(size_t N, T* input, T* dinput)
{
  for (size_t i = 0; i < N; i++)
  {
    dinput[i] = input[i] * input[i];
  }
}

template <typename T>
__attribute__((noinline))
std::vector<Triple<T>>
jac_f(size_t N, T* input)
{
  std::vector<Triple<T>> triplets;
  for (size_t i = 0; i < N; i++)
  {
    T* output   = __enzyme_todense<T*>((void*) ident_load<T>, (void*) ident_store<T>, i);
    T* d_output = __enzyme_todense<T*>((void*) sparse_load<T>, (void*) sparse_store<T>, i, &triplets);

    __enzyme_fwddiff<void>((void*) f<T>,
                           enzyme_const,
                           N,
                           enzyme_dup,
                           input,
                           output,
                           enzyme_dupnoneed,
                           (T*) 0x1,
                           d_output);
  }
  return triplets;
}

int main()
{
  size_t N = 5;

  double* x   = (double*) malloc(sizeof(double) * N);
  double  val = 0.0;
  for (int i = 0; i < N; ++i)
  {
    x[i]  = val;
    val  += 1.0;
  }

  auto res = jac_f<double>(N, x);

  printf("Number of elements %ld\n", res.size());

  for (auto& tup : res)
  {
    printf("%ld, %ld = %f\n", tup.row, tup.col, tup.val);
  }

  delete[] x;

  return 0;
}
