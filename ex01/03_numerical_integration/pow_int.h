#ifndef HASC_POW_INT_HH
#define HASC_POW_INT_HH
// Taken from https://github.com/ampl/gsl/blob/master/sys/pow_int.c

namespace hasc_ex01
{
double pow_uint(double x, unsigned int n)
{
  double value = 1.0;

  /* repeated squaring method 
   * returns 0.0^0 = 1.0, so continuous in x
   */
  do {
     if(n & 1) value *= x;  /* for n odd */
     n >>= 1;
     x *= x;
  } while (n);

  return value;
}

double pow_int(double x, int n)
{
  unsigned int un;

  if(n < 0) {
    x = 1.0/x;
    un = -n;
  } else {
    un = n;
  }

  return pow_uint(x, un);
}

} // namespace hasc_ex01

#endif // HASC_POW_INT_HH