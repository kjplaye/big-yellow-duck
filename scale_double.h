#include <math.h>

#define D2ULL(x) (*((unsigned long long *)&(x)))
#define LOG_BASE 0x1ff

//double with an augmented exponent.  number = value * 2^(LOG_BASE * exponent)
//We make sure that the exponent of the double is between 0x400 - LOG_BASE .. 0x400 + LOG_BASE
class scale_double
{
 public:
  double value;
  int exponent;
  
  scale_double() {value = 0;exponent = 0;}
  scale_double(double x) {value = x;exponent = 0;}
  double get_double(void) 
    {
      unsigned long long ex = (D2ULL(value) >> 52) & 0x7ff;
      double t;

      if (exponent == -1) 
	{
	  t = value;
	  D2ULL(t) ^= (ex << 52);
	  ex -= LOG_BASE;
	  D2ULL(t) ^= (ex << 52);
	  return t;
	}
      else if (exponent == 1) 
	{
	  t = value;
	  D2ULL(t) ^= (ex << 52);
	  ex += LOG_BASE;
	  D2ULL(t) ^= (ex << 52);
	  return t;
	}
      else if (exponent < 0) return 0;
      else if (exponent > 0) return -log(0);
      else return value;
    }
  
  int adjust(void)
    {
      unsigned long long ex = (D2ULL(value) >> 52) & 0x7ff;
      unsigned long long man = D2ULL(value) & ((((unsigned long long)1) << 52) - 1);
     
      if ((ex == 0x7ff || ex == 0) && man == 0) return -1;
      if (ex >= 0x400 + LOG_BASE)
	{
	  D2ULL(value) ^= (ex << 52);
	  
	  ex -= LOG_BASE;
	  exponent++;

	  D2ULL(value) ^= (ex << 52);
	  return 1;
	}
      if (ex <= 0x400 - LOG_BASE)
	{
	  D2ULL(value) ^= (ex << 52);

	  ex += LOG_BASE;
	  exponent--;

	  D2ULL(value) ^= (ex << 52);
	  return 1;
	}
      return 0;
    }

  scale_double & operator += (const scale_double & a);
  scale_double & operator -= (const scale_double & a);
  scale_double & operator *= (const scale_double & a);
  scale_double & operator /= (const scale_double & a);

  operator const double() {return get_double();} 
};

scale_double operator - (const scale_double & x)
{
  scale_double that(-x.value);
  that.exponent = x.exponent;
  return that;
}


scale_double operator * (const scale_double & x,const scale_double & y)
{
  scale_double that(x.value * y.value);
  that.exponent = x.exponent + y.exponent;
  that.adjust();

  return that;
}

scale_double operator + (const scale_double & x,const scale_double & y)
{
  scale_double that;
  double temp;
  unsigned long long ex;

  if (x.value == 0) return y;
  if (y.value == 0) return x;

  if (x.exponent == y.exponent)
    {
      that.exponent = x.exponent;
      that.value = x.value + y.value;
    }
  else if (x.exponent == y.exponent + 1)
    {
      temp = y.value;
      ex = (D2ULL(temp) >> 52) & 0x7ff;
      
      D2ULL(temp) ^= (ex << 52);      
      ex -= LOG_BASE;      
      D2ULL(temp) ^= (ex << 52);

      that.exponent = x.exponent;
      that.value = x.value + temp;
    }
  else if (x.exponent + 1 == y.exponent)
    {
      temp = x.value;
      ex = (D2ULL(temp) >> 52) & 0x7ff;
      
      D2ULL(temp) ^= (ex << 52);      
      ex -= LOG_BASE;      
      D2ULL(temp) ^= (ex << 52);

      that.exponent = y.exponent;
      that.value = y.value + temp;
    }
  else if (x.exponent > y.exponent) {that = x;}
  else if (x.exponent < y.exponent) {that = y;}

  that.adjust();
  return that;
}

scale_double operator - (const scale_double & x,const scale_double & y)
{
  return x + (- y);
}

scale_double operator / (const scale_double & x,const scale_double & y)
{
  if (y.value == 0) {printf("scale_double: divide by zero\n");exit(1);}
  if (x.value == 0) return 0;
  scale_double temp(1/y.value);
  temp.exponent = -y.exponent;
  return x * temp;
}

scale_double & scale_double::operator += (const scale_double & x)
{
  *this = *this + x;
  return *this;
}  

scale_double & scale_double::operator -= (const scale_double & x)
{
  *this = *this - x;
  return *this;
}  

scale_double & scale_double::operator *= (const scale_double & x)
{
  *this = *this * x;
  return *this;
}  

scale_double & scale_double::operator /= (const scale_double & x)
{
  *this = *this / x;
  return *this;
}  

scale_double operator == (const scale_double & x,const scale_double & y)
{
  return (x.exponent == y.exponent) && (x.value == y.value);
}

scale_double operator != (const scale_double & x,const scale_double & y)
{
  return !((x.exponent == y.exponent) && (x.value == y.value));
}


scale_double operator <= (const scale_double & x,const scale_double & y)
{
  if (x.value == 0) {return 0 <= y.value;}
  if (y.value == 0) {return x.value <= 0;}
  if (x.value < 0 && y.value > 0) {return 1;}
  if (x.value > 0 && y.value < 0) {return 0;}
  if (x.value > 0 && y.value > 0) 
    {
      if (x.exponent < y.exponent) return 1;
      if (x.exponent > y.exponent) return 0;
      return x.value <= y.value;
    }
  if (x.value < 0 && y.value < 0) 
    {
      if (x.exponent < y.exponent) return 1;
      if (x.exponent > y.exponent) return 0;
      return x.value <= y.value;
    }
  printf("scale_double operator <=: How did we get here?\n");
  return -1;
}

scale_double operator >= (const scale_double & x,const scale_double & y)
{
  return y <= x;
}

scale_double operator < (const scale_double & x,const scale_double & y)
{
  return (x <= y) && (x!=y);
}

scale_double operator > (const scale_double & x,const scale_double & y)
{
  return (x >= y) && (x!=y);
}

double log(const scale_double & x)
{
  return log(x.value) + log(2.0) * LOG_BASE * x.exponent;
}

double log2(const scale_double & x)
{
  return log2(x.value) + LOG_BASE * x.exponent;
}
