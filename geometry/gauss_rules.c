// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/gauss_rules.h"

void get_gauss_points(int n, real_t* points, real_t* weights)
{
  // FIXME: Need to look this over again.
  ASSERT(n >= 0);
  ASSERT(points != NULL);
  ASSERT(weights != NULL);

  int m = 2*(n-1); // order
  real_t p1, p2, p3;
  real_t pp, z;
  for (int i = 1; i <= n; ++i)
  {
    z = cos(M_PI * (i - 0.25) / (m + 0.5));

    while(1)
    {
      p1 = 1;
      p2 = 0;
      for (int j = 1; j <= m; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
      }
      // p1 is Legendre polynomial

      pp = m * (z*p1-p2) / (z*z - 1);

      if (ABS(p1/pp) < 2e-16) break;

      z = z - p1/pp;
    }

    z = ((1 - z) + p1/pp)/2;

    points[i-1] = z;
    points[n-i] = 1 - z;
    weights[i-1] = 1./(4*z*(1 - z)*pp*pp);
  }
}

void get_gauss_legendre_points(int n, real_t* points, real_t* weights)
{
  ASSERT(n >= 0);
  ASSERT(points != NULL);
  ASSERT(weights != NULL);

  // Tables of points and weights of degree up to 16.
  static const real_t x[16][8] =
   {{0.000000000000000},
    {0.577350269189626},
    {0.000000000000000, 0.774596669241483},
    {0.339981043584856, 0.861136311594053},
    {0.000000000000000, 0.538469310105683, 0.906179845938664},
    {0.238619186083197, 0.661209386466265, 0.932469514203152},
    {0.000000000000000, 0.405845151377397, 0.741531185599394, 0.949107912342759},
    {0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536},
    {0.000000000000000, 0.324253423403809, 0.613371432700590, 0.836031107326636, 0.968160239507626},
    {0.148874338981631, 0.433395394129247, 0.679409568299024, 0.865063366688985, 0.973906528517172},
    {0.000000000000000, 0.269543155952345, 0.519096129110681, 0.730152005574049, 0.887062599768095, 0.978228658146057},
    {0.125333408511469, 0.367831498918180, 0.587317954286617, 0.769902674194305, 0.904117256370475, 0.981560634246719},
    {0.000000000000000, 0.230458315955135, 0.448492751036447, 0.642349339440340, 0.801578090733310, 0.917598399222978, 0.984183054718588},
    {0.108054948707344, 0.319112368927890, 0.515248636358154, 0.687292904811685, 0.827201315069765, 0.928434883663574, 0.986283808696812},
    {0.000000000000000, 0.201194093997435, 0.394151347077563, 0.570972172608539, 0.724417731360170, 0.848206583410427, 0.937273392400706, 0.987992518020485},
    {0.095012509837637, 0.281603550779259, 0.458016777657227, 0.617876244402644, 0.755404408355003, 0.865631202387832, 0.944575023073233, 0.989400934991650}};

  static const real_t w[16][8] =
   {{2.0},
    {1.0},
    {0.888888888888889, 0.555555555555556},
    {0.652145154862546, 0.347854845137454},
    {0.568888888888889, 0.478628670499366, 0.236926885056189},
    {0.467913934572691, 0.360761573048139, 0.171324492379170},
    {0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168870},
    {0.362683783378362, 0.313706645877887, 0.222381034453374, 0.101228536290376},
    {0.330239355001260, 0.312347077040003, 0.260610696402935, 0.180648160694857, 0.081274388361574},
    {0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688},
    {0.272925086777901, 0.262804544510247, 0.233193764591990, 0.186290210927734, 0.125580369464905, 0.055668567116174},
    {0.249147045813403, 0.233492536538355, 0.203167426723066, 0.160078328543346, 0.106939325995318, 0.047175336386512},
    {0.232551553230874, 0.226283180262897, 0.207816047536889, 0.178145980761946, 0.138873510219787, 0.092121499837728, 0.040484004765316},
    {0.215263853463158, 0.205198463721290, 0.185538397477938, 0.157203167158194, 0.121518570687903, 0.080158087159760, 0.035119460331752},
    {0.202578241925561, 0.198431485327111, 0.186161000015562, 0.166269205816994, 0.139570677926154, 0.107159220467172, 0.070366047488108, 0.030753241996117},
    {0.189450610455069, 0.182603415044924, 0.169156519395003, 0.149595988816577, 0.124628971255534, 0.095158511682493, 0.062253523938648, 0.027152459411754}};

  int degree = 2 * (n - 1);
  if ((degree % 2) == 1)
  {
    points[0] = x[degree-1][0];
    weights[0] = w[degree-1][0];
    for (int i = 1; i < n; ++i)
    {
      points[2*i-1]  = -x[degree-1][i];
      weights[2*i-1] = w[degree-1][i];
      points[2*i]    = x[degree-1][i];
      weights[2*i]   = w[degree-1][i];
    }
  }
  else
  {
    for (int i = 0; i < n; ++i)
    {
      points[2*i]    = -x[degree-1][i];
      weights[2*i]   = w[degree-1][i];
      points[2*i+1]  = x[degree-1][i];
      weights[2*i+1] = w[degree-1][i];
    }
  }
}

noreturn void get_gauss_radau_points(int n, real_t* points, real_t* weights)
{
  ASSERT(n >= 2);
  ASSERT(points != NULL);
  ASSERT(weights != NULL);
  POLYMEC_NOT_IMPLEMENTED;
}

void get_gauss_lobatto_points(int n, real_t* points, real_t* weights)
{
  // FIXME: Have to re-examine this.
  ASSERT(n >= 2);
  ASSERT(points != NULL);
  ASSERT(weights != NULL);

  int degree = 2*(n-1);
  if (degree == 0)
  {
    points[0] = 0.5;
    weights[0] = 0.5;
  }
  else
  {
    points[0] = 0.;
    weights[0] = 0.;
    points[n] = 1.;
    weights[n] = 1.;
    if (n == 1) return;

    int m = (n - 1)/2, odd_n = n%2;

    if (!odd_n)
    {
      points[m+1] = 0.5;
      weights[m+1] = 1.0;
    }
    for (int i = 0; i < m; )
    {
      real_t y, z, d0, s0;
      z = cos(M_PI*(i + 1)/n);

      int k = 0;
      while (true)
      {
        // compute d0, s0 -- P'_p(z), P"_p(z)
        // (n+1)*P_{n+1}(z) = (2*n+1)*z*P_n(z)-n*P_{n-1}(z)
        // P'_{n+1}(z) = (2*n+1)* P_n(z)+P'_{n-1}(z)
        // P"_{n+1}(z) = (2*n+1)*P'_n(z)+P"_{n-1}(z)
        {
          real_t p0, p1, p2, d1;
          p2 = 1.;
          p1 = z;
          d0 = odd_n;
          d1 = 1 - odd_n;
          s0 = 0;
          for (int nn = 1; true; nn++)
          {
            p0 = ((2*nn+1)*z*p1 - nn*p2)/(nn + 1);
            if (nn%2 != odd_n)
            {
              d0 += (2*nn+1)*p1;
              s0 += (2*nn+1)*d1;
            }
            else
            {
              d1 += (2*nn+1)*p1;
            }
            if (nn == n - 1) break;
            p2 = p1;
            p1 = p0;
          }
        }

        if (ABS(d0/s0) < 2e-16) break;
        ++k;
        ASSERT(k < 6);

        z = z - d0/s0;
      }

      y = ((1 - z) + d0/s0)/2;

      points[++i] = y;
      weights[++i] = y;
      points[n-i] = 1 - y;
      weights[n-i] = y;
    }
  }
}

