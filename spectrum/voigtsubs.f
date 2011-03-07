
c voigt profile subroutines, written by Ted Williams
c using algorithm from Humlicek and Varghese.
c Free to redistribute but please retain the credits.

      subroutine voigt (wave,a,si,vgt,dvda)

c *******************************************************************
c
c  This subroutine evaluates a Voigt emission profile and its 
c  derivatives at a specified wavelength.  The parameters are:
c
c          a(1)  continuum level
c          a(2)  total line strength (area under Voigt)
c          a(3)  center wavelength
c          a(4)  Gaussian sigma component
c          a(5)  Lorentzian sigma
c
c  The full Gaussian sigma is the sum in quadrature of a(4) and
c  the instrumental gaussian sigma, si.
c
c *******************************************************************

      real a(5),dvda(5)

      parameter (root2=1.4142136e0)
      parameter (rootpi=1.7724539e0)


      sgsq = a(4)**2 + si**2
      sg = sqrt(sgsq)
      r2sg = 1.0 / (root2 * sg)
      x = (wave - a(3)) * r2sg
      y = a(5) * r2sg

      call voi(x,y,v,dvdx,dvdy)

      rpsg = r2sg / rootpi

c  this fixes the intensity problem (the reported intensities were too high)
c
      a2sg = a(2) / sg
c  
c     a2sg = a(2) * rpsg

      vgt = a(1) + a2sg * v

      dvda(1) = 1.0
      dvda(2) = rpsg * v
      dvda(3) = - a2sg * r2sg * dvdx
      dvda(4) = - a2sg * a(4) * (v + x * dvdx + y * dvdy) / sgsq
      dvda(5) = a2sg * r2sg * dvdy

      return
      end


c voi subroutine modified to give correct value for intensity

      subroutine voi (x,y,v,dvdx,dvdy)

C***********************************************************************
C
C     This routine computes the normalized Voigt function and its partial 
C     derivatives.  The integral of this function is the value of the 
C     Gaussian sigma of the Voigt function.
C
C     The computation is valid for  y > 0, and the maximum relative error
C     for the Voigt function is < 2.E-6
C
C     Subroutine adapted from:
C     J. Humlicek, J. Quant. Spectrosc. Radiat. Transfer 21, 309 (1979)
C     by Philip L. Varghese        <820315.0112>
C
C     Modified to compute derivatives by Ted Williams 5/21/92
C
C***********************************************************************

      real t(6), c(6), s(6)

      data t/0.314240376E+00, 0.947788391E+00, 0.159768264E+01,
     *       0.227950708E+01, 0.302063703E+01, 0.388972490E+01/
      data c/0.403621095E+00,-0.299993210E+00, 0.500980825E-02,
     *       0.399820281E-02,-0.956712138E-04, 0.199809468E-07/
      data s/0.555821146E+00, 0.922164680E-01,-0.619762811E-01,
     *       0.248076921E-02, 0.366661062E-04,-0.250346637E-06/

c          coefficients for non-normalized Voigt:
c      data c/0.101172805E01,-0.75197147E00,0.12557727E-01,
c     *       0.100220082E-01,-0.242068135E-03,0.500848061E-06/
c      data s/0.1393237E01,0.231152406E00,-0.155351466E00,
c     *       0.621836624E-02,0.919082986E-04,-0.627525958E-06/

      parameter (root2pi=2.50662827e0)


      v = 0.0
      dvdx = 0.0
      dvdy = 0.0

      y1 = y + 1.5
      y2 = y1 * y1

C  branch to region I or II depending on values of X and Y.

      if ((y.gt.0.85).or.(abs(x).lt.(18.1*y+1.65))) then

C  calculations for region I

         do i=1,6
            rm = x - t(i)
            dm = 1.0 / (rm * rm + y2)
            d1 = y1 * dm
            d2 = rm * dm
            d1d2 = d1 * d2
            rp = x + t(i)
            dp = 1.0 / (rp * rp + y2)
            d3 = y1 * dp
            d4 = rp * dp
            d3d4 = d3 * d4
            v  = v + c(i) * (d1 + d3) - s(i) * (d2 - d4)
            dvdx = dvdx - 2.0 *(c(i) * (d1d2 + d3d4)
     $             + s(i) * (0.5*(dm-dp) - d2*d2 + d4*d4))
            dvdy = dvdy + 2.0 * (s(i) * (d1d2 - d3d4)
     $             + c(i) * (0.5*(dm+dp) - d1*d1 - d3*d3))
         end do
      else

C  calculations for region II

         if (abs(x).lt.12.0) then
            v = exp(-x*x) / root2pi
            dvdx = -2.0 * x * v / root2pi
         end if

         y3 = y + 3.0
         do i=1,6
            rm = x - t(i)
            sm = rm * rm
            dm = 1.0 / (sm + y2)
            d1 = y1 * dm
            d2 = rm * dm
            tm = 1.0 / (sm + 2.25)
            d2tm = d2 * tm
            dmtm = dm * tm
            fm =  tm * (s(i)*y3*d2 + c(i)*(rm*d2 - 1.5*d1))
            rp = x + t(i)
            sp = rp * rp
            dp = 1.0 / (sp + y2)
            d3 = y1 * dp
            d4 = rp * dp
            tp = 1.0 / (sp + 2.25)
            d4tp = d4 * tp
            dptp = dp * tp
            fp = tp * (-s(i)*y3*d4 + c(i)*(rp*d4 - 1.5*d3))
            v = v + y * (fm + fp)
            dvdx = dvdx + y * (s(i)*y3*(dmtm - dptp)
     $           + 2.0 * (c(i)*(d2tm + d4tp)
     $           - fm*(d2 + rm*tm) - fp*(d4 + rp*tp)))
            dvdy = dvdy + fm + fp + y * (s(i)*(d2tm - d4tp)
     $           - 1.5*c(i)*(dmtm + dptp)
     $           - 2.0*(d1*fm + d3*fp))
         end do
      end if

      return
      end

