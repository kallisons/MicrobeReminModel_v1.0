alr<-function(t,s,p,bad=-9)
#
#-----------------------------------------------------------------------
#     this subroutine calculates adiabatic lapse(g) rate in the ocean
#     given in-situ temperature(t), salinity(s) and pressure(p).
#     the procedure used is bryden's(1973) polynomial formula.
#     **units**
#      temperature: degrees centigrade
#      salinity: permil
#      pressure: decibar
#      lapse rate: degrees per decibar
#     **reference**
#      h.l.bryden(1973); new polynomials for thermal expansion, adia-
#      batic temperature gradient and potential temperature of sea water.
#      deep-sea research, vol.20, pp 401-408.
#------------------------------------------------------------------------
#
	{
      t[t==bad]<-NA
      s[s==bad]<-NA
      p[p==bad]<-NA
      t2<-t*t
      t3<-t2*t
      sm35<-s-35.0
      gamma <-   0.35803e-1 + 
            0.85258e-2*t-0.68360e-4*t2+0.66228e-6*t3 +
            (0.18932e-2-0.42393e-4*t)*sm35 +
            (0.18741e-4-0.67795e-6*t+0.87330e-8*t2-0.54481e-10*t3)*p +
            (-0.11351e-6+0.27759e-8*t)*p*sm35 +
            (-0.46206e-9+0.18676e-10*t-0.21687e-12*t2)*p**2
#
      gamma*1.0e-3
#
	}
