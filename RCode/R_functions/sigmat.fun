sigmat<-function(t,s,p,pread=0.,dotheta=T,bad=-9.99)
     {
#
# t is temperature
# s is salinity
# p is sample pressure in decibars
# pread is reference pressure in decibars
#
      t[t==bad]<-NA
# change t to theta at pread if dotheta=T
      if (dotheta)
	      t<-thetap(t,s,p,pread,ni=3)
      s[s==bad]<-NA
      p<-pread*0.1
      t2<-t*t
      t3<-t2*t
      t4<-t3*t
      s2<-s*s
      s3half<-s**1.5
      p2<-p*p
#
      rw<-999.842594 + 6.793952e-2*t - 9.095290e-3*t2 +
              1.001685e-4*t3 - 1.120083e-6*t4 + 6.536332e-9*t**5
#
      rsto<-rw + (0.824493 - 4.0899e-3*t + 7.6438e-5*t2 -
              8.2467e-7*t3 + 5.3875e-9*t4) * s +
             (-5.72466e-3 + 1.0227e-4*t - 1.6546e-6*t2) * s3half +
             4.8314e-4 * s2
#
      kw<-19652.21 + 148.4206*t - 2.327105*t2 + 1.360477e-2*t3 -
              5.155288e-5*t4
#
      ksto<-kw + (54.6746 - 0.603459*t + 1.09987e-2*t2 -
              6.1670e-5*t3) * s +
             (7.944e-2 + 1.6483e-2*t - 5.3009e-4*t2) * s3half
#
      kstp<-ksto + (3.239908 + 1.43713e-3*t + 1.16092e-4*t2 -
              5.77905e-7*t3) * p +
             (2.2838e-3 - 1.0981e-5*t - 1.6078e-6*t2) * p * s +
             1.91075e-4 * p * s3half +
             (8.50935e-5 - 6.12293e-6*t + 5.2787e-8*t2) * p2 +
             (-9.9348e-7 + 2.0816e-8*t + 9.1697e-10*t2) * p2 * s
#
      density<-(rsto/(1.0-p/kstp))-1000.
      density
#
	}
