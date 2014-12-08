thetap<-function(t,s,p,pr=0.0,ni=3,bad=-9)
#
#-----------------------------------------------------------------------
#     this subroutine calculates potential temperature for an arbitrary
#     reference pressure level, using fofonoff's(1977) fourth-order
#     runge-kutta integration algorithm. the number of integration steps
#     (ni) is to be specified, although for most purposes ni=1 should
#     give an accurate enough result.
#     ** units**
#       in-situ temperature(t): degrees centigrade
#       salinity(s): permil
#       pressure(p): decibars
#       reference pressure(pr): decibars
#       potential temperature(pt): degrees centigrade
#     ** reference**
#       n.p.fofonoff(1977); computation of potential temperature of sea-
#       water for an arbitrary reference pressure. deep-sea research,
#       vol.24, pp.489-491.
#    Functions needed
#	alr - adiabatic lapse rate
#-----------------------------------------------------------------------
#
 {
      t[t==bad]<-NA
      s[s==bad]<-NA
      p[p==bad]<-NA
      dp<-(pr-p)/ni
      t0<-t
for(n in 1:ni)
	{
	pc<-p+dp*(n-1)
	dt1<-dp*alr(t0,s,pc)
	t1<-t0+0.5*dt1
	dt2<-dp*alr(t1,s,pc+0.5*dp)
	q1<-dt1
	t2<-t1+0.292893219*(dt2-q1)
	dt3<-dp*alr(t2,s,pc+0.5*dp)
	q2<-0.58578644*dt2 + .121320344*q1
	t3<-t2+1.707106781*(dt3-q2)
	dt4<-dp*alr(t3,s,pc+dp)
	q3<-3.414213562*dt3-4.121320344*q2
	t4<-t3+(dt4-2.0*q3)/6.0
	t0<-t4
	}
      t4
 }

