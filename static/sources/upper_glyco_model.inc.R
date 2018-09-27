# 2018-09-25 sokol@insa-toulouse.fr
# Upper glycolysis model taken from "CSB Lecture flux balance analysis" which cites E. Klipp, Systems Biology in Practice, 2005
upglyc=function(t, co, p) {
  # provide derivative vector for ODE solving

  # flux as function of conc and kin. parameters
  nu1=p["Vmax1"]*co["ATP"]*p["Glucose"]/(1+co["ATP"]/p["Katp1"]+p["Glucose"]/p["Kglucose1"]+co["ATP"]*p["Glucose"]/(p["Katp1"]*p["Kglucose1"]))
  nu2=p["k2"]*co["ATP"]*co["Gluc6P"]
  nu3=((p["Vfmax3"]/p["Kgluc6p3"])*co["Gluc6P"]-(p["Vrmax3"]/p["Kfruc6p3"])*co["Fruc6P"])/
    (1+co["Gluc6P"]/p["Kgluc6p3"]+co["Fruc6P"]/p["Kfruc6p3"])
  nu4=p["Vmax4"]*co["Fruc6P"]**2/(p["Kfruc6p4"]*(1+p["ka"]*(co["ATP"]/co["AMP"])**2)+co["Fruc6P"]**2)
  nu5=p["k5"]*co["Fruc16P2"]
  nu6=p["k6"]*co["ADP"]
  nu7=p["k7"]*co["ATP"]
  nu8=p["k8f"]*co["ATP"]*co["AMP"]-p["k8r"]*co["ADP"]**2

  # first derivatives in time of concentration vector co (named)
  Gluc6P=nu1-nu2-nu3
  Fruc6P=nu3-nu4
  Fruc16P2=nu4-nu5
  ATP=-nu1-nu2-nu4+nu6-nu7-nu8
  ADP=nu1+nu2+nu4-nu6+nu7+2*nu8
  AMP=-nu8
  list(c(Gluc6P, Fruc6P, Fruc16P2, ATP, ADP, AMP))
}
p=double(17)
names(p)=c(
 "Glucose",
 "Vmax1",
 "Katp1",
 "Kglucose1",
 "k2",
 "Vfmax3",
 "Vrmax3",
 "Kgluc6p3",
 "Kfruc6p3",
 "Vmax4",
 "Kfruc6p4",
 "ka",
 "k5",
 "k6",
 "k7",
 "k8f",
 "k8r")
p["Glucose"]=12.8174 # mM
p["Vmax1"]=1398 # mM/min ## 50.2747 mM/min
p["Katp1"]=0.1 # mM
p["Kglucose1"]=0.37 # mM
p["k2"]=2.26 # 1/mM/min
p["Vfmax3"]=140.282 # mM/min
p["Vrmax3"]=140.282 # mM/min
p["Kgluc6p3"]=0.80 # mM
p["Kfruc6p3"]=0.15 # mM
p["Vmax4"]=44.7287 # mM/min
p["Kfruc6p4"]=0.021 # mM^2
p["ka"]=0.15
p["k5"]=6.04622 # 1/min
p["k6"]=68.48 # 1/min
p["k7"]=3.21 # 1/min
p["k8f"]=432.9 # 1/mM/min
p["k8r"]=133.33 # 1/mM/min

co=c(Gluc6P=1, Fruc6P=1, Fruc16P2=1, ATP=1, ADP=1, AMP=1)
set.seed(7)
co=co*runif(co, 0.1, 5)

times=seq(0, 0.5, length.out = 101)
library(deSolve)
out=ode(y=co, times=times, func=upglyc, parms=p)
plot(out, xlab = "time", ylab = "-")


