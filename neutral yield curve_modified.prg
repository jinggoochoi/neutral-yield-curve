'written by choi(1510339)

close q

'!doml=1
'!iter=5000
'!ls=0.175162243'(estimate lambda_LW 참조)
!ls=0.05413984'(estimate lambda_LW 참조)

'if !doml=1 then
smpl 2001.1 2018.4

coef(4) a
coef(3) b
coef(3) g
coef(3) h
coef(7) s
'coef(3) e

genr y1=y(-1)
genr ystar1=ystar(-1)
genr deltaystar=(ystar-ystar(-1))-(ystar(-1)-ystar(-2))
genr lt1=lt(-1)
genr st1=st(-1)
genr ct1=ct(-1)

lt.hpf(100) lttrend
st.hpf(100) sttrend
ct.hpf(100) cttrend

genr sv1 = lt
genr sv2 = st
genr sv3 = ct
genr sv4 = lttrend
genr sv5 = sttrend
genr sv6 = cttrend

system m3

m3.append y=ystar+a(1)*(y1-ystar1)+(b(1)*lt(-1)+b(2)*st(-1)+b(3)*ct(-1))-(b(1)*(lttrend(-1)+h(1)*deltaystar)-b(2)*(sttrend(-1)+h(2)*deltaystar)-b(3)*(cttrend(-1)+h(3)*deltaystar))
'm3.append lt=a(2)*lt(-1)+(1-a(2))*lttrend
'm3.append st=a(3)*st(-1)+(1-a(3))*sttrend
'm3.append ct=a(4)*ct(-1)+(1-a(4))*cttrend
m3.append lt=a(2)*lt(-1)+(1-a(2))*(lttrend(-1)+h(1)*deltaystar)
m3.append st=a(3)*st(-1)+(1-a(3))*(sttrend(-1)+h(2)*deltaystar)
m3.append ct=a(4)*ct(-1)+(1-a(4))*(cttrend(-1)+h(3)*deltaystar)

m3.sur

show m3

system p3
p3.append lttrend=lttrend(-1)+h(1)*deltaystar
p3.append sttrend=sttrend(-1)+h(2)*deltaystar
p3.append cttrend=cttrend(-1)+h(3)*deltaystar
p3.sur
show p3

'set prior via data
vector(6) mp
mp(1) = lttrend(76)
mp(2) = sttrend(76)
mp(3) = cttrend(76)
mp(4) = lttrend(76)
mp(5) = sttrend(76)
mp(6) = cttrend(76)
               
'variance prior : variance of each trend

sym(6,6) vp
for !i=1 to 6
vp(!i,!i)=3
next

'set this value for pile-up problem
!sigmaL = 1.599036 'update depeneds on data
!sigmaS = 1.612594 'update depeneds on data
!sigmaC = 2.599014 'update depeneds on data
!sigmaY = 18.72893 'update depeneds on data

!lambdal=(m3.@se(2)/m3.@se(1))-!sigmaL/!sigmaY 
!lambdas=(m3.@se(3)/m3.@se(1))-!sigmaS/!sigmaY 
!lambdac=(m3.@se(4)/m3.@se(1))-!sigmaC/!sigmaY 

!lambdas1=p3.@se(2)/p3.@se(1)
!lambdac1=p3.@se(3)/p3.@se(1)

!s1=m3.@se(1)
!s2=((m3.@se(2))-!lambdal*m3.@se(1))^0.5
!s3=((m3.@se(3))-!lambdas*m3.@se(1))^0.5
!s4=((m3.@se(4))-!lambdac*m3.@se(1))^0.5
!s5=p3.@se(1)
!s6=p3.@se(2)
!s7=p3.@se(3)

'set initial parameter, by previous SUR estimation

param a(1) 0.649609 a(2) 0.373677 a(3) 0.516521 a(4) 0.152149 b(1) 0.111236 b(2) -0.083223 b(3) 0.100491  h(1) -3.742114 h(2) 3.997544 h(3) 7.984479
param s(1) !s1 s(2) !s2 s(3) !s3 s(4) !s4 s(5) !s5 s(6) !s6 s(7) !s7 'g(1) !lambdal g(2) !lambdas g(3) !lambdac

sspace e3

e3.append @signal y=ystar+a(1)*(y1-ystar1)+(b(1)*lt1+b(2)*st1+b(3)*ct1)-(b(1)*(sv4+h(1)*deltaystar)-b(2)*(sv5+h(2)*deltaystar)-b(3)*(sv6+h(3)*deltaystar))+[var=(s(1))^2] 
e3.append @signal lt=a(2)*lt1+(1-a(2))*(sv4+h(1)*deltaystar)+[var=(!lambdal*s(1)+s(2))^2]
e3.append @signal st=a(3)*st1+(1-a(3))*(sv5+h(2)*deltaystar)+[var=(!lambdas*s(1)+s(3))^2]
e3.append @signal ct=a(4)*ct1+(1-a(4))*(sv6+h(3)*deltaystar)+[var=(!lambdac*s(1)+s(4))^2]

e3.append @state sv1=sv1(-1)+h(1)*deltaystar+[var=(s(5))^2]
e3.append @state sv2=sv2(-1)+h(2)*deltaystar+[var=(!ls*s(5))^2]
e3.append @state sv3=sv3(-1)+h(3)*deltaystar+[var=(!lambdac1*s(5))^2]
e3.append @state sv4=sv1(-1)
e3.append @state sv5=sv2(-1)
e3.append @state sv6=sv3(-1)

'e3.append @evar cov(e1,e2) = (!ls*s(5))^d2
'e3.append @evar cov(e1,e3) = (!lambdac1*s(5))^2
'e3.append @evar cov(e2,e3) = (!ls*!lambdac1*s(5))^2

e3.append @mprior mp
e3.append @vprior vp
e3.ml(showopts,m=5000)

e3.makestates(t=smooth) sv1 sv2 sv3 sv4 sv5 sv6
e3.makestates(t=smoothse) se1 se2 se3 se4 se5 se6

show e3

'endif


