NB. Define mesh

load 'mesh.dat'

p =: px ,: py
t =: (t1 ,: t2) , t3

NB. Define local basis functions.

phi1 =: 3 : 0
'a b' =. y
1-(a+b)
)

phi2 =: 3 : 0
'a b' =. y
a
)

phi3 =: 3 : 0
'a b' =. y
b
)

bfuns =: phi1`phi2`phi3

NB. Quadrature points

X =: 0 1r2 1r2 ,: 1r2 0 1r2
W =: 1r6 1r6 1r6

NB. Local to global mapping

x1 =: (0{t){"1 p
x2 =: (1{t){"1 p
x3 =: (2{t){"1 p

A11 =: (0{x2)-(0{x1)
A12 =: (0{x3)-(0{x1)
A21 =: (1{x2)-(1{x1)
A22 =: (1{x3)-(1{x1)

detA =: (A11*A22)-(A12*A21)

invAt11 =: A22%detA
invAt21 =: (-A12)%detA
invAt12 =: (-A21)%detA
invAt22 =: A11%detA

evalbfun =: 3 : '(y{bfuns)`:6 X'
evalbfungrads =: 3 : '_2 ((y{bfuns)`:6 D. 1) \ , X'

nelems =: #0{t
nqp =: #0{X

mprod =: +/ . *

abs =: * *

assemble =: 3 : 0
ds=:0$1
is=:0$1
js=:0$1
for_jtr. (0,1,2) do.
 u =: (nelems,nqp) $ evalbfun jtr

 uhatx =: 0 { |: evalbfungrads jtr
 uhaty =: 1 { |: evalbfungrads jtr

 ux =: (invAt11 */ uhatx) + (invAt12 */ uhaty)
 uy =: (invAt21 */ uhatx) + (invAt22 */ uhaty)

 for_itr. (0,1,2) do.
  v =: (nelems,nqp) $ evalbfun itr

  vhatx =: 0 { |: evalbfungrads itr
  vhaty =: 1 { |: evalbfungrads itr

  vx =: (invAt11 */ vhatx) + (invAt12 */ vhaty)
  vy =: (invAt21 */ vhatx) + (invAt22 */ vhaty)
  
  ds=:ds,((((ux*vx)+(uy*vy)) mprod W) * abs(detA))
  is=:is,(itr{t)
  js=:js,(jtr{t)
 end.
end.
(": ds) fwrite 'output.dat'
(": is) fwrite 'outputi.dat'
(": js) fwrite 'outputj.dat'
)

assemble 1

exit''
