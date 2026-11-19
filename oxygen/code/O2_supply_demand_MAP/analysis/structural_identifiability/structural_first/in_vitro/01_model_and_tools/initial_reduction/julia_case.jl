# A bounded full five-live-class rational model, not a two-state substitute.
using StructuralIdentifiability,Random,Logging
const NM=StructuralIdentifiability.Nemo
out=ARGS[1];conditional=ARGS[2]=="conditional_k";mkpath(out)
names=vcat(["x$i" for i=1:7],["y$i" for i=1:5],["bmax","p0","a","k","w","alpha","dH","K","s2","s4"])
R,variables=NM.polynomial_ring(NM.QQ,names);F=NM.fraction_field(R);T=typeof(variables[1]);x=variables[1:7];y=variables[8:12];pars=variables[13:22]
bmax,p0,a,k,w,alpha,dH,K,s2,s4=F.(pars)
if conditional
 bmax=F(3)/4;p0=F(1)/10000;a=F(1)/100;w=F(1)/1000;alpha=F(2)/5;dH=F(1)/50;K=F(4)/5;s2=F(3)/10;s4=F(4)/5
end
function bounded(expr)
 count=max(length(numerator(expr)),length(denominator(expr)))
 if count>200000
  open(joinpath(out,"construction_limit.tsv"),"w") do f;println(f,"status\tterms\nterm_budget_exceeded\t",count);end
  error("Registered polynomial term allowance exceeded")
 end
 expr
end
nodes=[22,44,66,88,110];A=fill(F(0),7,7);M=copy(A);h=K/(K+F(1)/5)
destination(N)=N<=32 ? 1 : N<=54 ? 2 : N<=76 ? 3 : N<=98 ? 4 : 5
for (j,N) in enumerate(nodes)
 println("BUILD_PARENT ",N);flush(stdout)
 b=bounded(bmax/(1+alpha*h*F(N)/44));d=bounded(dH*h*F(N)/44);p=bounded(p0+a*d/(k+d));t=max(0,min(1,(N-44)//44));s=(1-t)*s2+t*s4
 B=fill(F(0),5);live=F(0)
 for m=0:N
  mass=bounded(F(binomial(big(N),m))*bounded(p^m)*bounded((1-p)^(N-m))*s^m)
  for child in (N-m,N+m)
   if child>0
    i=destination(child);B[i]=bounded(B[i]+mass);live=bounded(live+mass)
   end
  end
 end
 for i=1:5;A[i,j]=bounded(b*(1-w)*B[i]);end
 A[destination(2*N),j]+=b*w;A[j,j]-=b+d;A[6,j]=d;A[7,j]=bounded(b*(1-w)*(2-live))
 for i=1:5;M[i,j]=A[i,j];end
 M[j,j]+=d;M[7,j]=A[7,j]
end
A[6,6]=-F(1)/10000;A[7,7]=-F(1)/10000
L=sum(x[1:5]);capacity=F(big(10)^12)
eq=Dict{T,typeof(F(0))}(x[i]=>bounded(sum(A[i,j]*x[j] for j=1:7)-L/capacity*sum(M[i,j]*x[j] for j=1:7)) for i=1:7)
outputs=Dict{T,typeof(F(0))}(y[i]=>F(x[i]) for i=1:5)
ode=StructuralIdentifiability.ODE{T}(x,y,eq,outputs,T[])
if conditional
 # Remove syntactically unused fixed-coordinate symbols from the ring. Without
 # this operation the constructor would list them as parameters despite no use.
 fixed=Dict(pars[i]=>NM.QQ([3//4,1//10000,1//100,1//50,1//1000,2//5,1//50,4//5,3//10,4//5][i]) for i in eachindex(pars) if i!=4)
 ode=StructuralIdentifiability.set_parameter_values(ode,fixed)
end
open(joinpath(out,"ODE_READY"),"w") do f;println(f,"Rational seven-state model constructed");end
open(joinpath(out,"constructed_equations.txt"),"w") do f;show(f,ode);end
Random.seed!(2026111342)
result=assess_identifiability(ode;funcs_to_check=ode.parameters,prob_threshold=.999,loglevel=Logging.Info)
open(joinpath(out,"tool_parameter_results.tsv"),"w") do f
 println(f,"parameter\tclassification\tcondition")
 for (p,status) in result;println(f,p,'\t',status,'\t',conditional ? "Other nine reduced coordinates fixed" : "Ten reduced coordinates free");end
end
open(joinpath(out,"GLOBAL_COMPLETE"),"w") do f;println(f,"Tool result requires case/assumption review, including any Wronskian warning.");end
