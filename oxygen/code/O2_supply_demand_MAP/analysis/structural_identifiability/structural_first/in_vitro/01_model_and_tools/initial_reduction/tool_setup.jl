using Pkg
Pkg.add(PackageSpec(name="StructuralIdentifiability",version="0.5.36"))
Pkg.precompile()
using StructuralIdentifiability,Random,Logging
out=ARGS[1];mkpath(out)
Random.seed!(2026111341)
identified=@ODEmodel(x'(t)=-a*x(t),y(t)=x(t))
r1=assess_identifiability(identified;prob_threshold=0.999,loglevel=Logging.Info)
unidentified=@ODEmodel(x'(t)=-(a+b)*x(t),y(t)=x(t))
r2=assess_identifiability(unidentified;prob_threshold=0.999,loglevel=Logging.Info)
combo=assess_identifiability(unidentified;funcs_to_check=[sum(unidentified.parameters)],prob_threshold=0.999,loglevel=Logging.Info)
@assert any(string(k)=="a" && v==:globally for (k,v) in r1)
@assert all(v==:nonidentifiable for (k,v) in r2 if string(k) in ("a","b"))
@assert all(v==:globally for (k,v) in combo)
open(joinpath(out,"tool_smoke.tsv"),"w") do f
 println(f,"case\tquantity\tclassification")
 for (name,result) in [("identifiable_decay",r1),("unidentifiable_sum",r2),("identifiable_sum_combination",combo)],(k,v) in result
  println(f,name,'\t',k,'\t',v)
 end
end
open(joinpath(out,"tool_versions.txt"),"w") do f
 println(f,"Julia ",VERSION);Pkg.status(;io=f,mode=Pkg.PKGMODE_MANIFEST)
end
open(joinpath(out,"PASS"),"w") do f;println(f,"Known identifiable and nonidentifiable examples passed; no project claim.");end
