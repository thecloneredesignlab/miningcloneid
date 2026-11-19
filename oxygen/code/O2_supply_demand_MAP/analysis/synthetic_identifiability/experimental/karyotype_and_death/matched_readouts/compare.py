"""Summarize the replayed intervals and sampled two-parameter likelihood regions."""
from itertools import combinations
import math

import numpy as np
from scipy.spatial import Delaunay

from modal import HERE,TARGETS,case,read,stable,write
from profile import PAIR_LR,base

MODES=("stock_baseline","stock_death","stock_a7k")
TITLES=("Common baseline","+ Death proportions","+ 1% A7K karyotypes")


def pair_rows(mode,pa,pb):
    r,x=base(mode)
    pts={(float(r[f"x_{pa}"]),float(r[f"x_{pb}"])):(0.,"reference")}
    pair=read(HERE/f"pair_grid_{mode}.tsv")
    pairs=[t for t in pair if (t["parameter_a"],t["parameter_b"])==(pa,pb)]
    assert pairs
    for t in pairs+read(HERE/f"profile_points_{mode}.tsv"):
        if float(t["donor_margin"])<0 or not np.isfinite(float(t["loss"])):
            continue
        lr=float(t["loss"])-float(r["loss"])
        if not stable(t) and lr>PAIR_LR:
            continue
        xy=(round(float(t[f"x_{pa}"]),11),round(float(t[f"x_{pb}"]),11))
        if xy not in pts or pts[xy][0]>lr:
            pts[xy]=(lr,"joint_mesh" if t in pairs else "replayed_feasible")
    xy=np.array(list(pts))
    lr=np.array([a[0] for a in pts.values()])
    selected=lr<=PAIR_LR
    return xy,lr,selected,pairs


def area(tri,xy,lr):
    ans=0.
    for j in tri.simplices:
        vertices=zip(xy[j],lr[j])
        polygon=[]
        walk=list(vertices)
        for (x,a),(y,b) in zip(walk,walk[1:]+walk[:1]):
            if a<=PAIR_LR:
                polygon.append(x)
            if (a<=PAIR_LR)!=(b<=PAIR_LR):
                polygon.append(x+(y-x)*(PAIR_LR-a)/(b-a))
        if len(polygon)>=3:
            p=np.array(polygon)
            ans+=abs(np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1)))/2
    return ans


def joint_summary():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    options=list(combinations(TARGETS,2))
    fig,axs=plt.subplots(3,3,figsize=(12,10),layout="constrained")
    result=[]
    for i,(pa,pb) in enumerate(options):
        for j,mode in enumerate(MODES):
            xy,lr,supported,grid=pair_rows(mode,pa,pb)
            tri=Delaunay(xy)
            ax=axs[i,j]
            masked=mtri.Triangulation(xy[:,0],xy[:,1],tri.simplices)
            ax.tricontourf(masked,lr,levels=[min(-.05,lr.min()-.005),PAIR_LR],
                           colors=["#b6d8df"],alpha=.9)
            ax.tricontour(masked,lr,levels=[PAIR_LR],colors=["#176371"],linewidths=1.2)
            ax.scatter(xy[~supported,0],xy[~supported,1],s=7,c="#aeb7bf",alpha=.65)
            ax.scatter(xy[supported,0],xy[supported,1],s=9,c="#247c89",alpha=.8)
            r,_=base(mode)
            ax.scatter(float(r[f"x_{pa}"]),float(r[f"x_{pb}"]),marker="D",s=42,
                       c="#973e26",zorder=9,edgecolors="white",linewidths=.4)
            if i==0:
                ax.set_title(TITLES[j],fontsize=12)
            ax.set_xlabel(pa,fontsize=10)
            ax.set_ylabel(pb,fontsize=10)
            ax.set_xlim(-.02,1.02); ax.set_ylim(-.02,1.02)
            ax.set_xticks([0,.25,.5,.75,1]); ax.set_yticks([0,.25,.5,.75,1])
            ax.grid(alpha=.18,lw=.4)
            good=xy[supported]
            result.append(dict(mode=mode,parameter_a=pa,parameter_b=pb,mesh_fits=len(grid),
                    mesh_feasible=sum(float(t["donor_margin"])>=0. for t in grid),
                    mesh_stationary=sum(stable(t) for t in grid),
                    supported_vertices=len(good),min_a=good[:,0].min(),max_a=good[:,0].max(),
                    min_b=good[:,1].min(),max_b=good[:,1].max(),
                    interpolated_normalized_area=area(tri,xy,lr)))
    fig.suptitle("Pairwise 95% profile-likelihood support",fontweight="bold")
    fig.supxlabel("Registered fitted coordinates. Shading interpolates evaluated feasible points; diamond = best fit.",
                  fontsize=10)
    fig.savefig(HERE/"pairwise_support.png",dpi=165)
    plt.close(fig)
    write(HERE/"pairwise_support_summary.tsv",result)


def one_summary():
    baseline={r["parameter"]:r for r in read(HERE/f"intervals_{MODES[0]}.tsv")}
    report=[]
    for m in MODES:
        metrics=HERE/f"recovery_metrics_{m}.tsv"
        recovery=read(metrics) if metrics.exists() else []
        for row in read(HERE/f"intervals_{m}.tsv"):
            s=row["parameter"]
            old=baseline[s]
            width=float(row["x_upper"])-float(row["x_lower"])
            first=float(old["x_upper"])-float(old["x_lower"])
            recovered=[r for r in recovery if r["parameter"]==s]
            if row["precision_gate"]=="outside_supported":
                decision="outside_profile"
            elif len(recovered)==2 and all(int(r["allocated"])==100 for r in recovered):
                decision=("precision_pass" if all(float(r["rate"])>=.9 for r in recovered)
                           else "recovery_below_90_percent")
            else:
                decision=row["precision_gate"]
            report.append(dict(row,normalized_width=width,width_change_fraction=width/first-1.,
                recovery_rates="|".join(r["rate"] for r in recovered),verdict=decision))
    write(HERE/"profile_comparison.tsv",report)


if __name__=="__main__":
    one_summary()
    joint_summary()
