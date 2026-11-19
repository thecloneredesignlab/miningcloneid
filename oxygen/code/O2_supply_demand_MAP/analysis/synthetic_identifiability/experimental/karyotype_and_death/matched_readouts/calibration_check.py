"""Check classification bound and dated-stock allocation sensitivities."""
from concurrent.futures import as_completed
from fit import run_one,take
from modal import HERE,case,pool,read,write,stable


def run():
    tasks=[]
    seeds={"stock_death_logwide":["stock_death"],
           "stock_death_total":["stock_death","stock_baseline"],
           "a7k":["stock_a7k","baseline"]}
    for m,labels in seeds.items():
        k,c,y=case(m)
        for label in labels:
            z=c.anchor(read(HERE/f"best_{label}.tsv")[0])
            tasks.append((m,label,z,{},130,32.))
        if m=="stock_death_logwide":
            znew=z.copy(); znew[-4]=.9; znew[-2]=max(0.,z[-2]-.3)
            tasks.append((m,"wider_initial_offset",znew,{},130,32.))
    with pool(2) as ex:
        rows=[f.result() for f in as_completed(ex.submit(run_one,t) for t in tasks)]
    for m in seeds:
        r=[z for z in rows if z["mode"]==m]
        exact=take(r)
        write(HERE/f"classification_check_{m}.tsv",r)
        print(m,exact["loss"],exact["stationarity"],stable(exact),flush=True)


if __name__=="__main__":
    run()
