import re,sys
for f in sys.argv[1:]:
    src=open(f,encoding='latin-1').read().split('\n')
    func='?'; depth=0; stack=[]; out=[]
    brace=0; loopdepth=[]
    for i,l in enumerate(src,1):
        m=re.match(r'^[A-Za-z][\w\s\*&:<>,]*?(\w+::)?(\w+)\s*\(',l)
        if m and not l.strip().endswith(';') and brace==0: func=(m.group(1) or '')+m.group(2)
        code=re.sub(r'//.*','',l)
        if re.search(r'\b(for|while)\s*\(',code):
            loopdepth=[d for d in loopdepth if d<=brace]
            loopdepth.append(brace)
            if len(loopdepth)>=2: out.append((func,i,len(loopdepth)))
        brace+=code.count('{')-code.count('}')
        loopdepth=[d for d in loopdepth if d<=brace]
    agg={}
    for fn,i,d in out:
        a=agg.setdefault(fn,[i,0]); a[1]=max(a[1],d)
    if agg: print(f, '  '.join(f"{k}(l{v[0]},prof{v[1]})" for k,v in agg.items()))
    else: print(f,'(sem lacos aninhados)')
