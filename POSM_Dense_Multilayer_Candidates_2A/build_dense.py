"""Dense 0.30/0.30 mm multilayer MOT coils, series routed through B.Cu.
Requires numpy scipy shapely matplotlib sexpdata. Run after search_dense.py.
Same PCB, multiple copper layers; no separate stacked PCB assemblies.
"""
from pathlib import Path
import json, uuid, math, re, csv, collections
import numpy as np
from scipy.optimize import brentq, root
from shapely.geometry import LineString,Point,Polygon
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sexpdata

ROOT=Path(__file__).resolve().parent; MU=4*np.pi*1e-7; RHO=1.724e-8
W=.3; P=.6008;T=.07; V=.6;DRILL=.3;AR=1.85;BW=.6
def uid():return str(uuid.uuid4())
def xy(q):return f'{100+q[0]:.6f} {100-q[1]:.6f}'
def polar(r,a):return np.array([r*np.cos(a),r*np.sin(a)])
def sp(a,n,hand,phase,step=.5):
    t=np.linspace(0,n*2*np.pi,round(n*360/step)+1)
    r=a+P*t/(2*np.pi);phi=phase+hand*t
    return np.round(np.c_[r*np.cos(phi),r*np.sin(phi)],6)
def path3(points,z):return np.c_[points,np.full(len(points),z)]
def field(paths,q):
    B=np.zeros(3)
    for v in paths:
        v=np.asarray(v)*1e-3;a=q-v[:-1];b=q-v[1:]
        aa=np.linalg.norm(a,axis=1);bb=np.linalg.norm(b,axis=1)
        f=(aa+bb)/(aa*bb*(aa*bb+np.sum(a*b,axis=1)))
        B+=MU/(4*np.pi)*np.sum(np.cross(a,b)*f[:,None],axis=0)
    return B
def tensor(paths,q):
    d=1e-6
    return np.column_stack([(field(paths,q+np.eye(3)[j]*d)-field(paths,q-np.eye(3)[j]*d))/(2*d) for j in range(3)])
def length(p):return float(np.linalg.norm(np.diff(p,axis=0),axis=1).sum())

configs=json.loads((ROOT/'dense_parameters.json').read_text());allresults=[]
fig,axes=plt.subplots(2,2,figsize=(12,11),layout='constrained')
for ci,c in enumerate(configs):
    h=c['target_height_mm'];K=c['coil_layers'];L=c['copper_layers'];zs=np.array(c['layer_z_mm']);zb=-(c['board_thickness_mm']-T)
    layernames=['F.Cu']+[f'In{i}.Cu' for i in range(1,L-1)]+['B.Cu']
    physical_z=-np.arange(L)*(c['board_thickness_mm']-T)/(L-1)
    half=max(c['body_size_mm']/2,math.ceil(2*(c['outer_end_radius_mm']+2.5))/2)
    plus=np.array([half+2.5,1.6]);minus=np.array([half+2.5,-1.6]);phases=np.pi/K+np.arange(K)*2*np.pi/K
    # Backside links rotate by one spoke as they travel inward. Unlike straight
    # radial links, they do not pass through the preceding layer's inner via.
    def layout(a2,step=.5):
        coils=[];necks=[];outs=[];bridges=[];avs=[];dvs=[];paths=[]
        rout=a2+c['turns_outer']*P;rv=rout+.9
        for j,phase in enumerate(phases):
            inner=sp(c['inner_start_radius_mm'],c['turns_inner'],1,phase,step)
            outer=sp(a2,c['turns_outer'],-1,phase,step)
            av=np.round(polar(AR,phase),6);dv=np.round(polar(rv,phase),6)
            neck=np.array([av,inner[0]]);bridge=np.array([inner[-1],outer[0]]);output=np.array([outer[-1],dv])
            coils.append([inner,outer]);necks.append(neck);bridges.append(bridge);outs.append(output);avs.append(av);dvs.append(dv)
        links=[]
        for j in range(K-1):
            rr=np.linspace(rv,AR,241);aa=phases[j]+(2*np.pi/K)*(rv-rr)/(rv-AR)
            pp=np.round(np.c_[rr*np.cos(aa),rr*np.sin(aa)],6);pp[0]=dvs[j];pp[-1]=avs[j+1];links.append(pp)
        # IN travels through the free upper half of the unoccupied final sector;
        # OUT leaves through the lower half. No other backside link uses it.
        anchor_in=np.array([half-1.2,1.6]);rp=np.linalg.norm(anchor_in);ap=np.arctan2(anchor_in[1],anchor_in[0]);rr=np.linspace(rp,AR,241);aa=ap+(phases[0]-ap)*(rp-rr)/(rp-AR)
        incoming=np.round(np.c_[rr*np.cos(aa),rr*np.sin(aa)],6);incoming[0]=anchor_in;incoming[-1]=avs[0];incoming=np.vstack([plus,incoming])
        anchor_out=np.array([half-1.2,-1.6]);end_angle=phases[-1]-2*np.pi;rm=np.linalg.norm(anchor_out);am=np.arctan2(anchor_out[1],anchor_out[0]);rr=np.linspace(rv,rm,161);aa=np.linspace(end_angle,am,len(rr))
        outgoing=np.round(np.c_[rr*np.cos(aa),rr*np.sin(aa)],6);outgoing[0]=dvs[-1];outgoing[-1]=anchor_out;outgoing=np.vstack([outgoing,minus])
        paths.append(path3(incoming,zb))
        for j in range(K):
            paths += [np.array([[*avs[j],zb],[*avs[j],zs[j]]]),path3(necks[j],zs[j]),path3(coils[j][0],zs[j]),path3(bridges[j],zs[j]),path3(coils[j][1],zs[j]),path3(outs[j],zs[j]),np.array([[*dvs[j],zs[j]],[*dvs[j],zb]])]
            if j<K-1:paths.append(path3(links[j],zb))
        paths += [path3(outgoing,zb),path3(np.array([minus,plus]),zb)]
        return dict(coils=coils,necks=necks,outs=outs,bridges=bridges,avs=avs,dvs=dvs,links=links,incoming=incoming,outgoing=outgoing,paths=paths)
    target=np.array([0.,0.,h*1e-3]);a20=c['outer_start_radius_mm']
    # Refine outer starting radius including every internal connection and a
    # representative short external wire closure; not only isolated windings.
    def bz(a2):return field(layout(a2,step=1.)['paths'],target)[2]
    lo=max(c['inner_start_radius_mm']+c['turns_inner']*P+W+1.,a20-2.);hi=a20+2.
    assert bz(lo)*bz(hi)<0,(h,bz(lo),bz(hi))
    a2=brentq(bz,lo,hi,xtol=1e-9)
    d=layout(a2);J=tensor(d['paths'],target);gu=abs(J[2,2])*100
    I20=20/gu;I10=I20/2;assert I20<2.
    bias=-field(d['paths'],target)*I20*1e4
    null=root(lambda q:field(d['paths'],q),target,tol=1e-10);assert null.success
    # Fabrication geometry checks on all layers, including all through vias.
    copper_polys=[];all_vias=d['avs']+d['dvs']; min_gap=1e6;min_via_clear=1e6
    per=720
    for j in range(K):
        inner,outer=d['coils'][j]
        cp=[LineString(q).buffer(W/2,quad_segs=12) for q in [inner,outer,d['necks'][j],d['bridges'][j],d['outs'][j]]]
        copper_polys.append(cp)
        for pts in [inner,outer]:
            assert LineString(pts).is_simple
            for k in range(0,len(pts)-per-1,per//4):
                aa=LineString(pts[k:min(k+per//4+1,len(pts)-per)])
                bb=LineString(pts[max(0,k+per-per//8):min(k+per+per//4+per//8+1,len(pts))])
                min_gap=min(min_gap,aa.distance(bb)-W)
        for vj,v in enumerate(all_vias):
            if vj in [j,j+K]:continue
            disc=Point(v).buffer(V/2)
            for poly in cp:min_via_clear=min(min_via_clear,disc.distance(poly))
    assert min_gap>=.3-2e-6,(h,min_gap)
    assert min_via_clear>=.3-2e-6,(h,min_via_clear)
    backpaths=[d['incoming'],*d['links'],d['outgoing']];back_polys=[LineString(q).buffer(BW/2) for q in backpaths]
    min_back_clear=1e6
    for j,poly in enumerate(back_polys):
        assert poly.is_valid and poly.geom_type=='Polygon'
        for other in back_polys[j+1:]:min_back_clear=min(min_back_clear,poly.distance(other))
        assert poly.distance(Point(0,0).buffer(1))>=.5-1e-6
    assert min_back_clear>=.3-2e-6,(h,min_back_clear)
    min_back_via_clear=1e6
    for pts,poly in zip(backpaths,back_polys):
        for v in all_vias:
            if min(np.linalg.norm(v-pts[0]),np.linalg.norm(v-pts[-1]))<1e-5:continue
            min_back_via_clear=min(min_back_via_clear,poly.distance(Point(v).buffer(V/2)))
    assert min_back_via_clear>=.3-2e-6,(h,min_back_via_clear)
    outline=np.array([[-half,-half],[half,-half],[half,-5],[half+8,-5],[half+8,5],[half,5],[half,half],[-half,half],[-half,-half]])
    allowed=Polygon(outline).buffer(-.5).difference(Point(0,0).buffer(1.5))
    for poly in [q for layer in copper_polys for q in layer]+back_polys+[Point(v).buffer(V/2) for v in all_vias]+[Point(v).buffer(1.1) for v in [plus,minus]]:
        assert allowed.buffer(2e-6).covers(poly),(h,'copper leaves 0.5mm inset of board or cutout')
    # DC resistance includes all active copper and each used barrel length.
    coil_len=sum(length(q) for pair in d['coils'] for q in pair)
    link_len=sum(length(q) for q in d['necks']+d['bridges']+d['outs'])
    back_len=sum(length(q) for q in backpaths)
    Rcoil=RHO*coil_len*1e-3/(W*1e-3*T*1e-3)
    Rlinks=RHO*1e-3/(T*1e-3)*(link_len/(W*1e-3)+back_len/(BW*1e-3))
    Rvia=RHO*sum(2*abs(zb-z)*1e-3 for z in zs)/(np.pi*DRILL*1e-3*25e-6)
    R=Rcoil+Rlinks+Rvia
    name=f'DENSE_{int(h):02d}mm_{L}Cu_{K}WINDING_LAYERS';folder=ROOT/name;folder.mkdir(exist_ok=True)
    text=['(kicad_pcb (version 20221018) (generator posm_dense_multilayer)',f'(general (thickness {c["board_thickness_mm"]}))','(paper "A4")','(layers']
    for li,ln in enumerate(layernames):text.append(f'({31 if li==L-1 else li} "{ln}" signal)')
    text+=['(36 "B.SilkS" user "b.silkscreen") (37 "F.SilkS" user "f.silkscreen") (38 "B.Mask" user) (39 "F.Mask" user) (40 "Dwgs.User" user "user.drawings") (44 "Edge.Cuts" user) (49 "F.Fab" user))','(setup (stackup']
    dielectric=(c['board_thickness_mm']-L*T)/(L-1)
    for li,ln in enumerate(layernames):
        text.append(f'(layer "{ln}" (type "copper") (thickness {T}))')
        if li<L-1:text.append(f'(layer "dielectric {li+1}" (type "core") (thickness {dielectric:.9f}) (material "TO_BE_CONFIRMED") (epsilon_r 4.5) (loss_tangent 0.02))')
    text+=[') (pad_to_mask_clearance 0.05))','(net 0 "")','(net 1 "MOT_DISTRIBUTED_WINDING")']
    graph=[];n_arc=0;n_seg=0;arc_error=0.
    def group(label,ids):text.append(f'(group "{label}" (id {uid()}) (members '+ ' '.join('"'+u+'"' for u in ids)+'))')
    def segments(label,pts,width,layer,z):
        global n_seg
        ids=[]
        for a,b in zip(pts[:-1],pts[1:]):
            u=uid();ids.append(u);text.append(f'(segment (start {xy(a)}) (end {xy(b)}) (width {width}) (layer "{layer}") (net 1) (tstamp {u}))');graph.append((tuple([*a,z]),tuple([*b,z])));n_seg+=1
        group(label,ids)
    def arcs(label,pts,layer,z):
        global n_arc,arc_error
        ids=[]
        # Every 10 samples = 5 degrees. The middle sample sets the arc exactly.
        for k in range(0,len(pts)-1,10):
            a,m,b=pts[k],pts[k+5],pts[k+10];u=uid();ids.append(u)
            text.append(f'(arc (start {xy(a)}) (mid {xy(m)}) (end {xy(b)}) (width {W}) (layer "{layer}") (net 1) (tstamp {u}))');graph.append((tuple([*a,z]),tuple([*b,z])));n_arc+=1
            aa=2*(m-a);bb=2*(b-a);ctr=np.linalg.solve(np.vstack([aa,bb]),[np.dot(m,m)-np.dot(a,a),np.dot(b,b)-np.dot(a,a)])
            rad=np.linalg.norm(a-ctr);err=np.max(abs(np.linalg.norm(pts[k:k+11]-ctr,axis=1)-rad));arc_error=max(arc_error,float(err))
        group(label,ids)
    for j in range(K):
        ln=layernames[j];z=float(zs[j])
        arcs(f'L{j+1}_INNER_{c["turns_inner"]}T',d['coils'][j][0],ln,z)
        arcs(f'L{j+1}_OUTER_{c["turns_outer"]}T',d['coils'][j][1],ln,z)
        segments(f'L{j+1}_A_NECK',d['necks'][j],W,ln,z);segments(f'L{j+1}_B_C_BRIDGE',d['bridges'][j],W,ln,z);segments(f'L{j+1}_D_ESCAPE',d['outs'][j],W,ln,z)
        for v in [d['avs'][j],d['dvs'][j]]:
            text.append(f'(via (at {xy(v)}) (size {V}) (drill {DRILL}) (layers "F.Cu" "B.Cu") (net 1) (tstamp {uid()}))');graph.append((tuple([*v,z]),tuple([*v,zb])))
    for j,pts in enumerate(backpaths):segments(f'BACK_SERIES_PATH_{j}',pts,BW,'B.Cu',zb)
    # Same accepted wire terminal and strain-relief arrangement as before.
    text += [f'(footprint "POSM:Wire_Terminals_with_Strain_Relief" (layer "F.Cu") (tstamp {uid()}) (at {100+half+2.5:.6f} 100)', '(attr through_hole exclude_from_pos_files exclude_from_bom)','(fp_text reference "J1" (at 0 0) (layer "F.Fab") (effects (font (size 0.8 0.8) (thickness 0.12))))','(fp_text value "MOT_CURRENT" (at 0 4.2) (layer "F.Fab") (effects (font (size 0.8 0.8) (thickness 0.12))))','(pad "1" thru_hole circle (at 0 -1.6) (size 2.2 2.2) (drill 1) (layers "*.Cu" "*.Mask") (net 1 "MOT_DISTRIBUTED_WINDING"))','(pad "2" thru_hole circle (at 0 1.6) (size 2.2 2.2) (drill 1) (layers "*.Cu" "*.Mask") (net 1 "MOT_DISTRIBUTED_WINDING"))','(pad "" np_thru_hole circle (at 3 -2) (size 2 2) (drill 2) (layers "*.Cu" "*.Mask"))','(pad "" np_thru_hole circle (at 3 2) (size 2 2) (drill 2) (layers "*.Cu" "*.Mask"))',')']
    for label,y in [('IN',3.5),('OUT',-3.5)]:text.append(f'(gr_text "{label}" (at {xy([half+2.5,y])}) (layer "F.SilkS") (tstamp {uid()}) (effects (font (size 0.8 0.8) (thickness 0.12))))')
    outline=np.array([[-half,-half],[half,-half],[half,-5],[half+8,-5],[half+8,5],[half,5],[half,half],[-half,half],[-half,-half]])
    for a,b in zip(outline[:-1],outline[1:]):text.append(f'(gr_line (start {xy(a)}) (end {xy(b)}) (layer "Edge.Cuts") (width 0.05) (tstamp {uid()}))')
    text.append(f'(gr_circle (center 100 100) (end 101 100) (layer "Edge.Cuts") (width 0.05) (fill none) (tstamp {uid()}))')
    text.append(f'(gr_text "{name} | 0.30mm trace / >=0.30mm clear gap" (at 100 {100+half-2}) (layer "Dwgs.User") (tstamp {uid()}) (effects (font (size 0.8 0.8) (thickness 0.12))))');text.append(')')
    board=folder/(name+'.kicad_pcb');board.write_text('\n'.join(text))
    parsed=sexpdata.loads(board.read_text());assert len([s for s in parsed if isinstance(s,list) and str(s[0])=='arc'])==n_arc
    assert arc_error<.0001,(h,arc_error)
    deg=collections.Counter();adj=collections.defaultdict(list)
    for a,b in graph:
        a=tuple(np.round(a,6));b=tuple(np.round(b,6));deg[a]+=1;deg[b]+=1;adj[a].append(b);adj[b].append(a)
    assert sum(x==1 for x in deg.values())==2 and all(x in [1,2] for x in deg.values())
    seen=set();todo=[next(iter(adj))]
    while todo:
        v=todo.pop()
        if v not in seen:seen.add(v);todo.extend(adj[v])
    assert len(seen)==len(deg)
    assert abs(np.trace(J))<1e-7 and np.max(abs(J-J.T))<1e-7
    result={**c,'name':name,'outer_start_radius_mm':a2,'outer_end_radius_mm':a2+c['turns_outer']*P,'body_size_mm':2*half,'overall_board_size_mm':[2*half+8,2*half],'current_10G_A':I10,'current_20G_A':I20,'gradient_per_A_routed_G_per_cm':gu,'gradient_at_2A_G_per_cm':2*gu,'total_R20_ohm':R,'P_at10G_W':I10**2*R,'P_at20G_W':I20**2*R,'V_at20G_V':I20*R,'compensation_bias_at20G_G':bias.tolist(),'uncompensated_null_mm':(null.x*1000).tolist(),'gradient_tensor_at20G_G_per_cm':(J*I20*100).tolist(),'native_arc_count':n_arc,'native_segment_count':n_seg,'max_arc_fit_error_mm':arc_error,'minimum_turn_gap_mm':min_gap,'minimum_unrelated_via_to_active_layer_copper_mm':min_via_clear,'minimum_backside_path_clearance_mm':min_back_clear,'via_diameter_mm':V,'via_drill_mm':DRILL,'inner_via_radius_mm':AR,'back_routing_width_mm':BW,'backplane_z_mm':zb,'layer_names':layernames,'physical_layer_z_mm':physical_z.tolist(),'total_turns_per_group':[K*c['turns_inner'],K*c['turns_outer']],'series_connections':[{ 'from_layer':layernames[j],'outer_via':d['dvs'][j].tolist(),'to_layer':layernames[j+1],'inner_via':d['avs'][j+1].tolist()} for j in range(K-1)],'field_model_note':'closed centerline-filament model, all active layers, barrel lengths, backside series links, and straight local external terminal closure; ambient fields and actual external wires excluded'}
    (folder/'parameters_and_checks.json').write_text(json.dumps(result,indent=2)+'\n');allresults.append(result)
    # Explicit stackup proposal; all copper layers are 70um and evenly spaced
    # in midplane for this model. Fabricator must confirm a realizable stackup.
    with (folder/'proposed_stackup.csv').open('w') as f:
        wr=csv.writer(f);wr.writerow(['layer','copper_midplane_z_mm','copper_thickness_mm','purpose'])
        for j,ln in enumerate(layernames):wr.writerow([ln,float(physical_z[j]),T,'winding pair' if j<K else ('series routing' if ln=='B.Cu' else 'empty')])
    # First winding layer and complete backside routing shown separately.
    def draw(ax,back=False):
        ax.plot(outline[:,0],outline[:,1],color='#444',lw=1)
        if not back:
            for q,col in zip(d['coils'][0],['#b75a1b','#167d91']):
                poly=LineString(q).buffer(W/2);ax.add_patch(plt.Polygon(np.asarray(poly.exterior.coords),facecolor=col,edgecolor='none'))
            for pts in [d['necks'][0],d['bridges'][0],d['outs'][0]]:
                poly=LineString(pts).buffer(W/2);ax.add_patch(plt.Polygon(np.asarray(poly.exterior.coords),facecolor='#b75a1b',edgecolor='none'))
        else:
            for poly in back_polys:ax.add_patch(plt.Polygon(np.asarray(poly.exterior.coords),facecolor='#7256a4',edgecolor='none'))
        for v in all_vias:
            ax.add_patch(plt.Circle(v,V/2,color='#bca44b'));ax.add_patch(plt.Circle(v,DRILL/2,color='white'))
        for v in [plus,minus]:ax.add_patch(plt.Circle(v,1.1,color='#bca44b'));ax.add_patch(plt.Circle(v,.5,color='white'))
        for yy in [-2,2]:ax.add_patch(plt.Circle([half+5.5,yy],1,facecolor='white',edgecolor='#444'))
        ax.add_patch(plt.Circle([0,0],1,facecolor='white',edgecolor='#444'))
        ax.set(aspect='equal',xlim=(-half-1,half+9),ylim=(-half-1,half+1),xlabel='x (mm)',ylabel='y (mm)');ax.grid(alpha=.1)
    draw(axes.flat[ci]);axes.flat[ci].set_title(f'{h:g} mm | {c["turns_inner"]}+{c["turns_outer"]} turns x {K} layers\n10-20 G/cm: {I10:.2f}-{I20:.2f} A | 0.30 / 0.30 mm',fontsize=11)
    f,aa=plt.subplots(1,2,figsize=(13,6),layout='constrained');draw(aa[0]);draw(aa[1],True);aa[0].set_title('F.Cu winding pair (other winding layers rotated)');aa[1].set_title('B.Cu: all interlayer series connections');f.savefig(folder/'layout_front_back.png',dpi=170);plt.close(f)
    print(name,json.dumps({q:result[q] for q in ['current_10G_A','current_20G_A','P_at20G_W','body_size_mm','minimum_turn_gap_mm','minimum_unrelated_via_to_active_layer_copper_mm','minimum_backside_path_clearance_mm']}),flush=True)

fig.savefig(ROOT/'dense_comparison.png',dpi=170);plt.close(fig)
(ROOT/'dense_results.json').write_text(json.dumps(allresults,indent=2)+'\n')
with (ROOT/'comparison.csv').open('w') as f:
    wr=csv.writer(f);wr.writerow(['height_mm','N1_per_layer','N2_per_layer','winding_layers','total_copper_layers','I10_A','I20_A','P10_W','P20_W','V20_V','body_side_mm','overall_width_mm'])
    for d in allresults:wr.writerow([d['target_height_mm'],d['turns_inner'],d['turns_outer'],d['coil_layers'],d['copper_layers'],d['current_10G_A'],d['current_20G_A'],d['P_at10G_W'],d['P_at20G_W'],d['V_at20G_V'],d['body_size_mm'],d['overall_board_size_mm'][0]])
# Include the accepted standalone terminal footprint and per-board tables.
source=(ROOT/allresults[0]['name']/(allresults[0]['name']+'.kicad_pcb')).read_text();start=source.index('(footprint ');end=source.index('\n)\n(gr_text',start)+2;fp=source[start:end]
fp=fp.replace('(footprint "POSM:Wire_Terminals_with_Strain_Relief"','(footprint "Wire_Terminals_with_Strain_Relief" (version 20221018) (generator posm_wire_terminals)',1)
fp=re.sub(r' \(tstamp [^)]*\)','',fp,count=1);fp=re.sub(r' \(at [\d.-]+ [\d.-]+\)','',fp,count=1);fp=fp.replace('(fp_text reference "J1"','(fp_text reference "REF**"').replace(' (net 1 "MOT_DISTRIBUTED_WINDING")','')
(ROOT/'POSM.pretty').mkdir(exist_ok=True);(ROOT/'POSM.pretty/Wire_Terminals_with_Strain_Relief.kicad_mod').write_text(fp+'\n')
for d in allresults:
    folder=ROOT/d['name']
    (folder/'fp-lib-table').write_text('(fp_lib_table (lib (name "POSM") (type "KiCad") (uri "${KIPRJMOD}/../POSM.pretty") (options "") (descr "Included terminal footprint")))\n')
    project={'meta':{'filename':d['name']+'.kicad_pro','version':1},'board':{'design_settings':{'rules':{'min_copper_edge_clearance':0.5,'min_clearance':0.3,'min_track_width':0.3,'min_via_annular_width':0.15,'min_via_diameter':0.6,'min_through_hole_diameter':0.3}}}}
    (folder/(d['name']+'.kicad_pro')).write_text(json.dumps(project,indent=2)+'\n')
