import os,sys,subprocess
V=sys.argv[1]
W='/dartfs-hpc/scratch/f007058/audit_ngc6791/pyeep_work'
VARF=f'/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var/varfeh+031.{V}'
vals=[l.split()[0] for l in open(VARF) if l.strip()]
feh,yprim,afe=float(vals[0]),float(vals[1]),float(vals[2])
z=float(subprocess.run(['/dartfs-hpc/rc/home/8/f007058/mcdsep/run/gs98getz',str(feh),str(afe)],
                       capture_output=True,text=True).stdout.strip())
tracks=sorted(os.listdir(f'{W}/mesa/{V}'))
os.makedirs(f'{W}/eeps/{V}',exist_ok=True); os.makedirs(f'{W}/isos/{V}',exist_ok=True)
lines=['#version string, max 8 characters','pyeep',
 '#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)',
 '   %.4f %.5e   %.2f        %.2f     0.0'%(yprim,z,feh,afe),
 '#data directories: 1) history files, 2) eeps, 3) isochrones',
 f'{W}/mesa/{V}', f'{W}/eeps/{V}', f'{W}/isos/{V}',
 '#read history_columns','my_history_columns_basic.list',
 '#specify tracks',str(len(tracks))]+tracks+[
 '#specify isochrones',f'isochrones_{V}.txt','min_max','linear','81','4.0e9','12.0e9']
open(f'{W}/run/in_{V}','w').write('\n'.join(lines)+'\n')
print(f'var {V}: Y={yprim:.4f} Z={z:.5e} [Fe/H]={feh:.3f} [a/Fe]={afe:.2f}  tracks={len(tracks)}')
