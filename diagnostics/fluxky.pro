pro fluxky,numsteps
ndim=50000L
scale = 1e4
tstart=0.0
tend=200.0

ip = 4
ipem = ip+6                  ; +6 for EM contribution
iphist = ip                ; for fluxa.ps

ipx = 7      ;5, 7, 9 for r/a=0.8, 0.85, 0.9 respectively

istart=numsteps-9000
iend=numsteps-1

ntor=11
nsubd=12

dt = 1.0
nskip = 10

rmaj=1069.65
a=361.46
cs=sqrt(1.42/2.0)      ;T=1.42 at r/a=0.8, 1.245 at 0.85, 1.03 at 0.9
rhoi=2*cs/(1.0*1.0)             ;mi=2, q=1, B=1
ly=98.61

csova=cs/a
csovrmaj=csova*a/rmaj

;Te(nr2)=1.06, a=361.46, ne(nr2)=0.99, cn0e=0.666,
;eflxgb=ne*cn0e*Te(nr2) *sqrt(Te(nr2)/mi) * (rhoi/a)**2
;rhoi/a=mi*sqrt(Te/mi)/q(1)/a = 4.028e4, sqrt(Te/mi)=0.728
;nu*vu*Tu = 1e20*3.095e5*1e3*1.6e-19 = 4.88e9
;pflxgb=eflxgb/Te(nr2)=7.79e-6
;nu*vu = 3.095e25

eflxgb = 4.03e4                 ;w/m^2
pflxgb = 2.41e20                ;1/m^2s
print,'eflxgb,pflxgb=',eflxgb,pflxgb
flxu=fltarr(12)
flxu(0:2) = pflxgb
flxu(3:5) = eflxgb
flxu(6:8) = pflxgb
flxu(9:11) = eflxgb

channel=strarr(6)
channel(0)='e particle flux in 1/(m^2*s)'
channel(1)='D particle flux in 1/(m^2*s)'
channel(2)='C particle flux in 1/(m^2*s)'
channel(3)='e heat flux in w/m^2'
channel(4)='D heat flux in w/m^2'
channel(5)='C heat flux in w/m^2'

kyrhoi = fltarr(ntor)
for j=0,ntor-1 do begin
   kyrhoi(j) = (j+1)*2*3.1415926/ly *rhoi
endfor

t0 = fltarr(ndim)
ytmp = fltarr(ntor,12,nsubd,ndim)
yav = fltarr(ntor,12,nsubd)
ytot = fltarr(12,nsubd,ndim)

y5 = fltarr(ntor,12,nsubd,ndim)
y5av = fltarr(ntor,12,nsubd)
y5tot = fltarr(12,nsubd,ndim)

x=fltarr(ntor*nsubd)
n5=numsteps
ntmp=n5
infl='fluxa'
openr, 1, infl
for k=0L,ntmp-1L do begin
   for l=0,11 do begin
      readf,1,t,x
      t0(k) = t
      for j = 0,ntor-1 do begin
         for i = 0,nsubd-1 do begin
            m = j*nsubd+i
            ytmp(j,l,i,k) = x(m)
         endfor
      endfor
   endfor
endfor
t0 = t0*csovrmaj
tend=min([tend,t0(ntmp-1)])
close, 1

for k=0,ntmp-1 do begin
;   yav(ip,k) = (ytmp(ip,2,k)+ytmp(ip,3,k)+ytmp(ip,4,k)+ytmp(ip,5,k))/4
;   yav(ip,k) = (ytmp(ip,0,k)+ytmp(ip,1,k))/2
endfor


print,'istart,iend=', istart,iend

for j = 0,ntor-1 do begin
   for i = 0,nsubd-1 do begin 
      dum = 0.
      dum1 = 0.
      for k=istart,iend do begin
         dum = dum+ytmp(j,ip,i,k)
         dum1 = dum1+ytmp(j,ip+6,i,k)         
      endfor
      yav(j,ip,i) = dum/(iend-istart+1)
      yav(j,ip+6,i) = dum1/(iend-istart+1)      
   endfor
endfor

for i = 0,nsubd-1 do begin 
   for k=0,ntmp-1 do begin
      dum = 0.
      dum1 = 0.      
      for j = 0,ntor-1 do begin
         dum = dum+ytmp(j,ip,i,k)
         dum1 = dum1+ytmp(j,ip+6,i,k)         
      endfor
      ytot(ip,i,k) = dum
      ytot(ip+6,i,k) = dum1      
   endfor
endfor

;print,'fluxa average=', dum
y5 = ytmp
y5av = yav
y5tot = ytot

set_plot, 'ps'

;Specify the red component of each color:
RED = [0, 1, 1, 0, 0, 1]
;Specify the green component of each color:
GREEN = [0, 1, 0, 1, 0, 1]
;Specify the blue component of each color:
BLUE = [0, 1, 0, 0, 1, 0]
;Load the first six elements of the color table:
TVLCT, 255 * RED, 255 * GREEN, 255 * BLUE

RED =   [0,  1,   1,   0,   0,   0,   .54,  .75,  1.,  1,   .49, .93, .98, .66,   .2,   .5,   0.2,  0.7]
GREEN = [0,  1,   0,   1,   0,   1,   .17,  .75,  .41, 0,   .99, .51, .5,  .66,  .2,   .5,   0.7,   0.2]
BLUE =  [0,  0,   0,   0,   1,   1,   .89,  .75,  .71, 1,   .0,  .93, .45, .66,   .5,   0,   1,   1]

;Load the first six elements of the color table:
TVLCT, 255 * RED, 255 * GREEN, 255 * BLUE

huang=1
hong=2
lu=3
lan=4
cyan=5
blueviolet=6
silver=7
hotpink=8
magenta=9
lawngreen=10
violet=11
salmon=12
darkgrey=13

ldsq=1.3
xsq=[ldsq,-ldsq,-ldsq,ldsq,ldsq]
ysq=[ldsq,ldsq,-ldsq,-ldsq,ldsq]

ldtr=1
xtr=[0,-1,1,0]
ytr=[1,-0.5,-0.5,1]

A=findgen(31)*!pi*2/30
xcir=1.0*cos(A)
ycir=1.0*sin(A)

!p.multi=[0,1,1,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.0
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='fluxa.ps'
device, file=flnm, /color,/encapsulated, xsize=21, ysize=21

xmax = max([max(t0)])
ymax = max([   max(y5tot(iphist,1:6,0:n5-1))   ])
ymin = min([   min(y5tot(iphist,1:6,0:n5-1))   ])
plot, t0(0:n5-1),y5tot(iphist,4,0:n5-1), title='', $
  xstyle=1, ystyle=1,xtitle='tC!ls!N/R',$
  ytitle='Q!li',$ 
  xrange=[0,xmax*1.1], $
  yrange=[ymin,ymax*1.2],linestyle=0,color=0,$
  charsize=2,charthick=4,position=[0.3,0.3,0.9,0.9]


oplot, t0(0:n5-1), y5tot(iphist,5,0:n5-1),linestyle=0, color=cyan                ;g1000a, nzgyro=20
oplot, t0(0:n5-1), y5tot(iphist,6,0:n5-1),linestyle=0, color=lu                  ;g1000b, nzgyro=100
oplot, t0(0:n5-1), y5tot(iphist,7,0:n5-1),linestyle=0, color=lan                 ;g1000c. nzgyro=3
oplot, t0(0:n5-1), y5tot(iphist,8,0:n5-1),linestyle=0, color=hong                ;g1000c. nzgyro=3
oplot, t0(0:n5-1), y5tot(iphist,9,0:n5-1),linestyle=0, color=magenta                 ;g1000c. nzgyro=3

lxa=0.3
rina=0.65
rova=fltarr(12)
for i=0,11 do begin
   rova(i)=lxa/nsubd*(i+0.5)+rina
endfor
openw,1,'flxhis.dat'
printf,1,'column 1 is time, column 2 to 9 for fluxes at r/a='
printf,1, format='(8G10.4)',  rova(3),rova(4),rova(5),rova(6),rova(7),rova(8),rova(9),rova(10)
printf,1,'respectively '
printf,1,' '
for k = 0,n5-1 do begin
   printf,1,format='(10G12.5)',t0(k),y5tot(iphist,3,k),y5tot(iphist,4,k),y5tot(iphist,5,k),y5tot(iphist,6,k),y5tot(iphist,7,k), $
          y5tot(iphist,8,k),y5tot(iphist,9,k),y5tot(iphist,10,k)
endfor

xst=2
yst = 55
yoff=5
;xyouts, xst,yst, '100', color=lu
;xyouts, xst,yst-yoff, '20', color=cyan
;xyouts, xst,yst-yoff*2, '5', color=salmon
;xyouts, xst,yst-yoff*3, '4', color=0
;xyouts, xst,yst-yoff*4, '3', color=lan

;xyouts, 500,yst-yoff*4, '1', color=hong
;xyouts, 500,yst-yoff*5, '05', color=13
;xyouts, 500,yst-yoff*6, '6', color=blueviolet
;xyouts, 500,yst-yoff*7, '8', color=magenta

device, /close




flnm='fluxky.ps'
device, file=flnm, /color,/encapsulated, xsize=21, ysize=21

xpos=0.38
ypos=6
yoff=1.7
xskp=0.03
yskp=-0.4
ymax=max(y5av(0:ntor-1,ip,ipx))
ymin=min(y5av(0:ntor-1,ip,ipx))

if(ymax gt 0.0 and abs(ymax) ge abs(ymin))then begin
   ymin=min(-0.1*ymax,ymin)
endif


print,' '
print,'ip,ipx= ',ip,ipx
print,channel(ip)
dumtot = 0.
dum = 0.
dum1 = 0.
for j = 0,ntor-1 do begin
   dum = (y5av(j,ip,ipx)+y5av(j,ip,ipx+1))/2.0
   dum1 = (y5av(j,ipem,ipx)+y5av(j,ipem,ipx+1))/2.0
   dumtot = dumtot+dum+dum1
   print,j,kyrhoi(j),(dum+dum1) ;*flxu(ip)
endfor
print,'total=', dumtot*flxu(ip)

usersym,xcir,ycir,color=0, thick=5
plot, kyrhoi(0:ntor-1),y5av(0:ntor-1,ip,ipx), title='', $
  xstyle=1, ystyle=1,xtitle='ky*rhoi',xrange=[0,1.5],yrange=[ymin*1.2,ymax*1.2], $
  ytitle='Q!li',psym=8,linestyle=8,color=0, $
  charsize=2,charthick=4,position=[0.3,0.3,0.9,0.9]
;oplot,[xpos],[ypos],psym=8,linestyle=8,color=0

return
end
