pro rz

k=1L
nrgem=201
nthgem=201

rdata=dblarr(nrgem,nthgem)
zdata=dblarr(nrgem,nthgem)
x=dblarr(nthgem)

infl='rz.dat'
dumchar=' '
openr, 1, infl
z=0.d
for k = 0,nrgem-1 do begin
   readf,1,x
   for j = 0,nthgem-1 do begin
      rdata(k,j) = x(j)
   endfor

   readf,1,x
   for j = 0,nthgem-1 do begin
      zdata(k,j) = x(j)
   endfor
endfor

close, 1

set_plot, 'ps'

;dum=   [0   1    2    3    4    5     6     7    8    9    10   11   12   13   14   15    16    17] 
RED =   [0,  1,   1,   0,   0,   .9,   1.,  .2,   .4,  .5,   .2, .0,  .5,  0,   .2,   .5,   0.2,  0.7]
GREEN = [0,  1,   0,   1,   0,   .2,   .6,  .7,   .4,  .7,   .4, .3,  0,   .5,  .2,   .5,   0.7,   0.2]
BLUE =  [0,  0,   0,   0,   1,   .9,   .1,  1.,   .4,  .6,  .75,  .3,  0,   0,   .5,   0,   1,   1]

;Load the first six elements of the color table:
TVLCT, 255 * RED, 255 * GREEN, 255 * BLUE

hong=2
lu=3
lan=4
huang=5

!p.multi=[0,1,1,0]
    !P.Charthick=5.0
    !P.Charsize=1.5
    !X.Charsize=1.0
    !Y.Charsize=1.3
    !P.Thick=1.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='rz.ps'
device, file=flnm, /color, xsize=15, ysize=25,xoffset=4,yoffset=2 

xmin=min(rdata)
xmax=max(rdata)
ymin=min(zdata)
ymax=max(zdata)
print,xmin,xmax,ymin,ymax

rskip=1
nrpts=nrgem/rskip

tskip=1
ntpts=nthgem/tskip

xdata=fltarr(ntpts+1)
ydata=fltarr(ntpts+1)


print,'rskip,nrpts, nrlast=', rskip,nrpts, nrpts*rskip
print,'tskip,ntpts, ntlast=', tskip,ntpts, ntpts*tskip

plot, [0],[0], linestyle=0, xrange=[xmin,xmax],yrange=[ymin,ymax] ;,xticks=2,xtickv=[0.25,0.45,0.65]
for i = 0,nrpts-1 do begin
   for j = 0,ntpts-1 do begin
      xdata(j) = rdata(i*rskip,j*tskip)
      ydata(j) = zdata(i*rskip,j*tskip)
   endfor
   oplot,xdata(0:ntpts-2),ydata(0:ntpts-2),psym=3
endfor

device, /close



return
end
