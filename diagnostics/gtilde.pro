pro gtilde
  nxsrc=100L
  nesrc=200L

  ildu=0
  vcut=15.d0

  g=dblarr(nxsrc,nesrc)
  ga=dblarr(nxsrc,nesrc)  
  phasvol=dblarr(nxsrc,nesrc)
  infl='gtilde'
  openr, 1, infl
  for i = 0,nxsrc-1 do begin
     for j = 0,nesrc-1 do begin
        readf,1,k,k,x1,x2,x3,x4
        g(i,j)=x3
        ga(i,j)=x4        
        phasvol(i,j)=x2
     endfor
  endfor
  close,1

  goto, endread
  g0=dblarr(nxsrc,nesrc)
  infl='gtildi50000'
  openr, 1, infl
  for i = 0,nxsrc-1 do begin
     for j = 0,nesrc-1 do begin
        readf,1,k,k,x1,x2,x3,x4
        g0(i,j)=x3
     endfor
  endfor
  close,1

  endread:
  
;radial plot
  r=dblarr(nxsrc)
  tr=dblarr(nxsrc)       ;T(r)
  r0a=0.45
  lxa=0.8
  rin=0.05
  rout=0.85
  kappat=6.96
  wt=0.3
  a=0.36d0
  lref=1.0d0

  for i=0,nxsrc-1 do begin
     r(i)=rin+(rout-rin)/double(nxsrc)*(double(i)+0.5d0)         ;r/a
     tr(i)= exp(-kappat*wt*a/lref*tanh((r(i)-r0a)/wt))  
; OB    print,i,tr(i)
  endfor

  e=dblarr(nesrc)
  ecut=15.0
  for i=0,nesrc-1 do begin
     e(i)=(i+0.5)*ecut/nesrc
  endfor

  
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


!p.multi=[0,3,3,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.3
    !P.Thick=10.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='gx.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=1 
for k=0,8 do begin
   ieplt=1+k*20
   ymax=max(g(0:nxsrc-1,ieplt))

;   gx(0:nxsrc-1,ieplt)=gx(0:nxsrc-1,ieplt)/gxmax*ymax

   plot, r(0:nxsrc-1), g(0:nxsrc-1,ieplt),title=String(ieplt,format='(I3)'), xtitle='r/a',ytitle='',linestyle=0
   oplot,r(0:nxsrc-1),ga(0:nxsrc-1,ieplt),linestyle=0,color=lu
;   oplot,r(0:nxsrc-1),g0(0:nxsrc-1,ieplt),linestyle=1,color=hong
endfor
device,/close


flnm='geps.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=1 

i=50
j=100
for k=0,8 do begin
   ixplt=5+k*11
   ymax=max([   max(g(ixplt,i:j)),  max(ga(ixplt,i:j))   ])

;   gx(0:nxsrc-1,ieplt)=gx(0:nxsrc-1,ieplt)/gxmax*ymax

   plot, e(0:j), g(ixplt,0:j),title=String(ixplt,format='(I3)'), xtitle='e/T',ytitle='',linestyle=0,yrange=[0,ymax]
   oplot,e(0:j),ga(ixplt,0:j),linestyle=0,color=lu
;   oplot,e(0:j),g0(ixplt,0:j),linestyle=1,color=hong
endfor
device,/close

!p.multi=[0,1,1,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.3
    !P.Thick=10.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='gxone.ps'
device, file=flnm, /color, xsize=21, ysize=15,xoffset=0,yoffset=1 
   ieplt=70 ;1+k*20
   ymax=max(g(0:nxsrc-1,ieplt))

;   gx(0:nxsrc-1,ieplt)=gx(0:nxsrc-1,ieplt)/gxmax*ymax
   s=6.28^1.5
   plot, r(0:nxsrc-1), s*g(0:nxsrc-1,ieplt), xtitle='r/a',ytitle='g',linestyle=0
   oplot,r(0:nxsrc-1),s*ga(0:nxsrc-1,ieplt),linestyle=0,color=lu
;   oplot,r(0:nxsrc-1),g0(0:nxsrc-1,ieplt),linestyle=1,color=hong
device,/close

return
end
