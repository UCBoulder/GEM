pro zs,n
  istart=0
  iend=n ;1000
  iskip=floor(float(iend-1-istart)/4.)

  jp1=istart
  jp2=jp1+iskip
  jp3=jp1+iskip*2
  jp4=jp1+iskip*3
  jp5=jp1+iskip*4
  print,'iskip,jp5=', iskip,jp5

k=1L
;n=8
ipg=n-1
rlk=0.45

ndim=n

imx=128L
s=fltarr(imx)
t=fltarr(ndim)
phi=fltarr(2,ndim,imx)
apar=fltarr(2,ndim,imx)
x=fltarr(imx)

phimx=fltarr(2,ndim)
phimn=fltarr(2,ndim)
apamx=fltarr(2,ndim)
apamn=fltarr(2,ndim)
cut=fltarr(ndim)

for k = 0,ndim-1 do begin
   cut(k) = 1.0
endfor

dt = 1.0
nskip = 10
rmaj=706.7
a=254.4
csova=0.00278
csovrmaj=csova*a/rmaj

lx=152.65
dx=lx/imx
rin=0.05
rout=0.85
ds=(rout-rin)/imx
ir=(rlk-rin)/(rout-rin)*imx
for k = 0,imx-1 do begin
   s(k) = rin+k*ds
endfor

infl='zprof'
dumchar=' '
openr, 1, infl
z=0.d
for k = 0,n-1 do begin
   readf,1,dum,x
   t(k) = dum*csovrmaj
   for i = 0,imx-1 do begin
      phi(0,k,i) = x(i)/cut(k)
   endfor

   readf,1,dum,x
;   t(k) = dum
   for i = 0,imx-1 do begin
      apar(0,k,i) = x(i)/cut(k)
   endfor
endfor

close, 1

;relaxation
deltp=fltarr(n,imx)
deltpa=fltarr(imx)
delta=fltarr(imx)
dtn=fltarr(imx)
iskp=10
for k = 0,n-1 do begin
   for i = iskp,imx-iskp-1 do begin
      deltp(k,i)=(apar(0,k,i+iskp)-apar(0,k,i-iskp))/(2*iskp*dx)
   endfor
endfor
istart=n-2
iend=n-1
ndum=iend-istart
for i = 0,imx-1 do begin
   dum = 0.
   dum1 = 0.
   dum2 = 0.
   for k = istart,iend-1 do begin
      dum = dum+deltp(k,i)
      dum1 = dum1+apar(0,k,i)
      dum2 = dum2+phi(0,k,i)            
   endfor
   dum = dum/ndum
   dum1 = dum1/ndum
   dum2 = dum2/ndum   
   deltpa(i) = dum
   delta(i) = dum1
   dtn(i) = dum2   
endfor
print, 'ir, rlk, delta T prime   ', ir, rlk, deltpa(ir)


for k = 0,n-1 do begin
   dum = 0.0
   for i = 1,imx-1 do begin
      dum = dum+phi(0,k,i)*phi(0,k,i)
   endfor
   dum = sqrt(dum/double(imx))

   phimx(0,k) = dum               ;max(phi(k,0:imx-1))
   phimn(0,k) = min(phi(0,k,0:imx-1))

   dum = 0.0
   for i = 1,imx-1 do begin
      dum = dum+apar(0,k,i)*apar(0,k,i)
   endfor
   dum = sqrt(dum/double(imx))

   apamx(0,k) = dum 
   apamn(0,k) = min(apar(0,k,0:imx-1))
endfor



set_plot, 'ps'

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

!p.multi=[0,1,1,0]  ;[0, ncolumn, nrow,0]
    !P.Charthick=5.0
    !P.Charsize=1.5
    !X.Charsize=1.0
    !Y.Charsize=1.3
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='zs.ps'
device, file=flnm, /color, xsize=20, ysize=40 ;,xoffset=6,yoffset=2 

ymin=min(phi(0,ipg,0:imx-1))
ymax=max(phi(0,ipg,0:imx-1))
;plot, s, phi(0,ipg,0:imx-1), linestyle=0,xtitle='r/a',ytitle='!7d!3n',position=[0.3,0.6,0.9,0.9], $
;      yrange=[ymin,ymax],ystyle=1,yticks=2    ;,ytickv=[-0.02,0.05,0.1]

ymin=min(delta(0:imx-1))
ymax=max(delta(0:imx-1))
;plot, s, delta(0:imx-1), linestyle=0,xtitle='r/a',ytitle='<!7d!3T>',position=[0.3,0.6,0.9,0.9], $
;      yrange=[ymin,ymax*1.2],ystyle=1,yticks=2    ;,ytickv=[-0.02,0.05,0.1]

ymin=min(phi)
ymax=max(phi(0,jp1,0:imx-1))
plot, s, phi(0,jp1,0:imx-1), linestyle=0,xtitle='r/a',ytitle='<!7d!3n>',position=[0.3,0.1,0.9,0.4], $
      yrange=[ymin,ymax*1.2],ystyle=1,yticks=2                 ;,ytickv=[-0.1,0.0,0.1]
oplot,s,phi(0,jp2,0:imx-1),linestyle=0,color=cyan
oplot,s,phi(0,jp3,0:imx-1),linestyle=0,color=huang
oplot,s,phi(0,jp4,0:imx-1),linestyle=0,color=lan
oplot,s,phi(0,jp5,0:imx-1),linestyle=0,color=lu

;oplot,s,phi(0,n-1,0:imx-1),linestyle=0,color=magenta

device, /close


!p.multi=[0,1,1,0]
    !P.Charthick=5.0
    !P.Charsize=1.5
    !X.Charsize=1.0
    !Y.Charsize=1.3
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='deltpa.ps'
device, file=flnm, /color, xsize=20, ysize=20 ;,xoffset=6,yoffset=2 

ymin=min(deltpa(0:imx-1))
ymax=max(deltpa(0:imx-1))
plot, s(iskp:imx-iskp-1), deltpa(iskp:imx-iskp-1), linestyle=0,xtitle='r/a',ytitle='d<!7d!3T>/dr',position=[0.3,0.1,0.9,0.4], $
      xrange=[rin,rout],yrange=[ymin,ymax],ystyle=1,yticks=2                 ;,ytickv=[-0.1,0.0,0.1]

device, /close



!p.multi=[0,1,1,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.3
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='zsphit.ps'
device, file=flnm, /color, xsize=21, ysize=21,xoffset=0,yoffset=1 
kdum=0

plot, t, phimx(0,0:n-1),xtitle='t',ytitle='rms(dni)',linestyle=0

close

flnm='zsapat.ps'
device, file=flnm, /color, xsize=21, ysize=21,xoffset=0,yoffset=1 
kdum=0

plot, t, apamx(0,0:n-1),xtitle='t',ytitle='rms(dti)',linestyle=0

close




return
end
