pro gamma, n0
ndim=30000L
noff=10
ip=6
dt=3.0
xnplot=10
print,alog(2.71828)
rmaj=1058.60
a=385.62
cs=sqrt(1.0/2.0)
csova=cs/a
csovrmaj=csova*a/rmaj
print,'R/cs=', 1.0/csovrmaj

t0 = fltarr(ndim)
ya = dblarr(10,ndim)
ga = dblarr(ndim)

x1=1d0
x2=1d0
x3=1d0
x4=1d0
x5=1d0
x6=1d0
x7=1d0
x8=1d0
x9=1d0
;n0=2300
infl='flux'
openr, 1, infl
for k=0L,n0-1L do begin
readf,1, x0,x1,x2,x3,x4,x5,x6,x7,x8,x9
t0(k) = x0
ya(1,k) = x1
ya(2,k) = x2
ya(3,k) = x3
ya(4,k) = x4
ya(5,k) = x5
ya(6,k) = x6
ya(7,k) = x7
endfor
close, 1

for i=noff+1,n0-noff-1 do begin
   ga(i)=alog(abs(ya(ip,i+noff)/ya(ip,i-noff)))/(xnplot*dt*2*noff)/2.0
endfor


jump:
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

!p.multi=[0,1,1,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.0
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='gamma.ps'
device, file=flnm, /color,/encapsulated;, xsize=21, ysize=21

t0 = t0*csovrmaj
xmax = max([max(t0)])
ymax = max(ga(noff+1:n0-noff-1)) *0.2
ymin = min(ga(noff+1:n0-noff-1))
plot, t0(noff+1:n0-noff-1),ga(noff+1:n0-noff-1), title='', $
  xstyle=1, ystyle=1,xtitle='tC!ls!N/R',$
  ytitle='gamma',$ 
  xrange=[0,xmax*1.2], $
  yrange=[-ymax*0.5,ymax*0.5],linestyle=0,$
  charsize=2,charthick=4,color=0,position=[0.3,0.3,0.9,0.9]

device, /close

return
end
