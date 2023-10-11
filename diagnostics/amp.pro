pro amp,n5
scale = 1e4
ip = 1
ndim=20000L

rmaj=1058.60
a=385.62
cs=sqrt(1.0/2.0)
csova=cs/a
csovrmaj=csova*a/rmaj

t0 = fltarr(ndim)
y0 = fltarr(10,ndim)
y1 = fltarr(10,ndim)
y2 = fltarr(10,ndim)
y3 = fltarr(10,ndim)
y4 = fltarr(10,ndim)
y5 = fltarr(10,ndim)
y5im = fltarr(10,ndim)

y6 = fltarr(10,ndim)
y8 = fltarr(10,ndim)
y10 = fltarr(10,ndim)


omeg = 1./sqrt(0.00124)/sqrt(2)/(2.*3.33*1048.2)
print,'omeg= ', omeg
k=1L

;goto, jump1
;n5=1450
infl='./yyre'
openr, 1, infl
for k=0L,n5-1L do begin
readf,1, x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
t0(k) = x0
y5(1,k) = x1
y5im(1,k) = x2
y5(2,k) = x3
y5im(2,k) = x4
y5(3,k) = x5
y5im(3,k) = x6
y5(4,k) = x7
y5im(4,k) = x8
y5(5,k) = x9
y5im(5,k) = x10
endfor
close, 1




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

flnm='amp.ps'
device, file=flnm, /color,/encapsulated, xsize=21, ysize=21

t0 = t0*csovrmaj
xmax = max([max(t0)])
ymax = max([   max(y5(ip,0:n5-1)),max(y5im(ip,0:n5-1))   ])
ymin = min([   min(y5(ip,0:n5-1)),min(y5im(ip,0:n5-1))   ])
plot, t0(0:n5-1),y5(ip,0:n5-1), title='', $
  xstyle=1, ystyle=1,xtitle='!3tC!ls!N/R',$
  ytitle='A',$ 
  xrange=[0,xmax*1.1], $
  yrange=[ymin*1.01,ymax*1.2],linestyle=0,color=0,$
  charsize=2,charthick=4,position=[0.3,0.3,0.9,0.9]

oplot,t0(0:n5-1),y5im(ip,0:n5-1),linestyle=0,color=lu

device, /close

return
end
