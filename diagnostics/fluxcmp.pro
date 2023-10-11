pro fluxcmp
scale = 1e4
ip = 3
ndim=20000L

rmaj=706.7
a=254.4
csova=0.00278
csovrmaj=csova*a/rmaj
nsubd=8

t0 = fltarr(ndim)
ytmp = fltarr(10,nsubd,ndim)
yav = fltarr(10,ndim)


y0 = fltarr(10,nsubd,ndim)
y0av = fltarr(10,ndim)
y1 = fltarr(10,nsubd,ndim)
y2 = fltarr(10,nsubd,ndim)
y3 = fltarr(10,nsubd,ndim)
y4 = fltarr(10,nsubd,ndim)

y5 = fltarr(10,nsubd,ndim)
y5av = fltarr(10,ndim)

y6 = fltarr(10,nsubd,ndim)
y6av = fltarr(10,ndim)

y8 = fltarr(10,nsubd,ndim)
y8av = fltarr(10,ndim)

y10 = fltarr(10,nsubd,ndim)


omeg = 1./sqrt(0.00124)/sqrt(2)/(2.*3.33*1048.2)
print,'omeg= ', omeg
k=1L

;goto, jump1
n5=15000
ntmp=n5
infl='fluxa'
openr, 1, infl
for k=0L,ntmp-1L do begin
   for i=0,3 do begin
      readf,1,t,x0,x1,x2,x3,x4,x5,x6,x7
      t0(k) = t
      ytmp(i,0,k) = x0
      ytmp(i,1,k) = x1
      ytmp(i,2,k) = x2
      ytmp(i,3,k) = x3
      ytmp(i,4,k) = x4
      ytmp(i,5,k) = x5
      ytmp(i,6,k) = x6
      ytmp(i,7,k) = x7      
   endfor
endfor
close, 1

for k=0,ntmp-1 do begin
   yav(ip,k) = (ytmp(ip,2,k)+ytmp(ip,3,k)+ytmp(ip,4,k)+ytmp(ip,5,k))/4
endfor
y5 = ytmp
y5av = yav

n6=n5
ntmp=n6
infl='../x10-4x/fluxa'
openr, 1, infl
for k=0L,ntmp-1L do begin
   for i=0,3 do begin
      readf,1,t,x0,x1,x2,x3,x4,x5,x6,x7
      t0(k) = t
      ytmp(i,0,k) = x0
      ytmp(i,1,k) = x1
      ytmp(i,2,k) = x2
      ytmp(i,3,k) = x3
      ytmp(i,4,k) = x4
      ytmp(i,5,k) = x5
      ytmp(i,6,k) = x6
      ytmp(i,7,k) = x7      
   endfor
endfor
close, 1

for k=0,ntmp-1 do begin
   yav(ip,k) = (ytmp(ip,6,k)+ytmp(ip,3,k)+ytmp(ip,4,k)+ytmp(ip,5,k))/4
endfor
y6 = ytmp
y6av = yav

n8=n5
ntmp=n8
infl='../x10-2n/fluxa'
openr, 1, infl
for k=0L,ntmp-1L do begin
   for i=0,3 do begin
      readf,1,t,x0,x1,x2,x3,x4,x5,x6,x7
      t0(k) = t
      ytmp(i,0,k) = x0
      ytmp(i,1,k) = x1
      ytmp(i,2,k) = x2
      ytmp(i,3,k) = x3
      ytmp(i,4,k) = x4
      ytmp(i,5,k) = x5
      ytmp(i,6,k) = x6
      ytmp(i,7,k) = x7      
   endfor
endfor
close, 1

for k=0,ntmp-1 do begin
   yav(ip,k) = (ytmp(ip,6,k)+ytmp(ip,3,k)+ytmp(ip,4,k)+ytmp(ip,5,k))/4
endfor
y8 = ytmp
y8av = yav




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


flnm='fluxcmp.ps'
device, file=flnm, /color,/encapsulated, xsize=21, ysize=21

t0 = t0*csovrmaj
xmax = max([max(t0)])
ymax = max([   max(y5(0:3,0:nsubd-1,0:n5-1))   ])
ymin = min([   min(y5(0:3,0:nsubd-1,0:n5-1))   ])

plot, t0(0:n5-1),y5av(ip,0:n5-1), title='', $
  xstyle=1, ystyle=1,xtitle='tC!ls!N/R',$
  ytitle='Q!li',$ 
  xrange=[0,xmax*1.1], $
  yrange=[0,ymax*1.2],linestyle=0,color=0,$
  charsize=2,charthick=4,position=[0.3,0.3,0.9,0.9]

oplot, t0(0:n6-1), y6av(ip,0:n6-1),linestyle=0, color=lan
oplot, t0(0:n8-1), y8av(ip,0:n8-1),linestyle=0, color=lu

;xyouts,21,24,'!3GTC',charsize=1.4,color=lan
;xyouts,61,23,'!3ORB5',charsize=1.4,color=0
;xyouts,96,20,'!3TOY',charsize=1.4,color=lu

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

return
end
