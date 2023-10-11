pro wt
scale = 1e4
tstart=40.0
tend=100.0
ip = 0
ndim=20000L

dt = 1.0
nskip = 10
rmaj=706.7
a=254.4
csova=0.00278
csovrmaj=csova*a/rmaj
nsubd=10

t0 = fltarr(ndim)
ytmp = fltarr(10,nsubd,ndim)
yav = fltarr(10,ndim)


y0 = fltarr(10,nsubd,ndim)
y1 = fltarr(10,nsubd,ndim)
y2 = fltarr(10,nsubd,ndim)
y3 = fltarr(10,nsubd,ndim)
y4 = fltarr(10,nsubd,ndim)

y5 = fltarr(10,nsubd,ndim)
y5av = fltarr(10,ndim)

y6 = fltarr(10,nsubd,ndim)
y6av = fltarr(10,ndim)

y8 = fltarr(10,nsubd,ndim)
y10 = fltarr(10,nsubd,ndim)


;goto, jump1
n5=5000
ntmp=n5
infl='plot1'
openr, 1, infl
for k=0L,ntmp-1L do begin
   for i=0,0 do begin
      readf,1,t,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9
      t0(k) = t
      ytmp(i,0,k) = x0
      ytmp(i,1,k) = x1
      ytmp(i,2,k) = x2
      ytmp(i,3,k) = x3
      ytmp(i,4,k) = x4
      ytmp(i,5,k) = x5
      ytmp(i,6,k) = x6
      ytmp(i,7,k) = x7
      ytmp(i,8,k) = x8
      ytmp(i,9,k) = x9      
   endfor
endfor
t0 = t0*csovrmaj
close, 1

y5 = ytmp


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

flnm='wt.ps'
device, file=flnm, /color,/encapsulated, xsize=21, ysize=21

xmax = max([max(t0)])
plot, t0(0:n5-1),y5(ip,6,0:n5-1), title='', $
  xstyle=1, ystyle=1,xtitle='!3tC!ls!N/R',$
  ytitle='<|w|>',$ 
  xrange=[0,xmax*1.1], $
  yrange=[0,0.5],linestyle=0,color=0,$
  charsize=2,charthick=4,position=[0.3,0.3,0.9,0.9]


oplot, t0(0:n5-1), y5(ip,8,0:n5-1),linestyle=0, color=lu                ;g1000a, nzgyro=20


device, /close

return
end
