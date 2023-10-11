pro qidata

  npts = 10
  gammah=[0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
  qi=[0.107, 0.33, 0.47, 0.81, 1.31, 0.96, 1.02, 1.30, 0.76, 0.71]
  var=[0.0086, 0.067, 0.147, 0.39, 0.48, 0.26, 0.19, 0.47, 0.35, 0.35]




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

flnm='qidata.ps'
device, file=flnm, /color,/encapsulated, xsize=21, ysize=21

plot, gammah,qi, title='', $
  xstyle=1, ystyle=1,xtitle='gamma_H',xrange=[0,12],yrange=[0,1.8], $
  ytitle='Q!li',linestyle=0,color=0, $
  charsize=2,charthick=4,position=[0.3,0.3,0.9,0.9]

for k = 0,npts-1 do begin
   oplot,[gammah(k),gammah(k)],[qi(k)-var(k)/2, qi(k)+var(k)/2],linestyle=0
endfor

device, /close

return
end
