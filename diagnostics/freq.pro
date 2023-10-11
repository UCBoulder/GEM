pro freq
irplot=3

kmx=64
nfreq=400

s=strarr(20)

k=1L
nr=7
f=fltarr(nfreq)
pw=dblarr(kmx,nr,nfreq)

infl='freq'
dumchar=' '
openr, 1, infl
idum=1L
x1=1.0D
x2=1.0D
for k = 0,kmx-1 do begin
   gclr = k
   for i = 0,6 do begin
      readf,1,idum,x1,x2
      readf,1,idum,x1,x2
      readf,1,idum,x1,x2
      readf,1,idum,x1,x2
      readf,1,idum,x1,x2
      readf,1,dumchar

      ir = i
      for j = 0,nfreq-1 do begin
         readf,1,idum,x1,x2
         f(j)=x1
         pw(gclr,ir,j)=x2
      endfor
      readf,1,dumchar
    endfor
endfor

close, 1

for gclr=0,kmx-1 do begin
   for ir=0,nr-1 do begin
;      print,gclr,ir,pw(gclr,ir,100)
   endfor
endfor

set_plot, 'ps'

;Specify the red component of each color:
RED = [0, 1, 1, 0, 0, 1]
;Specify the green component of each color:
GREEN = [0, 1, 0, 1, 0, 1]
;Specify the blue component of each color:
BLUE = [0, 1, 0, 0, 1, 0]
;Load the first six elements of the color table:
TVLCT, 255 * RED, 255 * GREEN, 255 * BLUE

hong=2
lu=3
lan=4
huang=5

!p.multi=[0,2,4,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.0
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='freq.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=2 
;/encapsulated doesn't wrok with yoffset. offset and size in cm

pwmax = max(pw(0:kmx-1,0:6,0:nfreq-1))
for i=0,7 do begin
   gclr=i*8+4 ;kmx/2-5+i
    plot, f(0:nfreq-1),pw(gclr,irplot,0:nfreq-1)/pwmax, title='gclr='+String(gclr,format='(I3)'), $
      xstyle=1, ystyle=1,xtitle='omega', linestyle=0
endfor

device, /close

!p.multi=[0,2,4,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.0
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='freqr.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=2 
;/encapsulated doesn't wrok with yoffset. offset and size in cm

gclr=32                         ;kmx/2-5+i
pwmax = max(pw(gclr,0:6,0:nfreq-1))
for i=0,6 do begin
   plot, f(0:nfreq-1),pw(gclr,i,0:nfreq-1)/pwmax, title='r='+String(i,format='(I1)'), $
      xstyle=1, ystyle=1,xtitle='omega', linestyle=0
endfor

device, /close

return
end
