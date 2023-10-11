                i=int(xt/dx)
                j=int(yt/dy)
                k=0 !int(z3e(m)/dz)-gclr*kcnt

                wx0=float(i+1)-xt/dx 
                wx1=1.-wx0
                wy0=float(j+1)-yt/dy
                wy1=1.-wy0
                wz0=float(gclr*kcnt+k+1)-z3e(m)/dz
                wz1=1.-wz0
                x000=wx0*wy0*wz0
                x001=wx0*wy0*wz1
                x010=wx0*wy1*wz0
                x011=wx0*wy1*wz1
                x100=wx1*wy0*wz0
                x101=wx1*wy0*wz1
                x110=wx1*wy1*wz0
                x111=wx1*wy1*wz1     
                
                      
            aparhp = x000*apar(i,j,k) + x100*apar(i+1,j,k)  &
           + x010*apar(i,j+1,k) + x110*apar(i+1,j+1,k) + &
           x001*apar(i,j,k+1) + x101*apar(i+1,j,k+1) +   &
           x011*apar(i,j+1,k+1) + x111*apar(i+1,j+1,k+1)
          
