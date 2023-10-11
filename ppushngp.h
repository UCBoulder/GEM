                i=int(xt/dx+0.5)
                j=int(yt/dy+0.5)
                k=int(z2(m,ns)/dz+0.5)-gclr*kcnt

            exp1=exp1 + ex(i,j,k)
            eyp=eyp + ey(i,j,k)
            ezp =ezp + ez(i,j,k)
            delbxp = delbxp+delbx(i,j,k)     
            delbyp = delbyp+delby(i,j,k)
            dpdzp = dpdzp+dpdz(i,j,k)     
            dadzp = dadzp+dadz(i,j,k)     
            aparp = aparp+apar(i,j,k)
            if (ipbm == 1) then
              delbxhp = delbxhp+delbxh(i,j,k)
              delbyhp = delbyhp+delbyh(i,j,k)
              dahdzp = dahdzp+dahdz(i,j,k)
              aparhp = aparhp+aparh(i,j,k)
            end if 

